classdef ConverterModel
% CONVERTERMODEL  Model of a four-bed adiabatic SO2 converter.
%
%   150 MTPD sulfuric acid plant, double contact / double absorption.
%   Sulfur burning, V2O5 on silica gel, four adiabatic beds with interbed
%   cooling and interpass absorption of SO3 between beds 3 and 4.
%
%   All methods are static. Call them as ConverterModel.<name>(...).
%
%
%   GOVERNING EQUATIONS
%   -------------------
%   Each bed is treated as an adiabatic, pseudo-homogeneous, isobaric plug
%   flow reactor. Over a slice of thickness dz and cross-section Ac:
%
%       SO2 balance     Ac dz rhoB rm = -d(F_SO2)
%                       => (CA0 u0 / rhoB) df/dz = rm
%
%       energy balance  G cp dT = rm rhoB (-dH) dz
%
%   Eliminating dz between the two and integrating with the linear heat of
%   reaction dH(T) = -dHa + dHb T gives
%
%       T = (T0 + a (f - f0)) / (1 + b (f - f0))
%       a = CA0 u0 dHa / (G cp)        b = CA0 u0 dHb / (G cp)
%
%   This relation holds in KELVIN and in no other scale, because dH is
%   linear in absolute temperature.
%
%   Bed depth then follows from
%
%       dz = alpha df,      alpha = CA0 u0 / (rhoB rm)
%
%   integrated by the trapezoidal rule on a uniform grid in conversion.
%
%
%   KINETICS
%   --------
%       rm = (k1 pSO2 pO2 - k2 pSO3 sqrt(pO2)) / sqrt(pSO2)
%       ln k1 = 12.07 - 31000/(Rc T)
%       ln k2 = 22.75 -    E2/(Rc T)          Rc = 1.987 cal/(mol K)
%
%   k1 and k2 are referenced to atmospheres, so partial pressures are
%   formed on the same basis.
%
%
%   CEILING CONVERSION
%   ------------------
%   As conversion proceeds adiabatically the gas heats up, which lowers the
%   attainable conversion, while conversion itself raises pSO3. Together
%   these drive rm to zero at a finite conversion fstar. Since alpha
%   diverges there, no bed depth can reach or exceed fstar. kineticCeiling
%   locates it from the rate expression itself rather than from a separate
%   equilibrium correlation: the two do not coincide exactly, and a target
%   placed using the correlation can fall where the rate expression is
%   already negative, which makes the depth integral meaningless.

methods (Static)

%% ===================== BASIS ==========================================
function b = basis()
%BASIS  Design basis. Single definition, read by everything else.

% --- operating ---------------------------------------------------------
b.P    = [1.20 1.20 1.10 1.10];      % bar (abs), one per bed
b.Tin  = [417  427  437  417 ];      % degC, bed inlet temperature
b.ff   = [0.70 0.94 0.975 0.95];     % conversion target per bed
                                     %   beds 1-3: cumulative
                                     %   bed 4   : local, on residual SO2
b.etaAbs = 0.96;                     % interpass absorber SO3 removal

% --- geometry and catalyst ---------------------------------------------
b.D    = 3.0;                        % m,      bed diameter
b.rhoB = 550;                        % kg/m3,  bulk catalyst density
b.rhoP = 1350;                       % kg/m3,  particle density
b.dp   = 7.2e-3;                     % m,      surface-volume dia 6 Vp/Ap
b.Ac   = 0.25*pi*b.D^2;              % m2
b.eps  = 1 - b.rhoB/b.rhoP;          % voidage

% --- gas ---------------------------------------------------------------
b.cp   = 0.25*4.184;                 % J/(g K)
b.mu   = 3.7e-5;                     % Pa s, near 750 K

% --- kinetics ----------------------------------------------------------
b.E2   = 53600;                      % cal/mol, reverse activation energy
b.Rc   = 1.987;                      % cal/(mol K)

% --- thermochemistry ---------------------------------------------------
b.dHa  = 24.60*4184;                 % J/mol,      dH(T) = -dHa + dHb T
b.dHb  = 1.99e-3*4184;               % J/(mol K)
b.Tref = 273;                        % K, degC -> K offset

% --- limits, used by the optimisation ----------------------------------
b.Tmax  = 630;                       % degC, catalyst thermal limit
b.Xspec = 0.99875;                   % required overall SO2 conversion
b.dPmax = 15e3;                      % Pa, allowable total bed dP
b.Lmin  = 0.50;                      % m, minimum bed depth
% A shallow bed over a 3 m diameter distributes gas poorly: at L/D below
% about 0.15 the support grid, the quartz hold-down layer and the entrance
% region occupy a large fraction of the bed, and channelling makes the plug
% flow assumption untenable. Without this floor the optimiser drives bed 1
% to roughly a quarter of a metre, which is not a buildable bed.

% --- feed to bed 1, mol/s ----------------------------------------------
b.feed = struct('SO2',17.396,'O2',19.136,'N2',137.439,'CO2',0.002,'SO3',0);

% --- molar masses, g/mol -----------------------------------------------
b.MW = struct('SO2',64.06,'O2',32.00,'N2',28.01,'CO2',44.01,'SO3',80.06);

b.Npts = 1000;                        % grid points per bed
end


%% ===================== STREAM ========================================
function m = streamMass(s, MW)
%STREAMMASS  Total mass flow, g/s, from the molar flow vector.
m = s.SO2*MW.SO2 + s.O2*MW.O2 + s.N2*MW.N2 + s.CO2*MW.CO2 + s.SO3*MW.SO3;
end


function s = absorberOutlet(X3, b)
%ABSORBEROUTLET  Bed 4 inlet stream after interpass absorption.
%   SO3 formed up to cumulative conversion X3 is removed with efficiency
%   b.etaAbs, lowering pSO3 and reopening the driving force for the SO2
%   that remains.
s.SO2 = b.feed.SO2*(1 - X3);
s.O2  = b.feed.O2 - b.feed.SO2*X3/2;
s.SO3 = b.feed.SO2*X3*(1 - b.etaAbs);
s.N2  = b.feed.N2;
s.CO2 = b.feed.CO2;
end


%% ===================== BED PHYSICS ===================================
function [T, a, bcoef] = adiabaticT(f, f0, T0, s, mdot, b)
%ADIABATICT  Gas temperature along an adiabatic bed, K.
%   T0 and T are in kelvin; see the class header for the derivation.
G     = mdot/b.Ac;                   % g/(m2 s)
CA0u0 = s.SO2/b.Ac;                  % mol/(m2 s)
a     = (CA0u0/(G*b.cp))*b.dHa;      % K
bcoef = (CA0u0/(G*b.cp))*b.dHb;      % dimensionless
T     = (T0 + a*(f-f0))./(1 + bcoef*(f-f0));
end


function [rm, pSO2, pO2, pSO3] = rateSO2(f, T, P, s, b)
%RATESO2  Net oxidation rate of SO2, mol/(g catalyst . s).
%   f is conversion on the basis of s.SO2, the SO2 entering the bed.
%   The sign of rm is returned intact: a negative value means the gas has
%   passed the composition at which forward and reverse steps balance, and
%   the caller must decide what to do rather than have it hidden.
nSO2 = s.SO2*(1-f);
nO2  = s.O2 - s.SO2*f/2;
nSO3 = s.SO3 + s.SO2*f;
nTot = s.N2 + s.CO2 + nSO2 + nO2 + nSO3;

pSO2 = P*nSO2./nTot;                 % atm
pO2  = P*nO2 ./nTot;
pSO3 = P*nSO3./nTot;

k1 = exp(12.07 - 31000./(b.Rc*T));   % mol/(g s atm^1.5)
k2 = exp(22.75 - b.E2 ./(b.Rc*T));   % mol/(g s atm)

rm = (k1.*pSO2.*pO2 - k2.*pSO3.*sqrt(pO2))./sqrt(pSO2);
end


function fstar = kineticCeiling(P, T0, f0, s, mdot, b)
%KINETICCEILING  Conversion at which the net rate reaches zero along the
%   adiabatic path from (f0, T0). Returns NaN if there is no forward
%   driving force at the bed inlet.
    function r = rateAt(x)
        T = ConverterModel.adiabaticT(x, f0, T0, s, mdot, b);
        r = ConverterModel.rateSO2(x, T, P, s, b);
    end
if rateAt(f0) <= 0
    fstar = NaN; return
end
hi = 1 - 1e-9;
if rateAt(hi) > 0
    fstar = hi;                      % rate positive to complete conversion
else
    fstar = fzero(@rateAt, [f0 + 1e-12, hi]);
end
end


function [f, T, z, W] = bedMarch(P, T0, f0, ff, s, mdot, b, N)
%BEDMARCH  One adiabatic bed, marched in conversion.
%
%   Inputs   P bar(abs) | T0 K | f0,ff conversion on the basis of s.SO2
%            s stream mol/s | mdot g/s | b basis | N grid points
%   Outputs  f conversion | T K | z m cumulative | W kg catalyst
if nargin < 8, N = b.Npts; end

f  = linspace(f0, ff, N);
T  = ConverterModel.adiabaticT(f, f0, T0, s, mdot, b);
rm = ConverterModel.rateSO2(f, T, P, s, b);          % mol/(g s)

if any(rm <= 0)
    k = find(rm <= 0, 1);
    error('ConverterModel:pastCeiling', ...
        ['Net rate reaches zero at f = %.5f (T = %.1f C). The target ' ...
         '%.5f is at or beyond the ceiling for an inlet of %.1f C at ' ...
         '%.2f bar and is unreachable at any bed depth.'], ...
         f(k), T(k)-b.Tref, ff, T0-b.Tref, P);
end

alpha = (s.SO2/b.Ac)./(1000*b.rhoB*rm);              % m (1000 converts kg->g)
dz    = (alpha(1:end-1) + alpha(2:end))*(f(2)-f(1))/2;
z     = [0 cumsum(dz)];                              % m
W     = b.rhoB*b.Ac*z(end);                          % kg
end


function dP = ergun(L, Tavg, P, mdot, s, b)
%ERGUN  Pressure drop across a packed bed, Pa.
%
%       dP/L = 150 mu (1-e)^2 u/(e^3 dp^2) + 1.75 (1-e) rho u^2/(e^3 dp)
%
%   dp is the surface-volume diameter 6 Vp/Ap, which already carries the
%   particle shape and is not multiplied by sphericity again. Gas density
%   is taken at the mean bed temperature from the ideal gas law.
nTot  = s.SO2 + s.O2 + s.N2 + s.CO2 + s.SO3;         % mol/s
Mavg  = (mdot/1000)/nTot;                            % kg/mol
rho_g = P*1e5*Mavg/(8.314*Tavg);                     % kg/m3
u     = (mdot/1000)/(rho_g*b.Ac);                    % m/s, superficial
e     = b.eps;
dP = L*( 150*b.mu*(1-e)^2*u/(e^3*b.dp^2) ...
       + 1.75*(1-e)*rho_g*u^2/(e^3*b.dp) );
end


%% ===================== CONVERTER =====================================
function r = sizeAtTargets(ff, TinC, b)
%SIZEATTARGETS  Size the four beds to absolute conversion targets.
%
%   ff(1:3) cumulative, ff(4) local to bed 4. This is the direct sizing
%   calculation: the conversions are given and the depths come out.
%
%   Output r: L, W, Tout, Xcum, fstar, dP, profiles, Ltot, Wtot, Xtot, dPtot

r.L = zeros(1,4); r.W = zeros(1,4); r.Tout = zeros(1,4);
r.Xcum = zeros(1,4); r.fstar = zeros(1,4); r.dP = zeros(1,4);
r.z = cell(1,4); r.T = cell(1,4); r.f = cell(1,4);

s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
f0   = 0;  X3 = 0;

for j = 1:4
    if j == 4
        X3   = r.Xcum(3);
        s    = ConverterModel.absorberOutlet(X3, b);
        mdot = ConverterModel.streamMass(s, b.MW);
        f0   = 0;                                    % local basis
    end
    T0 = TinC(j) + b.Tref;
    r.fstar(j) = ConverterModel.kineticCeiling(b.P(j), T0, f0, s, mdot, b);
    [f, T, z, W] = ConverterModel.bedMarch(b.P(j), T0, f0, ff(j), s, mdot, b);

    r.L(j)    = z(end);
    r.W(j)    = W;
    r.Tout(j) = T(end) - b.Tref;
    r.dP(j)   = ConverterModel.ergun(z(end), 0.5*(T(1)+T(end)), b.P(j), mdot, s, b);
    r.z{j} = z;  r.T{j} = T - b.Tref;

    if j == 4
        r.f{j}    = X3 + (1 - X3)*f;                 % map onto cumulative
        r.Xcum(j) = X3 + (1 - X3)*ff(j);
    else
        r.f{j}    = f;
        r.Xcum(j) = ff(j);
        f0        = ff(j);
    end
end

r.stream4 = s;  r.mdot4 = mdot;
r.Ltot = sum(r.L);  r.Wtot = sum(r.W);
r.Xtot = r.Xcum(4); r.dPtot = sum(r.dP);
r.feasible = true;
end


function r = runByApproach(phi, TinC, b)
%RUNBYAPPROACH  Size the four beds from approach fractions.
%
%   Each bed outlet is specified as a fraction of the approach available
%   to that bed:
%
%       f_out = f_in + phi (fstar - f_in),      0 < phi < 1
%
%   Two properties follow, and both matter inside an optimiser:
%     - every phi in (0,1) gives a target strictly below the ceiling, so
%       the depth integral is finite and positive;
%     - f_out > f_in for every phi > 0, so conversion is monotone along
%       the train without ordering constraints.
%
%   r.feasible is false if any bed has no forward driving force at its
%   inlet; r.driving then holds the inlet rate for each bed.

r.L = nan(1,4); r.W = nan(1,4); r.Tout = nan(1,4);
r.Xcum = nan(1,4); r.fstar = nan(1,4); r.driving = nan(1,4); r.dP = nan(1,4);
r.z = cell(1,4); r.T = cell(1,4); r.f = cell(1,4);
r.feasible = true;

s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
X    = 0;  X3 = 0;

for j = 1:4
    if j == 4
        X3   = X;
        s    = ConverterModel.absorberOutlet(X3, b);
        mdot = ConverterModel.streamMass(s, b.MW);
        f0   = 0;
    else
        f0 = X;
    end

    T0  = TinC(j) + b.Tref;
    Tf0 = ConverterModel.adiabaticT(f0, f0, T0, s, mdot, b);
    r.driving(j) = ConverterModel.rateSO2(f0, Tf0, b.P(j), s, b);

    fstar = ConverterModel.kineticCeiling(b.P(j), T0, f0, s, mdot, b);
    r.fstar(j) = fstar;
    if isnan(fstar)
        r.feasible = false;
        return
    end

    ftgt = f0 + phi(j)*(fstar - f0);
    [f, T, z, W] = ConverterModel.bedMarch(b.P(j), T0, f0, ftgt, s, mdot, b);

    r.L(j)    = z(end);
    r.W(j)    = W;
    r.Tout(j) = T(end) - b.Tref;
    r.dP(j)   = ConverterModel.ergun(z(end), 0.5*(T(1)+T(end)), b.P(j), mdot, s, b);
    r.z{j} = z;  r.T{j} = T - b.Tref;

    if j == 4
        r.f{j} = X3 + (1 - X3)*f;
        X      = X3 + (1 - X3)*ftgt;
    else
        r.f{j} = f;
        X      = ftgt;
    end
    r.Xcum(j) = X;
end

r.stream4 = s;  r.mdot4 = mdot;
r.Ltot = sum(r.L);  r.Wtot = sum(r.W);
r.Xtot = X;         r.dPtot = sum(r.dP);
end


function phi = phiFromTargets(ff, TinC, b)
%PHIFROMTARGETS  Approach fractions equivalent to a set of conversion
%   targets, so a design specified by conversion can seed the optimiser.
phi  = nan(1,4);
s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
X    = 0;
for j = 1:4
    if j == 4
        s    = ConverterModel.absorberOutlet(X, b);
        mdot = ConverterModel.streamMass(s, b.MW);
        f0   = 0;
    else
        f0 = X;
    end
    T0     = TinC(j) + b.Tref;
    fstar  = ConverterModel.kineticCeiling(b.P(j), T0, f0, s, mdot, b);
    phi(j) = (ff(j) - f0)/(fstar - f0);
    X      = ff(j);
end
end

end % methods
end % classdef
