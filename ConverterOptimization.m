function out = ConverterOptimization(b, opts)
% CONVERTEROPTIMIZATION  Minimum catalyst inventory for the four-bed SO2
%   converter at a fixed overall conversion.
%
%   out = ConverterOptimization(b)
%   out = ConverterOptimization(b, opts)
%
%   b     basis struct from ConverterModel.basis
%   opts  optional struct, fields
%           phi0   1 x 4 starting approach fractions (default: from b.ff)
%           Tin0   1 x 4 starting inlet temperatures, degC (default b.Tin)
%           lbT    scalar or 1 x 4 lower bound on inlet temperature, degC
%           ubT    scalar or 1 x 4 upper bound on inlet temperature, degC
%                  (default 480; the beds settle near 460, so a bound at
%                  460 clips the solution by about half a degree)
%           phiMax scalar cap on the approach fraction (default 0.999);
%                  a numerical guard only, and it should not be active at
%                  the solution
%           display fmincon Display setting (default 'final'; use 'iter'
%                  to print every iteration, 'off' for silence). The
%                  iteration history is retained in out.history either way.
%
%   out   struct with fields x, fval, phi, Tin, result, start, exitflag,
%         output, lambda, history, Tmargin, Tlimited, active,
%         usedBestFeasible
%
%
%   PROBLEM
%   -------
%     variables    x = [phi1 phi2 phi3 phi4  Tin1 Tin2 Tin3 Tin4]
%     minimise     total catalyst volume,  Ac sum(L)
%     subject to   overall conversion       >= b.Xspec
%                  bed outlet temperature   <= b.Tmax
%                  every bed depth          >= b.Lmin
%                  total bed pressure drop  <= b.dPmax
%                  net rate at every bed inlet > 0
%                  bounds on phi and on the inlet temperatures
%
%
%   WHY APPROACH FRACTIONS AND NOT CONVERSIONS
%   ------------------------------------------
%   Bed depth follows from dz = (CA0 u0 / rhoB rm) df, and rm falls to zero
%   at a finite ceiling conversion on the adiabatic path. If the bed outlet
%   conversions are taken as free variables, the search can request targets
%   beyond that ceiling, where rm < 0 and the depth integral turns negative.
%   The objective then decreases as the design becomes more infeasible, and
%   the solver is drawn into a region with no physical meaning.
%
%   Specifying each bed as a fraction of its available approach,
%
%       f_out = f_in + phi (fstar - f_in),      0 < phi < 1,
%
%   confines every candidate to rm > 0. The ceiling fstar is recomputed for
%   each candidate, since it depends on that bed's inlet temperature and on
%   the composition reaching it.
%
%
%   ON THE BOUNDS
%   -------------
%   phiMax exists only to keep the depth integral finite, since depth
%   diverges as the ceiling is approached. It is a numerical guard and must
%   not be the thing that decides the answer. A tight value binds on several
%   beds at once and the optimum then reports a property of that number
%   rather than of the process, so the default is loose. checkBounds tests
%   its multiplier and out.active.phiAtMax records the outcome.
%
%   The minimum bed depth b.Lmin is the constraint that shapes the shallow
%   end of the design. It has a physical basis and costs little: at 0.5 m it
%   adds a few per cent to the catalyst inventory and holds the peak
%   temperature clear of the catalyst limit.
%
%   TWO REMARKS ON THE SOLUTION
%   ---------------------------
%   The temperature limit tends towards active on bed 1: a hotter inlet
%   raises the rate and shortens the bed, and the catalyst thermal limit is
%   what stops it. Check the margin reported at the end. A design sitting
%   exactly on the limit has no operating headroom, and sintering rates
%   climb steeply above it.
%
%   Catalyst volume alone is an incomplete objective. It omits the interbed
%   exchanger area implied by the cooling duties, and hotter inlets raise
%   those duties while shortening the beds. An annualised-cost objective
%   would move the optimum towards cooler inlets and longer beds.

if nargin < 2, opts = struct(); end
if ~isfield(opts,'Tin0'),    opts.Tin0    = b.Tin; end
if ~isfield(opts,'phi0'),    opts.phi0    = ConverterModel.phiFromTargets(b.ff, opts.Tin0, b); end
if ~isfield(opts,'lbT'),     opts.lbT     = 400; end
if ~isfield(opts,'ubT'),     opts.ubT     = 480; end
if ~isfield(opts,'phiMax'),  opts.phiMax  = 0.999; end
if ~isfield(opts,'display'), opts.display = 'final'; end

x0 = [opts.phi0(:).' opts.Tin0(:).'];

% phi is capped short of 1 because the depth integral diverges as the
% ceiling is approached. The cap is deliberately loose: it is a numerical
% guard, not a design constraint, and it should not be active at the
% solution. checkBounds below reports it if it is.
lb = [0.30 0.30 0.30 0.30, opts.lbT.*ones(1,4)];
ub = [opts.phiMax*ones(1,4), opts.ubT.*ones(1,4)];

% Central differences with a step well above the tolerance of the internal
% root find in kineticCeiling, so gradients are not dominated by its noise.
fopts = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'Display', opts.display, ...
    'FiniteDifferenceType','central', ...
    'FiniteDifferenceStepSize',1e-6, ...
    'OptimalityTolerance',1e-6, ...
    'ConstraintTolerance',1e-8, ...
    'MaxFunctionEvaluations',4000, ...
    'OutputFcn',@recordHistory);

% OptimalityTolerance is set near the resolution of a finite-difference
% gradient on this objective. Central differences with a step of 1e-6 on a
% quantity of order 10 leave a roundoff floor of about eps*|f|/h in each
% gradient component, and the constraint scaling below multiplies into the
% Lagrangian gradient. A tolerance far below that floor cannot be met, and
% the solver would terminate on step size with the same answer.

cachedRun([], []);                   % clear the evaluation cache

recordHistory([], [], 'reset');
[x, fval, exitflag, output, lambda] = fmincon(@(x) objective(x,b), x0, ...
    [], [], [], [], lb, ub, @(x) constraints(x,b), fopts);
out.history = recordHistory([], [], 'fetch');

% fmincon can terminate at a point that is not the best feasible point it
% visited, and says so through output.bestfeasible. The returned x is the
% last iterate, not necessarily the best one. Take the better of the two.
out.usedBestFeasible = false;
if isfield(output,'bestfeasible') && ~isempty(output.bestfeasible)
    bf = output.bestfeasible;
    if bf.fval < fval
        x    = bf.x;
        fval = bf.fval;
        out.usedBestFeasible = true;
    end
end

out.x        = x;
out.fval     = fval;
out.phi      = x(1:4);
out.Tin      = x(5:8);
out.result   = ConverterModel.runByApproach(out.phi, out.Tin, b);
out.start    = ConverterModel.runByApproach(opts.phi0, opts.Tin0, b);
out.exitflag = exitflag;
out.output   = output;
out.Tlimited = abs(max(out.result.Tout) - b.Tmax) < 0.5;
out.Tmargin  = b.Tmax - max(out.result.Tout);
out.lambda   = lambda;
out.active   = checkBounds(lambda, out.result, b);
end


function a = checkBounds(lambda, r, b)
%CHECKBOUNDS  Report which bounds and constraints hold the solution.
%
%   A bound is active when its Lagrange multiplier is non-zero, which is
%   what it means for that bound to be carrying load. Testing the distance
%   from the bound instead cannot separate a variable held there from one
%   whose optimum happens to lie nearby, and the two call for different
%   responses: the first means the bound decided the answer, the second
%   means it did not.
tol = 1e-8;
a.phiAtMax   = find(lambda.upper(1:4) > tol);
a.phiAtMin   = find(lambda.lower(1:4) > tol);
a.TinAtMax   = find(lambda.upper(5:8) > tol);
a.TinAtMin   = find(lambda.lower(5:8) > tol);
a.depthAtMin = find(r.L <= b.Lmin + 1e-4);
a.Tlimit     = find(r.Tout >= b.Tmax - 0.5);
a.anyPhiCap  = ~isempty(a.phiAtMax);
end


function stop = recordHistory(~, optimValues, state)
%RECORDHISTORY  Retain the iteration history so convergence can be examined
%   after the solve without printing every iteration to the console.
%
%   Called by fmincon through OutputFcn. Two extra states are accepted from
%   the caller: 'reset' to clear the store before a solve and 'fetch' to
%   return it afterwards.
persistent H
stop = false;
switch state
    case 'reset'
        H = struct('iteration',[],'fval',[],'feasibility',[],'optimality',[], ...
                   'stepsize',[]);
        return
    case 'fetch'
        stop = H;
        return
    case 'init'
        if isempty(H)
            H = struct('iteration',[],'fval',[],'feasibility',[], ...
                       'optimality',[],'stepsize',[]);
        end
    case 'iter'
        H.iteration(end+1)   = optimValues.iteration;
        H.fval(end+1)        = optimValues.fval;
        H.feasibility(end+1) = optimValues.constrviolation;
        H.optimality(end+1)  = optimValues.firstorderopt;
        if isfield(optimValues,'stepsize') && ~isempty(optimValues.stepsize)
            H.stepsize(end+1) = optimValues.stepsize;
        else
            H.stepsize(end+1) = NaN;
        end
end
end


%% ======================= LOCAL FUNCTIONS ==============================

function V = objective(x, b)
%OBJECTIVE  Total catalyst volume, m3.
r = cachedRun(x, b);
if ~r.feasible
    V = 1e3;                         % finite and large: keeps the search bounded
    return
end
V = b.Ac*r.Ltot;
end


function [c, ceq] = constraints(x, b)
%CONSTRAINTS  Inequality constraints, scaled to comparable magnitude.
%
%   Scaling matters. Conversion residuals are of order 1e-3 while
%   temperature residuals are of order 1e2. Left unscaled, the solver's
%   constraint tolerance is satisfied by the conversion constraint long
%   before that constraint means anything.

r = cachedRun(x, b);

if ~r.feasible
    % No forward driving force somewhere in the train. Report it through
    % the inlet-rate constraint so the solver is pushed back out.
    d = r.driving(:);
    d(isnan(d)) = -1;
    c   = [1; -1e6*d];
    ceq = [];
    return
end

c = [ 1e3*(b.Xspec - r.Xtot);        % overall conversion specification
      (r.Tout(:) - b.Tmax)/10;       % catalyst thermal limit
      10*(b.Lmin - r.L(:));          % minimum buildable bed depth
      (r.dPtot - b.dPmax)/1e3;       % hydraulic limit
      -1e6*r.driving(:) ];           % forward driving force at each inlet
ceq = [];
end


function r = cachedRun(x, b)
%CACHEDRUN  fmincon evaluates the objective and the constraints separately
%   at the same point; the bed march is the expensive part, so the last
%   evaluation is retained. Call with empty arguments to clear.
persistent xLast rLast
if isempty(x)
    xLast = []; rLast = []; r = [];
    return
end
if ~isempty(xLast) && isequal(x, xLast)
    r = rLast; return
end
r = ConverterModel.runByApproach(x(1:4), x(5:8), b);
xLast = x; rLast = r;
end
