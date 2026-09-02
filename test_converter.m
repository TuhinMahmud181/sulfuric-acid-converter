function tests = test_converter
%TEST_CONVERTER  Regression and invariant tests for the converter model.
%
%   Run with:  results = runtests('test_converter.m')
%
%   The tests guard the physics against refactoring. Two kinds are present:
%   invariants, which must hold for any correct implementation, and
%   regressions, which pin the sizing at the design basis so that a change
%   in the model cannot move the design without a test failing.
tests = functiontests(localfunctions);
end


%% ======================= REGRESSION ====================================

function testDesignSizing(tc)
% The design basis must size to the depths the plant is built on.
b = ConverterModel.basis();
r = ConverterModel.sizeAtTargets(b.ff, b.Tin, b);
tc.verifyEqual(r.L, [0.513 1.214 1.298 1.499], 'AbsTol', 0.005);
tc.verifyEqual(r.Ltot, 4.524,  'RelTol', 0.005);
tc.verifyEqual(r.Wtot, 17588,  'RelTol', 0.005);
tc.verifyEqual(r.Xtot, 0.99875,'AbsTol', 1e-5);
end


function testDesignTemperatures(tc)
b = ConverterModel.basis();
r = ConverterModel.sizeAtTargets(b.ff, b.Tin, b);
tc.verifyEqual(r.Tout, [616.4 496.1 447.1 426.0], 'AbsTol', 1.0);
end


%% ======================= ENERGY AND MASS ===============================

function testAdiabaticRelationClosesEnergyBalance(tc)
% The temperature returned by adiabaticT must satisfy the energy balance it
% is derived from, with the heat of reaction taken at the outlet
% temperature:
%
%     mdot cp (T - T0) = nSO2 (-dH(T)) (f - f0),   dH(T) = -dHa + dHb T
%
% Residual form, checked at several conversions.
b    = ConverterModel.basis();
s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
T0   = b.Tin(1) + b.Tref;

for fend = [0.10 0.35 0.60]
    T        = ConverterModel.adiabaticT(fend, 0, T0, s, mdot, b);
    sensible = mdot*b.cp*(T - T0);
    released = s.SO2*(b.dHa - b.dHb*T)*fend;
    tc.verifyEqual(sensible, released, 'RelTol', 1e-10);
end
end


function testAdiabaticSlopeAtZeroConversion(tc)
% At the bed inlet the balance reduces to dT/df = (a - b T0), which fixes
% the initial slope of the temperature profile independently of the
% closed-form rearrangement.
b    = ConverterModel.basis();
s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
T0   = b.Tin(1) + b.Tref;

h       = 1e-7;
slopeFD = (ConverterModel.adiabaticT(h, 0, T0, s, mdot, b) - T0)/h;
a       = s.SO2*b.dHa/(mdot*b.cp);
bc      = s.SO2*b.dHb/(mdot*b.cp);
tc.verifyEqual(slopeFD, a - bc*T0, 'RelTol', 1e-5);
end


function testStreamMassMatchesMoles(tc)
b = ConverterModel.basis();
m = ConverterModel.streamMass(b.feed, b.MW);
tc.verifyEqual(m, 5576.5, 'RelTol', 1e-3);
end


function testAbsorberConservesSulfurAndOxygen(tc)
% SO3 removal must not alter the SO2 or O2 carried forward.
b  = ConverterModel.basis();
X3 = 0.975;
s  = ConverterModel.absorberOutlet(X3, b);
tc.verifyEqual(s.SO2, b.feed.SO2*(1-X3),               'RelTol', 1e-12);
tc.verifyEqual(s.O2,  b.feed.O2 - b.feed.SO2*X3/2,     'RelTol', 1e-12);
tc.verifyEqual(s.SO3, b.feed.SO2*X3*(1-b.etaAbs),      'RelTol', 1e-12);
tc.verifyEqual(s.N2,  b.feed.N2,                       'RelTol', 1e-12);
end


%% ======================= CEILING AND MONOTONICITY ======================

function testRateVanishesAtCeiling(tc)
% The rate must be zero at fstar and positive just below it.
b     = ConverterModel.basis();
s     = b.feed;
mdot  = ConverterModel.streamMass(s, b.MW);
T0    = b.Tin(1) + b.Tref;
fstar = ConverterModel.kineticCeiling(b.P(1), T0, 0, s, mdot, b);

Tstar = ConverterModel.adiabaticT(fstar, 0, T0, s, mdot, b);
tc.verifyEqual(ConverterModel.rateSO2(fstar, Tstar, b.P(1), s, b), 0, 'AbsTol', 1e-12);

fb = fstar - 1e-4;
Tb = ConverterModel.adiabaticT(fb, 0, T0, s, mdot, b);
tc.verifyGreaterThan(ConverterModel.rateSO2(fb, Tb, b.P(1), s, b), 0);
end


function testTargetBeyondCeilingIsRejected(tc)
% A target at or beyond the ceiling has no solution and must raise an
% error rather than return a depth.
b    = ConverterModel.basis();
s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
T0   = b.Tin(1) + b.Tref;
tc.verifyError(@() ConverterModel.bedMarch(b.P(1), T0, 0, 0.80, s, mdot, b), ...
               'ConverterModel:pastCeiling');
end


function testDepthIncreasesWithTarget(tc)
% Depth must rise monotonically with the conversion target.
b    = ConverterModel.basis();
s    = b.feed;
mdot = ConverterModel.streamMass(s, b.MW);
T0   = b.Tin(1) + b.Tref;
L = arrayfun(@(ff) lastZ(b.P(1), T0, 0, ff, s, mdot, b), [0.55 0.62 0.66 0.69]);
tc.verifyTrue(all(diff(L) > 0));
end


function testDepthIncreasesWithApproachFraction(tc)
% The same monotonicity in the variable the optimisation uses.
b = ConverterModel.basis();
L = arrayfun(@(p) ConverterModel.runByApproach([p p p p], b.Tin, b).Ltot, ...
             [0.80 0.90 0.95 0.98]);
tc.verifyTrue(all(diff(L) > 0));
end


%% ======================= NUMERICS ======================================

function testGridConverged(tc)
% The conversion grid must be fine enough that the sizing does not move
% appreciably on refinement.
b  = ConverterModel.basis();
r1 = ConverterModel.sizeAtTargets(b.ff, b.Tin, b);
b2 = b; b2.Npts = 4*b.Npts;
r2 = ConverterModel.sizeAtTargets(b2.ff, b2.Tin, b2);
tc.verifyEqual(r2.Ltot, r1.Ltot, 'RelTol', 0.002);
end


function testPressureDropSmallEnoughForIsobaricBeds(tc)
% The isobaric assumption within a bed requires the drop to be a small
% fraction of the absolute pressure.
b = ConverterModel.basis();
r = ConverterModel.sizeAtTargets(b.ff, b.Tin, b);
tc.verifyLessThan(r.dPtot/(b.P(1)*1e5), 0.05);
end


function testApproachEquivalenceRoundTrip(tc)
% Converting conversion targets to approach fractions and back must
% reproduce the same bed depths.
b   = ConverterModel.basis();
phi = ConverterModel.phiFromTargets(b.ff, b.Tin, b);
rA  = ConverterModel.runByApproach(phi, b.Tin, b);
rB  = ConverterModel.sizeAtTargets(b.ff, b.Tin, b);
tc.verifyEqual(rA.L, rB.L, 'RelTol', 1e-6);
end


%% ======================= HELPERS =======================================

function z = lastZ(P, T0, f0, ff, s, mdot, b)
[~, ~, zz] = ConverterModel.bedMarch(P, T0, f0, ff, s, mdot, b);
z = zz(end);
end
