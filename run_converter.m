% RUN_CONVERTER  Driver for the four-bed SO2 converter calculation.
%
%   Sizes the beds at the design basis, then optionally minimises the
%   catalyst inventory at the same overall conversion and reports the two
%   side by side.
%
%   Requires ConverterModel.m and ConverterOptimization.m on the path.
%   The sizing needs no toolboxes; the optimisation needs fmincon.

clc; clear; close all;

RUN_OPTIMISATION = true;             % set false for sizing only

b = ConverterModel.basis();

%% ======================= SIZING AT THE DESIGN BASIS ====================

des = ConverterModel.sizeAtTargets(b.ff, b.Tin, b);

mdot0 = ConverterModel.streamMass(b.feed, b.MW);
fprintf('DESIGN BASIS\n');
fprintf('  Ac = %.3f m2,  feed %.1f g/s,  G = %.1f g/(m2 s)\n', ...
        b.Ac, mdot0, mdot0/b.Ac);
fprintf('  Absorber removes %.1f%% of SO3 between beds 3 and 4.\n', 100*b.etaAbs);
fprintf('  Bed 4 inlet: SO2 %.4f, O2 %.4f, SO3 %.4f mol/s, %.1f g/s\n', ...
        des.stream4.SO2, des.stream4.O2, des.stream4.SO3, des.mdot4);
if des.stream4.SO3 > des.stream4.SO2
    fprintf(2, ['  WARNING: SO3 exceeds SO2 at the bed 4 inlet. The residual\n' ...
                '           product suppresses the equilibrium shift the\n' ...
                '           second pass exists to obtain.\n']);
end

printBeds('SIZING AT THE DESIGN BASIS', des, b, b.ff);
plotProfiles(des, 'Design basis');

%% ======================= OPTIMISATION ==================================

if ~RUN_OPTIMISATION
    return
end

fprintf('\n');
opt = ConverterOptimization(b);
r   = opt.result;

fprintf('\nexitflag %d after %d iterations\n', opt.exitflag, opt.output.iterations);
if opt.usedBestFeasible
    fprintf(['  The final iterate was not the best feasible point visited;\n' ...
             '  the better one has been taken.\n']);
end
fprintf('\n bed    phi     Tin[C]   Tout[C]   fstar     X_cum      L[m]    dP[kPa]\n');
for j = 1:4
    fprintf('  %d    %.4f   %6.1f   %6.1f   %.5f   %.5f   %6.3f   %6.2f\n', ...
        j, opt.phi(j), opt.Tin(j), r.Tout(j), r.fstar(j), ...
        r.Xcum(j), r.L(j), r.dP(j)/1000);
end

fprintf('\n                      design     optimum\n');
fprintf(' total depth   [m]   %7.3f    %7.3f\n',  des.Ltot,  r.Ltot);
fprintf(' catalyst      [kg]  %7.0f    %7.0f    (%+.1f%%)\n', ...
        des.Wtot, r.Wtot, 100*(r.Wtot - des.Wtot)/des.Wtot);
fprintf(' conversion          %7.5f    %7.5f\n', des.Xtot,  r.Xtot);
fprintf(' peak T        [C]   %7.1f    %7.1f    (limit %.0f)\n', ...
        max(des.Tout), max(r.Tout), b.Tmax);
fprintf(' bed dP        [kPa] %7.2f    %7.2f\n', des.dPtot/1e3, r.dPtot/1e3);

fprintf(' min bed depth [m]   %7.3f    %7.3f    (floor %.2f)\n', ...
        min(des.L), min(r.L), b.Lmin);

plotConvergence(opt.history);

plotProfiles(r, 'Optimum');

figure('Name','Catalyst distribution');
bar([des.L(:) r.L(:)]);
set(gca,'XTickLabel',{'Bed 1','Bed 2','Bed 3','Bed 4'});
ylabel('Bed depth [m]');
legend('design','optimum','Location','northwest');
grid on; title('Catalyst distribution');


%% ======================= LOCAL FUNCTIONS ==============================

function plotConvergence(h)
%PLOTCONVERGENCE  First-order optimality and constraint violation against
%   iteration. Shows whether the solve settled or merely ran out of steps.
if isempty(h) || isempty(h.iteration), return, end
figure('Name','Convergence');
subplot(2,1,1);
semilogy(h.iteration, max(h.optimality,eps), '-o', 'LineWidth',1.3,'MarkerSize',3);
ylabel('first-order optimality'); grid on; title('Convergence');
subplot(2,1,2);
semilogy(h.iteration, max(h.feasibility,eps), '-o', 'LineWidth',1.3,'MarkerSize',3);
ylabel('constraint violation'); xlabel('iteration'); grid on;
end


function printBeds(header, r, b, ff)
%PRINTBEDS  One row per bed. f/fstar shows how close each bed runs to the
%   conversion at which its net rate reaches zero; a bed close to 1 there
%   has a depth that is weakly determined and sensitive to the kinetic
%   parameters.
fprintf('\n%s\n', header);
fprintf(' bed   L[m]    Tin[C]   Tout[C]    X_cum      fstar     f/fstar   dP[kPa]\n');
for j = 1:4
    fprintf('  %d   %6.3f   %6.1f   %6.1f   %8.5f   %7.5f   %.4f   %6.2f\n', ...
        j, r.L(j), b.Tin(j), r.Tout(j), r.Xcum(j), r.fstar(j), ...
        ff(j)/r.fstar(j), r.dP(j)/1000);
end
fprintf('\n Total bed depth     = %.4f m\n', r.Ltot);
fprintf(' Catalyst volume     = %.2f m3   (%.0f L per MTPD)\n', ...
        b.Ac*r.Ltot, 1000*b.Ac*r.Ltot/150);
fprintf(' Catalyst weight     = %.0f kg\n', r.Wtot);
fprintf(' Overall conversion  = %.5f\n', r.Xtot);
fprintf(' Bed pressure drop   = %.2f kPa  (%.2f%% of inlet pressure)\n', ...
        r.dPtot/1000, 100*r.dPtot/(b.P(1)*1e5));
end


function plotProfiles(r, tag)
%PLOTPROFILES  Temperature and cumulative conversion against bed depth.
figure('Name',['Bed profiles - ' tag]);
for j = 1:4
    subplot(2,2,j);
    yyaxis left;  plot(r.z{j}, r.T{j}, '-',  'LineWidth', 1.4);
    ylabel('T [\circC]');
    yyaxis right; plot(r.z{j}, 100*r.f{j}, '--','LineWidth', 1.4);
    ylabel('Cumulative conversion [%]');
    xlabel('Bed depth z [m]');
    title(sprintf('Bed %d', j)); grid on;
end
sgtitle(tag);
end
