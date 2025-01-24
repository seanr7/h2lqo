%% RUNME_ADVECDIFF_IRKA
% Script file to run all experiments on the advection diffusion problem
% with a quadratic cost function as the quantity of interest.

%
% This file is part of the archive Code, Data, and Results for Numerical 
% Experiments in "...".
% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

clc;
clear all;
close all;

% Get and set all paths
[rootpath, filename, ~] = fileparts(mfilename('fullpath'));
loadname                = [rootpath filesep() ...
    'data' filesep() filename];
savename                = [rootpath filesep() ...
    'results' filesep() filename];

% Add paths to drivers and data.
addpath([rootpath, '/drivers'])
addpath([rootpath, '/data'])
addpath([rootpath, '/results'])

% Write .log file, put in `out' folder
if exist([savename '.log'], 'file') == 2
    delete([savename '.log']);
end
outname = [savename '.log']';

diary(outname)
diary on; 

fprintf(1, ['SCRIPT: ' upper(filename) '\n']);
fprintf(1, ['========' repmat('=', 1, length(filename)) '\n']);
fprintf(1, '\n');


%% Problem data.
fprintf(1, 'Load problem data.\n');
fprintf(1, '------------------\n');

% Advection-diffusion problem taken from [Diaz et al., 2024]
% nx   = 3000;  No. of spatial grid points
% diff = 1;     Diffusion parameter
% adv  = 1;     Advection parameter

load('data/AdvecDiff_n3000.mat')
fprintf(1, '\n');

% Input, output dimensions.
[n, m] = size(B);
[p, ~] = size(Clin);
nx     = 3000;

%% Compute reduced order models.
r  = 30; % Order of approximation.
E  = speye(n, n);

% Boolean; set `true' to recompute reduced models.
recompute = true;
% recompute = false;
if recompute
    fprintf(1, '1. Computing reduced model using LQO-IRKA (default interpolation data).\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % Input opts; default interpolation data.
    optsIRKA             = struct();
    optsIRKA.poles       = -4*logspace(0, 4, r)';
    optsIRKA.tol         = 10e-8;
    optsIRKA.maxIter     = 200; 
    optsIRKA.checkInterp = true;
    morStart             = tic;
    
    [ErIRKA, ArIRKA, BrIRKA, ClinrIRKA, MquadrIRKA, infoIRKA] = miso_lqoirka(...
        E, A, B, Clin.', Mquad, r, optsIRKA);
    save('results/AdvecDiff3000_IRKA_r30.mat', 'ErIRKA', 'ArIRKA', 'BrIRKA', ...
        'ClinrIRKA', 'MquadrIRKA', 'infoIRKA')

    fprintf(1, '1. ROM COMPUTED VIA LQO-IRKA (default interpolation data) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % Write convergence data.
    write = 1;
    if write
        poleChange       = infoIRKA.poleChange';
        [numIterates, ~] = size(poleChange);
        conv             = [(1:numIterates)', poleChange];
        dlmwrite('results/AdvecDiff3000_IRKA_r30conv.dat', conv, 'delimiter', ...
            '\t', 'precision', 8);
    end

    fprintf(1, '2. Computing reduced model using LQO-BT.\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');

    morStart = tic;

    [ErBT, ArBT, BrBT, ClinrBT, MquadrBT, infoBT] = miso_lqobt(E, A, B, Clin.', ...
        Mquad, r);
    save('results/AdvecDiff3000_BT_r30.mat', 'ArBT', 'BrBT', 'ClinrBT', ...
        'MquadrBT', 'infoBT')

    fprintf(1, '2. ROM COMPUTED VIA LQO-BT IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % TODO: Update these computations ... Better strategy for points and
    % directions?
     
    % Methods of tangential interpolation.
   
    % Shifts and tangent directions.
    % shifts         = 1i*logspace(1, 4, r/2);
    % shifts         = [shifts, conj(shifts)];
    % shifts         = cplxpair(shifts);
    optsTangInterp                = struct();
    optsTangInterp.shifts         = -1*logspace(0, 6, r);
    optsTangInterp.rightTangents  = 10*rand(m, r);
    optsTangInterp.loLeftTangents = 10*rand(r, 1);
    tmp                           = 10*rand(r, r);                     
    optsTangInterp.qoLeftTangents = (tmp + tmp')/2;
    optsTangInterp.interpolation  = 'standard';

    fprintf(1, '3. Computing reduced-order model by tangential interpolation (standard).\n')
    fprintf(1, '-----------------------------------------------------------------------------\n');
    morStart = tic;

    [ErTangInterp, ArTangInterp, BrTangInterp, ClinrTangInterp, MquadrTangInterp] = ...
        miso_lqointerp(E, A, B, Clin.', Mquad, r, optsTangInterp);

    save('results/AdvecDiff3000_TangInterp_r30.mat', 'ArTangInterp', 'BrTangInterp', ...
        'ClinrTangInterp', 'MquadrTangInterp')

    fprintf(1, '3. ROM COMPUTED VIA TANGENTIAL INTERPOLATION (standard) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    fprintf(1, '4. Computing reduced-order model by tangential interpolation (mixed).\n')
    fprintf(1, '-----------------------------------------------------------------------------\n');
    morStart = tic;
    optsMixedInterp                = struct();
    optsMixedInterp.shifts         = -1*logspace(0, 6, r);
    optsMixedInterp.rightTangents  = 10*rand(m, r);
    optsMixedInterp.loLeftTangents = 10*rand(r, 1);
    tmp                            = 10*rand(r, r);                     
    optsMixedInterp.qoLeftTangents = (tmp + tmp')/2;
    optsMixedInterp.interpolation  = 'mixed';

    [ErMixedInterp, ArMixedInterp, BrMixedInterp, ClinrMixedInterp, MquadrMixedInterp] = ...
        miso_lqointerp(E, A, B, Clin.', Mquad, r, optsMixedInterp);

    save('results/AdvecDiff3000_MixedTangInterp_r30.mat', 'ArMixedInterp', 'BrMixedInterp', ...
        'ClinrMixedInterp', 'MquadrMixedInterp')

    fprintf(1, '4. ROM COMPUTED VIA TANGENTIAL INTERPOLATION (mixed) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

else 
    fprintf(1, 'Loading reduced-order models.\n');
    fprintf(1, '-----------------------------\n');
    load('results/AdvecDiff3000_IRKA_r30.mat')
    load('results/AdvecDiff3000_BT_r30.mat')   
    load('results/AdvecDiff3000_TangInterp_r30.mat') 
    load('results/AdvecDiff3000_MixedTangInterp_r30.mat') 
end

fprintf(1, '\n');


%% Simulate full and reduced-order outputs; sinusoidal input.
fprintf(1, 'Solving full-order AdvecDiff problem via ode15\n');
fprintf(1, '----------------------------------------------\n');

Tfin = 10;             % Final time for simulation
nt   = round(Tfin*64); % Number of time steps

fprintf(1, 'First, sinusoidal input u0(t) = 0.5*(cos(pi*t)+1)\n');
fprintf(1, '-------------------------------------------------\n');

% Set system inputs 
u0 = @(t) 0.5*(cos(pi*t)+1); % Dirichlet input on left boundary
u1 = @(t) zeros(size(t));    % Neumann input on right boundary

ode_rtol = 1e-6; 
tsteps   = linspace(0, Tfin, nt+1); % Time-steps
v0       = zeros(nx, 1);            % Initial condition 
options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(nx, nx), 'Jacobian', A, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
[t, v]   = ode15s(fAdvDiff, tsteps, v0, options);
% Note, v is nt \times nx
% Compute quadratic cost by 1/(2*n) * ||x(t) - \bf{1}||_2^2

fprintf(1, 'Simulate full output.\n');
fprintf(1, '--------------------.\n');

% Output has mixed linear and quadratic terms
y = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, length(t));

% Colors for plotting
ColMat(1,:) = [0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];

figure(1)
set(gca, 'fontsize', 10)
fs = 12; % Fontsize
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
hold on
grid on
xlabel('$t$','fontsize',fs,'interpreter','latex'); 
ylabel('$y(t)$','fontsize',fs,'interpreter','latex')


fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15\n');
fprintf(1, '--------------------------------------------------\n');

% 1. For LQOIRKA reduced model.
ode_rtol        = 1e-6; 
nxr             = r;
vr0             = zeros(nxr, 1); % Initial condition 
options         = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArIRKA, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffredIRKA = @(tr, vrIRKA)(ArIRKA*vrIRKA + BrIRKA*[u0(tr);u1(tr)]);
[tr, vrIRKA]    = ode15s(fAdvDiffredIRKA, tsteps, vr0, options);
% Note, vr is nt \times r

fprintf(1, 'Simulate IRKA reduced output\n');
fprintf(1, '----------------------------\n');
yrIRKA = ClinrIRKA.'*vrIRKA.';
for tt = 1:length(tr)
    yrIRKA(:,tt) = yrIRKA(:,tt) + vrIRKA(tt, :)*MquadrIRKA*vrIRKA(tt, :)';
end
yrIRKA = yrIRKA + (1/2) * ones(1, length(tr));

% 2. For LQOBT reduced model.
vr0           = zeros(nxr, 1); % Initial condition 
options       = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArBT, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffredBT = @(tr,vrBT)(ArBT*vrBT + BrBT*[u0(tr);u1(tr)]);
[tr, vrBT]    = ode15s(fAdvDiffredBT, tsteps, vr0, options);
% Note, vr is nt \times nxr

fprintf(1, 'Simulate bt reduced output\n');
fprintf(1, '--------------------------\n');
yrBT = ClinrBT.'*vrBT.';
for tt = 1:length(tr)
    yrBT(:,tt) = yrBT(:,tt) + vrBT(tt, :)*MquadrBT*vrBT(tt, :)';
end
yrBT = yrBT + (1/2) * ones(1, length(tr));

% 3. For tangential interpolation.
vr0                   = zeros(nxr, 1); % Initial condition 
options               = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                        'Mass', eye(r,r), 'Jacobian', ArIRKA, ...
                        'MStateDependence', 'none', 'Stats','on');
fAdvDiffredTangInterp = @(tr, vrTangInterp)(ArTangInterp*vrTangInterp + BrTangInterp*[u0(tr);u1(tr)]);
[tr, vrTangInterp]    = ode15s(fAdvDiffredTangInterp, tsteps, vr0, options);

fprintf(1, 'Simulate tang. reduced output\n');
fprintf(1, '----------------------------\n');
yrTangInterp = CLinrTangInterp*vrTangInterp.';
for tt = 1:length(tr)
    yrTangInterp(:,tt) = yrTangInterp(:,tt) + vrTangInterp(tt, :)*MQuadrTangInterp*vrTangInterp(tt, :)';
end
yrTangInterp = yrTangInterp + (1/2) * ones(1, length(tr));

% 3. For mixed tangential interpolation.
vr0                    = zeros(nxr, 1); % Initial condition 
options                = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                         'Mass', eye(r,r), 'Jacobian', ArIRKA, ...
                         'MStateDependence', 'none', 'Stats','on');
fAdvDiffredMixedInterp = @(tr, vrMixedInterp)(ArMixedInterp*vrMixedInterp + BrMixedInterp*[u0(tr);u1(tr)]);
[tr, vrMixedInterp]    = ode15s(fAdvDiffredMixedInterp, tsteps, vr0, options);

fprintf(1, 'Simulate tang. reduced output\n');
fprintf(1, '----------------------------\n');
yrMixedInterp = CLinrMixedInterp*vrMixedInterp.';
for tt = 1:length(tr)
    yrMixedInterp(:,tt) = yrMixedInterp(:,tt) + vrMixedInterp(tt, :)*MQuadrMixedInterp*vrMixedInterp(tt, :)';
end
yrMixedInterp = yrMixedInterp + (1/2) * ones(1, length(tr));

plotBool = true;
if plotBool
    plot(t, yrIRKA, '--', 'color',ColMat(2,:), LineWidth=1.5); hold on
    plot(t, yrBT, '-.', 'color',ColMat(3,:), LineWidth=1.5); 
    % plot(t, yrTangInterp, '-.', 'color',ColMat(4,:), LineWidth=1.5); 
    % plot(tr, yrMixedInterp, '-.', 'color',ColMat(5,:), LineWidth=1.5); 
    lgd = legend('$y(t)$', '$y_r(t)$', '$y_{r,bt}(t)$', '$y_{tang}(t)$', '$y_{mixed}(t)$', 'interpreter','latex', ...
        'FontName', 'Arial', 'location', 'northeast');
    fontsize(lgd, 10, 'points')
    
    subplot(2,1,2)
    semilogy(tr, abs(y - yrIRKA)./abs(y), '--','color',ColMat(2,:),LineWidth=1.5); 
    hold on;
    semilogy(tr, abs(y - yrBT)./abs(y), '-.','color',ColMat(3,:),LineWidth=1.5); 
    % semilogy(tr, abs(y - yrTangInterp)./abs(y), '-.','color',ColMat(4,:),LineWidth=1.5); 
    % semilogy(tr, abs(y - yrMixedInterp)./abs(y), '-.','color',ColMat(5,:),LineWidth=1.5); 
    xlabel('$t$','interpreter','latex'); 
    ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
        LineWidth=1.5)
    
    % Overwrite figure
    % print -depsc2 results/AdvecDiff1200_sinusoidal_r20_OutputPlot
    
    % Write data
    write = true;
    if write
        outputs = [t, y', yrIRKA', yrBT', yrTangInterp', yrMixedInterp.'];
        dlmwrite('results/AdvecDiff1200_sinusoidal_r20_Outputs.dat', outputs, 'delimiter', ...
            '\t', 'precision', 8);
        outputerrors = [t, (abs(y-yrIRKA)./abs(y))', (abs(y-yrBT)./abs(y))', (abs(y-yrTangInterp)./abs(y))', ...
            (abs(y-yrMixedInterp)./abs(y))'];
        dlmwrite('results/AdvecDiff1200_sinusoidal_r20_OutputErrors.dat', outputerrors, ...
            'delimiter', '\t', 'precision', 8);
    end
end


% Max errors.
data             = load('results/AdvecDiff1200_sinusoidal_r20_OutputErrors.dat');
errorIRKA        = data(:, 2);
errorBT          = data(:, 3);
errorTangInterp  = data(:, 4);
errorMixedInterp = data(:, 5);
fprintf(1, 'Max error due to LQO-IRKA           : %.10f \n', max(errorIRKA))
fprintf(1, 'Max error due to LQO-BT             : %.10f \n', max(errorBT))
fprintf(1, 'Max error due to interp (tangential): %.10f \n', max(errorTangInterp))
fprintf(1, 'Max error due to interp (mixed      : %.10f \n', max(errorMixedInterp))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, '\n');


%% Simulate full and reduced-order outputs; exponentially damped input.
% fprintf(1, 'Solving full-order AdvecDiff problem via ode15\n');
% fprintf(1, '----------------------------------------------\n');
% 
% fprintf(1, 'Second, exponentially damped input u0(t) = exp(-t/5)*t^2\n');
% fprintf(1, '--------------------------------------------------------\n');
% 
% Tfin   = 30;                      % Final time for simulation
% nt     = round(Tfin*10);          % Number of time steps
% tsteps = linspace(0, Tfin, nt+1); % Time-steps
% 
% % Reset relevant input
% u0 = @(t) exp(-t/5)*t^2; % Exponentially damped input on left boundary
% 
% options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
%                   'Mass', eye(nx, nx), 'Jacobian', A, ...
%                   'MStateDependence', 'none', 'Stats','on');
% 
% % Redefine problem
% fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
% [t, v]   = ode15s(fAdvDiff, tsteps, v0, options);
% % Note, v is nt \times nx
% % Compute quadratic cost by 1/(2*n) * ||x(t) - \bf{1}||_2^2
% 
% fprintf(1, 'Simulate full output.\n');
% fprintf(1, '--------------------.\n');
% 
% % Output has mixed linear and quadratic terms
% y = Clin*v';
% for tt = 1:length(t)
%     y(:, tt) = y(:, tt) + v(tt, :)*Mquad*v(tt, :)';
% end
% y = y + (1/2) * ones(1, length(t));
% 
% % Colors for plotting
% ColMat(1,:) = [0.8500    0.3250    0.0980];
% ColMat(2,:) = [0.3010    0.7450    0.9330];
% ColMat(3,:) = [0.9290    0.6940    0.1250];
% ColMat(4,:) = [0.4660    0.6740    0.1880];
% ColMat(5,:) = [0.4940    0.1840    0.5560];
% 
% figure(2)
% set(gca, 'fontsize', 10)
% fs = 12; % Fontsize
% subplot(2,1,1)
% plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
% hold on
% grid on
% xlabel('$t$','fontsize',fs,'interpreter','latex'); 
% ylabel('$y(t)$','fontsize',fs,'interpreter','latex')
% 
% 
% fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15\n');
% fprintf(1, '--------------------------------------------------\n');
% 
% % 1. For LQOIRKA reduced model.
% ode_rtol        = 1e-6; 
% nxr             = r;
% vr0             = zeros(nxr, 1); % Initial condition 
% options         = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
%                   'Mass', eye(r,r), 'Jacobian', ArIRKA, ...
%                   'MStateDependence', 'none', 'Stats','on');
% fAdvDiffredIRKA = @(tr, vrIRKA)(ArIRKA*vrIRKA + BrIRKA*[u0(tr);u1(tr)]);
% [tr, vrIRKA]    = ode15s(fAdvDiffredIRKA, tsteps, vr0, options);
% % Note, vr is nt \times r
% 
% fprintf(1, 'Simulate IRKA reduced output\n');
% fprintf(1, '----------------------------\n');
% yrIRKA = CLinrIRKA.'*vrIRKA.';
% for tt = 1:length(tr)
%     yrIRKA(:,tt) = yrIRKA(:,tt) + vrIRKA(tt, :)*MQuadrIRKA*vrIRKA(tt, :)';
% end
% yrIRKA = yrIRKA + (1/2) * ones(1, length(tr));

% % 2. For LQOBT reduced model.
% vr0           = zeros(nxr, 1); % Initial condition 
% options       = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
%                   'Mass', eye(r,r), 'Jacobian', ArBT, ...
%                   'MStateDependence', 'none', 'Stats','on');
% fAdvDiffredBT = @(tr,vrBT)(ArBT*vrBT + BrBT*[u0(tr);u1(tr)]);
% [tr, vrBT]    = ode15s(fAdvDiffredBT, tsteps, vr0, options);
% % Note, vr is nt \times nxr
% 
% fprintf(1, 'Simulate bt reduced output\n');
% fprintf(1, '--------------------------\n');
% yrBT = CLinrBT.'*vrBT.';
% for tt = 1:length(tr)
%     yrBT(:,tt) = yrBT(:,tt) + vrBT(tt, :)*MQuadrBT*vrBT(tt, :)';
% end
% yrBT = yrBT + (1/2) * ones(1, length(tr));
% 
% % 3. For tangential interpolation.
% vr0                   = zeros(nxr, 1); % Initial condition 
% options               = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
%                         'Mass', eye(r,r), 'Jacobian', ArIRKA, ...
%                         'MStateDependence', 'none', 'Stats','on');
% fAdvDiffredTangInterp = @(tr, vrTangInterp)(ArTangInterp*vrTangInterp + BrTangInterp*[u0(tr);u1(tr)]);
% [tr, vrTangInterp]    = ode15s(fAdvDiffredTangInterp, tsteps, vr0, options);
% 
% fprintf(1, 'Simulate tang. reduced output\n');
% fprintf(1, '----------------------------\n');
% yrTangInterp = CLinrTangInterp*vrTangInterp.';
% for tt = 1:length(tr)
%     yrTangInterp(:,tt) = yrTangInterp(:,tt) + vrTangInterp(tt, :)*MQuadrTangInterp*vrTangInterp(tt, :)';
% end
% yrTangInterp = yrTangInterp + (1/2) * ones(1, length(tr));
% 
% % 3. For mixed tangential interpolation.
% vr0                    = zeros(nxr, 1); % Initial condition 
% options                = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
%                          'Mass', eye(r,r), 'Jacobian', ArIRKA, ...
%                          'MStateDependence', 'none', 'Stats','on');
% fAdvDiffredMixedInterp = @(tr, vrMixedInterp)(ArMixedInterp*vrMixedInterp + BrMixedInterp*[u0(tr);u1(tr)]);
% [tr, vrMixedInterp]    = ode15s(fAdvDiffredMixedInterp, tsteps, vr0, options);
% 
% fprintf(1, 'Simulate tang. reduced output\n');
% fprintf(1, '----------------------------\n');
% yrMixedInterp = CLinrMixedInterp*vrMixedInterp.';
% for tt = 1:length(tr)
%     yrMixedInterp(:,tt) = yrMixedInterp(:,tt) + vrMixedInterp(tt, :)*MQuadrMixedInterp*vrMixedInterp(tt, :)';
% end
% yrMixedInterp = yrMixedInterp + (1/2) * ones(1, length(tr));
% 
% plot(t, yrIRKA, '--', 'color',ColMat(2,:), LineWidth=1.5); hold on
% plot(t, yrBT, '-.', 'color',ColMat(3,:), LineWidth=1.5); 
% plot(t, yrTangInterp, '-.', 'color',ColMat(4,:), LineWidth=1.5); 
% plot(t, yrMixedInterp, '-.', 'color',ColMat(5,:), LineWidth=1.5); 
% lgd = legend('$y(t)$', '$y_r(t)$', '$y_{r,bt}(t)$', '$y_{tang}(t)$', '$y_{mixed}(t)$', 'interpreter','latex', ...
%     'FontName', 'Arial', 'location', 'northeast');
% fontsize(lgd, 10, 'points')
% 
% subplot(2,1,2)
% semilogy(tr, abs(y - yrIRKA)./abs(y), '--','color',ColMat(2,:),LineWidth=1.5); 
% hold on;
% semilogy(tr, abs(y - yrBT)./abs(y), '-.','color',ColMat(3,:),LineWidth=1.5); 
% semilogy(tr, abs(y - yrTangInterp)./abs(y), '-.','color',ColMat(4,:),LineWidth=1.5); 
% semilogy(tr, abs(y - yrMixedInterp)./abs(y), '-.','color',ColMat(5,:),LineWidth=1.5); 
% xlabel('$t$','interpreter','latex'); 
% ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
%     LineWidth=1.5)
% 
% % Overwrite figure
% print -depsc2 results/AdvecDiff_exponential_r20_OutputPlot
% 
% % Write data
% write = 1;
% if write
%     outputs = [t, y', yrIRKA', yrBT', yrTangInterp', yrMixedInterp.'];
%     dlmwrite('results/AdvecDiff_exponential_r20_Outputs.dat', outputs, 'delimiter', ...
%         '\t', 'precision', 8);
%     outputerrors = [t, (abs(y-yrIRKA)./abs(y))', (abs(y-yrBT)./abs(y))', (abs(y-yrTangInterp)./abs(y))', ...
%         (abs(y-yrMixedInterp)./abs(y))'];
%     dlmwrite('results/AdvecDiff_exponential_r20_OutputErrors.dat', outputerrors, ...
%         'delimiter', '\t', 'precision', 8);
% end
% 
% fprintf(1, '\n');
% 

%% Hierarchy of reduced models.
fprintf(1, 'Computing hiearchy of reduced models using LQO-TSIA and LQO-BT\n')
fprintf(1, '--------------------------------------------------------------\n')

rmax         = 20;
relErrorIRKA = zeros(rmax/2, 1);
relErrorBT   = zeros(rmax/2, 1);
% Precompute H2 norm of full-order model for relative error.
P       = lyap(A, B*B');
Q       = lyap(A', Clin'*Clin + Mquad*P*Mquad);
fomNorm = sqrt(abs(trace(B'*Q*B)));

% Toggle input opts.
for r = 2:2:rmax  
    fprintf(1, 'Current reduced model is hierarchy; r=%d\n', r)
    fprintf(1, '----------------------------------------\n')
     % LQO-IRKA.
    optsIRKA             = struct();
    optsIRKA.poles       = -6*logspace(0, 5, r)';
    optsIRKA.tol         = 10e-6;
    optsIRKA.maxIter     = 500; 
    optsIRKA.checkInterp = false;
    optsIRKA.plotConv    = false;
    [ArIRKA, BrIRKA, ClinrIRKA, MquadrIRKA, infoIRKA] = miso_lqoirka(...
        A, B, Clin.', Mquad, r, optsIRKA);

    AerrIRKA           = [A, zeros(nx, r); zeros(r, nx), ArIRKA]; 
    BerrIRKA           = [B; BrIRKA];
    CerrIRKA           = [Clin, - ClinrIRKA']; 
    MerrIRKA           = [Mquad, zeros(nx, r); zeros(r, nx), -MquadrIRKA];
    PerrIRKA           = lyap(AerrIRKA, BerrIRKA*BerrIRKA');                  
    QerrIRKA           = lyap(AerrIRKA', CerrIRKA'*CerrIRKA + MerrIRKA*PerrIRKA*MerrIRKA);
    k                  = r/2;  
    relErrorIRKA(k, :) = sqrt(abs(trace(BerrIRKA'*QerrIRKA*BerrIRKA)))/fomNorm;

    fprintf(1, 'Relative H2 error of LQO-TSIA reduced model is ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relErrorIRKA(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

    % LQO-BT
    [ArBT, BrBT, ClinrBT, MquadrBT, infoBT] = miso_lqobt(A, B, Clin.', ...
        Mquad, r);

    % Build error realization, and compute H2 error using Gramian
    % formulation
    AerrBT        = [A, zeros(nx, r); zeros(r, nx), ArBT]; 
    BerrBT        = [B; BrBT];
    CerrBT        = [Clin, - ClinrBT']; 
    MerrBT        = [Mquad, zeros(nx, r); zeros(r, nx), -MquadrBT];
    PerrBT        = lyap(AerrBT, BerrBT*BerrBT');                  
    QerrBT        = lyap(AerrBT', CerrBT'*CerrBT + MerrBT*PerrBT*MerrBT);
    k             = r/2;  
    relErrorBT(k, :) = sqrt(abs(trace(BerrBT'*QerrBT*BerrBT)))/fomNorm;
    fprintf(1, 'Relative H2 error of LQO-BT reduced model is ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relErrorBT(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

end

write = 1;
if write
    relh2errors = [(2:2:rmax)', relErrorIRKA, relErrorBT];
    dlmwrite('results/AdvecDiff1200_IRKA_BT_H2errors.dat', relh2errors, 'delimiter', '\t', 'precision', ...
        8);
end


%% Plot H2 errors.
fprintf(1, 'Plotting H2 errors\n');
fprintf(1, '------------------\n');

figure(4)
set(gca, 'fontsize', 10)
semilogy(2:2:rmax, relErrorIRKA(1:rmax/2), '-o', 'color', ColMat(2,:), LineWidth=1.5)
hold on;
semilogy(2:2:rmax, relErrorBT(1:rmax/2), '-*', 'color', ColMat(3,:), LineWidth=1.5)
xlim([2,rmax])
lgd = legend('lqo-tsia', 'lqo-bt', 'interpreter', 'latex', 'FontName', 'Arial',...
    'location', 'northeast');
xlabel('$r$', 'interpreter', 'latex')
ylabel('$||\mathcal{G}-\mathcal{G}_r||_{\mathcal{H}_2}/||\mathcal{G}||_{\mathcal{H}_2}$',...
    'interpreter', 'latex')

print -depsc2 results/AdvecDiff1200_H2errors_Plot


%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
