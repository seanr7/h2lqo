%% RUNME_ADVECDIFF
% Script file to run all experiments on the advection diffusion problem
% with a quadratic cost function as the quantity of interest.

%
% This file is part of the archive Code, Data, and Results for Numerical 
% Experiments in "$\mathcal{H}_2$ optimal model reduction of linear 
% systems with multiple quadratic outputs".
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

% Add paths to drivers and data
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


%% Load base data.
fprintf(1, 'Load problem data\n');
fprintf(1, '-----------------\n');

% Advection-diffusion problem taken from [Diaz et al., 2024]
% nx   = 3000; % No. of spatial grid points
% adv  = 1;    % Advection parameter
% diff = 1;   % Diffusion parameter

load('data/AdvecDiff_n3000.mat')
fprintf(1, '\n');

% Input, output dimensions.
[n, m] = size(B);
[p, ~] = size(Clin);

%% Compute reduced-order models.

% Order of approximation.
r = 30; 

% Boolean; set `true' to recompute reduced models.
recompute = true;
% recompute = false;
if recompute
    fprintf(1, '1. Computing reduced model using LQO-TSIA (defailt, diagonal initialization).\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % Input opts.
    opts_stdInit                = struct();
    opts_stdInit.tol            = 10e-14;
    opts_stdInit.maxIter        = 200; 
    opts_stdInit.convMonitoring = 'approximate';

    % Initial model values are vanilla; diagonal matrix spread among
    % full-order model poles, and identity input, output matrices.
    opts_stdInit.Ar = -diag(logspace(0, 4, r));
    opts_stdInit.Br = eye(r, m);
    opts_stdInit.Cr = eye(p, r);
    opts_stdInit.Mr = eye(r, r);
    morStart        = tic;

    [ArTSIA_stdInit, BrTSIA_stdInit, ClinrTSIA_stdInit, MquadrTSIA_stdInit, infoTSIA_stdInit] ...
        = mimo_lqotsia(A, B, Clin, Mquad, r, opts_stdInit);
    save('results/AdvecDiff3000_TSIA_stdInit_r30.mat', 'ArTSIA_stdInit', 'BrTSIA_stdInit', ...
        'ClinrTSIA_stdInit', 'MquadrTSIA_stdInit', 'opts_stdInit', 'infoTSIA_stdInit')

    fprintf(1, '1. ROM COMPUTED VIA LQO-TSIA (default, diagonal initialization) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------\n');

    fprintf(1, '2. Computing reduced model using LQO-TSIA (truncated initialization).\n');
    fprintf(1, '---------------------------------------------------------------------\n');

    % Input opts.
    opts_truncInit                = struct();
    opts_truncInit.tol            = 10e-14;
    opts_truncInit.maxIter        = 200; 
    opts_truncInit.convMonitoring = 'approximate';

    % Initial model values are obtained by truncating full-order matrices.
    opts_truncInit.Ar = A(1:r, 1:r);
    opts_truncInit.Br = B(1:r, :);
    opts_truncInit.Cr = Clin(:, 1:r);
    opts_truncInit.Mr = Mquad(1:r, 1:r);
    morStart          = tic;

    [ArTSIA_truncInit, BrTSIA_truncInit, ClinrTSIA_truncInit, MquadrTSIA_truncInit, infoTSIA_truncInit] ...
        = mimo_lqotsia(A, B, Clin, Mquad, r, opts_truncInit);
    save('results/AdvecDiff3000_TSIA_truncInit_r30.mat', 'ArTSIA_truncInit', 'BrTSIA_truncInit', ...
        'ClinrTSIA_truncInit', 'MquadrTSIA_truncInit', 'opts_truncInit', 'infoTSIA_truncInit')

    fprintf(1, '2. ROM COMPUTED VIA LQO-TSIA (truncated initialization) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------\n');
 
    fprintf(1, '3. Computing reduced model using LQO-TSIA (eigs initialization).\n');
    fprintf(1, '---------------------------------------------------------------------\n');

    % Input opts.
    opts_eigsInit                = struct();
    opts_eigsInit.tol            = 10e-14;
    opts_eigsInit.maxIter        = 200; 
    opts_eigsInit.convMonitoring = 'approximate';

    % Initial model values are obtained from projecting the dominant
    % eigenvectors.
    eigsOpts     = struct();
    eigsOpts.p   = 100;
    eigsOpts.tol = 1e-10;
    [V, D, ~]    = eigs(A, r, 'bothendsreal', eigsOpts);
    [V, ~]       = qr(V, 'econ');

    iopts_eigsInit.Ar = (V.'*A*V);
    opts_eigsInit.Br = (V.'*B);
    opts_eigsInit.Cr = Clin*V;
    opts_eigsInit.Mr = V.'*Mquad*V;
    morStart         = tic;

    [ArTSIA_eigsInit, BrTSIA_eigsInit, ClinrTSIA_eigsInit, MquadrTSIA_eigsInit, infoTSIA_eigsInit] ...
        = mimo_lqotsia(A, B, Clin, Mquad, r, opts_eigsInit);
    save('results/AdvecDiff3000_TSIA_eigsInit_r30.mat', 'ArTSIA_eigsInit', 'BrTSIA_eigsInit', ...
        'ClinrTSIA_eigsInit', 'MquadrTSIA_eigsInit', 'opts_eigsInit', 'infoTSIA_eigsInit')

    fprintf(1, '3. ROM COMPUTED VIA LQO-TSIA (eigs initialization) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------\n');

    fprintf(1, '4. Computing reduced model using LQO-BT.\n');
    fprintf(1, '----------------------------------------\n');
    
    morStart                                = tic;
    [ArBT, BrBT, ClinrBT, MquadrBT, infoBT] = mimo_lqobt(A, B, Clin, Mquad, r);
    save('results/AdvecDiff3000_BT_r30.mat', 'ArBT', 'BrBT', 'ClinrBT', ...
        'MquadrBT', 'infoBT')

    fprintf(1, '4. ROM COMPUTED VIA LQO-BT IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------\n');
else 
    fprintf(1, 'Loading reduced-order models.\n');
    fprintf(1, '-----------------------------\n');
    load('results/AdvecDiff3000_TSIA_stdInit_r30.mat')
    load('results/AdvecDiff3000_TSIA_truncInit_r30.mat')
    load('results/AdvecDiff3000_TSIA_eigsInit_r30.mat')  
    load('results/AdvecDiff3000_BT_r30.mat')
end

fprintf(1, '\n');


%% Full-order simulation: sinusoidal input.
fprintf(1, 'Solving full-order AdvecDiff problem via ode15.\n');
fprintf(1, '-----------------------------------------------\n');

Tfin = 10;             % Final time for simulation
nt   = round(Tfin*64); % Number of time steps

fprintf(1, '1. Sinusoidal input: u0(t) = 0.5*(cos(pi*t)+1).\n');
fprintf(1, '--------------------------------------------------\n');

% Set system inputs.
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
% Note, v is nt \times nx.
% Compute quadratic cost by 1/(2*n) * ||x(t) - \bf{1}||_2^2.

fprintf(1, 'Simulate full output.\n');
fprintf(1, '--------------------.\n');

% Output has mixed linear and quadratic terms.
y = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + v(tt, :)*Mquad*v(tt, :)';
end
% Recover output.
y = y + (1/(2*nx)) * ones(1, length(t));

% Colors for plotting.
ColMat(1,:) = [0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];

% Plot full-order output.
figure(1)
set(gca, 'fontsize', 10)
fs = 12; % Fontsize
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
hold on
grid on
xlabel('$t$','fontsize',fs,'interpreter','latex'); 
ylabel('$y(t)$','fontsize',fs,'interpreter','latex')

%% Reduced-order simulations: sinusoidal input.
fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15.\n');
fprintf(1, '---------------------------------------------------\n');

fprintf(1, '1. Simulate LQO-TSIA (diagonal initialization) reduced output.\n');
fprintf(1, '--------------------------------------------------------------\n');

% Opts.
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArTSIA_stdInit, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffred          = @(tr, vrTSIA_stdInit)(ArTSIA_stdInit*vrTSIA_stdInit + BrTSIA_stdInit*[u0(tr);u1(tr)]);
[tr, vrTSIA_stdInit] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r.

yrTSIA_stdInit = ClinrTSIA_stdInit*vrTSIA_stdInit';
for tt = 1:length(tr)
    yrTSIA_stdInit(:,tt) = yrTSIA_stdInit(:,tt) + ...
        vrTSIA_stdInit(tt, :)*MquadrTSIA_stdInit*vrTSIA_stdInit(tt, :)';
end
yrTSIA_stdInit = yrTSIA_stdInit + (1/(2*nx))*ones(1, length(tr));

fprintf(1, '2. Simulate LQO-TSIA (truncated initialization) reduced output.\n');
fprintf(1, '--------------------------------------------------------------\n');

% Opts.
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArTSIA_truncInit, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffred            = @(tr, vrTSIA_truncInit)(ArTSIA_truncInit*vrTSIA_truncInit + BrTSIA_truncInit*[u0(tr);u1(tr)]);
[tr, vrTSIA_truncInit] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r.

yrTSIA_truncInit = ClinrTSIA_truncInit*vrTSIA_truncInit';
for tt = 1:length(tr)
    yrTSIA_truncInit(:,tt) = yrTSIA_truncInit(:,tt) + ...
        vrTSIA_truncInit(tt, :)*MquadrTSIA_truncInit*vrTSIA_truncInit(tt, :)';
end
yrTSIA_truncInit = yrTSIA_truncInit + (1/(2*nx))*ones(1, length(tr));

fprintf(1, '3. Simulate LQO-TSIA (eigs initialization) reduced output.\n');
fprintf(1, '--------------------------------------------------------\n');

% Opts.
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArTSIA_eigsInit, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffred         = @(tr, vrTSIA_eigsInit)(ArTSIA_eigsInit*vrTSIA_eigsInit + BrTSIA_eigsInit*[u0(tr);u1(tr)]);
[tr, vrTSIA_eigsInit] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r.

yrTSIA_eigsInit = ClinrTSIA_eigsInit*vrTSIA_eigsInit';
for tt = 1:length(tr)
    yrTSIA_eigsInit(:,tt) = yrTSIA_eigsInit(:,tt) + ...
        vrTSIA_eigsInit(tt, :)*MquadrTSIA_eigsInit*vrTSIA_eigsInit(tt, :)';
end
yrTSIA_eigsInit = yrTSIA_eigsInit + (1/(2*nx))*ones(1, length(tr));

fprintf(1, '4. Simulate LQO-BT reduced output.\n');
fprintf(1, '----------------------------------\n');

% Opts.
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArBT, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffredbt = @(tr,vrBT)(ArBT*vrBT + BrBT*[u0(tr);u1(tr)]);
[tr, vrBT]    = ode15s(fAdvDiffredbt, tsteps, vr0, options);
% Note, vr is nt \times nxr

yrBT = ClinrBT*vrBT';
for tt = 1:length(tr)
    yrBT(:,tt) = yrBT(:,tt) + vrBT(tt, :)*MquadrBT*vrBT(tt, :)';
end
yrBT = yrBT + (1/(2*nx))*ones(1, length(tr));

% Plots.
plot(t, yrTSIA_stdInit,   '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
plot(t, yrTSIA_truncInit, '-.', 'color', ColMat(3,:), LineWidth=1.5); 
plot(t, yrTSIA_eigsInit,  '--', 'color', ColMat(4,:), LineWidth=1.5); 
plot(t, yrBT,             '-.', 'color', ColMat(5,:), LineWidth=1.5); 
lgd = legend('$y(t)$', '$y_{r,TSIAstd}(t)$', '$y_{r,TSIAtrunc}(t)$', '$y_{r,TSIAeigs}(t)$', ...
    '$y_{r,BT}(t)$', 'interpreter','latex', 'FontName', 'Arial', 'location', ...
    'northeast');
fontsize(lgd, 10, 'points')

subplot(2,1,2)
semilogy(tr, abs(y - yrTSIA_stdInit)./abs(y),   '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
semilogy(tr, abs(y - yrTSIA_truncInit)./abs(y), '-.', 'color', ColMat(3,:), LineWidth=1.5); 
semilogy(tr, abs(y - yrTSIA_eigsInit)./abs(y),  '--', 'color', ColMat(4,:), LineWidth=1.5); 
semilogy(tr, abs(y - yrBT)./abs(y),             '-.', 'color', ColMat(5,:), LineWidth=1.5); 
xlabel('$t$','interpreter','latex'); 
ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
    LineWidth=1.5)

% Print errors.
fprintf(1, 'Order r = %d.\n', r)
fprintf(1, '--------------\n')
fprintf(1, 'Sinusoidal input: Relative L-infty error due to LQO-TSIA (standard initialization) : %.16f \n', max(abs(y - yrTSIA_stdInit)./abs(y)))
fprintf(1, 'Sinusoidal input: Relative L-infty error due to LQO-TSIA (truncated initialization): %.16f \n', max(abs(y - yrTSIA_truncInit)./abs(y)))
fprintf(1, 'Sinusoidal input: Relative L-infty error due to LQO-TSIA (eigs initialization)     : %.16f \n', max(abs(y - yrTSIA_eigsInit)./abs(y)))
fprintf(1, 'Sinusoidal input: Relative L-infty error due to LQO-BT                             : %.16f \n', max(abs(y - yrBT)./abs(y)))
fprintf(1, '------------------------------------------------------------\n')
fprintf(1, 'Sinusoidal input: Relative L-2 error due to LQO-TSIA (standard initialization) : %.16f \n', sum(abs(y - yrTSIA_stdInit))/sum(abs(y)))
fprintf(1, 'Sinusoidal input: Relative L-2 error due to LQO-TSIA (truncated initialization): %.16f \n', sum(abs(y - yrTSIA_truncInit))/sum(abs(y)))
fprintf(1, 'Sinusoidal input: Relative L-2 error due to LQO-TSIA (eigs initialization)     : %.16f \n', sum(abs(y - yrTSIA_eigsInit))/sum(abs(y)))
fprintf(1, 'Sinusoidal input: Relative L-2 error due to LQO-BT                             : %.16f \n', sum(abs(y - yrBT))/sum(abs(y)))
fprintf(1, '------------------------------------------------------------\n')

% Write data.
write = 1;
if write
    % Overwrite figure.
    print -depsc2 results/AdvecDiff3000_sinusoidal_r30_OutputPlot

    outputs = [t, y', yrTSIA_stdInit', yrTSIA_truncInit', yrTSIA_eigsInit', yrBT'];
    dlmwrite('results/AdvecDiff3000_sinusoidal_r30_Outputs.dat', outputs, 'delimiter', ...
        '\t', 'precision', 8);
    outputerrors = [t, (abs(y-yrTSIA_stdInit)./abs(y))', (abs(y-yrTSIA_truncInit)./abs(y))', ...
        (abs(y-yrTSIA_eigsInit)./abs(y))', (abs(y-yrBT)./abs(y))'];
    dlmwrite('results/AdvecDiff3000_sinusoidal_r30_OutputErrors.dat', outputerrors, ...
        'delimiter', '\t', 'precision', 8);
end

fprintf(1, '\n');

%% Full-order simulation: exponential input.
fprintf(1, 'Solving full-order AdvecDiff problem via ode15.\n');
fprintf(1, '-----------------------------------------------\n');

Tfin   = 30;                      % Final time for simulation
nt     = round(Tfin*10);          % Number of time steps
% tsteps = linspace(0, Tfin, nt+1); % Time-steps

fprintf(1, 'Exponentially damped input: u0(t) = exp(-t/5)*t^2.\n');
fprintf(1, '--------------------------------------------------\n');

% Reset relevant input
u0 = @(t) exp(-t/5)*t^2; 
u1 = @(t) zeros(size(t));    % Neumann input on right boundary

ode_rtol = 1e-6; 
tsteps   = linspace(0, Tfin, nt+1); % Time-steps
v0       = zeros(nx, 1);            % Initial condition 
options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol', ode_rtol, ...
                  'Mass', eye(nx, nx), 'Jacobian', A, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
[t, v]   = ode15s(fAdvDiff, tsteps, v0, options);
% Note, v is nt \times nx.
% Compute quadratic cost by 1/(2*n) * ||x(t) - \bf{1}||_2^2.

fprintf(1, 'Simulate full output.\n');
fprintf(1, '--------------------.\n');

% Output has mixed linear and quadratic terms.
y = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + v(tt, :)*Mquad*v(tt, :)';
end
% Recover output.
y = y + (1/(2*nx)) * ones(1, length(t));

% Colors for plotting.
ColMat(1,:) = [0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];

% Plot full-order output.
figure(2)
set(gca, 'fontsize', 10)
fs = 12; % Fontsize
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
hold on
grid on
xlabel('$t$','fontsize',fs,'interpreter','latex'); 
ylabel('$y(t)$','fontsize',fs,'interpreter','latex')

%% Reduced-order simulations: exponential input.
fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15.\n');
fprintf(1, '---------------------------------------------------\n');

fprintf(1, '1. Simulate LQO-TSIA (diagonal initialization) reduced output.\n');
fprintf(1, '--------------------------------------------------------------\n');

% Opts.
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArTSIA_stdInit, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffred          = @(tr, vrTSIA_stdInit)(ArTSIA_stdInit*vrTSIA_stdInit + BrTSIA_stdInit*[u0(tr);u1(tr)]);
[tr, vrTSIA_stdInit] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r.

yrTSIA_stdInit = ClinrTSIA_stdInit*vrTSIA_stdInit';
for tt = 1:length(tr)
    yrTSIA_stdInit(:,tt) = yrTSIA_stdInit(:,tt) + ...
        vrTSIA_stdInit(tt, :)*MquadrTSIA_stdInit*vrTSIA_stdInit(tt, :)';
end
yrTSIA_stdInit = yrTSIA_stdInit + (1/(2*nx))*ones(1, length(tr));

fprintf(1, '2. Simulate LQO-TSIA (truncated initialization) reduced output.\n');
fprintf(1, '--------------------------------------------------------------\n');

% Opts.
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArTSIA_truncInit, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffred            = @(tr, vrTSIA_truncInit)(ArTSIA_truncInit*vrTSIA_truncInit + BrTSIA_truncInit*[u0(tr);u1(tr)]);
[tr, vrTSIA_truncInit] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r.

yrTSIA_truncInit = ClinrTSIA_truncInit*vrTSIA_truncInit';
for tt = 1:length(tr)
    yrTSIA_truncInit(:,tt) = yrTSIA_truncInit(:,tt) + ...
        vrTSIA_truncInit(tt, :)*MquadrTSIA_truncInit*vrTSIA_truncInit(tt, :)';
end
yrTSIA_truncInit = yrTSIA_truncInit + (1/(2*nx))*ones(1, length(tr));

fprintf(1, '3. Simulate LQO-TSIA (eigs initialization) reduced output.\n');
fprintf(1, '--------------------------------------------------------\n');

% Opts.
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArTSIA_eigsInit, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffred         = @(tr, vrTSIA_BTInit)(ArTSIA_eigsInit*vrTSIA_BTInit + BrTSIA_eigsInit*[u0(tr);u1(tr)]);
[tr, vrTSIA_eigsInit] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r.

yrTSIA_eigsInit = ClinrTSIA_eigsInit*vrTSIA_eigsInit';
for tt = 1:length(tr)
    yrTSIA_eigsInit(:,tt) = yrTSIA_eigsInit(:,tt) + ...
        vrTSIA_eigsInit(tt, :)*MquadrTSIA_eigsInit*vrTSIA_eigsInit(tt, :)';
end
yrTSIA_eigsInit = yrTSIA_eigsInit + (1/(2*nx))*ones(1, length(tr));

fprintf(1, '4. Simulate LQO-BT reduced output.\n');
fprintf(1, '----------------------------------\n');

% Opts.
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', ArBT, ...
                  'MStateDependence', 'none', 'Stats','on');
% Simulation.
fAdvDiffredbt = @(tr,vrBT)(ArBT*vrBT + BrBT*[u0(tr);u1(tr)]);
[tr, vrBT]    = ode15s(fAdvDiffredbt, tsteps, vr0, options);
% Note, vr is nt \times nxr

yrBT = ClinrBT*vrBT';
for tt = 1:length(tr)
    yrBT(:,tt) = yrBT(:,tt) + vrBT(tt, :)*MquadrBT*vrBT(tt, :)';
end
yrBT = yrBT + (1/(2*nx))*ones(1, length(tr));

% Plots.
plot(t, yrTSIA_stdInit,   '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
plot(t, yrTSIA_truncInit, '-.', 'color', ColMat(3,:), LineWidth=1.5); 
plot(t, yrTSIA_eigsInit,  '--', 'color', ColMat(4,:), LineWidth=1.5); 
plot(t, yrBT,             '-.', 'color', ColMat(5,:), LineWidth=1.5); 
lgd = legend('$y(t)$', '$y_{r,TSIAstd}(t)$', '$y_{r,TSIAtrunc}(t)$', '$y_{r,TSIAeigs}(t)$', ...
    '$y_{r,BT}(t)$', 'interpreter','latex', 'FontName', 'Arial', 'location', ...
    'northeast');
fontsize(lgd, 10, 'points')

subplot(2,1,2)
semilogy(tr, abs(y - yrTSIA_stdInit)./abs(y),   '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
semilogy(tr, abs(y - yrTSIA_truncInit)./abs(y), '-.', 'color', ColMat(3,:), LineWidth=1.5); 
semilogy(tr, abs(y - yrTSIA_eigsInit)./abs(y),  '--', 'color', ColMat(4,:), LineWidth=1.5); 
semilogy(tr, abs(y - yrBT)./abs(y),             '-.', 'color', ColMat(5,:), LineWidth=1.5); 
xlabel('$t$','interpreter','latex'); 
ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
    LineWidth=1.5)

% Print errors.
fprintf(1, 'Order r = %d.\n', r)
fprintf(1, '--------------\n')
fprintf(1, 'Exponential input: Relative L-infty error due to LQO-TSIA (standard initialization) : %.16f \n', max(abs(y - yrTSIA_stdInit)./abs(y)))
fprintf(1, 'Exponential input: Relative L-infty error due to LQO-TSIA (truncated initialization): %.16f \n', max(abs(y - yrTSIA_truncInit)./abs(y)))
fprintf(1, 'Exponential input: Relative L-infty error due to LQO-TSIA (eigs initialization)     : %.16f \n', max(abs(y - yrTSIA_eigsInit)./abs(y)))
fprintf(1, 'Exponential input: Relative L-infty error due to LQO-BT                             : %.16f \n', max(abs(y - yrBT)./abs(y)))
fprintf(1, '------------------------------------------------------------\n')
fprintf(1, 'Exponential input: Relative L-2 error due to LQO-TSIA (standard initialization) : %.16f \n', sum(abs(y - yrTSIA_stdInit))/sum(abs(y)))
fprintf(1, 'Exponential input: Relative L-2 error due to LQO-TSIA (truncated initialization): %.16f \n', sum(abs(y - yrTSIA_truncInit))/sum(abs(y)))
fprintf(1, 'Exponential input: Relative L-2 error due to LQO-TSIA (eigs initialization)     : %.16f \n', sum(abs(y - yrTSIA_eigsInit))/sum(abs(y)))
fprintf(1, 'Exponential input: Relative L-2 error due to LQO-BT                             : %.16f \n', sum(abs(y - yrBT))/sum(abs(y)))
fprintf(1, '------------------------------------------------------------\n')

% Write data.
write = 1;
if write
    % Overwrite figure.
    print -depsc2 results/AdvecDiff3000_exponential_r30_OutputPlot

    outputs = [t, y', yrTSIA_stdInit', yrTSIA_truncInit', yrTSIA_eigsInit', yrBT'];
    dlmwrite('results/AdvecDiff3000_exponential_r30_Outputs.dat', outputs, 'delimiter', ...
        '\t', 'precision', 8);
    outputerrors = [t, (abs(y-yrTSIA_stdInit)./abs(y))', (abs(y-yrTSIA_truncInit)./abs(y))', ...
        (abs(y-yrTSIA_eigsInit)./abs(y))', (abs(y-yrBT)./abs(y))'];
    dlmwrite('results/AdvecDiff3000_exponential_r30_OutputErrors.dat', outputerrors, ...
        'delimiter', '\t', 'precision', 8);
end

fprintf(1, '\n');


%% Convergence study.
fprintf(1, 'Convergence plotting.\n');
fprintf(1, '---------------------\n');

fprintf(1, 'PRECOMPUTING FOM H2-NORM FOR EXACT CONVERGENCE MONITORING.\n')
fprintf(1, '---------------------------------------------------------.\n')
P       = lyap(A, B*B');
Q       = lyap(A', Clin'*Clin + Mquad*P*Mquad);
fomNorm = abs(trace(B'*Q*B));

% 1. Standard (diagonal) initialization.
% Reset input opts.
opts_stdInit                = struct();
opts_stdInit.tol            = 10e-14;
opts_stdInit.maxIter        = 200; 
opts_stdInit.convMonitoring = 'exact';
opts_stdInit.fomNorm        = fomNorm;

% Initial model values are vanilla; diagonal matrix spread among
% full-order model poles, and identity input, output matrices.
opts_stdInit.Ar = -diag(logspace(0, 4, r));
opts_stdInit.Br = eye(r, m);
opts_stdInit.Cr = eye(p, r);
opts_stdInit.Mr = eye(r, r);

[~, ~, ~, ~, infoTSIA_stdInit] = mimo_lqotsia(A, B, Clin, Mquad, r, opts_stdInit);

% 2. Truncated initialization.

% Reset iinput opts.
opts_truncInit                = struct();
opts_truncInit.tol            = 10e-14;
opts_truncInit.maxIter        = 200; 
opts_truncInit.convMonitoring = 'exact';
opts_truncInit.fomNorm        = fomNorm;

% Initial model values are obtained by truncating full-order matrices.
opts_truncInit.Ar = A(1:r, 1:r);
opts_truncInit.Br = B(1:r, :);
opts_truncInit.Cr = Clin(:, 1:r);
opts_truncInit.Mr = Mquad(1:r, 1:r);

[~, ~, ~, ~, infoTSIA_truncInit] = mimo_lqotsia(A, B, Clin, Mquad, r, opts_truncInit);
 
% 3. eigs initialization.

% Reset input opts.
opts_eigsInit                = struct();
opts_eigsInit.tol            = 10e-14;
opts_eigsInit.maxIter        = 200; 
opts_eigsInit.convMonitoring = 'exact';
opts_eigsInit.fomNorm        = fomNorm;

% Initial model values are obtained from projecting the dominant
% eigenvectors.
[V, D, ~] = eigs(A, r, 'bothendsreal', eigsOpts);
[V, ~]    = qr(V, 'econ');

opts_eigsInit.Ar = (V.'*A*V);
opts_eigsInit.Br = (V.'*B);
opts_eigsInit.Cr = Clin*V;
opts_eigsInit.Mr = V.'*Mquad*V;

[~, ~, ~, ~, infoTSIA_eigsInit] = mimo_lqotsia(A, B, Clin, Mquad, r, opts_eigsInit);

% Standard initialization.
changeInErrors_stdInit = infoTSIA_stdInit.changeInErrors;
changeInTails_stdInit  = infoTSIA_stdInit.changeInTails;
maxIter_stdInit        = max(size(changeInTails_stdInit));

% Truncated initialization.
changeInErrors_truncInit = infoTSIA_truncInit.changeInErrors;
changeInTails_truncInit  = infoTSIA_truncInit.changeInTails;
maxIter_truncInit        = max(size(changeInTails_truncInit));

% eigs initialization.
changeInErrors_eigsInit = infoTSIA_eigsInit.changeInErrors;
changeInTails_eigsInit  = infoTSIA_eigsInit.changeInTails;
maxIter_eigsInit        = max(size(changeInTails_eigsInit));

figure(3)
set(gca, 'fontsize', 10)
semilogy(1:maxIter_stdInit,   changeInErrors_stdInit,   '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
semilogy(1:maxIter_stdInit,   changeInTails_stdInit,    '-*', 'color', ColMat(2,:), LineWidth=1.5);
semilogy(1:maxIter_truncInit, changeInErrors_truncInit, '-o', 'color', ColMat(3,:), LineWidth=1.5);
semilogy(1:maxIter_truncInit, changeInTails_truncInit,  '-*', 'color', ColMat(3,:), LineWidth=1.5);
semilogy(1:maxIter_eigsInit,  changeInErrors_eigsInit,  '-o', 'color', ColMat(4,:), LineWidth=1.5); 
semilogy(1:maxIter_eigsInit,  changeInTails_eigsInit,   '-*', 'color', ColMat(4,:), LineWidth=1.5);

xMax =  max([maxIter_stdInit, maxIter_truncInit, maxIter_eigsInit]);
xlim([2, xMax])
xlabel('iteration count, $k$', 'interpreter', 'latex')

lgd = legend('Std (errors)', 'Std (tails)', 'Truncated (errors)', 'Truncated (tails)', ...
    'eigs (errors)', 'eigs (tails)', 'interpreter','latex', 'FontName', ...
    'Arial', 'location', 'northeast');
fontsize(lgd, 10, 'points')

% Write data
write = 1;
if write
    print -depsc2 results/AdvecDiff3000_r30_conv

    convStdInit   = [(1:maxIter_stdInit)', changeInErrors_stdInit', changeInTails_stdInit'];
    dlmwrite('results/AdvecDiff3000_r30_conv_stdInit.dat', convStdInit, 'delimiter', '\t', 'precision', ...
        12);
    convTruncInit = [(1:maxIter_truncInit)', changeInErrors_truncInit', changeInTails_truncInit'];
    dlmwrite('results/AdvecDiff3000_r30_conv_truncInit.dat', convTruncInit, 'delimiter', '\t', 'precision', ...
        12);
    conveigsInit    = [(1:maxIter_eigsInit)', changeInErrors_eigsInit', changeInTails_eigsInit'];
    dlmwrite('results/AdvecDiff3000_r30_conv_eigsInit.dat', conveigsInit, 'delimiter', '\t', 'precision', ...
        12);
end

%% Hierarchy of reduced models.
fprintf(1, 'Computing hiearchy of reduced models using LQO-TSIA and LQO-BT.\n')
fprintf(1, '---------------------------------------------------------------\n')

rmax                 = 30;
errorsTSIA_stdInit   = zeros(rmax/2, 1);
errorsTSIA_truncInit = zeros(rmax/2, 1);
errorsTSIA_BTInit    = zeros(rmax/2, 1);
errorsBT             = zeros(rmax/2, 1);

% Full model H2 norm saved from previous section.
% Use Gramians to pre-compute BT-MOR bases.
try
    U = chol(P);   
    U = U';
catch
    % If Gramians are not SPD due to roundoff.
    U = chol(P + eye(n, n)*10e-10);   
    U = U';
end
try
    L = chol(Q);   
    L = L';
catch
    % If Gramians are not SPD due to roundoff.
    L = chol(Q + eye(n, n)*10e-10);    
    L = L';
end
[Z, S, Y] = svd(U'*L);

for r = 2:2:rmax  
    fprintf(1, 'Current reduced model in hierarchy: r = %d.\n', r)
    fprintf(1, '-------------------------------------------\n')
    k = r/2; % Iterate counter

    fprintf(1, '1. Computing reduced model using LQO-TSIA (standard initialization).\n');
    fprintf(1, '---------------------------------------------------------------------\n');

    % Input opts.
    opts_stdInit                = struct();
    opts_stdInit.tol            = 10e-14;
    opts_stdInit.maxIter        = 200; 
    opts_stdInit.convMonitoring = 'exact';
    opts_stdInit.fomNorm        = fomNorm;

    % Initial model values are vanilla; diagonal matrix spread among
    % full-order model poles, and identity input, output matrices.
    opts_stdInit.Ar = -diag(logspace(0,4,r));
    opts_stdInit.Br = eye(r, m);
    opts_stdInit.Cr = eye(p, r);
    opts_stdInit.Mr = eye(r, r);
    [ArTSIA_stdInit, BrTSIA_stdInit, ClinrTSIA_stdInit, MquadrTSIA_stdInit, infoTSIA_stdInit] ...
        = mimo_lqotsia(A, B, Clin, Mquad, r, opts_stdInit);

    % Squared H2 errors relative to squared H2 error of FOM are computed
    % throughout iteration.
    errorsTSIA_stdInit(k) = infoTSIA_stdInit.errors(end);
    fprintf(1, 'Squared H2 error due to LQOTSIA (std) relative to that of the FOM is: ||G - Gr||^2_H2/||G||^2_H2 = %.16f\n', ...
        errorsTSIA_stdInit(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n');

    fprintf(1, '2. Computing reduced model using LQO-TSIA (truncated initialization).\n');
    fprintf(1, '---------------------------------------------------------------------\n');

    % Input opts.
    opts_truncInit                = struct();
    opts_truncInit.tol            = 10e-14;
    opts_truncInit.maxIter        = 200; 
    opts_truncInit.convMonitoring = 'exact';
    opts_truncInit.fomNorm        = fomNorm;

    % Initial model values are obtained by truncating full-order matrices.
    opts_truncInit.Ar = A(1:r, 1:r);
    opts_truncInit.Br = B(1:r, :);
    opts_truncInit.Cr = Clin(:, 1:r);
    opts_truncInit.Mr = Mquad(1:r, 1:r);
    [ArTSIA_truncInit, BrTSIA_truncInit, ClinrTSIA_truncInit, MquadrTSIA_truncInit, infoTSIA_truncInit] ...
        = mimo_lqotsia(A, B, Clin, Mquad, r, opts_truncInit);

    % Squared H2 errors relative to squared H2 error of FOM are computed
    % throughout iteration.
    errorsTSIA_truncInit(k) = infoTSIA_truncInit.errors(end);
    fprintf(1, 'Squared H2 error due to LQOTSIA (truncated) relative to that of the FOM is: ||G - Gr||^2_H2/||G||^2_H2 = %.16f\n', ...
        errorsTSIA_truncInit(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n');

    fprintf(1, '3. Computing reduced model using LQO-TSIA (eigs initialization).\n');
    fprintf(1, '---------------------------------------------------------------------\n');

    % Input opts.
    opts_eigsInit                = struct();
    opts_eigsInit.tol            = 10e-14;
    opts_eigsInit.maxIter        = 200; 
    opts_eigsInit.convMonitoring = 'exact';
    opts_eigsInit.fomNorm        = fomNorm;

    % Initial model values are obtained from projecting the dominant
    % eigenvectors.
    [V, D, ~] = eigs(A, r, 'bothendsreal', eigsOpts);
    [V, ~]    = qr(V, 'econ');

    opts_eigsInit.Ar = (V.'*A*V);
    opts_eigsInit.Br = (V.'*B);
    opts_eigsInit.Cr = Clin*V;
    opts_eigsInit.Mr = V.'*Mquad*V;
    [ArTSIA_eigsInit, BrTSIA_eigsInit, ClinrTSIA_eigsInit, MquadrTSIA_eigsInit, infoTSIA_eigsInit] ...
        = mimo_lqotsia(A, B, Clin, Mquad, r, opts_eigsInit);
    
    % Squared H2 errors relative to squared H2 error of FOM are computed
    % throughout iteration.
    errorsTSIA_BTInit(k) = infoTSIA_eigsInit.errors(end);
    fprintf(1, 'Squared H2 error due to LQOTSIA (BT) relative to that of the FOM is: ||G - Gr||^2_H2/||G||^2_H2 = %.16f\n', ...
        errorsTSIA_BTInit(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n');

    fprintf(1, '4. Computing reduced model using LQO-BT.\n');
    fprintf(1, '----------------------------------------\n');
    
    % Compute projection matrices; pre-computed factors.
    V = U*Z(:, 1:r)*S(1:r, 1:r)^(-1/2); % Right
    W = L*Y(:, 1:r)*S(1:r, 1:r)^(-1/2); % Left
    
    % Compute reduced order model via projection.
    ArBT = W'*A*V;   BrBT = W'*B;  ClinrBT = Clin*V;   MquadrBT = V'*Mquad*V;

    % H2 error using trace formula.
    AerrBT         = [A, zeros(nx, r); zeros(r, nx), ArBT]; 
    BerrBT         = [B; BrBT];
    CerrBT         = [Clin, - ClinrBT]; 
    MerrBT         = [Mquad, zeros(nx, r); zeros(r, nx), -MquadrBT];
    PerrBT         = lyap(AerrBT, BerrBT*BerrBT');                  
    QerrBT         = lyap(AerrBT', CerrBT'*CerrBT + MerrBT*PerrBT*MerrBT);
    errorsBT(k, :) = abs(trace(BerrBT'*QerrBT*BerrBT))/fomNorm;
    fprintf(1, 'Squared H2 error due to LQOBT relative to that of the FOM is: ||G - Gr||^2_H2/||G||^2_H2 = %.16f\n', ...
        errorsBT(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')
end

write = 1;
if write
    relh2errors = [(2:2:rmax)', errorsTSIA_stdInit, errorsTSIA_truncInit, errorsTSIA_BTInit, errorsBT];
    dlmwrite('results/AdvecDiff3000_H2errors.dat', relh2errors, 'delimiter', '\t', 'precision', ...
        8);
end


%% Plot H2 errors.
fprintf(1, 'Plotting H2 errors\n');
fprintf(1, '------------------\n');

figure(4)
set(gca, 'fontsize', 10)
semilogy(2:2:rmax, errorsTSIA_stdInit(1:rmax/2),   '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
semilogy(2:2:rmax, errorsTSIA_truncInit(1:rmax/2), '-o', 'color', ColMat(3,:), LineWidth=1.5)
semilogy(2:2:rmax, errorsTSIA_BTInit(1:rmax/2),    '-o', 'color', ColMat(4,:), LineWidth=1.5)
semilogy(2:2:rmax, errorsBT(1:rmax/2),             '-*', 'color', ColMat(5,:), LineWidth=1.5)
xlim([2,rmax])
lgd = legend('LQO-TSIA (std)', 'LQO-TSIA (trunc)', 'LQO-TSIA (eigs)', 'LQO-BT', 'interpreter', 'latex', 'FontName', 'Arial',...
    'location', 'northeast');
xlabel('$r$', 'interpreter', 'latex')
ylabel('$||\mathcal{G}-\mathcal{G}_r||_{\mathcal{H}_2}^2/||\mathcal{G}||_{\mathcal{H}_2}^2$',...
    'interpreter', 'latex')

print -depsc2 results/AdvecDiff3000_r30_H2errors_Plot


%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
