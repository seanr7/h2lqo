%% RUNME_ADVECDIFF
% Script file to run all experiments on the advection diffusion problem
% with a quadratic cost function as the quantity of interest.

%
% This file is part of the archive Code, Data, and Results for Numerical 
% Experiments in "$\mathcal{H}_2$ optimal model reduction of linear 
% systems with multiple quadratic outputs".
% Copyright (c) 2024 Sean Reiter, ...
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
% nx   = 300;   No. of spatial grid points
% diff = 1e-2;  Diffusion parameter
% adv  = 1;     Advection parameter

load('results/AdvecDiff_n300.mat')

fprintf(1, '\n');


%% Compute reduced-order models.

r = 30; % Approximation order for both approacbes

% Set to true to recompute reduced order models
recompute = 1;
if recompute
    fprintf(1, 'Computing reduced model using LQO-TSIA.\n');
    fprintf(1, '---------------------------------------\n');
    % Input opts
    opts         = struct();
    opts.tol     = 10e-14;
    opts.maxiter = 200; 
    % Linear quadratic output two-sided iteration algoritm (LQO-TSIA)
    [Ar_tsia, Br_tsia, Clinr_tsia, Mquadr_tsia, info_tsia] = mimolqo_tsia(...
        A, B, Clin, Mquad, r, opts);
    save('results/AdvecDiff_lqotsia_r30.mat', 'Ar_tsia', 'Br_tsia', 'Clinr_tsia', ...
        'Mquadr_tsia', 'info_tsia')

    fprintf(1, 'Computing reduced model using LQO-BT.\n');
    fprintf(1, '-------------------------------------\n');
    
    % Linear quadratic output balanced truncation (LQO-BT)
    [Ar_bt, Br_bt, Clinr_bt, Mquadr_bt, info_bt] = mimolqo_bt(A, B, Clin, ...
        Mquad, r);
    save('results/AdvecDiff_lqobt_r30.mat', 'Ar_bt', 'Br_bt', 'Clinr_bt', ...
        'Mquadr_bt', 'info_bt')
else 
    fprintf(1, 'Loading reduced-order models.\n');
    fprintf(1, '-----------------------------\n');
    load('results/AdvecDiff_lqotsia_r30.mat') % Load LQO-TSIA reduced model
    load('results/AdvecDiff_lqobt_r30.mat')   % Load LQO-BT reduced model
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
                  'Mass', E, 'Jacobian', A, ...
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
    y(:, tt) = y(:, tt) + (h/2)*v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, length(t));

% Colors for plotting
ColMat(1,:) = [ 0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [  0.9290    0.6940    0.1250];
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

% For tsia reduced model
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', Ar_tsia, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffred   = @(tr, vr_tsia)(Ar_tsia*vr_tsia + Br_tsia*[u0(tr);u1(tr)]);
[tr, vr_tsia] = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times r

fprintf(1, 'Simulate tsia reduced output\n');
fprintf(1, '----------------------------\n');
yr_tsia = Clinr_tsia*vr_tsia';
for tt = 1:length(tr)
    yr_tsia(:,tt) = yr_tsia(:,tt) + (h/2)*vr_tsia(tt, :)*Mquadr_tsia*vr_tsia(tt, :)';
end
yr_tsia = yr_tsia + (1/2) * ones(1, length(tr));

% For bt reduced model
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', Ar_bt, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffredbt = @(tr,vrbt)(Ar_bt*vrbt + Br_bt*[u0(tr);u1(tr)]);
[tr, vr_bt] = ode15s(fAdvDiffredbt, tsteps, vr0, options);
% Note, vr is nt \times nxr

fprintf(1, 'Simulate bt reduced output\n');
fprintf(1, '--------------------------\n');
yr_bt = Clinr_bt*vr_bt';
for tt = 1:length(tr)
    yr_bt(:,tt) = yr_bt(:,tt) + (h/2)*vr_bt(tt, :)*Mquadr_bt*vr_bt(tt, :)';
end
yr_bt = yr_bt + (1/2) * ones(1, length(tr));

plot(t, yr_tsia, '--', 'color',ColMat(2,:), LineWidth=1.5); hold on
plot(t, yr_bt, '-.', 'color',ColMat(3,:), LineWidth=1.5); 
lgd = legend('$y(t)$', '$y_r(t)$', '$y_{r,bt}(t)$', 'interpreter','latex', ...
    'FontName', 'Arial', 'location', 'northeast');
fontsize(lgd, 10, 'points')

subplot(2,1,2)
semilogy(tr, abs(y - yr_tsia)./abs(y), '--','color',ColMat(2,:),LineWidth=1.5); 
hold on;
semilogy(tr, abs(y - yr_bt)./abs(y), '-.','color',ColMat(3,:),LineWidth=1.5); 
xlabel('$t$','interpreter','latex'); 
ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
    LineWidth=1.5)

% Overwrite figure
print -depsc2 results/AdvecDiff_sinusoidal_r30_OutputPlot

% Write data
write = 1;
if write
    outputs = [t, y', yr_tsia', yr_bt'];
    dlmwrite('results/AdvecDiff_sinusoidal_r30_Outputs.dat', outputs, 'delimiter', ...
        '\t', 'precision', 8);
    outputerrors = [t, (abs(y-yr_tsia)./abs(y))', (abs(y-yr_bt)./abs(y))'];
    dlmwrite('results/AdvecDiff_sinusoidal_r30_OutputErrors.dat', outputerrors, ...
        'delimiter', '\t', 'precision', 8);
end

fprintf(1, '\n');


%% Simulate full and reduced-order outputs; exponentially damped input.
fprintf(1, 'Solving full-order AdvecDiff problem via ode15\n');
fprintf(1, '----------------------------------------------\n');

fprintf(1, 'Second, exponentially damped input u0(t) = exp(-t/5)*t^2\n');
fprintf(1, '--------------------------------------------------------\n');

Tfin   = 30;                      % Final time for simulation
nt     = round(Tfin*10);          % Number of time steps
tsteps = linspace(0, Tfin, nt+1); % Time-steps

% Reset relevant input
u0 = @(t) exp(-t/5)*t^2; % Exponentially damped input on left boundary

options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
                  'Mass', E, 'Jacobian', A, ...
                  'MStateDependence', 'none', 'Stats','on');

% Redefine problem
fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
[t, v]   = ode15s(fAdvDiff, tsteps, v0, options);

fprintf(1, 'Simulate full output.\n');
fprintf(1, '--------------------.\n');

% Output has mixed linear and quadratic terms
y = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + (h/2)*v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, length(t));

figure(2)
set(gca, 'fontsize', 10)
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
hold on
grid on
xlabel('$t$','fontsize',fs,'interpreter','latex'); 
ylabel('$y(t)$','fontsize',fs,'interpreter','latex')


fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15\n');
fprintf(1, '--------------------------------------------------\n');

% For tsia reduced model
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', Ar_tsia, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffred   = @(tr, vr_tsia)(Ar_tsia*vr_tsia + Br_tsia*[u0(tr);u1(tr)]);
[tr, vr_tsia] = ode15s(fAdvDiffred, tsteps, vr0, options);

fprintf(1, 'Simulate tsia reduced output\n');
fprintf(1, '----------------------------\n');
yr_tsia = Clinr_tsia*vr_tsia';
for tt = 1:length(tr)
    yr_tsia(:,tt) = yr_tsia(:,tt) + (h/2)*vr_tsia(tt, :)*Mquadr_tsia*vr_tsia(tt, :)';
end
yr_tsia = yr_tsia + (1/2) * ones(1, length(tr));

% For bt reduced model
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', Ar_bt, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffredbt = @(tr,vrbt)(Ar_bt*vrbt + Br_bt*[u0(tr);u1(tr)]);
[tr, vr_bt] = ode15s(fAdvDiffredbt, tsteps, vr0, options);

fprintf(1, 'Simulate bt reduced output\n');
fprintf(1, '--------------------------\n');
yr_bt = Clinr_bt*vr_bt';
for tt = 1:length(tr)
    yr_bt(:,tt) = yr_bt(:,tt) + (h/2)*vr_bt(tt, :)*Mquadr_bt*vr_bt(tt, :)';
end
yr_bt = yr_bt + (1/2) * ones(1, length(tr));

plot(t, yr_tsia, '--', 'color',ColMat(2,:), LineWidth=1.5); hold on
plot(t, yr_bt, '-.', 'color',ColMat(3,:), LineWidth=1.5); 
lgd = legend('$y(t)$', '$y_r(t)$', '$y_{r,bt}(t)$', 'interpreter','latex', ...
    'FontName', 'Arial', 'location', 'northeast');
fontsize(lgd, 10, 'points')

subplot(2,1,2)
semilogy(tr, abs(y - yr_tsia)./abs(y), '--','color',ColMat(2,:),LineWidth=1.5); 
hold on;
semilogy(tr, abs(y - yr_bt)./abs(y), '-.','color',ColMat(3,:),LineWidth=1.5); 
xlabel('$t$','interpreter','latex'); 
ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
    LineWidth=1.5)

% Overwrite figure
print -depsc2 results/AdvecDiff_exponential_r30_OutputPlot

% Write data
write = 1;
if write
    outputs = [t, y', yr_tsia', yr_bt'];
    dlmwrite('results/AdvecDiff_exponential_r30_Outputs.dat', outputs, 'delimiter', ...
        '\t', 'precision', 8);
    outputerrors = [t, (abs(y-yr_tsia)./abs(y))', (abs(y-yr_bt)./abs(y))'];
    dlmwrite('results/AdvecDiff_exponential_r30_OutputErrors.dat', outputerrors, ...
        'delimiter', '\t', 'precision', 8);
end

fprintf(1, '\n');


%% Convergence analysis.
fprintf(1, 'Plotting convergence of the method\n');
fprintf(1, '----------------------------------\n');

% Plot convergence
errs    = info_tsia.errs;
tails   = info_tsia.tails;
maxiter = max(size(errs));

changeinerrs     = zeros(1, maxiter);
changeinerrs(1)  = 1;
changeintails    = zeros(1, maxiter);
changeintails(1) = 1;
for i = 2:maxiter
    changeintails(i) = abs(tails(i) - tails(i-1))/tails(1);
    changeinerrs(i)  = abs(errs(i) - errs(i-1))/errs(1);
end

figure(3)
set(gca, 'fontsize', 10)
semilogy(1:maxiter, errs, '-o', 'color', ColMat(1,:), LineWidth=1.5)
hold on;
semilogy(1:maxiter, changeintails, '-*', 'color', ColMat(2,:), LineWidth=1.5)
% semilogy(1:maxiter, change_intails, '-*', 'color', ColMat(2,:), LineWidth=1.5)
% semilogy(1:maxiter, gradientAr, '-+', 'color', ColMat(3,:), LineWidth=1.5)
xlim([2,maxiter])
xlabel('iteration count, $k$', 'interpreter', 'latex')

print -depsc2 results/AdvecDiff_r30_conv

% Write data
write = 1;
if write
    % conv = [(1:maxiter)', errs', (abs(tails))', gradientAr'];
    conv = [(1:maxiter)', errs', changeintails'];
    dlmwrite('results/AdvecDiff_r30_conv.dat', conv, 'delimiter', '\t', 'precision', ...
        12);
end

%% Hierarchy of reduced models.
fprintf(1, 'Computing hiearchy of reduced models using LQO-TSIA and LQO-BT\n')
fprintf(1, '--------------------------------------------------------------\n')

rmax         = 30;
relerrs_tsia = zeros(rmax/2, 1);
relerrs_bt   = zeros(rmax/2, 1);
% Precompute H2 norm of full-order model for relative error
P        = lyap(A, B*B');
Q        = lyap(A', Clin'*Clin + Mquad*P*Mquad);
fom_norm = sqrt(abs(trace(B'*Q*B)));

% Toggle input opts
opts         = struct();
opts.tol     = 10e-14;
opts.maxiter = 200;
for r = 2:2:rmax  
    fprintf(1, 'Current reduced model is hierarchy; r=%d\n', r)
    fprintf(1, '----------------------------------------\n')
    
    % LQO-TSIA
    [Ar_tsia, Br_tsia, Clinr_tsia, Mquadr_tsia, info_tsia] = mimolqo_tsia(A, B, Clin, Mquad, r, opts);
    k                 = r/2;  
    errs              = info_tsia.errs(end);
    relerrs_tsia(k, :) = sqrt(errs);
    fprintf(1, 'Relative H2 error of LQO-TSIA reduced model is ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relerrs_tsia(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

    % LQO-BT
    [Ar_bt, Br_bt, Clinr_bt, Mquadr_bt, info] = mimolqo_bt(A, B, Clin, Mquad, ...
        r);

    % Build error realization, and compute H2 error using Gramian
    % formulation
    Aerrbt        = [A, zeros(nx, r); zeros(r, nx), Ar_bt]; 
    Berrbt        = [B; Br_bt];
    Cerrbt        = [Clin, - Clinr_bt]; 
    Merrbt        = [Mquad, zeros(nx, r); zeros(r, nx), -Mquadr_bt];
    Perrbt        = lyap(Aerrbt, Berrbt*Berrbt');                  
    Qerrbt        = lyap(Aerrbt', Cerrbt'*Cerrbt + Merrbt*Perrbt*Merrbt);
    k             = r/2;  
    relerrs_bt(k, :) = sqrt(abs(trace(Berrbt'*Qerrbt*Berrbt)))/fom_norm;
    fprintf(1, 'Relative H2 error of LQO-BT reduced model is ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relerrs_bt(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

end

write = 1;
if write
    relh2errors = [(2:2:rmax)', relerrs_tsia, relerrs_bt];
    dlmwrite('results/AdvecDiff_H2errors.dat', relh2errors, 'delimiter', '\t', 'precision', ...
        8);
end


%% Plot H2 errors.
fprintf(1, 'Plotting H2 errors\n');
fprintf(1, '------------------\n');

figure(4)
set(gca, 'fontsize', 10)
semilogy(2:2:rmax, relerrs_tsia(1:rmax/2), '-o', 'color', ColMat(2,:), LineWidth=1.5)
hold on;
semilogy(2:2:rmax, relerrs_bt(1:rmax/2), '-*', 'color', ColMat(3,:), LineWidth=1.5)
xlim([2,rmax])
lgd = legend('lqo-tsia', 'lqo-bt', 'interpreter', 'latex', 'FontName', 'Arial',...
    'location', 'northeast');
xlabel('$r$', 'interpreter', 'latex')
ylabel('$||\mathcal{G}-\mathcal{G}_r||_{\mathcal{H}_2}/||\mathcal{G}||_{\mathcal{H}_2}$',...
    'interpreter', 'latex')

print -depsc2 results/AdvecDiff_r30_H2errors_Plot


%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
