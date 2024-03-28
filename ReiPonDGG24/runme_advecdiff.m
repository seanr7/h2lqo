%% RUNME_ADVECDIFF
% Script file to run all experiments on the advection diffusion problem
% with a quadratic cost function as the quantity of interest.

%
% This file is part of the archive Code, Data, and Results for Numerical 
% Experiments in "..."
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
addpath([rootpath, '/QBDynamicsQBOutput_MATLAB'])

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

% TODO: Check best practice for using someone else's code in code you want
% to publish ... 
nx   = 300;  % No. of spatial grid points
diff = 1e-2; % Diffusion parameter
adv  = 1;    % Advection parameter

% Advection-diffusion problem taken from [Diaz et al., 2024]
[E, A, B, ~] = AdvDiff_Probgen_1D(nx, adv, diff);  
E            = full(E);     % E = I in this example
A            = full(A);    
B            = full(B);   
Clin         = eye(nx, nx); % Linear output matrix
Mquad        = eye(nx, nx); % Quadratic output matrix

fprintf(1, '\n');


%% Run ode15 and simulate output.
fprintf(1, 'Solve AdvecDiff problem via ode15\n');
fprintf(1, '---------------------------------\n');

Tfin = 2;               % final time for simulation
nt   = round(Tfin*128); % number of time steps

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
h      = 1/nx;
output = (h/2) * sum((Clin*v' - ones(nx,1)*ones(1,length(t)) ).^2, 1);

fprintf(1, 'Simulate full output.\n');
fprintf(1, '--------------------.\n');

% Output has mixed linear and quadratic terms
Clin = -h*ones(1, nx);
y    = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + (h/2)*v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, nt + 1);

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


%% Run algorithm.
fprintf(1, 'Compute reduced-order models using lqo-tsia.\n');
fprintf(1, '-------------------------------------------\n');

% Set approximation order
r = 25;

% Run two-sided iteration
[Ar, Br, Clinr, Mquadr, pole_history] = mimolqo_tsia(A, B, Clin, Mquad, ...
    r);

% Benchmark against balanced truncation approach
[Arbt, Brbt, Clinrbt, Mquadrbt, info] = mimolqo_bt(A, B, Clin, Mquad, r);


%% Run ode15 and simulate reduced output.
fprintf(1, 'Solve reduced AdvecDiff problem via ode15\n');
fprintf(1, '-----------------------------------------\n');

% For tsia reduced model
ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', Ar, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffred = @(tr,yr)(Ar*yr + Br*[u0(tr);u1(tr)]);
[tr, vr]    = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times nxr

fprintf(1, 'Simulate tsia reduced output\n');
fprintf(1, '----------------------------\n');
yr = Clinr*vr';
for tt = 1:length(tr)
    yr(:,tt) = yr(:,tt) + (h/2)*vr(tt, :)*Mquadr*vr(tt, :)';
end
yr = yr + (1/2) * ones(1, nt+1);

% For bt reduced model
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r,r), 'Jacobian', Arbt, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffred = @(tr,yrbt)(Arbt*yrbt + Brbt*[u0(tr);u1(tr)]);
[tr, vrbt]    = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times nxr

fprintf(1, 'Simulate bt reduced output\n');
fprintf(1, '--------------------------\n');
yrbt = Clinrbt*vrbt';
for tt = 1:length(tr)
    yrbt(:,tt) = yrbt(:,tt) + (h/2)*vrbt(tt, :)*Mquadrbt*vrbt(tt, :)';
end
yrbt = yrbt + (1/2) * ones(1, nt+1);

plot(t, yr, '-', 'color',ColMat(3,:), LineWidth=1.5); hold on
plot(t, yrbt, '-.', 'color',ColMat(4,:), LineWidth=1.5); 
lgd = legend('$y(t)$', '$y_r(t)$', '$y_{r,bt}(t)$', 'interpreter','latex','FontName','Arial',...
    'location', 'northeast');
fontsize(lgd, 10, 'points')

subplot(2,1,2)
plot(t, abs(y - yr)./abs(y), '-.','color',ColMat(3,:),LineWidth=1.5); 
hold on;
plot(t, abs(y - yrbt)./abs(y), '--.','color',ColMat(4,:),LineWidth=1.5); 
xlabel('$t$','interpreter','latex'); 
ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
    LineWidth=1.5)

% Overwrite figure
% saveas(figure(1), 'results/advecdiff_plots.png')
print -depsc2 results/advecdiff_output_plots

% Write data
write = 1;
if write
    outputs = [t, y', yr', yrbt'];
    dlmwrite('results/r25_outputs.dat', outputs, 'delimiter', '\t', 'precision', ...
        8);
    outputerrors = [t, (abs(y-yr)./abs(y))', (abs(y-yrbt)./abs(y))'];
    dlmwrite('results/r25_outputerrors.dat', outputerrors, ...
        'delimiter', '\t', 'precision', 8);
end


%% Plot convergence of poles.
fprintf(1, 'Plotting convergence of the method\n');
fprintf(1, '----------------------------------\n');
maxiter = max(size(pole_history));  pole_change = zeros(maxiter, 1);
for k = 2:maxiter
    pole_change(k) = abs(max(pole_history(:, k) - pole_history(:,k-1)));
end

figure(2)
set(gca, 'fontsize', 10)
semilogy(1:maxiter, pole_change, '-o', 'color',ColMat(3,:), LineWidth=1.5)
xlim([2,maxiter])
% golden_ratio = (sqrt(5)+1)/2;
% axes('position', [.125 .15 .75 golden_ratio-1])
xlabel('$k$', 'interpreter', 'latex')

print -depsc2 results/advecdiff_conv_plots

% Write data
write = 1;
if write
    conv = [(1:maxiter)', pole_change];
    dlmwrite('results/r25_conv.dat', conv, 'delimiter', '\t', 'precision', ...
        8);
end

%% Compute H2 errors for hierarchy of reduced models.
fprintf(1, 'Computing hiearchy of reduced models using lqo-tsia\n')
fprintf(1, '---------------------------------------------------\n')

rmax       = 30;
H2errors   = zeros(rmax/2, 1);
H2errorsbt = zeros(rmax/2, 1);
% Precompute H2 norm of full-order model for relative error
P        = lyap(A, B*B');
Q        = lyap(A', Clin'*Clin + Mquad*P*Mquad);
fom_norm = sqrt(abs(trace(B'*Q*B)));

% Toggle input opts
options      = struct();
opts.tol     = 10e-8;
opts.maxiter = 200;
for r = 2:2:rmax  
    fprintf(1, 'Current reduced model is hierarchy; r=%d\n', r)
    fprintf(1, '----------------------------------------\n')

    % tsia reduced models
    [Ar, Br, Clinr, Mquadr, poles] = mimolqo_tsia(A, B, Clin, Mquad, r, opts);

    % Build error realization, and compute H2 error using Gramian
    % formulation
    Aerr        = [A, zeros(nx, r); zeros(r, nx), Ar]; 
    Berr        = [B; Br];
    Cerr        = [Clin, - Clinr]; 
    Merr        = [Mquad, zeros(nx, r); zeros(r, nx), -Mquadr];
    Perr        = lyap(Aerr, Berr*Berr');                   % Error reachability Gramian
    Qerr        = lyap(Aerr', Cerr'*Cerr + Merr*Perr*Merr); % Error observability Gramian
    k           = r/2;  
    H2errors(k) = sqrt(abs(trace(Berr'*Qerr*Berr)))/fom_norm;
    fprintf(1, 'H2 error of current tsua reduced model is ||G-Gr||_H2 = %.16f\n', ...
        H2errors(k));

    % bt reduced models
    [Arbt, Brbt, Clinrbt, Mquadrbt, info] = mimolqo_bt(A, B, Clin, Mquad, ...
        r);

    % Build error realization, and compute H2 error using Gramian
    % formulation
    Aerrbt        = [A, zeros(nx, r); zeros(r, nx), Arbt]; 
    Berrbt        = [B; Brbt];
    Cerrbt        = [Clin, - Clinrbt]; 
    Merrbt        = [Mquad, zeros(nx, r); zeros(r, nx), -Mquadrbt];
    Perrbt        = lyap(Aerrbt, Berrbt*Berrbt');                  
    Qerrbt        = lyap(Aerrbt', Cerrbt'*Cerrbt + Merrbt*Perrbt*Merrbt);
    k             = r/2;  
    H2errorsbt(k) = sqrt(abs(trace(Berrbt'*Qerrbt*Berrbt)))/fom_norm;
    fprintf(1, 'H2 error of current bt reduced model is ||G-Gr||_H2 = %.16f\n', ...
        H2errorsbt(k));
end

write = 1;
if write
    h2errors = [(2:2:rmax)', H2errors, H2errorsbt];
    dlmwrite('results/h2errors.dat', h2errors, 'delimiter', '\t', 'precision', ...
        8);
end


%% Plot H2 errors.
fprintf(1, 'Plotting H2 errors\n');
fprintf(1, '------------------\n');

figure(3)
set(gca, 'fontsize', 10)
% golden_ratio = (sqrt(5)+1)/2;
% axes('position', [.125 .15 .75 golden_ratio-1])
semilogy(2:2:rmax, H2errors(1:rmax/2), '-o', 'color', ColMat(3,:), LineWidth=1.5)
hold on;
semilogy(2:2:rmax, H2errorsbt(1:rmax/2), '-*', 'color', ColMat(4,:), LineWidth=1.5)
xlim([2,rmax])
lgd = legend('lqo-tsia', 'lqo-bt', 'interpreter', 'latex', 'FontName', 'Arial',...
    'location', 'northeast');
xlabel('$r$', 'interpreter', 'latex')
ylabel('$||\mathcal{G}-\mathcal{G}_r||_{\mathcal{H}_2}/||\mathcal{G}||_{\mathcal{H}_2}$',...
    'interpreter', 'latex')

print -depsc2 results/advecdiff_h2errors_plots


%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
