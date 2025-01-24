%% RUNME_SOBUTTERFLY_IRKA
% Script file to run all experiments involving the Butterfly Gyroscope
% model with a quadratic output function.
%

%
% This file is part of the archive Code, Data, and Results for Numerical 
% Experiments in "..."
% Copyright (c) 2025 Sean Reiter,
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

clc;
clear all;
close all;

% Get and set all paths
[rootpath, filename, ~] = fileparts(mfilename('fullpath'));
loadname            = [rootpath filesep() ...
    'data' filesep() filename];
savename            = [rootpath filesep() ...
    'results' filesep() filename];

% Add paths to drivers and data
addpath([rootpath, '/drivers'])
addpath([rootpath, '/data'])

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
fprintf(1, 'Loading butterfly gyroscope benchmark problem.\n')
fprintf(1, '----------------------------------------------\n')

% From: https://morwiki.mpi-magdeburg.mpg.de/morwiki/index.php/Butterfly_Gyroscope
% Rayleigh Damping: D = alpha*M + beta*K
%   alpha = 0;
%   beta  = 1e-6;
load('data/Butterfly.mat')

% Rename.
Mso = M;    Dso = D;    Kso = K;    Bso = B;    Cso = C;

% Input, output dimensions.
[n, ~] = size(M);

% Quadratic-output (QO) matrix Qfo = C'*C.
outputs       = [3, 6, 9, 12]; % Displacements of electrodes in z direction.
Qfo           = spalloc(2*n, 2*n, nnz(C(outputs, :)'*C(outputs, :)));
Qfo(1:n, 1:n) = C(outputs, :)'*C(outputs, :);

% % MSD FOR TESTING.
% % MMESS MATLAB toolbox function.
% n1              = 10; 
% alpha           = .002; 
% beta            = alpha; 
% v               = 5;
% [Mso, Dso, Kso] = triplechain_MSD(n1, alpha, beta, v);
% 
% [n, ~] = size(Mso);
% Cpso   = ones(1, n);
% Bso    = ones(n, 1);
% 
% % Quadratic-output matrix. 
% Qfo           = spalloc(2*n, 2*n, nnz(Cpso' * Cpso));
% Qfo(1:n, 1:n) = Cpso'*Cpso; 

%% Compute reduced order models.
r  = 10; % Order of approximation.

% Boolean; set `true' to recompute reduced models.
% recompute = true;
recompute = false;
if recompute
    fprintf(1, '1. Computing reduced model using LQO-IRKA (default interpolation data).\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');
      
    tmp       = 1i*logspace(4, 6, r/2)';
    initPoles = ([tmp; conj(flipud(tmp))]);
    initPoles = -1000*ones(r, 1) + initPoles;

    % Input opts; default interpolation data.
    optsIRKA             = struct();
    optsIRKA.poles       = initPoles;
    optsIRKA.tol         = 10e-8;
    optsIRKA.maxIter     = 200; 
    optsIRKA.checkInterp = false;
    morStart             = tic;
    
    [Er_IRKA, Ar_IRKA, Br_IRKA, Qr_IRKA, infoIRKA] = siso_so_lqoirka(...
       Mso, Dso, Kso, Bso, Qfo, r, optsIRKA);
    save('results/Butterfly_IRKA_r30.mat', 'Er_IRKA', 'Ar_IRKA', 'Br_IRKA', ...
            'Qr_IRKA', 'infoIRKA')

    fprintf(1, '1. ROM COMPUTED VIA LQO-IRKA (default interpolation data) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % Write convergence data.
    write = 1;
    if write
        poleChange       = infoIRKA.poleChange';
        [numIterates, ~] = size(poleChange);
        conv             = [(1:numIterates)', poleChange];
        dlmwrite('results/Butterfly_IRKA_r30conv.dat', conv, 'delimiter', ...
            '\t', 'precision', 8);
    end
else 
    fprintf(1, 'Loading reduced-order models.\n');
    fprintf(1, '-----------------------------\n');
    load('results/Butterfly_IRKA_r30.mat')
end

fprintf(1, '\n');


%% Simulate full and reduced-order outputs; sinusoidal input.
fprintf(1, 'Solving full-order problem via ode15\n');
fprintf(1, '----------------------------------------------\n');

% Build first-order realization for simulation.
Efo                   = spalloc(2*n, 2*n, nnz(Mso) + n);            % Descriptor matrix; Efo = [I, 0: 0, Mso]
Efo(1:n, 1:n)         = speye(n);                                   % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = Mso;                                        % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(Kso) + nnz(Dso) + n); % Afo = [0, I; -Kso, -Dso]
Afo(1:n, n+1:2*n)     = speye(n);                                   % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -Kso;                                       % (2, 1) block is -Kso
Afo(n+1:2*n, n+1:2*n) = -Dso;                                       % (2, 2) block is -Dso 

Bfo             = spalloc(2*n, 1, nnz(Bso));                        % Bfo = [0; Bso];
Bfo(n+1:2*n, :) = Bso; 

Tfin = 50;             % Final time for simulation
nt   = round(Tfin*64); % Number of time steps

nx = 2*n; % Twice second-order system dimension

fprintf(1, '1. Sinusoidal input u(t) = 0.5*(cos(pi*t)+1)\n');
fprintf(1, '-------------------------------------------------\n');

fullOrderSolveStart = tic;
% Set system inputs 
u = @(t) 0.5*(cos(pi*t)+1); 

ode_rtol = 1e-6; 
tsteps   = linspace(0, Tfin, nt+1); % Time-steps
v0       = zeros(nx, 1);            % Initial condition 
options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
                  'Mass', Efo, 'Jacobian', Afo, ...
                  'MStateDependence', 'none', 'Stats','on');
f      = @(t,y)((Afo)*y + (Bfo)*u(t));
[t, v] = ode15s(f, tsteps, v0, options);
% Note, v is nt \times nx

fprintf(1, '1. FULL-ORDER SIMULATION (sinusoidal input) FINISHED IN %d s\n', toc(fullOrderSolveStart));
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, 'Simulate full output.\n');
fprintf(1, '--------------------.\n');

% Output has mixed linear and quadratic terms
y = zeros(1, length(t));
for tt = 1:length(t)
    y(:, tt) = v(tt, :)*Qfo*v(tt, :)';
end

% Colors for plotting
ColMat(1,:) = [0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];

figure(2)
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
redOrderSolveStart = tic;

% 1. For LQOIRKA reduced model.
ode_rtol        = 1e-6; 
nxr             = r;
vr0             = zeros(nxr, 1); % Initial condition 
options         = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', Er_IRKA, 'Jacobian', Ar_IRKA, ...
                  'MStateDependence', 'none', 'Stats','on');
fredIRKA     = @(tr, vr_IRKA)((Ar_IRKA)*vr_IRKA + (Br_IRKA)*u(tr));
[tr, vr_IRKA] = ode15s(fredIRKA, tsteps, vr0, options);
% Note, vr is nt \times r

fprintf(1, '1. FULL-ORDER SIMULATION (sinusoidal input) FINISHED IN %d s\n', toc(redOrderSolveStart));
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, 'Simulate IRKA reduced output.\n');
fprintf(1, '-----------------------------\n');
yr_IRKA = zeros(1, length(tr));
for tt = 1:length(tr)
    yr_IRKA(:,tt) = vr_IRKA(tt, :)*Qr_IRKA*vr_IRKA(tt, :)';
end

plotBool = true;
if plotBool
    plot(t, yr_IRKA, '--', 'color',ColMat(2,:), LineWidth=1.5); hold on
    % plot(t, yrBT, '-.', 'color',ColMat(3,:), LineWidth=1.5); 
    % plot(t, yrTangInterp, '-.', 'color',ColMat(4,:), LineWidth=1.5); 
    % plot(tr, yrMixedInterp, '-.', 'color',ColMat(5,:), LineWidth=1.5); 
    lgd = legend('$y(t)$', '$y_r(t)$', '$y_{r,bt}(t)$', '$y_{tang}(t)$', '$y_{mixed}(t)$', 'interpreter','latex', ...
        'FontName', 'Arial', 'location', 'northeast');
    fontsize(lgd, 10, 'points')
    
    subplot(2,1,2)
    semilogy(tr, abs(y - yr_IRKA)./abs(y), '--','color',ColMat(2,:),LineWidth=1.5); 
    hold on;
    % semilogy(tr, abs(y - yrBT)./abs(y), '-.','color',ColMat(3,:),LineWidth=1.5); 
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
        outputs = [t, y', yr_IRKA'];
        dlmwrite('results/Butterfly_sinusoidal_r30_Outputs.dat', outputs, 'delimiter', ...
            '\t', 'precision', 8);
        outputerrors = [t, (abs(y-yr_IRKA)./abs(y))'];
        dlmwrite('results/Butterfly_sinusoidal_r30_OutputErrors.dat', outputerrors, ...
            'delimiter', '\t', 'precision', 8);
    end
end

% Max errors.
data             = load('results/Butterfly_sinusoidal_r30_OutputErrors.dat');
errorIRKA        = data(:, 2);
% errorBT          = data(:, 3);
% errorTangInterp  = data(:, 4);
% errorMixedInterp = data(:, 5);
fprintf(1, 'Max error due to LQO-IRKA           : %.10f \n', max(errorIRKA))
% fprintf(1, 'Max error due to LQO-BT             : %.10f \n', max(errorBT))
% fprintf(1, 'Max error due to interp (tangential): %.10f \n', max(errorTangInterp))
% fprintf(1, 'Max error due to interp (mixed      : %.10f \n', max(errorMixedInterp))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, '\n');

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off