%% RUNME_ADVECDIFF_IRKA
% Script file to run all experiments on the advection diffusion problem
% with a quadratic cost function as the quantity of interest.

%
% This file is part of the archive Code, Data, and Results for Numerical 
% Experiments in "$\mathcal{H}_2$-optimal model reduction of linear
% quadratic-output systems by multivariate rational interpolation".
% Copyright (c) 2025 Sean Reiter
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

% Advection-diffusion problem taken from [Diaz et al., 2024].
% nx   = 3000;  No. of spatial grid points
% diff = 1;     Diffusion parameter
% adv  = 1;     Advection parameter

load('data/AdvecDiff_n3000.mat')
fprintf(1, '\n');

% State, input, output dimensions.
[n, m] = size(B);
[p, ~] = size(Clin);
nx     = 3000;

%% Compute reduced order models.
r = 30; % Order of approximation.
E = speye(n, n);

% Boolean; set `true' to recompute reduced models.
recompute = true;
% recompute = false;
if recompute
    %%
    fprintf(1, '1a. Computing reduced model using LQO-IRKA (eigs initialization).\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % Pick interpolation data from projected eigenproblem using 
    % dominant eigenpairs computed via MATLAB's 'eigs'.
    eigsOpts                = struct();
    eigsOpts.MaxIterations  = 100;
    eigsOpts.Tolerance      = 1e-10;
    [V, D, flag]            = eigs(A, r, 'smallestabs', eigsOpts); % E is identity
    [V, ~]                  = qr(V, 'econ');

    % Projected eigenproblem.
    Ar                    = V.'*A*V; 
    Br                    = V.'*B;    
    Clinr                 = Clin*V; 
    Mquadr                = V.'*Mquad*V;
    [XrIRKA_eigsInit, Lr] = eig(Ar);
    poles_eigsInit        = diag(Lr);

    % Sort.
    sortedMirroredPoles = cplxpair(-poles_eigsInit, 1e-6); 
    [~, sortIndex]      = ismember(sortedMirroredPoles, -(poles_eigsInit));    
    poles_eigsInit      = sortedMirroredPoles;

    % Sort eigenvectors accordingly.
    XrIRKA_eigsInit = XrIRKA_eigsInit(:, sortIndex);
    
    % Residues.
    XrInverse               = XrIRKA_eigsInit\eye(r, r); 
    rightTangents_eigsInit  = XrInverse*Br;
    rightTangents_eigsInit  = rightTangents_eigsInit.';
    loLeftTangents_eigsInit = XrIRKA_eigsInit.'*Clinr.';  
    qoLeftTangents_eigsInit = XrIRKA_eigsInit.'*Mquadr*XrIRKA_eigsInit; 

    morStart = tic;

    % Set input opts.
    optsIRKA_eigsInit                = struct();
    optsIRKA_eigsInit.poles          = poles_eigsInit;
    optsIRKA_eigsInit.rightTangents  = rightTangents_eigsInit;
    optsIRKA_eigsInit.loLeftTangents = loLeftTangents_eigsInit;
    optsIRKA_eigsInit.qoLeftTangents = qoLeftTangents_eigsInit;
    optsIRKA_eigsInit.tol            = 10e-10;
    optsIRKA_eigsInit.maxIter        = 200; 
    optsIRKA_eigsInit.checkInterp    = false;
    optsIRKA_eigsInit.computeH2      = false;

    [ErIRKA_eigsInit, ArIRKA_eigsInit, BrIRKA_eigsInit, ClinrIRKA_eigsInit, ...
        MquadrIRKA_eigsInit, infoIRKA_eigsInit] = miso_lqoirka(E, A, B, Clin.', ...
            Mquad, r, optsIRKA_eigsInit);
    save('results/AdvecDiff3000_IRKA_eigsInit_r30.mat', 'ErIRKA_eigsInit', 'ArIRKA_eigsInit', 'BrIRKA_eigsInit', ...
        'ClinrIRKA_eigsInit', 'MquadrIRKA_eigsInit', 'optsIRKA_eigsInit', 'infoIRKA_eigsInit')

    fprintf(1, '1a. ROM COMPUTED VIA LQO-IRKA (eigs initialization) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    %%
    fprintf(1, '1b. Computing reduced model using LQO-IRKA (imag initialization).\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');

    % Pick interpolation data to be imaginary shifts, 
    % orthonormal residues.

    % Compute shifts.
    tmp            = 1i*logspace(0, 3, r/2)';
    poles_imagInit = ([tmp; conj(flipud(tmp))]);

    % Sort.
    sortedMirroredPoles = cplxpair(-poles_imagInit, 1e-6); 
    [~, sortIndex]      = ismember(sortedMirroredPoles, -poles_imagInit);
    poles_imagInit      = sortedMirroredPoles;

    morStart = tic;

    % Set input opts. (Other inputs come from default options for this
    % initialization.)
    optsIRKA_imagInit                = struct();
    optsIRKA_imagInit.poles          = poles_imagInit;
    optsIRKA_imagInit.tol            = 10e-10;
    optsIRKA_imagInit.maxIter        = 200; 
    optsIRKA_imagInit.checkInterp    = false;
    optsIRKA_imagInit.computeH2      = false;
    
    [ErIRKA_imagInit, ArIRKA_imagInit, BrIRKA_imagInit, ClinrIRKA_imagInit, ...
        MquadrIRKA_imagInit, infoIRKA_imagInit] = miso_lqoirka(E, A, B, Clin.', ...
            Mquad, r, optsIRKA_imagInit);
    save('results/AdvecDiff3000_IRKA_imagInit_r30.mat', 'ErIRKA_imagInit', 'ArIRKA_imagInit', 'BrIRKA_imagInit', ...
        'ClinrIRKA_imagInit', 'MquadrIRKA_imagInit', 'optsIRKA_imagInit', 'infoIRKA_imagInit')

    fprintf(1, '1b. ROM COMPUTED VIA LQO-IRKA (imag initialization) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    %%
    fprintf(1, '2. Computing reduced model using LQO-BT.\n');
    fprintf(1, '-----------------------------------------------------------------------------\n');

    morStart = tic;

    [ArBT, BrBT, ClinrBT, MquadrBT, infoBT] = miso_lqobt(A, B, Clin.', ...
        Mquad, r);
    ErBT                                    = eye(r, r);
    save('results/AdvecDiff3000_BT_r30.mat', 'ErBT', 'ArBT', 'BrBT', 'ClinrBT', ...
        'MquadrBT', 'infoBT')

    fprintf(1, '2. ROM COMPUTED VIA LQO-BT IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    %%
    fprintf(1, '3a. Computing reduced-order model by tangential interpolation (Petrov-Galerkin projection, mixed conditions, eigs initialization as data).\n')
    fprintf(1, '-----------------------------------------------------------------------------\n');
    morStart = tic;

    % Shifts and tangent directions; use data from eigs initialization.
    optsOneStepInterp_eigsInit                = struct();
    optsOneStepInterp_eigsInit.shifts         = -poles_eigsInit;
    optsOneStepInterp_eigsInit.rightTangents  = rightTangents_eigsInit;
    optsOneStepInterp_eigsInit.loLeftTangents = loLeftTangents_eigsInit;                 
    optsOneStepInterp_eigsInit.qoLeftTangents = qoLeftTangents_eigsInit;
    optsOneStepInterp_eigsInit.interpolation  = 'mixed';
    optsOneStepInterp_eigsInit.checkInterp    = false;

    [ErOneStepInterp_eigsInit, ArOneStepInterp_eigsInit, BrOneStepInterp_eigsInit, ...
        ClinrOneStepInterp_eigsInit, MquadrOneStepInterp_eigsInit] = miso_lqointerp(...
            E, A, B, Clin.', Mquad, r, optsOneStepInterp_eigsInit);

    save('results/AdvecDiff3000_OneStepInterp_eigsInit_r30.mat', 'ErOneStepInterp_eigsInit', 'ArOneStepInterp_eigsInit', ...
        'BrOneStepInterp_eigsInit', 'ClinrOneStepInterp_eigsInit', 'MquadrOneStepInterp_eigsInit')

    fprintf(1, '3a. ROM COMPUTED VIA TANGENTIAL INTERPOLATION (Petrov-Galerkin projection, mixed conditions, eigs initialization as data) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

    %%
    fprintf(1, '3b. Computing reduced-order model by tangential interpolation (Petrov-Galerkin projection, mixed conditions, imag. initialization as data).\n')
    fprintf(1, '-----------------------------------------------------------------------------\n');
    morStart = tic;

    % Shifts and tangent directions; use extremal eigenpairs of full-order
    % problem (initialization in IRKA iteration).

    % Sort shifts into complex conjugate pairs.
    poles_imagInit = cplxpair(poles_imagInit); 

    optsOneStepInterp_imagInit                = struct();
    optsOneStepInterp_imagInit.shifts         = -poles_imagInit;
    optsOneStepInterp_imagInit.rightTangents  = eye(m, r);
    optsOneStepInterp_imagInit.loLeftTangents = eye(r, 1);                 
    optsOneStepInterp_imagInit.qoLeftTangents = eye(r, r);
    optsOneStepInterp_imagInit.interpolation  = 'mixed';
    optsOneStepInterp_imagInit.checkInterp    = false;

    [ErOneStepInterp_imagInit, ArOneStepInterp_imagInit, BrOneStepInterp_imagInit, ...
        ClinrOneStepInterp_imagInit, MquadrOneStepInterp_imagInit] = miso_lqointerp( ...
            E, A, B, Clin.', Mquad, r, optsOneStepInterp_imagInit);

    save('results/AdvecDiff3000_OneStepInterp_imagInit_r30.mat', 'ErOneStepInterp_imagInit', 'ArOneStepInterp_imagInit', ...
        'BrOneStepInterp_imagInit', 'ClinrOneStepInterp_imagInit', 'MquadrOneStepInterp_imagInit')

    fprintf(1, '3b. ROM COMPUTED VIA TANGENTIAL INTERPOLATION (Petrov-Galerkin projection, mixed conditions, imag. initialization as data) IN %d s\n', toc(morStart));
    fprintf(1, '-----------------------------------------------------------------------------\n');

else 
    fprintf(1, 'Loading reduced-order models.\n');
    fprintf(1, '-----------------------------\n');
    % Results in paper.
    load('results/AdvecDiff3000_IRKA_eigsInit_r30.mat')
    load('results/AdvecDiff3000_IRKA_imagInit_r30.mat')
    load('results/AdvecDiff3000_BT_r30.mat')   
    load('results/AdvecDiff3000_OneStepInterp_eigsInit_r30.mat') 
    load('results/AdvecDiff3000_OneStepInterp_imagInit_r30.mat') 
end

%% Convergence monitoring.
% Pre-compute H2 norm of the full-order model for convergence monitoring.
normStart = tic;
fprintf(1, 'PRE-COMPUTING H2 NORM OF FULL-ORDER MODEL.\n');
fprintf(1, '--------------------------------------------------\n');
P   = lyap(A, B*B');
rhs = Clin'*Clin + Mquad*P*Mquad; 
Q   = lyap(A', rhs);

% Squared H2 norm of full-order model.
fomH2norm = sqrt(abs(trace(B'*Q*B))); 
fprintf(1, 'H2 NORM COMPUTED IN %d s\n', toc(normStart));
fprintf(1, '--------------------------------------------------\n');

% Re-run IRKA iterations while computing H2 errors. 
computeError = false;
if computeError
    % Input opts loaded with reduced models. 
    %%
    % Set to compute H2 errors throughout.
    optsIRKA_eigsInit.computeH2 = true;
    optsIRKA_eigsInit.Q         = Q;
    optsIRKA_eigsInit.fomH2norm = fomH2norm;
    [~, ~, ~, ~, ~, infoIRKA_eigsInit] = miso_lqoirka(E, A, B, Clin.', ...
            Mquad, r, optsIRKA_eigsInit);
    % Write convergence data.
    write = 1;
    if write
        H2errors         = infoIRKA_eigsInit.H2errors';
        [numIterates, ~] = size(H2errors);
        conv             = [(1:numIterates)', H2errors];
        dlmwrite('results/AdvecDiff3000_IRKA_eigsInit_r30conv.dat', conv, 'delimiter', ...
            '\t', 'precision', 8);
    end

    %%
    % Set to compute H2 errors throughout.
    optsIRKA_imagInit.computeH2 = true;
    optsIRKA_imagInit.Q         = Q;
    optsIRKA_imagInit.fomH2norm = fomH2norm;
    [~, ~, ~, ~, ~, infoIRKA_imagInit] = miso_lqoirka(E, A, B, Clin.', ...
            Mquad, r, optsIRKA_imagInit);
    % Write convergence data.
    write = 1;
    if write
        H2errors         = infoIRKA_imagInit.H2errors';
        [numIterates, ~] = size(H2errors);
        conv             = [(1:numIterates)', H2errors];
        dlmwrite('results/AdvecDiff3000_IRKA_imagInit_r30conv.dat', conv, 'delimiter', ...
            '\t', 'precision', 8);
    end
end
fprintf(1, '\n');

%% Plots (convergence).
% Colors for plotting.
ColMat(1,:) = [0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];

% Plot H2 error convergence (if computed).
if computeError
    figure(1)
    H2errorsIRKA_eigsInit  = infoIRKA_eigsInit.H2errors';
    H2errorsIRKA_imagInit  = infoIRKA_imagInit.H2errors';
    nIterIRKA_eigsInit     = length(H2errorsIRKA_eigsInit);
    nIterIRKA_imagInit     = length(H2errorsIRKA_imagInit);
    set(gca, 'fontsize', 10)

    semilogy(1:nIterIRKA_eigsInit,  H2errorsIRKA_eigsInit',  '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
    semilogy(1:nIterIRKA_imagInit,  H2errorsIRKA_imagInit',  '-x', 'color', ColMat(3,:), LineWidth=1.5);

    legend('LQO-IRKA (eigs init)', 'LQO-IRKA (imag. init)', ...
        'interpreter', 'latex', 'FontName', 'Arial', 'location', 'northeast');
    xlabel('Iteration count', 'interpreter', 'latex')
    ylabel('Change in $||\mathcal{G}-\mathcal{G}_r||_{\mathcal{H}_2}/||\mathcal{G}||_{\mathcal{H}_2}$',...
        'interpreter', 'latex')
end


%% Output simulation 1.
% Sinc input.
Tfin = 10;             % Final time for simulation
nt   = round(Tfin*64); % Number of time steps

fprintf(1, '1. Sinc input u(t) = 5*sinc(t).\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
fullOrderSolveStart = tic;

% Set system inputs. (Note: sinc function uses Signal Processing Toolbox)
u0 = @(t) 5*sinc(t);    
u1 = @(t) zeros(size(t)); % Neumann input on right boundary

fprintf(1, 'Solving full-order AdvecDiff problem via ode15.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');

ode_rtol = 1e-6; 
tsteps   = linspace(0, Tfin, nt+1); % Time-steps
v0       = zeros(nx, 1);            % Initial condition 
options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
                  'Mass', E, 'Jacobian', A, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
[t, v]   = ode15s(fAdvDiff, tsteps, v0, options);
% Note: v is nt \times nx.

fprintf(1, '1. FULL-ORDER SIMULATION (sinc input) FINISHED IN %d s\n', toc(fullOrderSolveStart));
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, 'Simulate full-order output.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
y = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, length(t));

figure(2)
set(gca, 'fontsize', 10)
fs = 12; % Fontsize
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
hold on
grid on
xlabel('$t$',   'fontsize', fs, 'interpreter', 'latex'); 
ylabel('$y(t)$','fontsize', fs, 'interpreter', 'latex')

%%
fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
reducedOrderSolvesStart = tic;

ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); 

% LQO-IRKA.
fprintf(1, 'Simulate IRKA reduced outputs.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
%% 
% 1a. LQO-IRKA (eigs init).
% Eigs init.
options                  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                            'Mass', eye(r,r), 'Jacobian', ErIRKA_eigsInit\ArIRKA_eigsInit, ...
                            'MStateDependence', 'none', 'Stats','on');
fAdvDiffredIRKA_eigsInit = @(tr, vrIRKA_eigsInit)((ErIRKA_eigsInit\ArIRKA_eigsInit)*vrIRKA_eigsInit ...
                             + (ErIRKA_eigsInit\BrIRKA_eigsInit)*[u0(tr);u1(tr)]);

% Solver.
[trIRKA_eigsInit, vrIRKA_eigsInit] = ode15s(fAdvDiffredIRKA_eigsInit, tsteps, vr0, options);

% Reduced output.
yrIRKA_eigsInit = ClinrIRKA_eigsInit.'*vrIRKA_eigsInit.';
for tt = 1:length(trIRKA_eigsInit)
    yrIRKA_eigsInit(:,tt) = yrIRKA_eigsInit(:,tt) + vrIRKA_eigsInit(tt, :)*MquadrIRKA_eigsInit*vrIRKA_eigsInit(tt, :)';
end
yrIRKA_eigsInit = yrIRKA_eigsInit + (1/2)*ones(1, length(trIRKA_eigsInit));

%%
% 1b. LQO-IRKA (imag init).
options                  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                            'Mass', eye(r,r), 'Jacobian', ErIRKA_imagInit\ArIRKA_imagInit, ...
                            'MStateDependence', 'none', 'Stats','on');
fAdvDiffredIRKA_imagInit = @(tr, vrIRKA_imagInit)((ErIRKA_imagInit\ArIRKA_imagInit)*vrIRKA_imagInit ...
                             + (ErIRKA_imagInit\BrIRKA_imagInit)*[u0(tr);u1(tr)]);

% Solver.
[trIRKA_imagInit, vrIRKA_imagInit] = ode15s(fAdvDiffredIRKA_imagInit, tsteps, vr0, options);

% Reduced output.
yrIRKA_imagInit = ClinrIRKA_imagInit.'*vrIRKA_imagInit.';
for tt = 1:length(trIRKA_imagInit)
    yrIRKA_imagInit(:,tt) = yrIRKA_imagInit(:,tt) + vrIRKA_imagInit(tt, :)*MquadrIRKA_imagInit*vrIRKA_imagInit(tt, :)';
end
yrIRKA_imagInit = yrIRKA_imagInit + (1/2)*ones(1, length(trIRKA_imagInit));

%%
% 2. LQO-BT.
fprintf(1, 'Simulate BT reduced output\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
options       = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r, r), 'Jacobian', ArBT, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffredBT = @(tr,vrBT)(ArBT*vrBT + BrBT*[u0(tr);u1(tr)]);

% Solver.
[trBT, vrBT] = ode15s(fAdvDiffredBT, tsteps, vr0, options);

% Reduced output.
yrBT = ClinrBT.'*vrBT.';
for tt = 1:length(trBT)
    yrBT(:,tt) = yrBT(:,tt) + vrBT(tt, :)*MquadrBT*vrBT(tt, :)';
end
yrBT = yrBT + (1/2)*ones(1, length(trBT));

%%
% One-step interpolation.
fprintf(1, 'Simulate tangential interpolation (Petrov-Galerkin, mixed) reduced outputs.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
%%
% 3a. One-step interpolation (eigs init).
options                           = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                                    'Mass', eye(r, r), 'Jacobian', ErOneStepInterp_eigsInit\ArOneStepInterp_eigsInit, ...
                                    'MStateDependence', 'none', 'Stats','on');
fAdvDiffredOneStepInterp_eigsInit = @(tr, vrOneStepInterp_eigsInit)((ErOneStepInterp_eigsInit\ArOneStepInterp_eigsInit)*vrOneStepInterp_eigsInit ...
                                    + (ErOneStepInterp_eigsInit\BrOneStepInterp_eigsInit)*[u0(tr);u1(tr)]);

% Solver.
[trOneStepInterp_eigsInit, vrOneStepInterp_eigsInit] = ode15s(fAdvDiffredOneStepInterp_eigsInit, tsteps, vr0, options);

% Reduced output.
yrOneStepInterp_eigsInit = ClinrOneStepInterp_eigsInit.'*vrOneStepInterp_eigsInit.';
for tt = 1:length(trOneStepInterp_eigsInit)
    yrOneStepInterp_eigsInit(:,tt) = yrOneStepInterp_eigsInit(:,tt) + vrOneStepInterp_eigsInit(tt, :)*MquadrOneStepInterp_eigsInit*vrOneStepInterp_eigsInit(tt, :)';
end
yrOneStepInterp_eigsInit = yrOneStepInterp_eigsInit + (1/2)*ones(1, length(trOneStepInterp_eigsInit));

%%
% 3b. One-step interpolation (imag init).
options                           = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                                    'Mass', eye(r, r), 'Jacobian', ErOneStepInterp_imagInit\ArOneStepInterp_imagInit, ...
                                    'MStateDependence', 'none', 'Stats','on');
fAdvDiffredOneStepInterp_imagInit = @(tr, vrOneStepInterp_imagInit)((ErOneStepInterp_imagInit\ArOneStepInterp_imagInit)*vrOneStepInterp_imagInit ...
                                    + (ErOneStepInterp_imagInit\BrOneStepInterp_imagInit)*[u0(tr);u1(tr)]);

% Solver.
[trOneStepInterp_imagInit, vrOneStepInterp_imagInit] = ode15s(fAdvDiffredOneStepInterp_imagInit, tsteps, vr0, options);

% Reduced output.
yrOneStepInterp_imagInit = ClinrOneStepInterp_imagInit.'*vrOneStepInterp_imagInit.';
for tt = 1:length(trOneStepInterp_imagInit)
    yrOneStepInterp_imagInit(:,tt) = yrOneStepInterp_imagInit(:,tt) + vrOneStepInterp_imagInit(tt, :)*MquadrOneStepInterp_imagInit*vrOneStepInterp_imagInit(tt, :)';
end
yrOneStepInterp_imagInit = yrOneStepInterp_imagInit + (1/2)*ones(1, length(trOneStepInterp_imagInit));

fprintf(1, '1. REDUCED-ORDER SIMULATIONS (sinc input) FINISHED IN %d s\n', toc(reducedOrderSolvesStart));
fprintf(1, '-----------------------------------------------------------------------------\n');

%% Plotting.
% Since input.
plotBool = true;
if plotBool
    plot(trIRKA_eigsInit,           yrIRKA_eigsInit,         '--', 'color', ColMat(2,:), LineWidth=1.5); hold on
    plot(trIRKA_imagInit,           yrIRKA_imagInit,         '-.', 'color', ColMat(2,:), LineWidth=1.5);
    plot(trBT,                      yrBT,                    '--', 'color', ColMat(3,:), LineWidth=1.5); 
    plot(trOneStepInterp_eigsInit,  yrOneStepInterp_eigsInit,  '--', 'color', ColMat(5,:), LineWidth=1.5); 
    plot(trOneStepInterp_imagInit,  yrOneStepInterp_imagInit,  '-.', 'color', ColMat(5,:), LineWidth=1.5); 
    lgd = legend('$y(t)$', '$y_{IRKA,eigs}(t)$', '$y_{IRKA,imag}(t)$', ...
            '$y_{BT}(t)$', '$y_{OneStepInterp,eigs}(t)$', ...
            '$y_{OneStepInterp,imag}(t)$', 'interpreter','latex', 'FontName', 'Arial', ...
            'location', 'northeast');
    fontsize(lgd, 10, 'points')
    
    subplot(2,1,2)  
    semilogy(trIRKA_eigsInit,           abs(y - yrIRKA_eigsInit)./abs(y),           '--', 'color', ColMat(2,:), LineWidth=1.5); hold on;
    semilogy(trIRKA_imagInit,           abs(y - yrIRKA_imagInit)./abs(y),           '-.', 'color', ColMat(2,:), LineWidth=1.5);
    semilogy(trBT,                      abs(y - yrBT)./abs(y),                      '-.', 'color', ColMat(3,:), LineWidth=1.5); 
    semilogy(trOneStepInterp_eigsInit,  abs(y - yrOneStepInterp_eigsInit)./abs(y),  '--', 'color', ColMat(5,:), LineWidth=1.5); 
    semilogy(trOneStepInterp_imagInit,  abs(y - yrOneStepInterp_imagInit)./abs(y),  '-.', 'color', ColMat(5,:), LineWidth=1.5); 
    xlabel('$t$','interpreter','latex'); 
    ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
        LineWidth=1.5)
    
    % Write data
    write = true;
    if write
        outputs = [t, y', yrIRKA_eigsInit', yrIRKA_imagInit', yrBT', ...
             yrOneStepInterp_eigsInit', yrOneStepInterp_imagInit'];
        dlmwrite('results/AdvecDiff3000_sinc_r30_Outputs.dat', outputs, 'delimiter', ...
            '\t', 'precision', 8);
        outputerrors = [t, (abs(y - yrIRKA_eigsInit)./abs(y))', (abs(y - yrIRKA_imagInit)./abs(y))', ...
            (abs(y - yrBT)./abs(y))', (abs(y - yrOneStepInterp_eigsInit)./abs(y))', ...
            (abs(y - yrOneStepInterp_imagInit)./abs(y))'];
        dlmwrite('results/AdvecDiff3000_sinc_r30_OutputErrors.dat', outputerrors, ...
            'delimiter', '\t', 'precision', 8);
    end
end

% Print errors.
fprintf(1, 'Order r = %d.\n', r)
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, 'sinc input: Relative L-infty error due to LQO-IRKA (eigs init)      : %.16f \n', max(abs(y - yrIRKA_eigsInit)./abs(y)))
fprintf(1, 'sinc input: Relative L-infty error due to LQO-IRKA (imag init)      : %.16f \n', max(abs(y - yrIRKA_imagInit)./abs(y)))
fprintf(1, 'sinc input: Relative L-infty error due to LQO-BT                    : %.16f \n', max(abs(y - yrBT)./abs(y)))
fprintf(1, 'sinc input: Relative L-infty error due to interp (mixed, eigs init) : %.16f \n', max(abs(y - yrOneStepInterp_eigsInit)./abs(y)))
fprintf(1, 'sinc input: Relative L-infty error due to interp (mixed, imag init) : %.16f \n', max(abs(y - yrOneStepInterp_imagInit)./abs(y)))
fprintf(1, '-----------------------------------------------------------------------------\n');
fprintf(1, 'sinc input: Relative L-2 error due to LQO-IRKA (eigs init)          : %.16f \n', sqrt(sum(abs(y - yrIRKA_eigsInit).^2)/sum(abs(y).^2)))
fprintf(1, 'sinc input: Relative L-2 error due to LQO-IRKA (imag init)          : %.16f \n', sqrt(sum(abs(y - yrIRKA_imagInit).^2)/sum(abs(y).^2)))
fprintf(1, 'sinc input: Relative L-2 error due to LQO-BT                        : %.16f \n', sqrt(sum(abs(y - yrBT).^2)/sum(abs(y))))
fprintf(1, 'sinc input: Relative L-2 error due to interp (mixed, eigs init)     : %.16f \n', sqrt(sum(abs(y - yrOneStepInterp_eigsInit).^2)/sum(abs(y).^2)))
fprintf(1, 'sinc input: Relative L-2 error due to interp (mixed, imag init)     : %.16f \n', sqrt(sum(abs(y - yrOneStepInterp_imagInit).^2)/sum(abs(y).^2)))
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, '\n');

%% Output simulation 2.
% Exponential input.
fprintf(1, 'Solving full-order AdvecDiff problem via ode15.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');

Tfin = 10;                      % Final time for simulation
nt   = round(Tfin*64);          % Number of time steps

fprintf(1, '2. exponential input u(t) = exp(-t/5)*sin(4*pi*t).\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
fullOrderSolveStart = tic;

% Reset relevant input.
u0 = @(t) exp(-t/5)*sin(4*pi*t); % Exponentially damped input on left boundary

ode_rtol = 1e-6; 
tsteps   = linspace(0, Tfin, nt+1); % Time-steps
v0       = zeros(nx, 1);            % Initial condition 
options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
                  'Mass', E, 'Jacobian', A, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
[t, v]   = ode15s(fAdvDiff, tsteps, v0, options);
% Note, v is nt \times nx
% Compute quadratic cost by 1/(2*n) * ||x(t) - \bf{1}||_2^2.

fprintf(1, '2. FULL-ORDER SIMULATION (exponential input) FINISHED IN %d s\n', toc(fullOrderSolveStart));
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, 'Simulate full-order output.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
y = Clin*v';
for tt = 1:length(t)
    y(:, tt) = y(:, tt) + v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, length(t));

figure(3)
set(gca, 'fontsize', 10)
fs = 12; % Fontsize
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5); hold on; grid on
xlabel('$t$',   'fontsize', fs, 'interpreter', 'latex'); 
ylabel('$y(t)$','fontsize', fs, 'interpreter', 'latex')

%%
fprintf(1, 'Solving reduced-order AdvecDiff problems via ode15.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
reducedOrderSolvesStart = tic;

ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 

% LQO-IRKA.
fprintf(1, 'Simulate IRKA reduced outputs.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
%% 
% 1a. LQO-IRKA (eigs init).
% Eigs init.
options                  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                            'Mass', eye(r,r), 'Jacobian', ErIRKA_eigsInit\ArIRKA_eigsInit, ...
                            'MStateDependence', 'none', 'Stats','on');
fAdvDiffredIRKA_eigsInit = @(tr, vrIRKA_eigsInit)((ErIRKA_eigsInit\ArIRKA_eigsInit)*vrIRKA_eigsInit ...
                             + (ErIRKA_eigsInit\BrIRKA_eigsInit)*[u0(tr);u1(tr)]);

% Solver.
[trIRKA_eigsInit, vrIRKA_eigsInit] = ode15s(fAdvDiffredIRKA_eigsInit, tsteps, vr0, options);

% Reduced output.
yrIRKA_eigsInit = ClinrIRKA_eigsInit.'*vrIRKA_eigsInit.';
for tt = 1:length(trIRKA_eigsInit)
    yrIRKA_eigsInit(:,tt) = yrIRKA_eigsInit(:,tt) + vrIRKA_eigsInit(tt, :)*MquadrIRKA_eigsInit*vrIRKA_eigsInit(tt, :)';
end
yrIRKA_eigsInit = yrIRKA_eigsInit + (1/2)*ones(1, length(trIRKA_eigsInit));

%%
% 1b. LQO-IRKA (imag init).
options                  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                            'Mass', eye(r,r), 'Jacobian', ErIRKA_imagInit\ArIRKA_imagInit, ...
                            'MStateDependence', 'none', 'Stats','on');
fAdvDiffredIRKA_imagInit = @(tr, vrIRKA_imagInit)((ErIRKA_imagInit\ArIRKA_imagInit)*vrIRKA_imagInit ...
                             + (ErIRKA_imagInit\BrIRKA_imagInit)*[u0(tr);u1(tr)]);

% Solver.
[trIRKA_imagInit, vrIRKA_imagInit] = ode15s(fAdvDiffredIRKA_imagInit, tsteps, vr0, options);

% Reduced output.
yrIRKA_imagInit = ClinrIRKA_imagInit.'*vrIRKA_imagInit.';
for tt = 1:length(trIRKA_imagInit)
    yrIRKA_imagInit(:,tt) = yrIRKA_imagInit(:,tt) + vrIRKA_imagInit(tt, :)*MquadrIRKA_imagInit*vrIRKA_imagInit(tt, :)';
end
yrIRKA_imagInit = yrIRKA_imagInit + (1/2)*ones(1, length(trIRKA_imagInit));

%%
% 2. LQO-BT.
fprintf(1, 'Simulate BT reduced output\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
options       = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', eye(r, r), 'Jacobian', ArBT, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffredBT = @(tr,vrBT)(ArBT*vrBT + BrBT*[u0(tr);u1(tr)]);

% Solver.
[trBT, vrBT] = ode15s(fAdvDiffredBT, tsteps, vr0, options);

% Reduced output.
yrBT = ClinrBT.'*vrBT.';
for tt = 1:length(trBT)
    yrBT(:,tt) = yrBT(:,tt) + vrBT(tt, :)*MquadrBT*vrBT(tt, :)';
end
yrBT = yrBT + (1/2)*ones(1, length(trBT));

%%
% One-step interpolation.
fprintf(1, 'Simulate tangential interpolation (Petrov-Galerkin, mixed) reduced outputs.\n');
fprintf(1, '-----------------------------------------------------------------------------\n');
%%
% 3a. One-step interpolation (eigs init).
options                            = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                                      'Mass', eye(r, r), 'Jacobian', ErOneStepInterp_eigsInit\ArOneStepInterp_eigsInit, ...
                                      'MStateDependence', 'none', 'Stats','on');
fAdvDiffredOneStepInterp_eigsInit = @(tr, vrOneStepInterp_eigsInit)((ErOneStepInterp_eigsInit\ArOneStepInterp_eigsInit)*vrOneStepInterp_eigsInit ...
                                    + (ErOneStepInterp_eigsInit\BrOneStepInterp_eigsInit)*[u0(tr);u1(tr)]);

% Solver.
[trOneStepInterp_eigsInit, vrOneStepInterp_eigsInit] = ode15s(fAdvDiffredOneStepInterp_eigsInit, tsteps, vr0, options);

% Reduced output.
yrOneStepInterp_eigsInit = ClinrOneStepInterp_eigsInit.'*vrOneStepInterp_eigsInit.';
for tt = 1:length(trOneStepInterp_eigsInit)
    yrOneStepInterp_eigsInit(:,tt) = yrOneStepInterp_eigsInit(:,tt) + vrOneStepInterp_eigsInit(tt, :)*MquadrOneStepInterp_eigsInit*vrOneStepInterp_eigsInit(tt, :)';
end
yrOneStepInterp_eigsInit = yrOneStepInterp_eigsInit + (1/2)*ones(1, length(trOneStepInterp_eigsInit));

%%
% 3b. One-step interpolation (imag init).
options                           = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                                    'Mass', eye(r, r), 'Jacobian', ErOneStepInterp_imagInit\ArOneStepInterp_imagInit, ...
                                    'MStateDependence', 'none', 'Stats','on');
fAdvDiffredOneStepInterp_imagInit = @(tr, vrOneStepInterp_imagInit)((ErOneStepInterp_imagInit\ArOneStepInterp_imagInit)*vrOneStepInterp_imagInit ...
                                    + (ErOneStepInterp_imagInit\BrOneStepInterp_imagInit)*[u0(tr);u1(tr)]);

% Solver.
[trOneStepInterp_imagInit, vrOneStepInterp_imagInit] = ode15s(fAdvDiffredOneStepInterp_imagInit, tsteps, vr0, options);

% Reduced output.
yrOneStepInterp_imagInit = ClinrOneStepInterp_imagInit.'*vrOneStepInterp_imagInit.';
for tt = 1:length(trOneStepInterp_imagInit)
    yrOneStepInterp_imagInit(:,tt) = yrOneStepInterp_imagInit(:,tt) + vrOneStepInterp_imagInit(tt, :)*MquadrOneStepInterp_imagInit*vrOneStepInterp_imagInit(tt, :)';
end
yrOneStepInterp_imagInit = yrOneStepInterp_imagInit + (1/2)*ones(1, length(trOneStepInterp_imagInit));

fprintf(1, '1. REDUCED-ORDER SIMULATIONS (sinc input) FINISHED IN %d s\n', toc(reducedOrderSolvesStart));
fprintf(1, '-----------------------------------------------------------------------------\n');

%% Plotting.
% Since input.
plotBool = true;
if plotBool
    plot(trIRKA_eigsInit,           yrIRKA_eigsInit,          '--', 'color', ColMat(2,:), LineWidth=1.5); hold on
    plot(trIRKA_imagInit,           yrIRKA_imagInit,          '-.', 'color', ColMat(2,:), LineWidth=1.5);
    plot(trBT,                      yrBT,                     '--', 'color', ColMat(3,:), LineWidth=1.5); 
    plot(trOneStepInterp_eigsInit,  yrOneStepInterp_eigsInit, '--', 'color', ColMat(5,:), LineWidth=1.5); 
    plot(trOneStepInterp_imagInit,  yrOneStepInterp_imagInit, '-.', 'color', ColMat(5,:), LineWidth=1.5); 
    lgd = legend('$y(t)$', '$y_{IRKA,eigs}(t)$', '$y_{IRKA,imag}(t)$', ...
            '$y_{BT}(t)$','$y_{OneStepInterp,eigs}(t)$', ...
            '$y_{OneStepInterp,imag}(t)$', 'interpreter','latex', 'FontName', 'Arial', ...
            'location', 'northeast');
    fontsize(lgd, 10, 'points')
    
    subplot(2,1,2)  
    semilogy(trIRKA_eigsInit,           abs(y - yrIRKA_eigsInit)./abs(y),           '--', 'color', ColMat(2,:), LineWidth=1.5); hold on;
    semilogy(trIRKA_imagInit,           abs(y - yrIRKA_imagInit)./abs(y),           '-.', 'color', ColMat(2,:), LineWidth=1.5);
    semilogy(trBT,                      abs(y - yrBT)./abs(y),                      '-.', 'color', ColMat(3,:), LineWidth=1.5); 
    semilogy(trOneStepInterp_eigsInit,  abs(y - yrOneStepInterp_eigsInit)./abs(y),  '--', 'color', ColMat(4,:), LineWidth=1.5); 
    semilogy(trOneStepInterp_imagInit,  abs(y - yrOneStepInterp_imagInit)./abs(y),  '-.', 'color', ColMat(4,:), LineWidth=1.5); 
    xlabel('$t$','interpreter','latex'); 
    ylabel('$|y(t) - y_r(t)|/|y(t)|$ ', 'fontsize', fs, 'interpreter', 'latex', ...
        LineWidth=1.5)
    
    % Write data
    write = true;
    if write
        outputs = [t, y', yrIRKA_eigsInit', yrIRKA_imagInit', yrBT', ...
             yrOneStepInterp_eigsInit', yrOneStepInterp_imagInit'];
        dlmwrite('results/AdvecDiff3000_exponential_r30_Outputs.dat', outputs, 'delimiter', ...
            '\t', 'precision', 8);
        outputerrors = [t, (abs(y - yrIRKA_eigsInit)./abs(y))', (abs(y - yrIRKA_imagInit)./abs(y))', ...
            (abs(y - yrBT)./abs(y))', (abs(y - yrOneStepInterp_eigsInit)./abs(y))', ...
            (abs(y - yrOneStepInterp_imagInit)./abs(y))'];
        dlmwrite('results/AdvecDiff3000_exponential_r30_OutputErrors.dat', outputerrors, ...
            'delimiter', '\t', 'precision', 8);
    end
end

% Print errors.
fprintf(1, 'Order r = %d.\n', r)
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, 'exponential input: Relative L-infty error due to LQO-IRKA (eigs init)      : %.16f \n', max(abs(y - yrIRKA_eigsInit)./abs(y)))
fprintf(1, 'exponential input: Relative L-infty error due to LQO-IRKA (imag init)      : %.16f \n', max(abs(y - yrIRKA_imagInit)./abs(y)))
fprintf(1, 'exponential input: Relative L-infty error due to LQO-BT                    : %.16f \n', max(abs(y - yrBT)./abs(y)))
fprintf(1, 'exponential input: Relative L-infty error due to interp (mixed, eigs init) : %.16f \n', max(abs(y - yrOneStepInterp_eigsInit)./abs(y)))
fprintf(1, 'exponential input: Relative L-infty error due to interp (mixed, imag init) : %.16f \n', max(abs(y - yrOneStepInterp_imagInit)./abs(y)))
fprintf(1, '-----------------------------------------------------------------------------\n');
fprintf(1, 'exponential input: Relative L-2 error due to LQO-IRKA (eigs init)          : %.16f \n', sqrt(sum(abs(y - yrIRKA_eigsInit).^2)/sum(abs(y).^2)))
fprintf(1, 'exponential input: Relative L-2 error due to LQO-IRKA (imag init)          : %.16f \n', sqrt(sum(abs(y - yrIRKA_imagInit).^2)/sum(abs(y).^2)))
fprintf(1, 'exponential input: Relative L-2 error due to LQO-BT                        : %.16f \n', sqrt(sum(abs(y - yrBT).^2)/sum(abs(y))))
fprintf(1, 'exponential input: Relative L-2 error due to interp (mixed, eigs init)     : %.16f \n', sqrt(sum(abs(y - yrOneStepInterp_eigsInit).^2)/sum(abs(y).^2)))
fprintf(1, 'exponential input: Relative L-2 error due to interp (mixed, imag init)     : %.16f \n', sqrt(sum(abs(y - yrOneStepInterp_imagInit).^2)/sum(abs(y).^2)))
fprintf(1, '-----------------------------------------------------------------------------\n');

fprintf(1, '\n');

%% Errors.
% Relative H2 errors.
H2errorIRKA_eigsInit           = compute_lqoH2_error(A, B, Clin', Mquad, ErIRKA_eigsInit\ArIRKA_eigsInit, ...
    ErIRKA_eigsInit\BrIRKA_eigsInit, ClinrIRKA_eigsInit, MquadrIRKA_eigsInit, fomH2norm);
H2errorIRKA_imagInit           = compute_lqoH2_error(A, B, Clin', Mquad, ErIRKA_imagInit\ArIRKA_imagInit, ...
    ErIRKA_imagInit\BrIRKA_imagInit, ClinrIRKA_imagInit, MquadrIRKA_imagInit, fomH2norm);
H2errorBT                      = compute_lqoH2_error(A, B, Clin', Mquad, ArBT, BrBT, ClinrBT, MquadrBT, fomH2norm);
H2errorOneSidedInterp_eigsInit = compute_lqoH2_error(A, B, Clin', Mquad, ErOneStepInterp_eigsInit\ArOneStepInterp_eigsInit, ...
    ErOneStepInterp_eigsInit\BrOneStepInterp_eigsInit, ClinrOneStepInterp_eigsInit, MquadrOneStepInterp_eigsInit, fomH2norm);
H2errorOneSidedInterp_imagInit = compute_lqoH2_error(A, B, Clin', Mquad, ErOneStepInterp_imagInit\ArOneStepInterp_imagInit, ...
    ErOneStepInterp_imagInit\BrOneStepInterp_imagInit, ClinrOneStepInterp_imagInit, MquadrOneStepInterp_imagInit, fomH2norm);

fprintf(1, '-----------------------------------------------------------------------------\n');
fprintf(1, 'Relative H-2 error due to LQO-IRKA (eigs init)      : %.16f \n', H2errorIRKA_imagInit)
fprintf(1, 'Relative H-2 error due to LQO-IRKA (imag init)      : %.16f \n', H2errorIRKA_eigsInit)
fprintf(1, 'Relative H-2 error due to LQO-BT                    : %.16f \n', H2errorBT)
fprintf(1, 'Relative H-2 error due to interp (mixed, eigs init) : %.16f \n', H2errorOneSidedInterp_eigsInit)
fprintf(1, 'Relative H-2 error due to interp (mixed, imag init) : %.16f \n', H2errorOneSidedInterp_imagInit)
fprintf(1, '-----------------------------------------------------------------------------\n');

%% Hierarchy of reduced models.
fprintf(1, 'Computing hiearchy of reduced models.\n')
fprintf(1, '--------------------------------------------------------------\n')

rmax                  = 30;
relErrorIRKA_eigsInit = zeros(rmax/2, 1);
relErrorIRKA_imagInit = zeros(rmax/2, 1);
relErrorBT            = zeros(rmax/2, 1);

% Gramians computed above.
% Precompute BT bases. 
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
    fprintf(1, 'CURRENT ORDER IN REDUCED MODEL HIERARCHY: r = %d\n', r)
    fprintf(1, '----------------------------------------\n')
  
    %%
    % 1a. LQO-IRKA (eigs initialization).

    % Input opts; default interpolation data from projected eigenproblem
    % using dominant eigenpairs computed via MATLAB's 'eigs'.
    eigsOpts                = struct();
    eigsOpts.MaxIterations  = 100;
    eigsOpts.Tolerance      = 1e-10;
    [V, D, flag]            = eigs(A, r, 'smallestabs', eigsOpts); % E is identity
    [V, ~]                  = qr(V, 'econ');

    % Get residues from projected eigenproblem.
    Ar     = V.'*A*V; 
    Br     = V.'*B;   
    Clinr  = Clin*V;
    Mquadr = V.'*Mquad*V;
    [XrIRKA_eigsInit, Lr] = eig(Ar); poles_eigsInit = diag(Lr);
    
    % Residues.
    XrInverse               = XrIRKA_eigsInit\eye(r, r); 
    rightTangents_eigsInit  = XrInverse*Br;
    rightTangents_eigsInit  = rightTangents_eigsInit.';
    loLeftTangents_eigsInit = XrIRKA_eigsInit.'*Clinr.';  
    qoLeftTangents_eigsInit = XrIRKA_eigsInit.'*Mquadr*XrIRKA_eigsInit; 

    % Set input opts.
    optsIRKA_eigsInit                = struct();
    optsIRKA_eigsInit.poles          = poles_eigsInit;
    optsIRKA_eigsInit.rightTangents  = rightTangents_eigsInit;
    optsIRKA_eigsInit.loLeftTangents = loLeftTangents_eigsInit;
    optsIRKA_eigsInit.qoLeftTangents = qoLeftTangents_eigsInit;
    optsIRKA_eigsInit.tol            = 10e-10;
    optsIRKA_eigsInit.maxIter        = 200; 
    optsIRKA_eigsInit.checkInterp    = false;
    optsIRKA_eigsInit.plotConv       = false;
    
    [ErIRKA, ArIRKA, BrIRKA, ClinrIRKA, MquadrIRKA, infoIRKA_eigsInit] = ...
        miso_lqoirka(E, A, B, Clin.', Mquad, r, optsIRKA_eigsInit);

    H2errorIRKA = compute_lqoH2_error(A, B, Clin', Mquad, (ErIRKA\ArIRKA), (ErIRKA\BrIRKA), ClinrIRKA, MquadrIRKA, fomH2norm);
    k           = r/2;  

    % Save error.
    relErrorIRKA_eigsInit(k, :) = H2errorIRKA;

    fprintf(1, 'RELATIVE H2 ERROR OF LQO-IRKA (eigs initialization) REDUCED MODEL IS: ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relErrorIRKA_eigsInit(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

    %%
    % 1b. LQO-IRKA (imag. initialization).

    % Input opts; pick interpolation data 'naively' to be imaginary shifts,
    % orthonormal residues.

    % Compute shifts.
    tmp            = 1i*logspace(0, 3, r/2)';
    poles_imagInit = ([tmp; conj(flipud(tmp))]);
    poles_imagInit = -1e3*ones(r, 1) + poles_imagInit;

    % Set input opts. (Other inputs come from default options for this
    % initialization.)
    optsIRKA_imagInit                = struct();
    optsIRKA_imagInit.poles          = poles_imagInit;
    optsIRKA_imagInit.tol            = 10e-10;
    optsIRKA_imagInit.maxIter        = 200; 
    optsIRKA_imagInit.checkInterp    = false;

    [ErIRKA, ArIRKA, BrIRKA, ClinrIRKA, MquadrIRKA, infoIRKA_imagInit] = ...
        miso_lqoirka(E, A, B, Clin.', Mquad, r, optsIRKA_imagInit);


    H2errorIRKA = compute_lqoH2_error(A, B, Clin', Mquad, (ErIRKA\ArIRKA), (ErIRKA\BrIRKA), ClinrIRKA, MquadrIRKA, fomH2norm);
    k           = r/2;  

    % Save error.
    relErrorIRKA_imagInit(k, :) = H2errorIRKA;

    fprintf(1, 'RELATIVE H2 ERROR OF LQO-IRKA (imag. initialization) REDUCED MODEL IS: ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relErrorIRKA_imagInit(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

    %%
    % 2. LQO-BT

    % Compute projection matrices; pre-computed factors.
    V = U*Z(:, 1:r)*S(1:r, 1:r)^(-1/2); % Right
    W = L*Y(:, 1:r)*S(1:r, 1:r)^(-1/2); % Left
    
    % Compute reduced order model via projection.
    ArBT = W'*A*V;   BrBT = W'*B;  ClinrBT = Clin*V;   MquadrBT = V'*Mquad*V;
    ErBT = eye(r, r);

    H2errorBT = compute_lqoH2_error(A, B, Clin', Mquad, ArBT, BrBT, ClinrBT', MquadrBT, fomH2norm);
    k         = r/2;  

    % Save error.
    relErrorBT(k, :) = H2errorBT;
    fprintf(1, 'RELATIVE H2 ERROR OF LQO-BT REDUCED MODEL IS: ||G - Gr||_H2/||G||_H2 = %.16f\n', ...
        relErrorBT(k, :));
    fprintf(1, '------------------------------------------------------------------------------\n')

end

% Write.
write = 1;
if write
    relH2Errors = [(2:2:rmax)', relErrorIRKA_eigsInit, relErrorIRKA_imagInit, ...
        relErrorBT];
    dlmwrite('results/AdvecDiff3000_H2errors.dat', relH2Errors, 'delimiter', '\t', 'precision', ...
        8);
end

%% Plot H2 errors.
fprintf(1, 'Plotting H2 errors\n');
fprintf(1, '------------------\n');

figure
set(gca, 'fontsize', 10)
semilogy(2:2:rmax, relErrorIRKA_eigsInit(1:rmax/2), '-o', 'color', ColMat(2,:), LineWidth=1.5); hold on;
semilogy(2:2:rmax, relErrorIRKA_imagInit(1:rmax/2), '-x', 'color', ColMat(2,:), LineWidth=1.5);
semilogy(2:2:rmax, relErrorBT(1:rmax/2),            '-o', 'color', ColMat(4,:), LineWidth=1.5)
xlim([2,rmax])
lgd = legend('LQO-IRKA (eigs init)', 'LQO-IRKA (imag. init)', 'LQO-BT', 'interpreter', ...
    'latex', 'FontName', 'Arial', 'location', 'northeast');
xlabel('$r$', 'interpreter', 'latex')
ylabel('$||\mathcal{G}-\mathcal{G}_r||_{\mathcal{H}_2}/||\mathcal{G}||_{\mathcal{H}_2}$',...
    'interpreter', 'latex')

print -depsc2 results/AdvecDiff3000_H2errors_Plot

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
