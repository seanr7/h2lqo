%% RUNME_LQOIRKA_MSD
% Test script to compute reduced models of a simple mass spring damper
% system using LQOIRKA.

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model reduction of dynamical systems with 
% root mean squared error"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

clc;
clear all;
close all;

% Get and set all paths.
[rootpath, filename, ~] = fileparts(mfilename('fullpath'));
loadname            = [rootpath filesep() ...
    'data' filesep() filename];
savename            = [rootpath filesep() ...
    'results' filesep() filename];

% Paths to drivers and data.
addpath([rootpath, '/drivers'])
addpath([rootpath, '/data'])

% Include intended order to not overwrite other log files.
if exist([savename 'r25' '.log'], 'file') == 2
    delete([savename 'r25' '.log']);
end
outname = [savename 'r25' '.log']';

diary(outname)
diary on; 

fprintf(1, ['SCRIPT: ' upper(filename) '\n']);
fprintf(1, ['========' repmat('=', 1, length(filename)) '\n']);
fprintf(1, '\n');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASE DATA.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'Creating toy model...\n')
fprintf(1, '-------------------------\n');

% MMESS MATLAB toolbox function.
n1              = 10; 
alpha           = .002; 
beta            = alpha; 
v               = 5;
[Mso, Dso, Kso] = triplechain_MSD(n1, alpha, beta, v);

[n, ~] = size(Mso);
Cso    = ones(1, n);
Bso    = ones(n, 1);

% First-order realization for interpolation checks.
Efo                   = spalloc(2*n, 2*n, nnz(Mso) + n); % Descriptor matrix; Efo = [I, 0: 0, Mso]
Efo(1:n, 1:n)         = speye(n);                        % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = Mso;                             % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(Kso) + nnz(Dso) + n); % Afo = [0, I; -Kso, -Dso]
Afo(1:n, n+1:2*n)     = speye(n);                                   % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -Kso;                                       % (2, 1) block is -Kso
Afo(n+1:2*n, n+1:2*n) = -Dso;                                       % (2, 2) block is -Dso 

Bfo             = spalloc(2*n, 1, nnz(Bso)); % bfo = [0; bso];
Bfo(n+1:2*n, :) = Bso;       

% Quadratic-output matrix. 
Qfo           = spalloc(2*n, 2*n, nnz(Cso' * Cso));
Qfo(1:n, 1:n) = Cso'*Cso; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REDUCED MODELS.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input opts.
opts          = struct();
r             = 4;
opts.maxiter  = 50; 
opts.tol      = 10e-10;
opts.plotConv = true;
opts.poles    = -1*logspace(-1, 4, r)';
opts.qoRes    = ones(r, r);

fprintf(1, 'Computing reduced model of order r = %d via LQO-IRKA.\n', r)
fprintf(1, '-----------------------------------------------------\n')

[Efo_r, Afo_r, Bfo_r, Qfo_r, info] = siso_lqoirka(Mso, Dso, Kso, Bso, Qfo, r, opts);

poleHistory = info.poleHistory;
qoRes       = info.qoRes; 

% Final poles.
poles = poleHistory(:, end - 1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATION CHECKS.                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define transfer functions.
HQuad    = @(s1, s2) ((s2*Efo - Afo)\Bfo).'*Qfo*((s1*Efo - Afo)\Bfo);
HQuadRed = @(s1, s2) ((s2*Efo_r - Afo_r)\Bfo_r).'*Qfo_r*((s1*Efo_r - Afo_r)\Bfo_r);

fprintf(1, '1. Quadratic, right tangential interpolation conditions.\n')
 fprintf(1, '-----------------------------------------------------------------\n')
for j = 1:r
    for k = 1:r
        fprintf(1, 'Relative error in condition (%d, %d): %d\n', j, k, ...
            norm(HQuad(-poles(j), -poles(k)) - HQuadRed(-poles(j), -poles(k)), 2)/ ...
            norm(HQuad(-poles(j), -poles(k)), 2))
        fprintf(1, '-----------------------------------------------------------------\n')
    end
end

HQuadPartials1    = @(s1, s2) -((s2*Efo - Afo)\Bfo).'*Qfo*((s1*Efo - Afo)\(Efo*((s1*Efo - Afo)\Bfo)));
HQuadRedPartials1 = @(s1, s2) -((s2*Efo_r - Afo_r)\Bfo_r).'*Qfo_r*((s1*Efo_r - Afo_r)\(Efo_r*((s1*Efo_r - Afo_r)\Bfo_r)));
HQuadPartials2    = @(s1, s2) -((s2*Efo - Afo)\(Efo*((s2*Efo - Afo)\Bfo))).'*Qfo*((s1*Efo - Afo)\Bfo);
HQuadRedPartials2 = @(s1, s2) -((s2*Efo_r - Afo_r)\(Efo_r*((s2*Efo_r - Afo_r)\Bfo_r))).'*Qfo_r*((s1*Efo_r - Afo_r)\Bfo_r);

fprintf(1, '2. Mixed, bi-tangential Hermite interpolation conditions.\n')
fprintf(1, '-----------------------------------------------------------------\n')
for k = 1:r
    % Quadratic component of mixed term
    quadTermFO = 0;
    quadTermRO = 0;
    for j = 1:r
        quadTermFO = quadTermFO + qoRes(k, j)*HQuadPartials1(-poles(k), -poles(j)) + ...
            qoRes(j, k)*HQuadPartials2(-poles(j), -poles(k));
        quadTermRO = quadTermRO + qoRes(k, j)*HQuadRedPartials1(-poles(k), -poles(j)) + ...
            qoRes(j, k)*HQuadRedPartials2(-poles(j), -poles(k));
    end
    fprintf(1, 'Relative error in condition %d: %d\n', k, norm(quadTermFO - quadTermRO, 2)/norm(quadTermFO, 2))
    fprintf(1, '-----------------------------------------------------------------\n')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINISHED.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
