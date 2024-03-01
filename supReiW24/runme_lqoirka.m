%% RUNME
% Script file to reduced models of the vibro-acoustic plate (plateTVA) data 
% set using the various benchmark interpolatory approaches.

%
% This file is part of the archive Code, Data and Results for Numerical 
% Experiments in "Interpolatory model order reduction of large-scale 
% dynamical systems with root mean squared error measures"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
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

%% Load base data.
% To compute root mean squared output, treat plateTVA model as a linear
% quadratic output system.

fprintf(1, 'Loading plateTVA model...\n')
fprintf(1, '-------------------------\n');
load('data/plateTVA_n201900m1q28278')
n_nodes = full(sum(sum(C)));

%% Convert plate model to first-order from second-order.
fprintf(1, 'Converting second-order realization to first-order linear quadratic output system\n')
tic
[n, ~] = size(M);

E_qo = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is (sparse) mass matrix

A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(E) + n);  % A_qo = [0, I; -K, -E]
A_qo(1:n, n+1:2*n) = speye(n); % (1, 2) block of A_qo
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B;

% Our quadratic output matrix is C' * C
Q_qo = spalloc(2*n, 2*n, nnz(C' * C));
Q_qo(1:n, 1:n) = C' * C; 
fprintf(1, 'First-order realization built in %.2f s\n',toc)
fprintf(1, '--------------------------------------------\n');

%% Computed reduced-models.
% Compute reduced-models using the iterative rational Krylov algorithm for
% linear quadratic output systems.
fprintf(1, 'Computing reduced-model via IRKA\n')
fprintf(1, '--------------------------------\n')

% Set order of reduction
r = 50;
% Initialization parameters
opts.maxiter  = 15;
opts.tol      = 10e-4;
opts.plotconv = false;
% Initial selection of poles and residues
poles      = -logspace(1, 3, r)';
opts.poles = poles;
tmp        = 10 *  rand(r, r); sores = (tmp+tmp')/2; 
opts.sores = sores;
opts.fores = [];

[E_qo_r, A_qo_r, B_qo_r, ~, Q_qo_r, info] = sisolqo_irka(E_qo, A_qo, B_qo, [], ...
    Q_qo, r, opts);

poles = info.pole_hist;
SOres = info.sores; 
filename = 'results/plateTVAlqo_r50_lqoirka.mat';
save(filename, 'E_qo_r', 'A_qo_r', 'B_qo_r', 'Q_qo_r', 'poles', 'SOres') 

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off