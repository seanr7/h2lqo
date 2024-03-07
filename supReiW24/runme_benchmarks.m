%% RUNME_BENCHMARKS
% Script file to reduced models of the vibro-acoustic plate (plateTVA) data 
% set using the various benchmark interpolatory approaches.

%
% This file is part of the archive Code and Results for Numerical 
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
load('data/plateTVA_n201900m1q28278.mat')
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
% Interpolatory reduced models are computed using the following
% approaches
%   1. Linfty + Galerkin projection
%   2. Linfty + Petrov-Galerkin projection
%   3. avg + Galerkin projection
%   4. avg + Petrov-Galerkin projection
% All methods are implemented in the function
% 'drivers/interpolatory_solves.m'

% Compute linear solves at the designated shifts
shifts = 1i*linspace(1,2*pi*251, 250);
% Set order of reduction
r = 50;
fprintf(1, '1. Linfty sampling and Galerkin projection\n')
fprintf(1, '------------------------------------------\n')
% Set input opts
opts.compress = 'Linfty';
opts.proj = 'g';
% Set to true to re-compute primitive bases, else use those saved in
% 'results/'
recomp_bases = false;
if recomp_bases
    Vprim = []; Wprim = [];
else
    load('results/prim_bases_g') 
end
recomp_evals = true;
if recomp_evals
    H_shifts = [];
else
    load('results/H_shifts') 
end
% Set input opts
opts.recomp_bases = recomp_bases;
opts.recomp_evals = recomp_evals;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
opts.H_shifts = H_shifts;

% Compute model reduction bases
[~, ~, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(E_qo, A_qo, ...
    B_qo, Q_qo, shifts, r, opts);
% Compute corresponding reduced model
E_qo_r50_Linfty_g = Worth'*E_qo*Vorth; A_qo_r50_Linfty_g = Worth'*A_qo*Vorth; 
Q_qo_r50_Linfty_g = Vorth'*Q_qo*Vorth; B_qo_r50_Linfty_g = Worth'*B_qo;

% Save shifts 
% save('H_shifts.mat', 'H_shifts');

filename = 'results/plateTVAlqo_r50_Linfty_g_redux.mat';
save(filename, 'E_qo_r50_Linfty_g', 'A_qo_r50_Linfty_g', 'B_qo_r50_Linfty_g', ...
    'Q_qo_r50_Linfty_g', 'pW', 'pV') 

%%
fprintf(1, '2. Linfty sampling and Petrov-Galerkin projection\n')
fprintf(1, '-------------------------------------------------\n')
% Set input opts
opts.compress = 'Linfty';
opts.proj = 'pg';
% Set to true to re-compute primitive bases, else use those saved in
% 'results/'
recomp_bases = false;
if recomp_bases
    Vprim = []; Wprim = [];
else
    load('results/prim_bases_pg') 
end
recomp_evals = true;
if recomp_evals
    H_shifts = [];
else
    load('results/H_shifts') 
end
% Set input opts
opts.recomp_bases = recomp_bases;
opts.recomp_evals = recomp_evals;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
opts.H_shifts = H_shifts;

% Compute model reduction bases
[~, ~, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(E_qo, A_qo, ...
    B_qo, Q_qo, shifts, r, opts);
% Compute corresponding reduced model
E_qo_r50_Linfty_pg = Worth'*E_qo*Vorth; A_qo_r50_Linfty_pg = Worth'*A_qo*Vorth; 
Q_qo_r50_Linfty_pg = Vorth'*Q_qo*Vorth; B_qo_r50_Linfty_pg = Worth'*B_qo;

filename = 'results/plateTVAlqo_r50_Linfty_pg_redux.mat';
save(filename, 'E_qo_r50_Linfty_pg', 'A_qo_r50_Linfty_pg', 'B_qo_r50_Linfty_pg', ...
    'Q_qo_r50_Linfty_pg', 'pW', 'pV') 

%%
fprintf(1, '3. Pivoted QR and Galerkin projection\n')
fprintf(1, '-------------------------------------\n')
% Set input opts
opts.compress = 'avg';
opts.proj = 'g';
% Set to true to re-compute primitive bases, else use those saved in
% 'results/'
% Note: At this point in the script, primitive bases for Galerkin and
% Petrov-Galerkin projection have been computed and saved, and can be
% recycled for faster computation.
recomp_bases = false;
if recomp_bases
    Vprim = []; Wprim = [];
else
    load('results/prim_bases_g') 
end
recomp_evals = true;
if recomp_evals
    H_shifts = [];
else
    load('results/H_shifts') 
end
% Set input opts
opts.recomp_bases = recomp_bases;
opts.recomp_evals = recomp_evals;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
opts.H_shifts = H_shifts;

% Compute model reduction bases
[~, ~, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(E_qo, A_qo, ...
    B_qo, Q_qo, shifts, r, opts);
% Compute corresponding reduced model
E_qo_r50_avg_g = Worth'*E_qo*Vorth; A_qo_r50_avg_g = Worth'*A_qo*Vorth; 
Q_qo_r50_avg_g = Vorth'*Q_qo*Vorth; B_qo_r50_avg_g = Worth'*B_qo;

filename = 'results/plateTVAlqo_r50_avg_g.mat';
save(filename, 'E_qo_r50_avg_g', 'A_qo_r50_avg_g', 'B_qo_r50_avg_g', ...
    'Q_qo_r50_avg_g', 'pW', 'pV') 

%%
fprintf(1, '4. Pivoted QR and Petrov-Galerkin projection\n')
fprintf(1, '--------------------------------------------\n')
% Set input opts
opts.compress = 'avg';
opts.proj = 'pg';
% Set to true to re-compute primitive bases, else use those saved in
% 'results/'
recomp_bases = false;
if recomp_bases
    Vprim = []; Wprim = [];
else
    load('results/prim_bases_pg') 
end
recomp_evals = true;
if recomp_evals
    H_shifts = [];
else
    load('results/H_shifts') 
end
% Set input opts
opts.recomp_bases = recomp_bases;
opts.recomp_evals = recomp_evals;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
opts.H_shifts = H_shifts;

% Compute model reduction bases
[~, ~, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(E_qo, A_qo, ...
    B_qo, Q_qo, shifts, r, opts);
% Compute corresponding reduced model
E_qo_r50_avg_pg = Worth'*E_qo*Vorth; A_qo_r50_avg_pg = Worth'*A_qo*Vorth; 
Q_qo_r50_avg_pg = Vorth'*Q_qo*Vorth; B_qo_r50_avg_pg = Worth'*B_qo;

filename = 'results/plateTVAlqo_r50_avg_pg.mat';
save(filename, 'E_qo_r50_avg_pg', 'A_qo_r50_avg_pg', 'B_qo_r50_avg_pg', ...
    'Q_qo_r50_avg_pg', 'pW', 'pV')

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off