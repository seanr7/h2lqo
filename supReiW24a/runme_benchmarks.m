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

% Rename damping matrix
D = E; 


%% Convert plate model to first-order from second-order.
% Used in computing reduced models via benchmark approaches
fprintf(1, 'Converting second-order realization to first-order linear quadratic output system\n')
tic

[n, ~] = size(M);

Efo                   = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
Efo(1:n, 1:n)         = speye(n);                      % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = M;                             % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(K) + nnz(D) + n); % A_qo = [0, I; -K, -D]
Afo(1:n, n+1:2*n)     = speye(n);                               % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -K;                                     % (2, 1) block is -D
Afo(n+1:2*n, n+1:2*n) = -D;                                     % (2, 2) block is -K

Bfo             = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
Bfo(n+1:2*n, :) = B; 

% Our quadratic output matrix is C' * C
Qfo           = spalloc(2*n, 2*n, nnz(C' * C));
Qfo(1:n, 1:n) = C' * C; 
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

%%
fprintf(1, '1. Linfty sampling and Galerkin projection\n')
fprintf(1, '------------------------------------------\n')
% Set input opts
opts              = struct();
opts.compress     = 'Linfty';
opts.proj         = 'g';
% Set to true to recompute primitive bases and tf evals, else use those 
% saved in 'results/'
opts.recomp_bases = false;
opts.recomp_evals = false;

if ~(opts.recomp_bases)
    load('results/prim_bases_g') 
    opts.Vprim = Vprim_g;
    opts.Wprim = Wprim_g;
end
if ~(opts.recomp_evals)
    load('results/H_shifts') 
    opts.H_shifts = H_shifts; 
end

r = 100;
% Compute interpolatory bases
[Wprim_g, Vprim_g, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(M, D, ...
    K, Qfo, B, shifts, r, opts);

% Save shifts and bases if first pass, else, can comment out
% save('results/H_shifts.mat', 'H_shifts');
% save('results/prim_bases_g', 'Vprim_g', 'Wprim_g');

% Now, compute reduced models for orders r = 25, 50, 75, 100
fprintf(1, '1a. Computing Linfty Galerkin reduced model, order r = 100\n')
fprintf(1, '----------------------------------------------------------\n')
% Compute corresponding reduced model
Efo_r100_Linfty_g = Worth'*Efo*Vorth; Afo_r100_Linfty_g = Worth'*Afo*Vorth; 
Qfo_r100_Linfty_g = Vorth'*Qfo*Vorth; Bfo_r100_Linfty_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r100_Linfty_g.mat';
save(filename, 'Efo_r100_Linfty_g', 'Afo_r100_Linfty_g', 'Bfo_r100_Linfty_g', ...
    'Qfo_r100_Linfty_g', 'pW', 'pV') 

% Redo input opts for r = 25, 50, 75 models
opts              = struct();
opts.compress     = 'Linfty';
opts.proj         = 'g';
% Bases should now be saved, can reuse
opts.recomp_bases = false;
opts.recomp_evals = false;
% Data saved from above
opts.Vprim        = Vprim_g;
opts.Wprim        = Wprim_g;
opts.H_shifts     = H_shifts;

fprintf(1, '1b. Computing Linfty Galerkin reduced model, order r = 75\n')
fprintf(1, '----------------------------------------------------------\n')
r = 75;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r75_Linfty_g = Worth'*Efo*Vorth; Afo_r75_Linfty_g = Worth'*Afo*Vorth; 
Qfo_r75_Linfty_g = Vorth'*Qfo*Vorth; Bfo_r75_Linfty_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r75_Linfty_g.mat';
save(filename, 'Efo_r75_Linfty_g', 'Afo_r75_Linfty_g', 'Bfo_r75_Linfty_g', ...
    'Qfo_r75_Linfty_g', 'pW', 'pV') 

fprintf(1, '1c. Computing Linfty Galerkin reduced model, order r = 50\n')
fprintf(1, '----------------------------------------------------------\n')
r = 50;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r50_Linfty_g = Worth'*Efo*Vorth; Afo_r50_Linfty_g = Worth'*Afo*Vorth; 
Qfo_r50_Linfty_g = Vorth'*Qfo*Vorth; Bfo_r50_Linfty_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r50_Linfty_g.mat';
save(filename, 'Efo_r50_Linfty_g', 'Afo_r50_Linfty_g', 'Bfo_r50_Linfty_g', ...
    'Qfo_r50_Linfty_g', 'pW', 'pV') 

fprintf(1, '1d. Computing Linfty Galerkin reduced model, order r = 25\n')
fprintf(1, '----------------------------------------------------------\n')
r = 25;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r25_Linfty_g = Worth'*Efo*Vorth; Afo_r25_Linfty_g = Worth'*Afo*Vorth; 
Qfo_r25_Linfty_g = Vorth'*Qfo*Vorth; Bfo_r25_Linfty_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r25_Linfty_g.mat';
save(filename, 'Efo_r25_Linfty_g', 'Afo_r25_Linfty_g', 'Bfo_r25_Linfty_g', ...
    'Qfo_r25_Linfty_g', 'pW', 'pV') 


%%
fprintf(1, '2. Linfty sampling and Petrov-Galerkin projection\n')
fprintf(1, '-------------------------------------------------\n')
% Set input opts
opts              = struct();
opts.compress     = 'Linfty';
opts.proj         = 'pg';
% Set to true to recompute primitive bases and tf evals, else use those 
% saved in 'results/'
opts.recomp_bases = true;
opts.recomp_evals = true;

if ~(opts.recomp_bases)
    load('results/prim_bases_pg') 
    opts.Vprim = Vprim_pg;
    opts.Wprim = Wprim_pg;
end
if ~(opts.recomp_evals)
    load('results/H_shifts') 
    opts.H_shifts = H_shifts;
end

r = 100;
% Compute interpolatory bases
[Wprim_pg, Vprim_pg, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(...
    M, D, K, Qfo, B, shifts, r, opts);

% Save shifts and bases if first pass, else, can comment out
save('results/H_shifts.mat', 'H_shifts');
save('results/prim_bases_pg', 'Vprim_pg', 'Wprim_pg');

% Now, compute reduced models for orders r = 25, 50, 75, 100
fprintf(1, '2a. Computing Linfty Petrov-Galerkin reduced model, order r = 100\n')
fprintf(1, '----------------------------------------------------------\n')
% Compute corresponding reduced model
Efo_r100_Linfty_pg = Worth'*Efo*Vorth; Afo_r100_Linfty_pg = Worth'*Afo*Vorth; 
Qfo_r100_Linfty_pg = Vorth'*Qfo*Vorth; Bfo_r100_Linfty_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r100_Linfty_pg.mat';
save(filename, 'Efo_r100_Linfty_pg', 'Afo_r100_Linfty_pg', 'Bfo_r100_Linfty_pg', ...
    'Qfo_r100_Linfty_pg', 'pW', 'pV') 

% Redo input opts for r = 25, 50, 75 models
opts              = struct();
opts.compress     = 'Linfty';
opts.proj         = 'pg';
% Bases should now be saved, can reuse
opts.recomp_bases = false;
opts.recomp_evals = false;
% Data saved from above
opts.Vprim        = Vprim_pg;
opts.Wprim        = Wprim_pg;
opts.H_shifts     = H_shifts;

fprintf(1, '2b. Computing Linfty Petrov-Galerkin reduced model, order r = 75\n')
fprintf(1, '----------------------------------------------------------\n')
r = 75;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r75_Linfty_pg = Worth'*Efo*Vorth; Afo_r75_Linfty_pg = Worth'*Afo*Vorth; 
Qfo_r75_Linfty_pg = Vorth'*Qfo*Vorth; Bfo_r75_Linfty_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r75_Linfty_pg.mat';
save(filename, 'Efo_r75_Linfty_pg', 'Afo_r75_Linfty_pg', 'Bfo_r75_Linfty_pg', ...
    'Qfo_r75_Linfty_pg', 'pW', 'pV') 

fprintf(1, '2c. Computing Linfty Petrov-Galerkin reduced model, order r = 50\n')
fprintf(1, '----------------------------------------------------------\n')
r = 50;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r50_Linfty_pg = Worth'*Efo*Vorth; Afo_r50_Linfty_pg = Worth'*Afo*Vorth; 
Qfo_r50_Linfty_pg = Vorth'*Qfo*Vorth; Bfo_r50_Linfty_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r50_Linfty_pg.mat';
save(filename, 'Efo_r50_Linfty_pg', 'Afo_r50_Linfty_pg', 'Bfo_r50_Linfty_pg', ...
    'Qfo_r50_Linfty_pg', 'pW', 'pV') 

fprintf(1, '2d. Computing Linfty Petrov-Galerkin reduced model, order r = 25\n')
fprintf(1, '----------------------------------------------------------\n')
r = 25;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r25_Linfty_pg = Worth'*Efo*Vorth; Afo_r25_Linfty_pg = Worth'*Afo*Vorth; 
Qfo_r25_Linfty_pg = Vorth'*Qfo*Vorth; Bfo_r25_Linfty_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r25_Linfty_pg.mat';
save(filename, 'Efo_r25_Linfty_pg', 'Afo_r25_Linfty_pg', 'Bfo_r25_Linfty_pg', ...
    'Qfo_r25_Linfty_pg', 'pW', 'pV') 



%%
fprintf(1, '3. Pivoted QR and Galerkin projection\n')
fprintf(1, '-------------------------------------\n')
% Set input opts
opts              = struct();
opts.compress     = 'avg';
opts.proj         = 'g';
% Set to true to recompute primitive bases and tf evals, else use those 
% saved in 'results/'
opts.recomp_bases = false;
opts.recomp_evals = false;

if ~(opts.recomp_bases)
    load('results/prim_bases_g') 
    opts.Vprim = Vprim_g;
    opts.Wprim = Wprim_g;
end
if ~(opts.recomp_evals)
    load('results/H_shifts') 
    opts.H_shifts = H_shifts;
end

r = 100;
% Compute interpolatory bases
[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, shifts, ...
    r, opts);

% Now, compute reduced models for orders r = 25, 50, 75, 100
fprintf(1, '3a. Computing Averaged Galerkin reduced model, order r = 100\n')
fprintf(1, '----------------------------------------------------------\n')
% Compute corresponding reduced model
Efo_r100_avg_g = Worth'*Efo*Vorth; Afo_r100_avg_g = Worth'*Afo*Vorth; 
Qfo_r100_avg_g = Vorth'*Qfo*Vorth; Bfo_r100_avg_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r100_avg_g.mat';
save(filename, 'Efo_r100_avg_g', 'Afo_r100_avg_g', 'Bfo_r100_avg_g', ...
    'Qfo_r100_avg_g', 'pW', 'pV') 

% Redo input opts for r = 25, 50, 75 models
opts              = struct();
opts.compress     = 'avg';
opts.proj         = 'g';
% Bases should now be saved, can reuse
opts.recomp_bases = false;
opts.recomp_evals = false;
% Data saved from above
opts.Vprim        = Vprim_g;
opts.Wprim        = Wprim_g;
opts.H_shifts     = H_shifts;


fprintf(1, '3b. Computing Averaged Galerkin reduced model, order r = 75\n')
fprintf(1, '----------------------------------------------------------\n')
r = 75;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r75_avg_g = Worth'*Efo*Vorth; Afo_r75_avg_g = Worth'*Afo*Vorth; 
Qfo_r75_avg_g = Vorth'*Qfo*Vorth; Bfo_r75_avg_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r75_avg_g.mat';
save(filename, 'Efo_r75_avg_g', 'Afo_r75_avg_g', 'Bfo_r75_avg_g', ...
    'Qfo_r75_avg_g', 'pW', 'pV') 

fprintf(1, '3c. Computing Averaged Galerkin reduced model, order r = 50\n')
fprintf(1, '----------------------------------------------------------\n')
r = 50;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r50_avg_g = Worth'*Efo*Vorth; Afo_r50_avg_g = Worth'*Afo*Vorth; 
Qfo_r50_avg_g = Vorth'*Qfo*Vorth; Bfo_r50_avg_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r50_avg_g.mat';
save(filename, 'Efo_r50_avg_g', 'Afo_r50_avg_g', 'Bfo_r50_avg_g', ...
    'Qfo_r50_avg_g', 'pW', 'pV') 

fprintf(1, '3d. Computing Averaged Galerkin reduced model, order r = 25\n')
fprintf(1, '----------------------------------------------------------\n')
r = 25;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r25_avg_g = Worth'*Efo*Vorth; Afo_r25_avg_g = Worth'*Afo*Vorth; 
Qfo_r25_avg_g = Vorth'*Qfo*Vorth; Bfo_r25_avg_g = Worth'*Bfo;

filename = 'results/plateTVAlqo_r25_avg_g.mat';
save(filename, 'Efo_r25_avg_g', 'Afo_r25_avg_g', 'Bfo_r25_avg_g', ...
    'Qfo_r25_avg_g', 'pW', 'pV') 


%%
fprintf(1, '4. Pivoted QR and Petrov-Galerkin projection\n')
fprintf(1, '--------------------------------------------\n')
% Set input opts
opts              = struct();
opts.compress     = 'avg';
opts.proj         = 'pg';
% Set to true to recompute primitive bases and tf evals, else use those 
% saved in 'results/'
opts.recomp_bases = false;
opts.recomp_evals = false;

if ~(opts.recomp_bases)
    load('results/prim_bases_pg') 
    opts.Vprim = Vprim_pg;
    opts.Wprim = Wprim_pg;
end
if ~(opts.recomp_evals)
    load('results/H_shifts') 
    opts.H_shifts = H_shifts;
end

r = 100;
% Compute interpolatory bases
[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, shifts, ...
    r, opts);

% Now, compute reduced models for orders r = 25, 50, 75, 100
fprintf(1, '4a. Computing Averaged Petrov-Galerkin reduced model, order r = 100\n')
fprintf(1, '----------------------------------------------------------\n')
% Compute corresponding reduced model
Efo_r100_avg_pg = Worth'*Efo*Vorth; Afo_r100_avg_pg = Worth'*Afo*Vorth; 
Qfo_r100_avg_pg = Vorth'*Qfo*Vorth; Bfo_r100_avg_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r100_avg_pg.mat';
save(filename, 'Efo_r100_avg_pg', 'Afo_r100_avg_pg', 'Bfo_r100_avg_pg', ...
    'Qfo_r100_avg_pg', 'pW', 'pV') 

% Redo input opts for r = 25, 50, 75 models
opts              = struct();
opts.compress     = 'avg';
opts.proj         = 'pg';
% Bases should now be saved, can reuse
opts.recomp_bases = false;
opts.recomp_evals = false;
% Data saved from above
opts.Vprim        = Vprim_pg;
opts.Wprim        = Wprim_pg;
opts.H_shifts     = H_shifts;


fprintf(1, '4b. Computing Averaged Petrov-Galerkin reduced model, order r = 75\n')
fprintf(1, '----------------------------------------------------------\n')
r = 75;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r75_avg_pg = Worth'*Efo*Vorth; Afo_r75_avg_pg = Worth'*Afo*Vorth; 
Qfo_r75_avg_pg = Vorth'*Qfo*Vorth; Bfo_r75_avg_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r75_avg_pg.mat';
save(filename, 'Efo_r75_avg_pg', 'Afo_r75_avg_pg', 'Bfo_r75_avg_pg', ...
    'Qfo_r75_avg_pg', 'pW', 'pV') 

fprintf(1, '4c. Computing Averaged Petrov-Galerkin reduced model, order r = 50\n')
fprintf(1, '----------------------------------------------------------\n')
r = 50;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r50_avg_pg = Worth'*Efo*Vorth; Afo_r50_avg_pg = Worth'*Afo*Vorth; 
Qfo_r50_avg_pg = Vorth'*Qfo*Vorth; Bfo_r50_avg_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r50_avg_pg.mat';
save(filename, 'Efo_r50_avg_pg', 'Afo_r50_avg_pg', 'Bfo_r50_avg_pg', ...
    'Qfo_r50_avg_pg', 'pW', 'pV') 

fprintf(1, '4d. Computing Averaged Petrov-Galerkin reduced model, order r = 25\n')
fprintf(1, '----------------------------------------------------------\n')
r = 25;

[~, ~, Worth, Vorth, ~, pW, pV] = interpolatory_solves(M, D, K, Qfo, B, ...
    shifts, r, opts);

% Compute corresponding reduced model
Efo_r25_avg_pg = Worth'*Efo*Vorth; Afo_r25_avg_pg = Worth'*Afo*Vorth; 
Qfo_r25_avg_pg = Vorth'*Qfo*Vorth; Bfo_r25_avg_pg = Worth'*Bfo;

filename = 'results/plateTVAlqo_r25_avg_pg.mat';
save(filename, 'Efo_r25_avg_pg', 'Afo_r25_avg_pg', 'Bfo_r25_avg_pg', ...
    'Qfo_r25_avg_pg', 'pW', 'pV') 


%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off