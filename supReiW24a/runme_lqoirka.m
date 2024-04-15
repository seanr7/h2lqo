%% RUNME_LQOIRKA
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
if exist([savename 'r100' '.log'], 'file') == 2
    delete([savename 'r100' '.log']);
end
outname = [savename 'r100' '.log']';

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
% Define quadratic output matrix Qfo = C' * C
Qfo           = spalloc(2*n, 2*n, nnz(C' * C));
Qfo(1:n, 1:n) = C' * C; 

%% Computed reduced-model using lqoirka.
fprintf(1, 'Computing reduced-model via lqoirka\n')
fprintf(1, '--------------------------------\n')

% Set order of reduction and input opts
opts          = struct();
r             = 100;
opts.maxiter  = 50;
opts.tol      = 10e-4;
opts.plotconv = false;

[Efo_r, Afo_r, Bfo_r, Qfo_r, info] = lqoirka(M, D, K, B, Qfo, r, opts);

poles    = info.pole_hist;
sores    = info.sores; 
filename = 'results/plateTVAlqo_r100_lqoirka.mat';
save(filename, 'Efo_r', 'Afo_r', 'Bfo_r', 'Qfo_r', 'poles', 'sores') 

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
