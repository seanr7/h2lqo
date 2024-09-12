%% RUNME_LQOIRKA
% Script to compute reduced models of the vibro-acoustic plate (plateTVA) 
% dataset using LQOIRKA.

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
fprintf(1, 'Loading plateTVA model...\n')
fprintf(1, '-------------------------\n');
load('data/plateTVA_n201900m1q28278.mat')

% Damping matrix.
D = E; 

% Quadratic-output (QO) matrix Qfo = C'*C.
[n, ~]        = size(M);
Qfo           = spalloc(2*n, 2*n, nnz(C' * C));
Qfo(1:n, 1:n) = C'*C; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REDUCED MODELS.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set input opts.
opts          = struct();
r             = 100;
opts.maxiter  = 50; 
opts.tol      = 10e-6;
opts.plotConv = false;
opts.poles    = -1*logspace(-1, 4, r);
opts.qoRes    = ones(r, r);

fprintf(1, 'Computing reduced model of order r = %d via LQO-IRKA.\n', r)
fprintf(1, '-----------------------------------------------------\n')

[Efo_r, Afo_r, Bfo_r, Qfo_r, info] = siso_lqoirka(M, D, K, B, Qfo, r, opts);

poleHistory = info.poleHistory;
qoRes       = info.qoRes; 
filename    = 'results/plateTVAROM_r100_LQOIRKA.mat';
save(filename, 'Efo_r', 'Afo_r', 'Bfo_r', 'Qfo_r', 'poleHistory', 'qoRes') 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINISHED.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
