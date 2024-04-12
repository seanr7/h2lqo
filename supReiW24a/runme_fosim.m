%% RUNME_FOSIM
% Script file to run full-order simulation of the vibro-acoustic plate 
% (plateTVA) data set.

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

%% Full-order simulation.
% Frequencies to sample at in range of 0 - 250 hz
s = 1i*linspace(0, 2*pi*250, 500);
s_hz = imag(s)/2/pi;
% recompute = true;
recompute = false;
if recompute == true
    fprintf(1, 'Beginning full order simulation\n')
    fprintf(1, '-------------------------------\n');
    overall_start = tic;
    res = zeros(1,length(s));
    for ii=1:length(s)
        fprintf(1, 'Frequency step %d, f=%.2f Hz ...\n',ii,imag(s(ii))/2/pi)
        current_iter = tic;
        tmp = (s(ii)*E_qo - A_qo) \ B_qo;
        res(ii) = sqrt((tmp'*Q_qo*tmp)/n_nodes); 
        fprintf(1, 'Current iteration of full order simulation finished in %.2f s\n',toc(current_iter))
        fprintf(1, '----------------------------------------------------------------------\n');
    end
    fprintf(1, 'Reduced order simulations finished in %.2f s\n', toc(overall_start))
    fprintf(1, '--------------------------------------------------\n');
    mag = 10*log10(abs(res)/1e-9);
    filename = 'data/fosim_data.mat';
    save(filename,'res','s_hz','mag')
else
    fprintf('Not re-running the full-order simulation; loading saved data from file fosim_data.mat\n')
    load('results/fosim_data.mat')
end

figure('name', 'Transfer function')
golden_ratio = (sqrt(5)+1)/2;
axes('position', [.125 .15 .75 golden_ratio-1])
plot(s_hz,mag,LineWidth=1)
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
