%% RUNME
% Script file to reduced models of the vibro-acoustic plate (plateTVA) data 
% set using the various benchmark interpolatory approaches.

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
fullpath = matlab.desktop.editor.getActiveFilename;
[rootpath, filename, ~] = fileparts(fullpath(3:end));
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
load('data/plateTVA_n201900m1q28278_fo')
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

%% 
% Simulate full-order model 
% 250 frequencies to sample at from 0hz - 250hz (s) are given in '*.mat' file 
% Instead, use 500 frequences from 0hz - 500hz

s = 1i*linspace(0,2*pi*500, 501);% Comment out to keep original freq range

% recompute = true;
recompute = false;
if recompute == true
    fprintf('Beginning full-order simulation. Estimated time of completion is %.2f s\n', 250*15.72)
    overall_start = tic;
    res = zeros(1,length(s));
    for ii=1:length(s)
        fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
        current_iter = tic;
        tmp = (s(ii) * E_qo - A_qo) \ B_qo;
        % res(ii) is H2(s(ii), s(ii)) = sqrt(((s(ii) * E_qo - A_qo)') \ Q_qo ((s(ii) * E_qo - A_qo) \ eye(n, n)))
        % i.e., sqrt() of quadratic-output transfer function
        res(ii) = sqrt((tmp' * Q_qo * tmp) / n_nodes); % Q: So really, we want sqrt(H_2(s(ii), s(ii))/n_nodes ? (Just to put it in my language..)
        fprintf('Iteration of full-order simulation finished in %.2f s\n',toc(current_iter))
    end
    fprintf('Full-order simulation finished; total time of completion is %.2f s\n', toc(overall_start))
    f = imag(s)/2/pi;
    mag = 10*log10(abs(res)/1e-9);
    frpintf('Saving full-order simulation data\n')
    filename = 'FOsim_data.mat';
    save(filename,'res','f','mag')
    movefile FOsim_data.mat data/FOsim_data.mat

else
    fprintf('Not re-running the full-order simulation; loading saved data from file FOSIM_data.mat\n')
    % load('FOSIM_data.mat')
end

% figure('name','Transfer function')
% plot(f,mag)
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')

%% 
% Now, compute LQO-ROM via 
%   - LQO-IRKA
%   - More to come ...


r = 50; % Order
% Initialization parameters
poles_prev = -logspace(1, 3, r)'; % Spread 
tmp = 10 *  rand(r, r);
SO_res_prev = (tmp+tmp')/2; 
itermax = 15;
tol = 10e-4;

fprintf('Beginning construction of order %d the LQO-ROM via LQO-IRKA\n', r)

[E_qo_r, A_qo_r, B_qo_r, ~, Q_qo_r, poles, ~, SO_res, pole_history] = lqo_irka(E_qo, A_qo, B_qo, ...
    [], Q_qo, poles_prev, [], SO_res_prev, itermax, tol, 1);

filename = 'plateTVA_r50_redux_lqoirka.mat';
save(filename, 'E_qo_r', 'A_qo_r', 'B_qo_r', 'Q_qo_r', 'pole_history', 'SO_res') 
movefile plateTVA_r50_redux_lqoirka.mat ~/data/plateTVA_r50_redux_lqoirka.mat

%%
% Now, do reduced-order simulations

fprintf('Beginning reduced-order simulation\n')
overall_start = tic;

res_ro = zeros(1,length(s));
for ii=1:length(s)
    fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    tmp = (s(ii) * E_qo_r - A_qo_r) \ B_qo_r;
    res_ro(ii) = sqrt((tmp' * Q_qo_r * tmp) / n_nodes);
    fprintf('Iteration of reduced-order simulation finished in %.2f s\n',toc(current_iter))
end
fprintf('Reduced-order simulation finished; time of completion is %.2f s/n', toc(overall_start))

f_ro = imag(s)/2/pi;
mag_ro = 10*log10(abs(res_ro)/1e-9);

fprintf('Saving reduced-order simulation data')
filename = 'ROsim_plateTVA_r50_redux_lqoirka.mat';
movefile ROsim_plateTVA_r50_redux_lqoirka.mat data/ROsim_plateTVA_r50_redux_lqoirka.mat
save(filename,'res_ro','f_ro','mag_ro')
