%%
% Script to test second order structured solves
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

%% Solve linear system via first order realization, second order realization

% First order solve; no structure taken into account
tic
s1      = -1 + 1i;
fosolve = (s1 * E_qo - A_qo)\B_qo;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc)
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
tic
tmp     = (s1*M + E)\B; 
sosolve = [(1/s1).*(tmp - ((s1^2).*M + s1.*E + K)\(K*tmp)); ...
    s1.*((s1^2).*M + s1.*E + K)\B];
fprintf(1, 'Single second-order structured solve finished in %.2f s\n',toc)
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, 'Relative solution error %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, 'Testing now for frequency closer to zero')
fprintf(1, '---------------------------------------------------------------\n');

% First order solve; no structure taken into account
tic
s1      = 10e-6*(-1 + 1i);
fosolve = (s1 * E_qo - A_qo)\B_qo;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc)
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
tic
tmp     = (s1*M + E)\B; 
sosolve = [(1/s1).*(tmp - ((s1^2).*M + s1.*E + K)\(K*tmp)); ...
    s1.*((s1^2).*M + s1.*E + K)\B];
fprintf(1, 'Single second-order structured solve finished in %.2f s\n',toc)
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, 'Relative solution error %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');