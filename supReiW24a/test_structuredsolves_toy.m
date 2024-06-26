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

% fprintf(1, 'Loading plateTVA model...\n')
% fprintf(1, '-------------------------\n');
% load('data/plateTVA_n201900m1q28278.mat')
% n_nodes = full(sum(sum(C)));

%%%% TOY TEST %%%%
fprintf(1, 'Loading toy model...\n')
fprintf(1, '-------------------------\n');
n1 = 10; alpha=.002; beta=alpha; v = 5;

[M, D, K]=triplechain_MSD(n1, alpha, beta, v);

M = full(M);
D = full(D);
K = full(K);
E = D;

C  = ones(1,size(K,1));
B  = ones(size(K,1),1);
B1 = B;
%%%% TOY TEST %%%%

%% Convert plate model to first-order from second-order.
fprintf(1, 'Converting second-order realization to first-order linear quadratic output system\n')
tic

[n, ~] = size(M);

E_qo                   = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n)         = speye(n);                      % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M;                             % (2, 2) block is (sparse) mass matrix

A_qo                   = spalloc(2*n, 2*n, nnz(K) + nnz(E) + n); % A_qo = [0, I; -K, -E]
A_qo(1:n, n+1:2*n)     = speye(n);                               % (1, 2) block of A_qo
A_qo(n+1:2*n, 1:n)     = -K;                                     % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E;                                     % (2, 2) block is -stiffness matrix

B_qo             = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B; 
B1               = B;           

% Our quadratic output matrix is C' * C
Q_qo           = spalloc(2*n, 2*n, nnz(C' * C));
Q_qo(1:n, 1:n) = C' * C; 
fprintf(1, 'First-order realization built in %.2f s\n',toc)
fprintf(1, '--------------------------------------------\n');

%% Test solves.
fprintf(1, '1. Purely real frequency.\n');
fprintf(1, '---------------------------------------------------------------\n');
s = -1000;

fprintf(1, 'v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B]\n');
fprintf(1, '---------------------------------------------------------------\n');
Bhat = B_qo;

tic
% First order solve; no structure taken into account
fosolve = (s * E_qo - A_qo)\Bhat;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc);
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
% rhs has form Bhat = [0; B], here
sosolve = so_structured_solve(M, E, K, B1, s, 0, 1);

fprintf(1, 'Relative solution error, v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B] %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error, v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B] %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, 'w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] \n');
fprintf(1, '---------------------------------------------------------------\n');
Bhat = Q_qo*fosolve;
tic
% First order solve; no structure taken into account
fosolve = ((s * E_qo - A_qo)')\Bhat;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc);
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
% rhs has form Bhat = [B; 0], here
B2      = Bhat(1:n, :);
sosolve = so_structured_solve(M, E, K, B2, s, 1, 1);

fprintf(1, 'Relative solution error, w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error, w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, '\n')

%% 
fprintf(1, '2. Purely imaginary frequency.\n');
fprintf(1, '---------------------------------------------------------------\n');
s = -1000*1i;

fprintf(1, 'v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B]\n');
fprintf(1, '---------------------------------------------------------------\n');
Bhat = B_qo;

tic
% First order solve; no structure taken into account
fosolve = (s * E_qo - A_qo)\Bhat;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc);
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
% rhs has form Bhat = [0; B], here
sosolve = so_structured_solve(M, E, K, B1, s, 0, 1);

fprintf(1, 'Relative solution error, v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B] %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error, v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B] %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, 'w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] \n');
fprintf(1, '---------------------------------------------------------------\n');
Bhat = Q_qo*fosolve;
tic
% First order solve; no structure taken into account
fosolve = ((s * E_qo - A_qo)')\Bhat;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc);
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
% rhs has form Bhat = [B; 0], here
B2      = Bhat(1:n, :);
sosolve = so_structured_solve(M, E, K, B2, s, 1, 1);

fprintf(1, 'Relative solution error, w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error, w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, '\n')

%%
fprintf(1, '3. Mixed frequency.\n');
fprintf(1, '---------------------------------------------------------------\n');
s = -10 + 1000*1i;

fprintf(1, 'v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B]\n');
fprintf(1, '---------------------------------------------------------------\n');
Bhat = B_qo;

tic
% First order solve; no structure taken into account
fosolve = (s * E_qo - A_qo)\Bhat;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc);
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
% rhs has form Bhat = [0; B], here
sosolve = so_structured_solve(M, E, K, B1, s, 0, 1);

fprintf(1, 'Relative solution error, v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B] %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error, v = (s*E - A)^{-1}*Bhat, rhs Bhat = [0; B] %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, 'w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] \n');
fprintf(1, '---------------------------------------------------------------\n');
Bhat = Q_qo*fosolve;
tic
% First order solve; no structure taken into account
fosolve = ((s * E_qo - A_qo)')\Bhat;
fprintf(1, 'Single first-order solve finished in %.2f s\n',toc);
fprintf(1, '---------------------------------------------------------------\n');

% Second order solve accounting for structure
% rhs has form Bhat = [B; 0], here
B2      = Bhat(1:n, :);
sosolve = so_structured_solve(M, E, K, B2, s, 1, 1);

fprintf(1, 'Relative solution error, w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] %.16f \n', norm(sosolve-fosolve)/norm(fosolve));
fprintf(1, '---------------------------------------------------------------\n');
fprintf(1, 'Absolute solution error, w = ((s*E - A)^H)^{-1}*Bhat, rhs Bhat = [B; 0] %.16f \n', norm(sosolve-fosolve));
fprintf(1, '---------------------------------------------------------------\n');

fprintf(1, '\n')

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off