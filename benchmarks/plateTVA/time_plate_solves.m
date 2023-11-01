%%
% Author: Sean Reiter (seanr7@vt.edu)
clear
close all
%% 
% Some intial timings...
load('plateTVA_n201900m1q28278_full')
n_nodes = full(sum(sum(C)));

%% Convert plate model to FO (first-order) from SO (second-order)
% Model is given in SO-form
% Necessarily, need to conver to FO to do LQO_IRKA for now
[n, ~] = size(M);

E_qo = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is mass matrix

A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(E) + n);  % A_qo = [0, I; -K, -E]
A_qo(1:n, 1:n) = speye(n); % (1, 1) block of A_qo
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
% No scalar output in this example; only QO

% Our `M' matrix (i.e., the quadratic output matrix) is C' * C
M_qo = C' * C; % Double check this...

%% Time one sparse linear solve + one IRKA iteration
si = s(1);

diary solveTimes
% Start timer
tic
solve = (-si * E_qo - A_qo)\B_qo; 
fprintf('Single linear solve (-si * E_qo - A_qo)\B_qo finished in %.2f s\n',toc)

diary off