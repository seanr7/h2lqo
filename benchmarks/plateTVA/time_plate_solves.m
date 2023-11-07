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
B_qo(n+1:2*n, :) = B;
% No scalar output in this example; only QO

% Our `M' matrix (i.e., the quadratic output matrix) is C' * C
M_qo = spalloc(2*n, 2*n, nnz(C' * C));
M_qo(1:n, 1:n) = C' * C; % Double check this...

%% Time one sparse linear solve + one IRKA iteration
si = s(1);

% Start timer
tic
solve = (-si * E_qo - A_qo) \ B_qo;
fprintf('Single linear solve finished in %.2f s\n', toc)

% How sparse is the original matrix + its resulting solution?
fprintf('Number of nonzero entries in the resolvent %d \n', nnz(si * E_qo - A_qo))
fprintf('So, resolvent is %f percent sparse \n', 100 * (n^2 - nnz(si * E_qo - A_qo)) / n^2)

fprintf('Nummber of nonzero entries in the solution is %d \n', nnz(solve))
fprintf('So, solution is %f percent sparse \n', 100 * (n - nnz(solve)) / n)

tic
solve_mat = (-conj(si) * E_qo' - A_qo') \ M_qo;
fprintf('Linear solve against matrix M_qo finished in %.2f s\n', toc)
fprintf('Number of nonzero entries in solution is %d \n', nnz(solve_mat))
fprintf('So, solution matrix is %f percent spares \n', 100 * (n^2 - nnz(solve_mat)) / n^2)

tic
mat_vec = solve_mat * solve;
fin = toc;
fprintf('Matvec of arising in LQO-IRK finished in %.2f s\n', fin)

tic 
pre_comp = (-conj(si) * E_qo' - A_qo') \ speye(n);
fprintf('Pre-compute action of (-conj(si) * E_qo - A_qo)^{-1} against identity, finished in %.2f s\n', toc)
fprintf('Number of nonzer entries in result %d\n', nnz(pre_comp))
fprintf('So, result is %f percent sparse', 100 * (n^2 - nnz(pre_comp)) / n^2)

