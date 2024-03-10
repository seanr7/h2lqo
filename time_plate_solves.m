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
A_qo(1:n, n+1:2*n) = speye(n); % (1, 2) block of A_qo
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B;
% No scalar output in this example; only QO

% Our `M' matrix (i.e., the quadratic output matrix) is C' * C
M_qo = spalloc(2*n, 2*n, nnz(C' * C));
M_qo(1:n, 1:n) = C' * C; % Double check this...

%% Time sparse-dense Sylvester solve

% Make it dense, to handle the `worst-case'
r = 250;
tmp = rand(r, r);
A_qo_r = (tmp\(-1*diag(logspace(0, 3, r))*tmp));
M_qo_r = eye(r, r);
B_qo_r = ones(r, 1);
E_qo_r = eye(r, r);

tic
fprintf('Beginning sparse-dense Sylvester solve for X\n')
[X, ~, ~, ~, ~] = mess_sylvester_sparse_dense(A_qo, 'N', A_qo_r, 'T', -B_qo*B_qo_r', E_qo, E_qo_r);
fprintf('Sparse-dense Sylvester solve of left projection matrix, X, finished in %.2f s\n', toc)
% And, Z \in \Rnr satisfies 
%   @math: A'*Z*Er + E'*Z*Ar - 2*M*X*Mr - C*Cr' = 0
rhs = -2*M_qo*X*M_qo_r';

fprintf('Beginning sparse-dense Sylvester solve for Z')
[Z, ~, ~, ~, ~] = mess_sylvester_sparse_dense(A_qo, 'T', A_qo_r, 'N', rhs, E_qo, E_qo_r);
fprintf('Sparse-dense Sylvester solve of right projection matrix, Z, finished in %.2f s\n', toc)

% Orthonormalize projection matrices

[Vr, ~] = qr(X, "econ");    

%% Get idea of spectrum of A
% k = 10;
% tic
% L_small_real = eigs(A_qo, E_qo, 10,'smallestreal');
% fprintf('Computation of %d eigenvalues of A, E with smallest real part finished in %.2f s\n', k, toc)
% fprintf('Eig-%d: %.8f\n', 1:k, L_small_real)
% tic
% L_small_abs = eigs(A_qo, E_qo, 10,'smallestabs');
% fprintf('Computation of %d eigenvalues of A, E with smallest magnitude finished in %.2f s\n', k, toc)
% fprintf('Eig-%d: %.8f\n', 1:k, L_small_abs)
% tic
% L_big_real = eigs(A_qo, E_qo, 10,'largestreal');
% fprintf('Computation of %d eigenvalues of A, E with largest real part finished in %.2f s\n', k, toc)
% fprintf('Eig-%d: %.8f\n', 1:k, L_big_real)
% tic
% L_big_abs = eigs(A_qo, E_qo, 10,'largestabs');
% fprintf('Computation of %d eigenvalues of A, E with largest magnitude finished in %.2f s\n', k, toc)
% fprintf('Eig-%d: %.8f\n', 1:k, L_big_abs)

%% Time one sparse linear solve + one IRKA iteration
% si = s(1);
% 
% % Start timer
% tic
% solve = (-si * E_qo - A_qo) \ B_qo;
% fprintf('Single linear solve finished in %.2f s\n', toc)
% 
% % How sparse is the original matrix + its resulting solution?
% fprintf('Number of nonzero entries in the resolvent %d \n', nnz(si * E_qo - A_qo))
% fprintf('So, resolvent is %f percent sparse \n', 100 * (n^2 - nnz(si * E_qo - A_qo)) / n^2)
% 
% fprintf('Nummber of nonzero entries in the solution is %d \n', nnz(solve))
% fprintf('So, solution is %f percent sparse \n', 100 * (n - nnz(solve)) / n)

%% Now, let's time an IRKA iteration, roughly

% r = 200;
% V_r = sparse(2*n, r);
% 
% % Start the clock
% overall_time = tic;
% fprintf('Beginning IRKA iteration')
% 
% % First, fill out columns of Vr; used in computing Wr
% k = 1;
% while k <= r
%     % 1. Construction of Vr enforces r + r^2 interpolation conditions:
%     %   @math: H1(poles(k)) = H1r(poles(k)), k = 1, ..., r
%     %   @math: H2(poles(i), poles(j)) = H2(poles(i), poles(j)), 
%     %           i, j = 1, ..., r
%     % Save time with spalloc?
%     this_iter_time = tic;
%     tmp = (-s(k) * E_qo - A_qo) \ B_qo; 
%     V_r(:, k) = spalloc(2*n, 1, nnz(tmp));
%     V_r(:, k) = tmp;
%     k = k + 1;
%     fprintf('Current iterate finished in %.2f s\n',toc(this_iter_time))
% end
% 
% tic
% fprintf('Timing QR pass')
% [Q_r, ~] = qr(V_r, "econ");   
% fprintf('QR finished in %.2f s\n', toc)
% % How sparse is the original matrix + its resulting solution?
% fprintf('Is Q even sparse? %f \n', issparse(Q_r))
% fprintf('Number of nonzero entries in orth basis %d \n', nnz(Q_r))
% 
% fprintf('Total elapsed time: %.2f s\n',toc(overall_time))


%%
% fprintf('Now, trying iteration with a different sparsity assignment scheme')
% 
% r = 200;
% 
% % Start the clock
% overall_time = tic;
% fprintf('Beginning IRKA iteration')
% 
% % First, fill out columns of Vr; used in computing Wr
% k = 1;
% while k <= r
%     % 1. Construction of Vr enforces r + r^2 interpolation conditions:
%     %   @math: H1(poles(k)) = H1r(poles(k)), k = 1, ..., r
%     %   @math: H2(poles(i), poles(j)) = H2(poles(i), poles(j)), 
%     %           i, j = 1, ..., r
%     % Save time with spalloc?
%     this_iter_time = tic;
%     tmp = (-s(k) * E_qo - A_qo) \ B_qo; 
%     [i, j, values] = find(tmp); % Get row, col, and values of nz elements of sparse matrix
%     if k == 1
%         V_r_alt = sparse(i, j, values, 2*n, 1);
%     else
%         V_r_alt = [V_r_alt, sparse(i, j, values, 2*n, 1)];
%     end
%     k = k + 1;
%     fprintf('Current iterate finished in %.2f s\n',toc(this_iter_time))
% end
% 
% fprintf('Are they the same projection matrices?: %.2f \n', norm(V_r - V_r_alt, 2))
% 
% tic
% fprintf('Timing QR pass')
% [Q_r, ~] = qr(V_r_alt, "econ");   
% fprintf('QR finished in %.2f s\n', toc)
% % How sparse is the original matrix + its resulting solution?
% fprintf('Is Q even sparse? %f \n', issparse(Q_r))
% fprintf('Number of nonzero entries in orth basis %d \n', nnz(Q_r))
% 
% fprintf('Total elapsed time: %.2f s\n',toc(overall_time))
