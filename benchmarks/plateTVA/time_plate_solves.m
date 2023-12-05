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

%% Now, let's time an IRKA iteration, roughly


V_r = sparse(2*n, r);
r = 200;

% Start the clock
overall_time = tic;
fprintf('Beginning IRKA iteration')

% First, fill out columns of Vr; used in computing Wr
k = 1;
while k <= r
    % 1. Construction of Vr enforces r + r^2 interpolation conditions:
    %   @math: H1(poles(k)) = H1r(poles(k)), k = 1, ..., r
    %   @math: H2(poles(i), poles(j)) = H2(poles(i), poles(j)), 
    %           i, j = 1, ..., r
    % Save time with spalloc?
    tic = this_iter_time;
    tmp = (-s(k) * E_qo - A_qo) \ B_qo; 
    V_r(:, k) = spalloc(2*n, 1, nnz(tmp));
    V_r(:, k) = tmp;
    k = k + 1;
    fprintf('Current iterate finished in %.2f s\n',toc(this_iter_time))
end

tic
fprintf('Timing QR pass')
[Q_r, ~] = qr(V_r, "econ");   
fprintf('QR finished in %.2f s\n', toc)
% How sparse is the original matrix + its resulting solution?
fprintf('Is Q even sparse? %f \n', issparse(Q_r))
fprintf('Number of nonzero entries in orth basis %d \n', nnz(Q_r))

fprintf('Total elapsed time: %.2f s\n',toc(overall_time))

%%
fprintf('Now, trying iteration with a different sparsity assignment scheme')

V_r_alt = sparse(2*n, r);
r = 200;

% Start the clock
overall_time = tic;
fprintf('Beginning IRKA iteration')

% First, fill out columns of Vr; used in computing Wr
k = 1;
while k <= r
    % 1. Construction of Vr enforces r + r^2 interpolation conditions:
    %   @math: H1(poles(k)) = H1r(poles(k)), k = 1, ..., r
    %   @math: H2(poles(i), poles(j)) = H2(poles(i), poles(j)), 
    %           i, j = 1, ..., r
    % Save time with spalloc?
    tic = this_iter_time;
    tmp = (-s(k) * E_qo - A_qo) \ B_qo; 
    [i, j, values] = find(tmp); % Get row, col, and values of nz elements of sparse matrix
    V_r_alt(:, k) = sparse(i, j, values, 2*n, 1);
    k = k + 1;
    fprintf('Current iterate finished in %.2f s\n',toc(this_iter_time))
end

fprintf('Are they the same projection matrices?: %.2f \n', norm(V_r - V_r_alt, 2))

tic
fprintf('Timing QR pass')
[Q_r, ~] = qr(V_r_alt, "econ");   
fprintf('QR finished in %.2f s\n', toc)
% How sparse is the original matrix + its resulting solution?
fprintf('Is Q even sparse? %f \n', issparse(Q_r))
fprintf('Number of nonzero entries in orth basis %d \n', nnz(Q_r))

fprintf('Total elapsed time: %.2f s\n',toc(overall_time))