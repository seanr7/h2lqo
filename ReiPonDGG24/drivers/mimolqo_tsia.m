function [Er, Ar, Br, Cr, Mr, poles] = mimolqo_tsia(E, A, B, C, M, r, opts)
% MIMOLQO_TSIA Two-sided iteration algorithm for model-order reudction of
% linear systems with multiple quadratic outputs

% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)

% Virginia Tech, Department of Mathematics
% Last editied: 2/27/2024

% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Er, Ar, Br, Cr, Mr)
%   using the two-sided iteration algorithm given in 
%   "$\mathcal{H}_2%-optimal model reduction of linear systems with 
%   multiple quadratic outputs"
%   At each iteration, two Sylvester equations
%
%       A*X*Er' + E*X*Ar' + B*Br' = 0                                   (1)
%       A'*Z*Er + E'*Z*Ar - 2*M1*X*M1r - ... - 2*Mp*X*Mpr - C*Cr' = 0   (2)
%
%   Are solved, and a linear quadratic output reduced model is obtained via
%   projeciton W = Z and V = X.
%   It is assumed that the eigenvalues of (s*E-A) lie in the open left
%   half-plane.

% INPUTS:
%   E    - invertible descriptor matrix with dimensions n x n in (1),
%          if empty set to eye(n, n)
%   A    - state matrix with dimensions n x n in (1)
%   B    - input matrix with dimensions n x m in (1)
%   C    - linear output matrix with dimensions p x n in (2)
%          if empty set to zeros(p, n)
%   M    - 3d-array of (symmetric) quadratic output matrices with 
%          dimensions p x n x n in (2)
%          if empty set to zeros(n, n)
%   r    - order of reduction
%   opts - structure, containing the following optional entries:
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | tol             | nonnegative scalar, tolerance for convergence     |
%   |                 | based on the poles of the reduced model           |
%   |                 | (default 10e-6)                                   |
%   +-----------------+---------------------------------------------------+
%   | maxiter         | positive integer, maximum number of iteration     |
%   |                 | steps                                             |
%   |                 | (default 100)                                     |
%   +-----------------+---------------------------------------------------+
%   | Er              | Initial descriptor matrix                         |
%   |                 | (default eye(r, r))                               |
%   +-----------------+---------------------------------------------------+
%   | Ar              | Initial state matrix                              |
%   |                 | (default -diag(logspace(-1,2,r))                  |
%   +-----------------+---------------------------------------------------+
%   | Br              | Initial input matrix                              |
%   |                 | (default eye(n, m))                               |
%   +-----------------+---------------------------------------------------+
%   | Cr              | Initial linear-output matrix                      |
%   |                 | (default eye(p, r))                               |
%   +-----------------+---------------------------------------------------+
%   | Mr              | Initial quadratic-output matrices                 |
%   |                 | (default repmat(eye(n,n), 1, 1, p)                |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Er    - reduced descriptor matrix with dimensions n x n in (1),
%   Ar    - reduced state matrix with dimensions r x r in (1)
%   Br    - reduced descriptor matrix with dimensions r x m in (1)
%   Cr    - reduced linear output matrix with dimensions p x r in (2)
%           If C is zero then Cr is zeros(p, r)
%   Mr    - 3d-array of reduced (symmetric) quadratic output matrices with 
%           dimensions p x r x r in (2) 
%           If M is zero then Mr is zeros(r, r)
%   poles - history of poles throughout the iteration 
%%
% Grab state, input, output dimensions
% Check input matrices.
n = size(A, 1);
m = size(B, 2);
p = size(C, 1);

if isempty(E)
    E = eye(n, n);
end

% Check and set inputs
if (nargin < 7) 
    opts = struct(); % Empty struct 
end

if ~isfield(opts, 'tol')
    opts.tol = 10e-6;
end
if ~isfield(opts, 'maxiter')
    opts.maxiter = 100;
end
if ~isfield(opts, 'Er')
    opts.Er = eye(r, r);
end
if ~isfield(opts, 'Ar')
    opts.Ar = -diag(logspace(-1, 2, r)); 
end
if ~isfield(opts, 'Br')
    opts.Br = eye(r, m);  
end
if ~isfield(opts, 'Cr')
    if isempty(C)
        opts.Cr = zeros(p, r);
    else
        opts.Cr = eye(p, r);
    end
end
if ~isfield(opts, 'Mr')
    if isempty(M)
        opts.Mr = repmat(zeros(r, r), 1, 1, p);
    else
        opts.Mr = repmat(eye(r, r), 1, 1, p);
    end
end

%% Begin algorithm.
overall_start = tic;
fprintf(1, 'Initializing algorithm\n')
fprintf(1, '----------------------\n');
Er = opts.Er;   Ar = opts.Ar;   Br = opts.Br;   Cr = opts.Cr;  
Mr = opts.Mr;  

% Set counter and tolerance to enter while
iter = 1;   err(iter) = eps + 1; 
while (err(iter) > opts.tol && iter <= opts.maxiter)
    iter_start = tic; % Start timer on this iteration
    fprintf(1, 'Current iterate is k = %d\n', iter)
    fprintf(1, '---------------------------------------\n');
    % Solve equation (1) for n x r right projection matrix X
    [X, ~, ~, ~, ~] = mess_sylvester_sparse_dense(A, 'N', Ar, 'T', B*Br', ...
        E, Er);
    % For some reason; Convergence is VERY poor if -B*Br' not present...

    % Solve equation (2) for n x r left projection matrix Z
    rhs = -C'*Cr;
    for i = 1:p
        rhs = rhs - 2*M(:, :, i)*X*Mr(:, :, i);
    end
    [Z, ~, ~, ~, ~] = mess_sylvester_sparse_dense(A, 'T', Ar, 'N', rhs, E, Er);
    
    % Orthonormalize projection matrices
    [V, ~] = qr(X, "econ");     [W, ~] = qr(Z, "econ");
    % Compute reduced model via projection
    Er = W'*E*V;   Ar = W'*A*V;   Br = W'*B;    Cr = C*V;  Mr = V'*M*V;

    % End the clock
    fprintf(1, 'Current iterate finished in %.2f s\n',toc(iter_start))
    fprintf(1, 'End of current iterate k = %d\n', iter)
    fprintf(1, '---------------------------------------\n');

    iter = iter + 1;    
    % Get eigenvalues of matrix pencil (s*Er-Ar) (reduced system poles)
    [~, Lr] = eig(Ar, Er);     poles(:, iter) = diag(Lr);
    err(iter) = max(abs(poles(:, iter) - poles(:, iter - 1)));
    fprintf('Change in poles is is %.2f \n', err(iter))
    fprintf(1, '---------------------------------------\n');
end

if iter == (opts.maxiter + 1)
    fprintf('Algorithm has terminated due to reaching the max no. of iterations; total time elapsed is %.2f s\n', toc(overall_start))
    fprintf(1, '---------------------------------------\n');
else
    fprintf('Algorithm has converged in %d iterations\n', iter)
    fprintf('Total time elapsed is %.2f s\n', toc(overall_start))
    fprintf(1, '---------------------------------------\n');
end
end