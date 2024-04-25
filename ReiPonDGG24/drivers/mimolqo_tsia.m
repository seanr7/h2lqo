function [Ar, Br, Cr, Mr, info] = mimolqo_tsia(A, B, C, M, r, opts)
% MIMOLQO_TSIA Two-sided iteration algorithm for model-order reudction of
% linear systems with multiple quadratic outputs
%
% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Ar, Br, Cr, Mr)
%   using the two-sided iteration algorithm given in 
%   "$\mathcal{H}_2%-optimal model reduction of linear systems with 
%   multiple quadratic outputs"
%   At each iteration, two Sylvester equations
%
%       A*X + X*Ar' + B*Br' = 0                                   (1)
%       A'*Z + Z*Ar - 2*M1*X*M1r - ... - 2*Mp*X*Mpr - C*Cr' = 0   (2)
%
%   Are solved, and a linear quadratic output reduced model is obtained via
%   projeciton W = Z and V = X.
%   It is assumed that the eigenvalues of (s*I-A) lie in the open left
%   half-plane.
%
% INPUTS:
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
%   |                 | based on the H2error of the reduced model         |
%   |                 | (default 10e-4)                                   |
%   +-----------------+---------------------------------------------------+
%   | maxiter         | positive integer, maximum number of iteration     |
%   |                 | steps                                             |
%   |                 | (default 100)                                     |
%   +-----------------+---------------------------------------------------+
%   | Ar              | Initial state matrix                              |
%   |                 | (default -diag(logspace(1,3,r))                   |
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
%   Ar   - reduced state matrix with dimensions r x r in (1)
%   Br   - reduced descriptor matrix with dimensions r x m in (1)
%   Cr   - reduced linear output matrix with dimensions p x r in (2)
%          If C is zero then Cr is zeros(p, r)
%   Mr   - 3d-array of reduced (symmetric) quadratic output matrices 
%          with dimensions p x r x r in (2) 
%          If M is zero then Mr is zeros(r, r)
%   info - structure, containing the following information for monitoring
%          convergence
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | errs            | history of squared relative H2 errors throughout  |
%   |                 | the iteration                                     |
%   +-----------------+---------------------------------------------------+
%   | tails           | variable parts of the relative H2 error           |
%   +-----------------+---------------------------------------------------+ 
%   | gradientAr      | absolute gradient w.r.t Ar of the squared H2 error|
%   +-----------------+---------------------------------------------------+ 
%

%
% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 4/23/2024
%

%%
% Grab state, input, output dimensions
% Check input matrices.
n = size(A, 1);
m = size(B, 2);
if isempty(C)
    try
        [~, ~, p] = size(M, 1);
    catch
        p = 1;
    end
else
    p = size(C, 1);
end

% Check and set inputs
if (nargin < 6) 
    opts = struct(); % Empty struct 
end

if ~isfield(opts, 'tol')
    opts.tol = 10e-6;
end
if ~isfield(opts, 'maxiter')
    opts.maxiter = 100;
end
if ~isfield(opts, 'Ar')
    opts.Ar = -diag(logspace(-3, -1, r)); 
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

% If no linear output matrix, set to zero
if isempty(C)
    C = zeros(p, n);
end

%% Begin algorithm.
overall_start = tic;
fprintf(1, 'Initializing algorithm\n')
fprintf(1, '----------------------\n');
Ar = opts.Ar;   Br = opts.Br;   Cr = opts.Cr;   Mr = opts.Mr;  

% Precompute the H2 norm of the full-order model for rel error tracking
norm_start = tic;
fprintf(1, 'Precomputing H2 norm of full order model\n')
fprintf(1, '----------------------------------------\n')
P   = lyap(A, B*B');
rhs = C'*C; 
for i = 1:p
    rhs = rhs + M(:, :, i)*P*M(:, :, i);
end
Q  = lyap(A', rhs);
h  = trace(B'*Q*B); % h = ||G||_{\mathcal{H}_2}^2
fprintf(1, 'Norm computed in %d s\n', toc(norm_start))
fprintf(1, '--------------------------\n') 

% Compute sqaured realtive H2 error due to initial reduced model
%   ||G - Gr||_{\mathcal{H}_2}^2/||G||_{\mathcal{H}_2}^2 = 
%       (h + tr(Br'*Qr*Br) + 2tr(B'*Y*Br))/h
% Here, Y denotes the off diagonal block of the quadratic output error 
% Gramian Q
% 1. Computed reduced quadratic output observability Gramian Qr and Y
Pr    = lyap(Ar, Br*Br'); 
X     = lyap(A, Ar', B*Br');
rhsQr = Cr'*Cr; 
rhsY  = -C'*Cr;
% For computing the gradient w.r.t Ar
Qrone = lyap(Ar', rhsQr);
Yone  = lyap(A', Ar, rhsY);
for i = 1:p
    rhsQr = rhsQr + Mr(:, :, i)*Pr*Mr(:, :, i);
    rhsY  = rhsY - M(:, :, i)*X*Mr(:, :, i);
end
Qr = lyap(Ar', rhsQr);
Y  = lyap(A', Ar, rhsY);

% 2. Compute squared relative H2 error from trace formula
iter             = 1;
tails(iter)      = trace(Br'*Qr*Br) + 2*trace(B'*Y*Br);               % `Tail' of H2 error
gradientAr(iter) = norm(2*((2*Qr - Qrone)*Pr + (2*Y' - Yone')*X), 2); % Gradient of squared H2 error w.r.t Ar
errs(iter)       = (h + tails(iter))/h;                               % Squared relative H2 error

% Set convergence criterion to enter while
changein_errs(iter) = opts.tol + 1;

while (changein_errs(iter) > opts.tol && iter <= opts.maxiter)
    iter_start = tic; % Start timer on this iteration
    fprintf(1, 'Current iterate is k = %d\n', iter)
    fprintf(1, '---------------------------------------\n')
    % Solve equation (1) for n x r right projection matrix X
    X = lyap(A, Ar', B*Br');

    % Solve equation (2) for n x r left projection matrix Z
    rhs = -C'*Cr;
    for i = 1:p
        rhs = rhs - 2*M(:, :, i)*X*Mr(:, :, i);
    end
    Z = lyap(A', Ar, rhs);
    
    % Orthonormalize projection matrices ...
    [V, ~] = qr(X, "econ");     [W, ~] = qr(Z, "econ");
    % ... and compute reduced model via projection
    obl = W'*V;
    Ar  = (obl)\(W'*A*V);   Br = (obl)\(W'*B);   Cr = C*V;   
    for i = 1:p
        Mr(:, :, i) = V'*M(:, :, i)*V;    
    end
  
    iter = iter + 1;  
    fprintf(1, 'Computing H2 error at the current iteration\n')
    fprintf(1, '-------------------------------------------\n')
    % Compute H2 error at this iteration
    %   ||G - Gr||_{\mathcal{H}_2}^2 = h + tr(Br'*Qr*Br) + 2tr(B'*Y*Br)
    % Here, Y denotes off diagonal block of error Gramian Q
    % 1. Computed reduced quadratic output observability Gramian Qr and Y
    Pr    = lyap(Ar, Br*Br'); 
    rhsQr = Cr'*Cr; 
    rhsY  = -C'*Cr;
    % For computing the gradient w.r.t Ar
    Qrone = lyap(Ar', rhsQr);
    Yone  = lyap(A', Ar, rhsY);
    for i = 1:p
        rhsQr = rhsQr + Mr(:, :, i)*Pr*Mr(:, :, i); % No factor of 2 for Y
        rhsY  = rhsY - M(:, :, i)*X*Mr(:, :, i);    % X was computed above
    end
    Qr = lyap(Ar', rhsQr);
    Y  = lyap(A', Ar, rhsY);
    % 2. Compute squared relative H2 error from trace formula
    tails(iter)      = trace(Br'*Qr*Br) + 2*trace(B'*Y*Br);               % `Tail' of H2 error
    gradientAr(iter) = norm(2*((2*Qr - Qrone)*Pr + (2*Y' - Yone')*X), 2); % Gradient of squared H2 error w.r.t Ar
    errs(iter)       = abs((h + tails(iter))/h);                          % Squared relative error
    fprintf(1, 'Squared relative H2 error of current tsia reduced model is ||G - Gr||_H2^2/||G||_H2^2 = %.12f\n', ...
        errs(iter));

    % Compare change in relative H2 error with previous iteration to track
    % convergence
    changein_errs(iter) = abs(errs(iter) - errs(iter-1))/errs(1);
    fprintf('Relative change in relative squared H2 errors is |e(j) - e(j-1)|/e(1) = %.12f \n', changein_errs(iter))
    fprintf(1, '----------------------------------------\n')

    % End the clock
    fprintf(1, 'Current iterate finished in %.2f s\n', toc(iter_start))
    fprintf(1, 'End of current iterate k = %d\n', iter)
    fprintf(1, '---------------------------------------\n')
end

% Save output info
info            = struct();
info.errs       = errs;    
info.tails      = tails;
info.gradientAr = gradientAr;

if iter == (opts.maxiter + 1)
    fprintf('Algorithm has terminated due to reaching the max no. of iterations; total time elapsed is %.2f s\n', ...
        toc(overall_start))
    fprintf(1, '---------------------------------------\n')
else
    fprintf('Algorithm has converged in %d iterations\n', iter)
    fprintf('Total time elapsed is %.2f s\n', toc(overall_start))
    fprintf(1, '---------------------------------------\n')
end
end