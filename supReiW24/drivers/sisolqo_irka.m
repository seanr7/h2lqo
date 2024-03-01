function [Er, Ar, br, cr, Mr, info] = sisolqo_irka(E, A, b, c, M, r, opts)
% SISOLQO_IRKA Iterative rational Krylov algorithm for model-order
% reudction of linear systems with a single quadratic output and single
% input
%
%
% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Er, Ar, br, cr, Mr)
%   via the iterative rational Krylov algorithm implemented in 
%   "Interpolatory model order reduction of large-scale dynamical systems
%   with root mean squared error measures"
%   At each iteration, construct n x r Petrov-Galerkin projection matrices
%   accroding to 
%
%       Vr(:, k) = ((-poles(k))*E - A)\b;                               (1)
%       Wr(:, k) = -fores(k)*((-poles(k))*E' - A')*c - sum_{j=1}^r 
%           sores(k, i)*(((-poles(k))*E' - A'))*M*((-poles(k))*E - A)\b (2)
%
%   Are solved, and a linear quadratic output reduced model is obtained via
%   projeciton. 
%   It is assumed that the eigenvalues of (s*E-A) lie in the open left
%   half-plane.

% INPUTS:
%   E    - invertible descriptor matrix with dimensions n x n in (1),
%          if empty set to eye(n, n)
%   A    - state matrix with dimensions n x n in (1)
%   b    - input matrix with dimensions n x 1 in (1)
%   c    - linear output matrix with dimensions n x 1 in (2)
%          if empty set to zeros(n, 1)
%   M    - symmetric quadratic output matrices with n x n in (2)
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
%   | plotconv        | bool, do we plot convergence of poles?            |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | poles           | Initial pole selection, r x 1 matrix              |
%   |                 | (default -(logspace(-1, 2, r)')                   |
%   +-----------------+---------------------------------------------------+
%   | fores           | Initial first-order residue selection, r x 1      |
%   |                 | matrix                                            |
%   |                 | (default if ~isempty(c) 10*rand(r,1)              |
%   |                 |          else zeros(r,1)                          |
%   +-----------------+---------------------------------------------------+
%   | sores           | Initial second-order residue selection, r x r     |
%   |                 | matrix                                            |
%   |                 | (default tmp = 10*rand(r, r);                     |
%   |                 |        sores = (tmp + tmp')/2)                    |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Er   - reduced descriptor matrix with dimensions n x n in (1),
%   Ar   - reduced state matrix with dimensions r x r in (1)
%   br   - reduced descriptor matrix with dimensions r x 1 in (1)
%   cr   - reduced linear output matrix with dimensions r x 1 in (2)
%           If c is zero then cr is zeros(1, r)
%   Mr   - reduced symmetric quadratic output matrix with dimensions r x r 
%           in (2)
%   info - structure, containing the following optional entries:
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | pole_hist       | Convergence history of reduced-order poles        |
%   +-----------------+---------------------------------------------------+
%   | fores           | Final first-order residues                        |
%   +-----------------+---------------------------------------------------+
%   | sores           | Final second-order residues                       |
%   +-----------------+---------------------------------------------------+

% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)

% Virginia Tech, Department of Mathematics
% Last editied: 2/29/2024

%%
% Grab state, input, output dimensions
% Check input matrices.
n = size(A, 1);

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
if ~isfield(opts, 'plotconv')
    opts.plotconv = 1;
end
if ~isfield(opts, 'poles')
    opts.poles = -(logspace(-1, 2, r)'); 
end
if ~isfield(opts, 'fores')
    if isempty(c)
        opts.fores = zeros(r, 1);
    else
        opts.fores = 10*rand(r, 1);  
    end
end 
if~isfield(opts, 'sores')
    tmp = 10*rand(r, r);
    opts.sores = (tmp + tmp')/2;
end

if isempty(E)
    E = eye(n, n);
end
pureqo = 0;
if isempty(c)
    % Set bool for output being purely quadratic
    pureqo = 1;
    cr = zeros(r, 1);   
    % Safety check, enforce zero resides
    opts.fores = zeros(r, 1);
end
%% Begin algorithm.
overall_start = tic;
fprintf(1, 'Initializing algorithm\n')
fprintf(1, '----------------------\n');
poles = opts.poles;   fores = opts.fores;   sores = opts.sores;  
pole_hist(:, 1) = poles;

% Set counter and tolerance to enter while
iter = 1;   err(iter) = opts.tol + 1; 
while (err(iter) > opts.tol && iter <= opts.maxiter)
    iter_start = tic; % Start timer on this iteration
    fprintf(1, 'Current iterate is k = %d\n', iter)
    fprintf(1, '---------------------------------------\n');
    % Space allocation for left and right projection matrices
    Vr = zeros(n, r);     Wr = zeros(n, r);

    % Fill out columns in Vr in (1); pre-computed for constructing Wr
    k = 1;
    while k <= r
        Vr(:, k) = ((-poles(k)) * E - A)\b; 
        k = k + 1;
    end
    % Fill out columns of Wr in (2)
    k = 1; 
    while k <= r
        tmp = zeros(n, 1); 
        i = 1;
        while i <= r 
            % Grab columns of Vr, sum sores(i,k)*(-poles(i) * E - A)\b over
            % i = 1, ..., r
            tmp = tmp + (sores(i, k) * Vr(:, i)); 
            i = i + 1;
        end
        % Premultiply by - 2 * M
        tmp = - 2 * M * tmp; 
        % If output is not purely quadratic, account for linear term
        if ~pureqo 
            tmp = tmp - fores(k) * c';
        end
        % Compute final linear solve against summed term
        Wr(:, k) = ((-poles(k) * E - A)'\tmp);
        k = k + 1;
    end
    
    % Orthonormalize projection matrices
    [Vr, ~] = qr(Vr, "econ");     [Wr, ~] = qr(Wr, "econ");
    % Compute reduced model via projection
    Wrt = Wr';  Vrt = Vr';
    Er = Wrt * E * Vr;   Ar = Wrt * A * Vr;   br = Wrt * b;
    if ~pureqo % Compute reduced linear output term
        cr = c * Vr;  
        fores_prev = fores;
    end
    Mr = Vrt * M * Vr;

    % Save poles and residues from last iteration
    poles_prev = poles; sores_prev = sores;
    
    % Solve r x r generalized eigenvalue problem s*Er - Ar to get new 
    % reduced-order poles and residues
    [Xr, Lr] = eig(Ar, Er);     poles = diag(Lr);
    try % Attempt to sort poles into complex conjugate pairs
        poles = cplxpair(poles, 10e-3); 
    catch % If model has become overtly complex, throw 
        warning('Warning! Reduced model has become overtly complex; cannot sort eigenvalues into complex conjugate pairs. Returning unsorted poles')
    end
    % Track poles
    pole_hist(:, iter+1) = poles;

    Xrinv = Xr\eye(r, r); mod_br = Er\br;
    if ~pureqo % Compute FO residues
        fores_prev = fores; 
        fores = (cr * Xr)'.*(Xrinv * (mod_br));    
    end
    sores = ((Xrinv * (mod_br))'.*((Xr' * Mr * Xr))).*(Xrinv * (mod_br));

    % End the clock
    fprintf(1, 'Current iterate finished in %.2f s\n',toc(iter_start))
    fprintf(1, 'End of current iterate k = %d\n', iter)
    fprintf(1, '---------------------------------------\n');

    % Track convergence of poles
    iter = iter + 1;    
    err(iter) = max(abs(poles - poles_prev));
    fprintf('Change in poles is is %.2f \n', err(iter))
    fprintf(1, '---------------------------------------\n');
end

% Once broken, poles and residues from the previous iteration are used to
% check interpolation. So, return these in info struct
info           = struct();
info.pole_hist = pole_hist;    
info.sores     = sores_prev;
if ~pureqo
    info.fores = fores_prev;
else
    info.fores = zeros(r, 1);
end

if iter == (opts.maxiter + 1)
    fprintf('Algorithm has terminated due to reaching the max no. of iterations; total time elapsed is %.2f s\n', toc(overall_start))
    fprintf(1, '---------------------------------------\n');
else
    fprintf('Algorithm has converged in %d iterations\n', iter)
    fprintf('Total time elapsed is %.2f s\n', toc(overall_start))
    fprintf(1, '---------------------------------------\n');
end

% Plot convergence of poles automatically if requested
if opts.plotconv == true
    semilogy(1:iter, err, '-o', LineWidth=1.5)
    xlim([1,iter])
    xlabel('Iteration count')
    ylabel('Magnitude of change in \lambda')
end
end