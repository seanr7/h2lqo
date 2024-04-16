function [Er, Ar, br, Qr, info] = lqoirka(Mso, Dso, Kso, bso, Qfo, r, opts)
% SISOLQO_IRKA Specialized implementation of the Iterative rational Krylov 
% algorithm for model-order reduction of linear quadratic output systems
% (lqoirka) with (1) A single, purely quadratic output (no linear
% component) (2) single input and (3) underlying second-order structure.
%
%
% DESCRIPTION:
%   Computes a linear quadratic output (lqo) reduced model (Er, Ar, br, Qr)
%   via the iterative rational Krylov algorithm (irka) implemented in 
%   "Interpolatory model order reduction of large-scale dynamical systems
%   with root mean squared error measures"
%   The full-order system matrices has the 2*n order first-order 
%   realization given by
%
%       Efo = [I  0;   0     Mso]; (1)
%       Afo = [I  0;  -Kso  -Dso]; (2)
%       bfo = [0; bso];            (3)
%
%   And Q is the first-order quadratic output matrix
%   At each iteration, construct n x r Petrov-Galerkin projection matrices
%   accroding to 
%
%       Vr(:, k) = ((-poles(k))*Efo - Afo)\bfo;                         (4)
%       Wr(:, k) = -2*sum_{j=1}^r sores(k, i)*
%         ((-poles(k))*Efo' - Afo')\(Qfo*(((-poles(k))*Efo - Afo)\bfo)) (5)
%
%   Are solved, and a linear quadratic output reduced model is obtained via
%   projeciton. 
%   The linear solves are not computed directly as in (4), (5). Instead, 
%   they are computed using the specialized companion function 
%   'so_structured_solve.m', which leverages the underlying second-order
%   structure of the first-order system matrices in (1) - (3). 
%   See 'help so_structured_solve' for details.
%   It is assumed that the eigenvalues of (s*Efo-Afo) lie in the open left
%   half-plane.
%
% INPUTS:
%   Mso  - sparse second order mass matrix with dimensions n x n in 
%          (1)
%   Dso  - sparse second order damping matrix with dimensions n x n 
%          in (2)
%   Kso  - sparse second order stiffness matrix with dimensions n x n 
%          in (2)
%   bso  - sparse second order input matrix with dimensions n x 1 in 
%          (3)
%   Qfo  - sparse symmetric quadratic output matrices with dimensions 
%          2*n x 2*n, such that Qfo(1:n, 1:n) is nonzero, Qfo is zero 
%          elsewhere
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
%   | sores           | Initial second-order residue selection, r x r     |
%   |                 | matrix                                            |
%   |                 | (default tmp = 10*rand(r, r);                     |
%   |                 |        sores = (tmp + tmp')/2)                    |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Er   - reduced descriptor matrix with dimensions r x r 
%   Ar   - reduced state matrix with dimensions r x r 
%   br   - reduced descriptor matrix with dimensions r x 1 
%   Qr   - reduced symmetric quadratic output matrix with dimensions r x r 
%   info - structure, containing the following entries:
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | pole_hist       | Convergence history of reduced-order poles        |
%   +-----------------+---------------------------------------------------+
%   | sores           | Final second-order residues                       |
%   +-----------------+---------------------------------------------------+

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model order reduction of large-scale 
% dynamical systems with root mean squared error measures"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

% Virginia Tech, Department of Mathematics
% Last editied: 4/15/2024

%%
% Grab state, input, output dimensions
% Check input matrices.
n = size(Mso, 1);

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
if~isfield(opts, 'sores')
    tmp = 10*rand(r, r);
    opts.sores = (tmp + tmp')/2;
end

% Build first-order realization for later projection
Efo                   = spalloc(2*n, 2*n, nnz(Mso) + n); % Descriptor matrix; Efo = [I, 0: 0, Mso]
Efo(1:n, 1:n)         = speye(n);                        % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = Mso;                             % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(Kso) + nnz(Dso) + n); % Afo = [0, I; -Kso, -Dso]
Afo(1:n, n+1:2*n)     = speye(n);                                   % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -Kso;                                       % (2, 1) block is -Kso
Afo(n+1:2*n, n+1:2*n) = -Dso;                                       % (2, 2) block is -Dso 

bfo             = spalloc(2*n, 1, nnz(bso)); % bfo = [0; bso];
bfo(n+1:2*n, :) = bso; 

%% Begin algorithm.
overall_start = tic;
fprintf(1, 'Initializing algorithm\n')
fprintf(1, '----------------------\n');
poles           = opts.poles;  
sores           = opts.sores;  
pole_hist(:, 1) = poles;

% Set counter and tolerance to enter while
iter      = 1;   
err(iter) = opts.tol + 1; 
while (err(iter) > opts.tol && iter <= opts.maxiter)
    iter_start = tic; % Start timer on this iteration
    fprintf(1, 'CURRENT ITERATE IS k = %d\n', iter)
    fprintf(1, '---------------------------------------\n');
    % Space allocation for left and right projection matrices
    Vr = zeros(2*n, r);     Wr = zeros(2*n, r);

    % Fill out columns in Vr in (4); pre-computed for constructing Wr
    for k = 1:r
        fprintf(1, 'Computing column i = %d of Vr\n', k)
        fprintf(1, '---------------------------------------\n');
        if abs(poles(k)) < 10e-4
            fprintf(1, 'WARNING! Current pole is small in absolute value; structured linear solve may be inaccurate!')
        end
        % Option 0 in 'so_structured_solve.m' implements (-poles(k) * E - A)\b)
        v        = so_structured_solve(Mso, Dso, Kso, bso, -poles(k), 0, 1);
        Vr(:, k) = v;
    end
    % Fill out columns of Wr in (5)
    for k = 1:r
        tmp = zeros(2*n, 1); 
        for i = 1:r
            % Grab columns of Vr, sum sores(i,k)*(-poles(i) * E - A)\b over
            % i = 1, ..., r
            tmp = tmp + (sores(i, k) * Vr(:, i)); 
        end
        % Premultiply by - 2 * Q
        tmp = - 2 * Qfo * tmp; 
        fprintf(1, 'Computing column i = %d of Wr\n', k)
        fprintf(1, '---------------------------------------\n');
        % Option 1 in 'so_structured_solve.m' implements ((-poles(k) * E - A)'\tmp)
        rhs      = tmp(1:n, :);
        w        = so_structured_solve(Mso, Dso, Kso, rhs, -poles(k), 1, 1);
        Wr(:, k) = w;
    end
    
    % Orthonormalize projection matrices
    [Vr, ~] = qr(Vr, "econ");     [Wr, ~] = qr(Wr, "econ");
    % Compute reduced model via projection
    Wrt = Wr';  Vrt = Vr';
    fprintf(1, 'Computing reduced model via projection!\n')
    fprintf(1, '---------------------------------------\n');
    % To project down, build first order realization from second order
    Er = Wrt * Efo * Vr;   Ar = Wrt * Afo * Vr;   br = Wrt * bfo;
    Qr = Vrt * Qfo * Vr;

    % Save poles and residues from last iteration
    poles_prev = poles; sores_prev = sores;
    
    % Solve r x r generalized eigenvalue problem s*Er - Ar to get new 
    % reduced-order poles and residues
    [Xr, Lr] = eig(Ar, Er);     poles = diag(Lr);
    try % Attempt to sort poles into complex conjugate pairs
        poles = cplxpair(poles, 10e-3); 
    catch % If model has become overtly complex, throw 
        warning('WARNING! Reduced model has become overtly complex; cannot sort eigenvalues into complex conjugate pairs. Returning unsorted poles')
    end
    % Track poles
    pole_hist(:, iter+1) = poles;

    Xrinv = Xr\eye(r, r); mod_br = Er\br;
    sores = ((Xrinv * (mod_br))'.*((Xr' * Qr * Xr))).*(Xrinv * (mod_br));

    % End the clock
    fprintf(1, 'CURRENT ITERATE FINISHED IN %.2f s\n', toc(iter_start))
    fprintf(1, 'END OF CURRENT ITERATE k = %d\n', iter)
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
