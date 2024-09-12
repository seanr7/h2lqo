function [Er, Ar, Br, Qr, info] = siso_lqoirka(Mso, Dso, Kso, Bso, Qfo, r, ...
                                                opts)
%SISO_LQOIRKA Specialized implementation of the iterative rational Krylov
% algorithm for single input, single (quadratic) output systems described
% by first-order ODEs with underlying second-order structure.
%
% SYNTAX:
%   [Er, Ar, Br, Qr, info] = siso_lqoirka(Mso, Dso, Kso, Bso, Qfo, r, opts)
%   [Er, Ar, Br, Qr]       = siso_lqoirka(Mso, Dso, Kso, Bso, Qfo, r)
%
% DESCRIPTION:
%   Function computes a LQO reduced-order model described by the state
%   quadruple (Er, Ar, Br, Qr) via LQO-IRKA. It is assumed that the
%   full-order matrix operators enjoy are a order 2*n linearization of a
%   second-order system, so that
%
%       Efo = [I  0;   0     Mso]; (1)
%       Afo = [I  0;  -Kso  -Dso]; (2)
%       Bfo = [0; Bso];            (3)
%
%   And Qfo is the first-order quadratic output matrix (e.g., Qfo =
%   Cpso'*Cpso).
%
%   At each iteration, an interpolatory projection is obtained using V and
%   W computed by:
%
%       V(:, k) = (mirroredPoles(k))*In - Afo)\Bfo;                    (4)
%       W(:, k) = ((mirroredPoles(k))*In - Afo.')\(2*(soRes(k, 1)*Qfo*V(:, i) + ... 
%                   + soRes(k, 1)*Qfo*V(:, r))                         (5) 
%
%   Are solved, and a linear quadratic output reduced model is obtained via
%   projeciton. 
%   The linear solves are not computed naively as in (4), (5). Instead, 
%   they are computed using the specialized companion function 
%   'so_structured_solve.m', which leverages the underlying second-order
%   structure of the first-order system matrices in (1) - (3). 
%   See 'help so_structured_solve' for details.
%   It is assumed that the eigenvalues of (s*Efo-Afo) lie in the open left
%   half-plane, and Qfo is symmetric.
%
% INPUTS:
%   Mso  - sparse second order mass matrix with dimensions n x n in 
%          (1)
%   Dso  - sparse second order damping matrix with dimensions n x n 
%          in (2)
%   Kso  - sparse second order stiffness matrix with dimensions n x n 
%          in (2)
%   Bso  - sparse second order input matrix with dimensions n x 1 in 
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
%   |                 | based on the change in reduced model poles        |
%   |                 | (default 10e-6)                                   |
%   +-----------------+---------------------------------------------------+
%   | maxIter         | positive integer, maximum number of iteration     |
%   |                 | steps                                             |
%   |                 | (default 100)                                     |
%   +-----------------+---------------------------------------------------+
%   | plotConv        | bool, do we plot convergence of poles?            |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | checkInterp     | bool, do we check interpolation conditions using  |
%   |                 | parameters from the previous iterate?             | 
%   |                 | ** not recommend for large problems **            |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | poles           | Initial pole selection, r x 1 matrix              |
%   |                 | (default -(logspace(-1, 5, r)')                   |
%   +-----------------+---------------------------------------------------+
%   | qoRes           | Initial r x r matrix of quadratic-output transfer |
%   |                 | function residues                                 |
%   |                 | (default tmp = 10*rand(r, r);                     |
%   |                 |        sores = (tmp + tmp')/2)                    |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Er   - reduced descriptor matrix with dimensions r x r 
%   Ar   - reduced state matrix with dimensions r x r 
%   Br   - reduced descriptor matrix with dimensions r x 1 
%   Qr   - reduced symmetric quadratic output matrix with dimensions r x r 
%   info - structure, containing the following entries:
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | poleHistory     | Convergence history of reduced-order poles        |
%   +-----------------+---------------------------------------------------+
%   | qoRes           | Final second-order residues                       |
%   +-----------------+---------------------------------------------------+
%   | poleChange      | Convergence history of max error between poles    |
%   +-----------------+---------------------------------------------------+
%

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model reduction of dynamical systems 
% with root mean squared error"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

% Virginia Tech, Department of Mathematics
% Last editied: 9/10/2024

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State dimension.
n = size(Mso, 1);

% Check and set inputs.
if (nargin < 7) 
    opts = struct(); 
end
if ~isfield(opts, 'tol')
    opts.tol = 10e-6;
end
if ~isfield(opts, 'maxIter')
    opts.maxIter = 100;
end
if ~isfield(opts, 'plotConv')
    opts.plotConv = 1;
end
if ~isfield(opts, 'checkInterp')
    opts.checkInterp = 1;
end
if ~isfield(opts, 'poles')
    opts.poles = -(logspace(-1, 5, r)'); 
end
if~isfield(opts, 'qoRes')
    tmpSum = 10*rand(r, r);
    opts.qoRes = (tmpSum + tmpSum')/2;
end

% First-order realization.
Efo                   = spalloc(2*n, 2*n, nnz(Mso) + n); % Descriptor matrix; Efo = [I, 0: 0, Mso]
Efo(1:n, 1:n)         = speye(n);                        % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = Mso;                             % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(Kso) + nnz(Dso) + n); % Afo = [0, I; -Kso, -Dso]
Afo(1:n, n+1:2*n)     = speye(n);                                   % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -Kso;                                       % (2, 1) block is -Kso
Afo(n+1:2*n, n+1:2*n) = -Dso;                                       % (2, 2) block is -Dso 

Bfo             = spalloc(2*n, 1, nnz(Bso)); % Bfo = [0; Bso];
Bfo(n+1:2*n, :) = Bso; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overallStart = tic;
fprintf(1, 'INITIALIZING ALGORITHM.\n')
    fprintf(1, '--------------------------------------------------\n');
mirroredPoles     = -opts.poles;
qoRes             = opts.qoRes;  
poleHistory(:, 1) = mirroredPoles;

% Set iteration count and tolerance to engage `while'.
iterate             = 1;   
poleChange(iterate) = opts.tol + 1; 

while (poleChange(iterate) > opts.tol && iterate <= opts.maxIter)
    thisIterStart = tic;
    fprintf(1, 'CURRENT ITERATE IS k = %d\n', iterate)
    fprintf(1, '---------------------------------------\n');

    % Space allocation for model reduction bases.
    VPrim = zeros(2*n, r);     WPrim = zeros(2*n, r);

    % Fill out columns in V in (4); pre-computed for constructing W.
    for k = 1:r
        fprintf(1, 'COMPUTING COLUMN k = %d of V\n', k)
        fprintf(1, '---------------------------------------\n');
        if abs(mirroredPoles(k)) < 10e-4
            fprintf(1, 'WARNING! Current pole is small in absolute value; structured linear solve may be inaccurate!')
        end
        % Option 0 in 'so_structured_solve.m' implements 
        % (mirroredPoles(k)*Efo - Afo)\Bfo).
        v           = so_structured_solve(Mso, Dso, Kso, Bso, mirroredPoles(k), ...
                                            0, 1);
        VPrim(:, k) = v;
    end

    % Fill out columns of Wr in (5).
    for k = 1:r
        tmpSum = zeros(2*n, 1); 
        for i = 1:r
            % Grab columns of V, form the sum:
            %   2*(qoRes(k, 1)*M*V(:, i) + ... 
            %       + qoRes(k, r)*M*V(:, r))
            % Add scalar*vector products, post-multiply by 2*Qfo.
            tmpSum = tmpSum + (qoRes(k, i) * VPrim(:, i)); 
        end
        tmpSum = 2*Qfo*tmpSum; 
        fprintf(1, 'COMPUTING COLUMN k = %d of W\n', k)
        fprintf(1, '---------------------------------------\n');
        % Option 1 in 'so_structured_solve.m' implements 
        % ((mirroredPoles(k)*Efo - Afo)'\tmpSum).
        rhs         = tmpSum(1:n, :);
        w           = so_structured_solve(Mso, Dso, Kso, rhs, mirroredPoles(k), ...
                                            1, 1);
        WPrim(:, k) = w;
    end
    
    % Step to ensure projection matrices are real-valued.
    % (It is assumed that poles, and tangents, are sorted into complex
    % conjugate pairs from the previous iteration.)
    V = zeros(2*n, r);    W = zeros(2*n, r);
    k = 1;
    while k <= r
        if (k < r && mirroredPoles(k) == conj(mirroredPoles(k + 1)))
            V(:, k)     = real(VPrim(:, k));
            V(:, k + 1) = imag(VPrim(:, k));
            W(:, k)     = real(WPrim(:, k));
            W(:, k + 1) = imag(WPrim(:, k));
            k           = k + 2;
        else
            V(:, k) = real(VPrim(:, k));
            W(:, k) = real(WPrim(:, k));
            k       = k + 1;
        end
    end
    
    % Orthonormalize.
    [V, ~] = qr(V, "econ");     [W, ~] = qr(W, "econ");

    % Projection.
    Wt = W';  Vt = V';
    fprintf(1, 'COMPUTING REDUCED-ORDER MODEL VIA PROJECTION.\n')
    fprintf(1, '----------------------------------------------\n');
    Er = Wt*Efo*V;   Ar = Wt*Afo*V;   Br = Wt*Bfo;
    Qr = Vt*Qfo*V;

    % Save poles and residues from last iteration for checks.
    mirroredPolesPrev = mirroredPoles; 
    qoResPrev         = qoRes;
    
    % Compute new poles and residues.
    [Xr, Lr] = eig(Ar, Er);     poles = diag(Lr);

    % Look for unstable poles.
    mirroredPoles = zeros(r, 1);
    for i = 1:r
        if real(poles(i)) < 0
            mirroredPoles(i) = -poles(i);
        else
            mirroredPoles(i) = poles(i);
        end
    end
    try % Attempt to sort poles into complex conjugate pairs.
        sortedMirroredPoles = cplxpair(mirroredPoles); 
        [~, sortIndex]      = ismember(sortedMirroredPoles, mirroredPoles);
        % Sort eigenvectors accordingly.
        Xr            = Xr(:, sortIndex);
        mirroredPoles = sortedMirroredPoles;
    catch % If model has become overtly complex, throw.
        warning('WARNING! Reduced model has become overtly complex; cannot sort eigenvalues into complex conjugate pairs. Returning unsorted poles.')
    end
    poleHistory(:, iterate + 1) = mirroredPoles;

    try % Attempt to sort poles into complex conjugate pairs
        poles = cplxpair(poles, 10e-3); 
    catch % If model has become overtly complex, throw 
        warning('WARNING! Reduced model has become overtly complex; cannot sort eigenvalues into complex conjugate pairs. Returning unsorted poles')
    end
    % Track poles.
    poleHistory(:, iterate + 1) = poles;

    % New residues.
    XrInv = Xr\eye(r, r); 
    modBr = Er\Br;
    qoRes = ((XrInv*(modBr))'.*((Xr' * Qr * Xr))).*(XrInv*(modBr));

    % Timings.
    fprintf(1, 'CURRENT ITERATE FINISHED IN %.2f s\n', toc(thisIterStart))
    fprintf(1, 'END OF CURRENT ITERATE k = %d\n', iterate)
    fprintf(1, '--------------------------------------------------\n');

    % Convergence tracking.
    iterate             = iterate + 1;    
    poleChange(iterate) = max(abs((mirroredPoles - mirroredPolesPrev)./mirroredPoles));
    fprintf('CHANGE IN POLES FROM PREVIOUS MODEL ITERATE IS %.12f \n', poleChange(iterate))
    fprintf(1, '--------------------------------------------------\n');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TERMINATION.                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return parameters from final iterate.
info                = struct();
info.poleHistory    = poleHistory;    
info.qoRes          = qoResPrev;
info.poleChange     = poleChange;

if iterate == (opts.maxIter + 1)
    fprintf('ALGORITHM HAS TERMINATED DUE TO REACHING MAX NO. OF ITERATIONS; TOTAL TIME ELAPSED IS %.2f s\n', toc(overallStart))
     fprintf(1, '--------------------------------------------------\n');
else
    fprintf('ALGORITHM HAS CONVERGED IN %d ITERATIONS.\n', iterate)
    fprintf('TOTAL TIME ELAPSED IS %.2f s\n', toc(overallStart))
    fprintf(1, '--------------------------------------------------\n');
end

% Plot convergence of poles.
if opts.plotConv == true
    semilogy(1:iterate, poleChange, '-o', LineWidth=1.5)
    xlim([1,iterate])
    xlabel('Iteration count')
    ylabel('Magnitude of change in \lambda')
end
end
