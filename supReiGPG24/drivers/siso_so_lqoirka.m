function [Efo_r, Afo_r, Bfo_r, Qfo_r, info] = siso_so_lqoirka(Mso, Dso, Kso, Bso, Qfo, r, opts)
% MISO_LQOIRKA Iterative Rational Krylov Algorithm (IRKA) for model-order 
% reduction of linear quadratic-output systems with a single input and a
% single output. Tailored for first-order linear systems obtained by
% lifting a second-order dynamical system, with a single quadratic output.
%
% SYNTAX:
%   [Efo_r, Afo_r, Bfo_r, Qfo_r, info] = siso_lqoirka(Mso, Dso, Kso, Bso, Qfo, r, opts)
%   [Efo_r, Afo_r, Bfo_r, Qfo_r]       = siso_lqoirka(Mso, Dso, Kso, Bso, Qfo, r)
%
% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Ar, Br, Cr, Qr)
%   using the Iterative Rational Krylov Algorithm presented in 
%   "..."
%   The full-order system matrices have the nfo = 2*n order first-order 
%   realization obtained from lifting a second-order system
%
%       Efo = [I  0;   0     Mso];   (0a)
%       Afo = [I  0;  -Kso  -Dso];   (0b)
%       Bfo = [0; Bso];              (0c)
%       Qfo = [Cpso'*Cpso, 0; 0, 0]; (0d)
%
%   At each iteration, an interpolatory reduced model is obtained using the
%   model reduction bases computed according to
%
%       V(:, k) = ((mirroredPoles(k))*Efo - Afo)\Bfo); (1)
%       W(:, k) = ((mirroredPoles(k))*Efo.' - Afo.')\(2*(qoLeftTangents(k, 1)*Qfo*V(:, i)) + ... 
%                   + (qoLeftTangents(k, r)*Qfo*V(:, r)));                   (2) 
%
%   `mirroredPoles',  are the mirror images of the poles of the previous 
%   model iterate. V and W are constructred to ensure realness of the 
%   reduced-order model, and then orthonormalized.
%
%   The linear solves are not computed naively as in (1), (2). Instead, 
%   they are computed using the specialized companion function 
%   'so_structured_solve.m', which leverages the underlying second-order
%   structure of the first-order system matrices in (0a) - (0d). 
%
%   It is assumed that the eigenvalues of (s*Efo - Afo) lie in the open left
%   half-plane, and that E is nonsingular.
%
% INPUTS:
%   Mso  - sparse second order mass matrix with dimensions n x n in 
%          (0a)
%   Dso  - sparse second order damping matrix with dimensions n x n 
%          in (0b)
%   Kso  - sparse second order stiffness matrix with dimensions n x n 
%          in (0b)
%   Bso  - sparse second order input matrix with dimensions n x 1 in 
%          (0c)
%   Qfo  - sparse symmetric quadratic output matrix with dimensions 
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
%   | poles           | initial pole selection, r x 1 matrix              |
%   |                 | (default -(logspace(0, 4, r)')                    |
%   +-----------------+---------------------------------------------------+
%   | qoLeftTangents  | initial quadratic-output left tangent directions  |
%   |                 | (scalars) r x r matrix                            |
%   |                 | (default tmp          = eye(r, r);                |
%   |                 |        soLeftTangents = (tmp + tmp')/2)           |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Efo_r   - reduced descriptor matrix with dimensions r x r
%   Afo_r   - reduced state matrix with dimensions r x r
%   Bfo_r   - reduced input matrix with dimensions r x 1
%   Qfo_r   - reduced (symmetric) quadratic output matrix with dimensions 
%             r x r
%   info    - structure, containing the following information for monitoring
%             convergence
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | poleHistory     | Convergence history of all reduced-order poles    |
%   +-----------------+---------------------------------------------------+
%   | qoLeftTangents  | converged quadratic-ouput left tangent directions |
%   +-----------------+---------------------------------------------------+
%   | poleChange      | Convergence history of max error between poles    |
%   +-----------------+---------------------------------------------------+
%

%
% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 1/24/2025
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem dimensions.
n          = size(Mso, 1);
poleChange = 1;

% Check and set input opts.
if (nargin < 7) 
    opts = struct(); % Empty struct.
end

if ~isfield(opts, 'tol')
    opts.tol = 10e-6;
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.tol\n')
    fprintf(1, '------------------------------------------------\n')
end
if ~isfield(opts, 'maxIter')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.maxIter\n')
    fprintf(1, '------------------------------------------------\n')
    opts.maxIter = 100;
end
if ~isfield(opts, 'poles')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.poles\n')
    fprintf(1, '------------------------------------------------\n')
    opts.poles = -(logspace(0, 4, r)'); 
end
if ~isfield(opts, 'qoLeftTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.qoLeftTangents\n')
    fprintf(1, '------------------------------------------------\n')
    tmp                 = eye(r, r);                     
    opts.qoLeftTangents = (tmp + tmp')/2;
end
if~isfield(opts, 'plotConv')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.plotConv\n')
    fprintf(1, '------------------------------------------------\n')
    opts.plotConv = true;
end
if~isfield(opts, 'checkInterp')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.checkInterp\n')
    fprintf(1, '------------------------------------------------\n')
    opts.checkInterp = true;
end

% Build first-order realization for later projection.
Efo                   = spalloc(2*n, 2*n, nnz(Mso) + n);            % Descriptor matrix; Efo = [I, 0: 0, Mso]
Efo(1:n, 1:n)         = speye(n);                                   % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = Mso;                                        % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(Kso) + nnz(Dso) + n); % Afo = [0, I; -Kso, -Dso]
Afo(1:n, n+1:2*n)     = speye(n);                                   % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -Kso;                                       % (2, 1) block is -Kso
Afo(n+1:2*n, n+1:2*n) = -Dso;                                       % (2, 2) block is -Dso 

Bfo             = spalloc(2*n, 1, nnz(Bso));                        % Bfo = [0; Bso];
Bfo(n+1:2*n, :) = Bso; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overallStart = tic;
fprintf(1, 'INITIALIZING ALGORITHM.\n')
fprintf(1, '--------------------------------------------------\n');
mirroredPoles  = -opts.poles;
qoLeftTangents = opts.qoLeftTangents;

% Set iteration count and tolerance to engage `while'.
iterate             = 1;   
poleChange(iterate) = opts.tol + 1; 

while (poleChange(iterate) > opts.tol && iterate <= opts.maxIter)
    thisIterStart = tic; 
    fprintf(1, 'CURRENT ITERATE IS k = %d\n', iterate)
    fprintf(1, '--------------------------------------------------\n');

    % Space allocation for model reduction bases.
    VPrim = zeros(2*n, r);     WPrim = zeros(2*n, r);

    % Fill out columns in V in (1); pre-computed for constructing W.
    for k = 1:r
        fprintf(1, 'COMPUTING COLUMN k = %d of V\n', k)
        fprintf(1, '--------------------------------------------------\n');
        % Option 0 implments v = ((mirroredPoles(k)*Efo - Afo)\Bfo);
        v           = so_structured_solve(Mso, Dso, Kso, Bso, mirroredPoles(k), ...
                        0, 1);
        VPrim(:, k) = v;
    end

    % Fill out columns of W in (2).
    for k = 1:r
        fprintf(1, 'COMPUTING COMLUMN k = %d of W\n', k)
        fprintf(1, '--------------------------------------------------\n');
        % Option 1 implements w = (mirroredPoles(k)*Efo.' - Afo.')\tmpSum);
        tmpSum = zeros(2*n, 1); 
        for i = 1:r
            % Grab columns of V, form the sum:
            %   2*(soLeftTangents(k, 1)*M*V(:, i) + ... 
            %       + soLeftTangents(k, r)*M*V(:, r))
            % Add scalar*vector products, post-multiply by 2*M.
            tmpSum = tmpSum + (qoLeftTangents(k, i)*VPrim(:, i)); 
        end
        tmpSum      = 2*Qfo*tmpSum; % Need -?
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
    
    % Orthonormalize projection matrices.
    [V, ~] = qr(V, "econ");     [W, ~] = qr(W, "econ");

    % Compute reduced-order model via projection.
    fprintf(1, 'COMPUTING REDUCED-ORDER MODEL VIA PROJECTION.\n')
    fprintf(1, '--------------------------------------------------\n');
    Efo_r = W.'*Efo*V;   Afo_r = W.'*Afo*V;    Bfo_r = W.'*Bfo;
    Qfo_r = V.'*Qfo*V;

    % Save poles and residues from last iteration for checks.
    mirroredPolesPrev  = mirroredPoles;  
    qoLeftTangentsPrev = qoLeftTangents;

    % Interpolation checks.
    if opts.checkInterp
        % Pre-compute.
        VrPrim   = zeros(r, r);
        QfoVec   = reshape(Qfo,  [(2*n)^2, 1]);
        Qfo_rVec = reshape(Qfo_r, [r^2, 1]);
        for k = 1 :r
            VrPrim(:, k) = ((mirroredPolesPrev(k)*Efo_r - Afo_r)\Bfo_r);
        end

        fprintf(1, '1. QUADRATIC, RIGHT TANGENTIAL CONDITIONS.\n')
        fprintf(1, '--------------------------------------------------\n');
        for i = 1:r
            for j = 1:r
                fprintf(1, 'ABSOLUTE ERROR IN CONDITION (%d, %d): %d\n', i, j, ...
                    norm(QfoVec.'*kron(VPrim(:, i), VPrim(:, j)) -  Qfo_rVec.'*kron(VrPrim(:, i), VrPrim(:, j)), 2)/norm(QfoVec.'*kron(VPrim(:, i), VPrim(:, j)), 2))
                fprintf(1, '--------------------------------------------------\n');
            end
        end

        fprintf(1, '2. MIXED, HERMITE TANGENTIAL CONDITIONS.\n')
        fprintf(1, '--------------------------------------------------\n');
        for k = 1:r
            roQuadTerm = 0;
            foQuadTerm = 0;
            for i = 1:r
                foQuadTerm = foQuadTerm + qoLeftTangentsPrev(k, i)*QfoVec.'*kron(-((mirroredPolesPrev(k)*Efo - Afo)\(Efo*VPrim(:, k))), VPrim(:, i)) + ...
                    qoLeftTangentsPrev(i, k)*QfoVec.'*kron(VPrim(:, i), -((mirroredPolesPrev(k)*Efo - Afo)\(Efo*VPrim(:, k))));

                roQuadTerm = roQuadTerm + qoLeftTangentsPrev(k, i)*Qfo_rVec.'*kron(-((mirroredPolesPrev(k)*Efo_r - Afo_r)\(Efo_r*VrPrim(:, k))), VrPrim(:, i)) + ...
                    qoLeftTangentsPrev(i, k)*Qfo_rVec.'*kron(VrPrim(:, i), -((mirroredPolesPrev(k)*Efo_r - Afo_r)\(Efo_r*VrPrim(:, k))));
            end

            fprintf(1, 'ABSOLUTE ERROR IN CONDITION %d: %d\n', k, ...
                norm(foQuadTerm - roQuadTerm, 2)/norm(foQuadTerm, 2))
            fprintf(1, '--------------------------------------------------\n');
        end
    end
    
    % Compute new reduced-order poles and residues.
    [Xr, Lr] = eig(Afo_r, Efo_r);     poles = diag(Lr);
    
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
    
    % New residues.
    % qoLeftTangents = Xr.'*Qfo_r*Xr; 
    XrInv          = Xr\eye(r, r); 
    modBr          = Efo_r\Bfo_r;
    qoLeftTangents = ((XrInv*(modBr)).'.*(Xr.'*Qfo_r*Xr)).*(XrInv*(modBr));

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
info.qoLeftTangents = qoLeftTangentsPrev;
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