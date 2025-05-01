function [Er, Ar, Br, cr, Mr, info] = miso_lqoirka(E, A, B, c, M, r, opts)
% MISO_LQOIRKA Iterative Rational Krylov Algorithm (IRKA) for model-order 
% reduction of linear quadratic-output (LQO) systems with multiple inputs, 
% and a single output.
%
% SYNTAX:
%   [Er, Ar, Br, cr, Mr, info] = miso_lqoirka(E, A, B, c, M, r, opts)
%   [Er, Ar, Br, cr, Mr]       = miso_lqoirka(E, A, B, c, M, r, opts)
%
% DESCRIPTION:
%   Computes a LQO reduced model (Ar, Br, cr, Mr) using the LQO-IRKA 
%   algorithm presented in: 
%   "$\mathcal{H}_2$-optimal model reduction of linear
%   quadratic-output systems by multivariate rational interpolation"
%   At each iteration, an interpolatory reduced model is obtained using the
%   model reduction bases computed according to
%
%       V(:, k) = ((mirroredPoles(k))*E - A)\(B*rightTangents(:, k))); (1)
%       W(:, k) = ((mirroredPoles(k))*E.' - A.')\(2*(qoLeftTangents(k, 1)*M*V(:, i)) + ... 
%                   + (qoLeftTangents(k, r)*M*V(:, r)) + ...
%                       c*loLeftTangents(k)))                           (2) 
%
%   `mirroredPoles',  are the mirror images of the poles of the previous 
%   model iterate.
%   V and W are constructred to ensure realness of the reduced-order
%   model, and then orthonormalized.
%   It is assumed that the eigenvalues of (s*E - A) lie in the open left
%   half-plane, and that E is nonsingular.
%
% INPUTS:
%   E    - descriptor matrix with dimensions n x n in (1)
%   A    - state matrix with dimensions n x n in (1)
%   B    - input matrix with dimensions n x m in (1)
%   c    - linear output matrix with dimensions n x 1 in (2)
%   M    - (symmetric) quadratic output matrix with dimensions n x n in (2)
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
%   |                 | (default 0)                                       |
%   +-----------------+---------------------------------------------------+
%   | checkInterp     | bool, do we check interpolation conditions using  |
%   |                 | parameters from the previous iterate?             | 
%   |                 | ** not recommend for large problems **            |
%   |                 | (default 0)                                       |
%   +-----------------+---------------------------------------------------+
%   | computeH2       | bool, do we compute the H2 error of the current   |   
%   |                 | model iterate?                                    |
%   |                 | ** not recommend for large problems **            |
%   |                 | (default 0)                                       |
%   +-----------------+---------------------------------------------------+
%   | fomH2norm       | H2 norm of the full-order model                   |
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%   | Q               | quadratic-output observability Gramian            |
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%   | poles           | initial pole selection, r x 1 matrix              |
%   |                 | (default -(logspace(0, 4, r)')                    |
%   +-----------------+---------------------------------------------------+
%   | rightTangents   | initial m x r matrix of right tangent directions  |
%   |                 | (default eye(m, r))                               |
%   +-----------------+---------------------------------------------------+
%   | loLeftTangents  | initial linear-output left tangent directions     |
%   |                 | (scalars) r x 1 matrix                            |
%   |                 | (default eye(r, 1))                               |
%   +-----------------+---------------------------------------------------+
%   | qoLeftTangents  | initial quadratic-output left tangent directions  |
%   |                 | (scalars) r x r matrix                            |
%   |                 | (default eye(r, r))                               |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Er   - reduced descriptor matrix with dimensions r x r in (1)
%   Ar   - reduced state matrix with dimensions r x r in (1)
%   Br   - reduced input matrix with dimensions r x m in (1)
%   cr   - reduced linear output matrix with dimensions r x 1 in (2)
%   Mr   - reduced (symmetric) quadratic output matrix with dimensions 
%          r x r in (2) 
%   info - structure, containing the following information for monitoring
%          convergence
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | poleHistory     | convergence history of all reduced-order poles    |
%   +-----------------+---------------------------------------------------+
%   | rightTangents   | converged right tangent directions                |
%   +-----------------+---------------------------------------------------+
%   | loLeftTangents  | converged linear-output left tangent directions   |
%   +-----------------+---------------------------------------------------+
%   | qoLeftTangents  | converged quadratic-ouput left tangent directions |
%   +-----------------+---------------------------------------------------+
%   | poleChange      | convergence history of max error between poles    |
%   +-----------------+---------------------------------------------------+
%   | H2errors        | convergence history of H2 errors                  |
%   +-----------------+---------------------------------------------------+
%

%
% Copyright (c) 2025 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 4/30/2025
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem dimensions.
n          = size(A, 1);
m          = size(B, 2);
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
if ~isfield(opts, 'rightTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.rightTangents\n')
    fprintf(1, '------------------------------------------------\n')
    opts.rightTangents = eye(m, r);
end
if ~isfield(opts, 'loLeftTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.loLeftTangents\n')
    fprintf(1, '------------------------------------------------\n')
    opts.loLeftTangents = eye(r, 1);
end
if ~isfield(opts, 'qoLeftTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.qoLeftTangents\n')
    fprintf(1, '------------------------------------------------\n')                  
    opts.qoLeftTangents = eye(r, r);   
end
if~isfield(opts, 'plotConv')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.plotConv\n')
    fprintf(1, '------------------------------------------------\n')
    opts.plotConv = false;
end
if~isfield(opts, 'checkInterp')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.checkInterp\n')
    fprintf(1, '------------------------------------------------\n')
    opts.checkInterp = false;
end
if~isfield(opts, 'computeH2')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.computeH2\n')
    fprintf(1, '------------------------------------------------\n')
    opts.computeH2 = false;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overallStart = tic;
fprintf(1, 'INITIALIZING ALGORITHM.\n')
fprintf(1, '--------------------------------------------------\n');
% Assumed to be pre-sorted into conjugate pairs.
mirroredPoles  = -opts.poles;
rightTangents  = opts.rightTangents;
loLeftTangents = opts.loLeftTangents;
qoLeftTangents = opts.qoLeftTangents;

% If tracking H2 errors, pre-compute H2 norm of the full-order model.
if opts.computeH2
    if~isfield(opts, 'fomH2norm')
        normStart = tic;
        fprintf(1, 'PRE-COMPUTING H2 NORM OF FULL-ORDER MODEL.\n');
        fprintf(1, '--------------------------------------------------\n');
        P   = lyap(A, B*B');
        rhs = c*c' + M*P*M; 
        Q   = lyap(A', rhs);
        
        % Squared H2 norm of full-order model.
        fomH2norm = sqrt(abs(trace(B'*Q*B))); 
        fprintf(1, 'H2 NORM COMPUTED IN %d s\n', toc(normStart));
        fprintf(1, '--------------------------------------------------\n');
    else
        fprintf(1, 'H2 NORM OF FULL-ORDER MODEL PASSED AS AN ARGUMENT.\n');
        fprintf(1, '--------------------------------------------------\n');
        fomH2norm = opts.fomH2norm;
        Q         = opts.Q;
    end
end

% Set iteration count and tolerance to engage `while'.
% If not computing H2 errors, below is all that is returned.
iterate             = 1;   
poleChange(iterate) = opts.tol + 1; 
H2errors(iterate)   = opts.tol + 1; 

while (poleChange(iterate) > opts.tol && iterate <= opts.maxIter)
    thisIterStart = tic; 
    fprintf(1, 'CURRENT ITERATE IS k = %d\n', iterate)
    fprintf(1, '--------------------------------------------------\n');

    % Space allocation for model reduction bases.
    VPrim = zeros(n, r);     WPrim = zeros(n, r);

    % Fill out columns in V in (1); pre-computed for constructing W.
    for k = 1:r
        % fprintf(1, 'COMPUTING COLUMN k = %d of V\n', k)
        % fprintf(1, '--------------------------------------------------\n');
        v           = ((mirroredPoles(k)*E - A)\(B*rightTangents(:, k)));
        VPrim(:, k) = v;
    end

    % Fill out columns of W in (2).
    for k = 1:r
        % fprintf(1, 'COMPUTING COMLUMN k = %d of W\n', k)
        % fprintf(1, '--------------------------------------------------\n');
        tmpSum = zeros(n, 1); 
        for i = 1:r
            % Grab columns of V, form the sum:
            %   2*(soLeftTangents(k, 1)*M*V(:, i) + ... 
            %       + soLeftTangents(k, r)*M*V(:, r))
            % Add scalar*vector products, post-multiply by 2*M.
            tmpSum = tmpSum + (qoLeftTangents(k, i)*VPrim(:, i)); 
        end
        tmpSum      = 2*M*tmpSum;
        w           = (mirroredPoles(k)*E.' - A.')\(c*loLeftTangents(k, :) + tmpSum);
        WPrim(:, k) = w;
    end

    % Step to ensure projection matrices are real-valued.
    % (It is assumed that poles, and tangents, are sorted into complex
    % conjugate pairs from the previous iteration.)
    V = zeros(n, r);    W = zeros(n, r);
    k = 1;
    while k <= r
        if (k < r && mirroredPoles(k) == conj(mirroredPoles(k + 1)))
            V(:, k)     = real(VPrim(:, k));
            V(:, k + 1) = imag(VPrim(:, k));
            W(:, k)     = real(WPrim(:, k));
            W(:, k + 1) = imag(WPrim(:, k));
            k           = k + 2;
        else
            % Take real parts to ignore away any complex noise.
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
    Er = W.'*E*V;   Ar = W.'*A*V;    Br = W.'*B;
    cr = V.'*c;     Mr = V.'*M*V;

    % Save poles and residues from last iteration for checks.
    mirroredPolesPrev  = mirroredPoles;  
    rightTangentsPrev  = rightTangents;
    loLeftTangentsPrev = loLeftTangents;
    qoLeftTangentsPrev = qoLeftTangents;

    % Interpolation checks.
    if opts.checkInterp
        % Pre-compute.
        VrPrim = zeros(r, r);
        MVec   = reshape(M,  [n^2, 1]);
        MrVec  = reshape(Mr, [r^2, 1]);
        for k = 1 :r
            VrPrim(:, k) = ((mirroredPolesPrev(k)*Er - Ar)\(Br*rightTangentsPrev(:, k)));
        end

        fprintf(1, '1. LINEAR, RIGHT TANGENTIAL CONDITIONS.\n')
        fprintf(1, '--------------------------------------------------\n');
        for k = 1:r
            fprintf(1, 'RELATIVE ERROR IN CONDITION %d: %d\n', k, ...
                norm(c.'*VPrim(:, k) - cr.'*VrPrim(:, k), 2)/norm(c.'*VPrim(:, k), 2))
            fprintf(1, '--------------------------------------------------\n');
        end

        fprintf(1, '2. QUADRATIC, RIGHT TANGENTIAL CONDITIONS.\n')
        fprintf(1, '--------------------------------------------------\n');
        for i = 1:r
            for j = 1:r
                fprintf(1, 'RELATIVE ERROR IN CONDITION (%d, %d): %d\n', i, j, ...
                    norm(MVec.'*kron(VPrim(:, i), VPrim(:, j)) -  MrVec.'*kron(VrPrim(:, i), VrPrim(:, j)), 2)/norm(MVec.'*kron(VPrim(:, i), VPrim(:, j)), 2))
                fprintf(1, '--------------------------------------------------\n');
            end
        end

        fprintf(1, '3. MIXED, HERMITE TANGENTIAL CONDITIONS.\n')
        fprintf(1, '--------------------------------------------------\n');
        for k = 1:r
            roQuadTerm = 0;
            foQuadTerm = 0;
            for i = 1:r
                foQuadTerm = foQuadTerm + qoLeftTangentsPrev(k, i)*MVec.'*kron(-((mirroredPolesPrev(k)*E - A)\(E*VPrim(:, k))), VPrim(:, i)) + ...
                    qoLeftTangentsPrev(i, k)*MVec.'*kron(VPrim(:, i), -((mirroredPolesPrev(k)*E - A)\(E*VPrim(:, k))));

                roQuadTerm = roQuadTerm + qoLeftTangentsPrev(k, i)*MrVec.'*kron(-((mirroredPolesPrev(k)*Er - Ar)\(Er*VrPrim(:, k))), VrPrim(:, i)) + ...
                    qoLeftTangentsPrev(i, k)*MrVec.'*kron(VrPrim(:, i), -((mirroredPolesPrev(k)*Er - Ar)\(Er*VrPrim(:, k))));
            end
            foLinTerm = -loLeftTangentsPrev(k)*c.'*((mirroredPolesPrev(k)*E - A)\(E*VPrim(:, k)));
            roLinTerm = -loLeftTangentsPrev(k)*cr.'*((mirroredPolesPrev(k)*Er - Ar)\(Er*VrPrim(:, k)));

            fprintf(1, 'RELATIVE ERROR IN CONDITION %d: %d\n', k, ...
                norm(foLinTerm + foQuadTerm - (roLinTerm + roQuadTerm), 2)/norm(foLinTerm + foQuadTerm, 2))
            fprintf(1, '--------------------------------------------------\n');
        end
    end
    
    % Compute new reduced-order poles and residues.
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
        sortedMirroredPoles = cplxpair(mirroredPoles, 1e-6); 
        [~, sortIndex]      = ismember(sortedMirroredPoles, mirroredPoles);
        % Sort eigenvectors accordingly.
        Xr            = Xr(:, sortIndex);
        mirroredPoles = sortedMirroredPoles;
    catch % If model has become overtly complex, throw.
        warning('WARNING! Reduced model has become overtly complex; cannot sort eigenvalues into complex conjugate pairs. Returning unsorted poles.')
    end
    poleHistory(:, iterate + 1) = poles;
    
    % New residues.
    XrInverse      = Xr\eye(r, r); 
    rightTangents  = XrInverse*(Er\Br);
    rightTangents  = rightTangents.';
    loLeftTangents = Xr.'*cr;  
    qoLeftTangents = Xr.'*Mr*Xr; 

    % Timings.
    fprintf(1, 'CURRENT ITERATE FINISHED IN %.2f s\n', toc(thisIterStart))
    fprintf(1, 'END OF CURRENT ITERATE k = %d\n', iterate)
    fprintf(1, '--------------------------------------------------\n');

    % Convergence tracking.
    iterate             = iterate + 1;    
    poleChange(iterate) = max(abs((mirroredPoles - mirroredPolesPrev)./mirroredPoles));
    fprintf('CHANGE IN POLES FROM PREVIOUS MODEL ITERATE IS %.12f \n', poleChange(iterate))
    fprintf(1, '--------------------------------------------------\n');

    % Compute current H2 error, if desired.
    if opts.computeH2
        errorStart = tic;
        fprintf(1, 'COMPUTING H2 ERROR OF CURRENT MODEL ITERATE.\n');
        fprintf(1, '--------------------------------------------------\n');
        % Reduced model error; use Gramian-based formula.
        Berr = [B; Er\Br];
        
        % Off-diagonal, (1,2) and (2,1), blocks of error Gramians.
        % Xr = lyap(A,  (Er\Ar).', B*(Er\Br).');
        Xr = mess_sylvester_sparse_dense(A, 'N', (Er\Ar), 'T', B*(Er\Br).', speye(n, n), eye(r, r));
        % Zr = lyap(A', (Er\Ar),  -c*cr.' - M*Xr*Mr);
        Zr = mess_sylvester_sparse_dense(A, 'T', (Er\Ar), 'N', (-c*cr.' - M*Xr*Mr), speye(n, n), eye(r, r));
    
        % Reduced Gramians; (2,2) block of error Gramian.
        Pr = lyap(Er\Ar,    (Er\Br)*(Er\Br)');                  
        Qr = lyap((Er\Ar)', cr*cr' + Mr*Pr*Mr);

        % Error Gramian.
        % Qerr = [Q, Zr; Zr.', Qr];

        % Save relative H2 error of current model iterate.
        % H2errors(:, iterate) = sqrt(abs(trace(Berr'*Qerr*Berr)))/fomH2norm;
        H2errors(:, iterate) = sqrt(abs(fomH2norm^2 + trace((Er\Br)'*Qr*(Er\Br)) + ...
            2*trace(B'*Zr*(Er\Br))))/fomH2norm;

        fprintf(1, 'H2 ERROR OF CURRENT MODEL ITERATE IS %.8f\n', H2errors(:, iterate))
        fprintf(1, '--------------------------------------------------\n');
        

        fprintf(1, 'H2 ERROR COMPUTED IN %d s\n', toc(errorStart));
        fprintf(1, '--------------------------------------------------\n');
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TERMINATION.                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Return parameters from final iterate.
info                = struct();
info.poleHistory    = poleHistory;    
info.rightTangents  = rightTangentsPrev;
info.foLeftTangents = loLeftTangentsPrev;
info.qoLeftTangents = qoLeftTangentsPrev;
info.poleChange     = poleChange;
info.H2errors       = H2errors;

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