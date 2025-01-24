function [Er, Ar, Br, cr, Mr] = miso_lqointerp(E, A, B, c, M, r, opts)
% MISO_LQOINTERP Driver function for model-order reduction of linear 
% quadratic-output systems with multiple inputs, and a single output using
% methods of tangential interpolation.
%
% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Er, Ar, Br, cr, Mr)
%   using the methods of tangential interpolation presented in.
%   "..."
%   For the given data, an interpolatory reduced model is computed using
%   one of the following strategies.
%   If opts.interpolation == 'standard', Galerkin projection is performed
%   with V = W, where V is computed according to
%
%       V(:, k) = ((shifts(k))*E - A)\(B*rightTangents(:, k))); (1)
%
%   If opts.interpolation == 'mixed', Petrov-Galerkin projection is
%   performed, where V is computed according to (1) and W is computed
%   according to
%
%       W(:, k) = ((shifts(k))*E.' - A.')\(2*(soRes(k, 1)*M*V(:, i)) + ... 
%                   + (soLeftTangents(k, r)*M*V(:, r)) + ...
%                       c*foLeftTangents(k)))                           (2) 
%
%   The matrices are constructred to ensure realness of the reduced-order
%   model, and then orthonormalized.
%   It is assumed that the eigenvalues of (s*E-A) lie in the open left
%   half-plane, and all interpolation data are sorted into complex
%   conjugate pairs.
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
%   | checkInterp     | bool, do we check interpolation conditions using  |
%   |                 | the given interpolation data?                     | 
%   |                 | ** not recommend for large problems **            |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | shifts          | interpolation shifts, r x 1 matrix                |
%   |                 | (default -(logspace(0, 4, r)')                    |
%   +-----------------+---------------------------------------------------+
%   | rightTangents   | m x r matrix of right tangent directions          |
%   |                 | (default 10*rand(m, r))                           |
%   +-----------------+---------------------------------------------------+
%   | loLeftTangents  | linear-output left tangent directions             |
%   |                 | (scalars) r x 1 matrix                            |
%   |                 | (default 10*rand(r, 1))                           |
%   +-----------------+---------------------------------------------------+
%   | qoLeftTangents  | quadratic-output left tangent directions          |
%   |                 | (scalars) r x r matrix                            |
%   |                 | (default tmp          = 10*rand(r, r);            |
%   |                 |        qoLeftTangents = (tmp + tmp')/2)           |
%   +-----------------+---------------------------------------------------+
%   | interpolation   | type of tangential interpolation to perform       |
%   |                 | (default 'standard')                              |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Er   - reduced descriptor matrix with dimensions r x r in (1)
%   Ar   - reduced state matrix with dimensions r x r in (1)
%   Br   - reduced input matrix with dimensions r x m in (1)
%   cr   - reduced linear output matrix with dimensions r x 1 in (2)
%   Mr   - reduced (symmetric) quadratic output matrix with dimensions 
%          r x r in (2) 
%

%
% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 1/20/2025
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem dimensions.
n = size(A, 1);
m = size(B, 2);

% Check and set input opts.
if (nargin < 7) 
    opts = struct(); % Empty struct.
end

if ~isfield(opts, 'shifts')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.shifts\n')
    fprintf(1, '------------------------------------------------\n')
    opts.shifts = -(logspace(0, 4, r)'); 
end
if ~isfield(opts, 'rightTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.rightTangents\n')
    fprintf(1, '------------------------------------------------\n')
    opts.rightTangents = 10*rand(m, r);
end
if ~isfield(opts, 'loLeftTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.loLeftTangents\n')
    fprintf(1, '------------------------------------------------\n')
    opts.loLeftTangents = 10*rand(r, 1);
end
if ~isfield(opts, 'qoLeftTangents')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.qoLeftTangents\n')
    fprintf(1, '------------------------------------------------\n')
    tmp                 = 10*rand(r, r);                     
    opts.qoLeftTangents = (tmp + tmp')/2;
end
if~isfield(opts, 'checkInterp')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.checkInterp\n')
    fprintf(1, '------------------------------------------------\n')
    opts.checkInterp = true;
end
if~isfield(opts, 'interpolation')
    fprintf(1, 'SETTING DEFAULT VALUE FOR opts.interpolation\n')
    fprintf(1, '------------------------------------------------\n')
    opts.interpolation = 'standard';
end

assert(strcmp(opts.interpolation, 'standard') || strcmp(opts.interpolation, 'mixed'), ...
    'Interpolation strategy not implemented!')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REDUCTION.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1, 'GETTING INTERPOLATION DATA.\n')
fprintf(1, '--------------------------------------------------\n');
shifts         = opts.shifts;
rightTangents  = opts.rightTangents;
loLeftTangents = opts.loLeftTangents;
qoLeftTangents = opts.qoLeftTangents;

% Space allocation for model reduction bases.
VPrim = zeros(n, r);     WPrim = zeros(n, r);
for k = 1:r
    % fprintf(1, 'COMPUTING COLUMN k = %d of V\n', k)
    % fprintf(1, '--------------------------------------------------\n');
    v           = ((shifts(k)*E - A)\(B*rightTangents(:, k)));
    VPrim(:, k) = v;
end

if strcmp(opts.interpolation, 'standard')
    % Compute reduced model using standard (no mixed) tangential
    % interpolation and Galerkin projection with V = W. 
    WPrim = VPrim; 
else
    % Compute reduced model using mixed tangential interpolation and
    % Petrov-Galerkin projection.
    % Fill out columns of W in (2).
    for k = 1:r
        % fprintf(1, 'COMPUTING COMLUMN k = %d of W\n', k)
        % fprintf(1, '--------------------------------------------------\n');
        tmpSum = zeros(n, 1); 
        for i = 1:r
            % Grab columns of V, form the sum:
            %   2*(qoLeftTangents(k, 1)*M*V(:, i) + ... 
            %       + qoLeftTangents(k, r)*M*V(:, r))
            % Add scalar*vector products, post-multiply by 2*M.
            tmpSum = tmpSum + (qoLeftTangents(k, i)*VPrim(:, i)); 
        end
        tmpSum      = 2*M*tmpSum;
        w           = (shifts(k)*E.' - A.')\(c*loLeftTangents(k, :) + tmpSum);
        WPrim(:, k) = w;
    end
end

% Step to ensure projection matrices are real-valued.
% (It is assumed that shifts, and tangents, are sorted into complex
% conjugate pairs.)
V = zeros(n, r);    W = zeros(n, r);
k = 1;
while k <= r
    if (k < r && shifts(k) == conj(shifts(k + 1)))
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
[V, ~] = qr(VPrim, "econ");     [W, ~] = qr(WPrim, "econ");

fprintf(1, 'COMPUTING REDUCED-ORDER MODEL VIA PROJECTION.\n')
fprintf(1, '--------------------------------------------------\n');
Er = W.'*E*V;   Ar = W.'*A*V;   Br = W.'*B;    
cr = V.'*c;     Mr = V.'*M*V;

% Interpolation checks.
if opts.checkInterp
    % Pre-compute.
    VrPrim = zeros(r, r);
    MVec   = reshape(M,  [n^2, 1]);
    MrVec  = reshape(Mr, [r^2, 1]);
    for k = 1 :r
        VrPrim(:, k) = ((shifts(k)*Er - Ar)\(Br*rightTangents(:, k)));
    end

    fprintf(1, '1. LINEAR, RIGHT TANGENTIAL CONDITIONS.\n')
    fprintf(1, '--------------------------------------------------\n');
    for k = 1:r
        fprintf(1, 'ABSOLUTE ERROR IN CONDITION %d: %d\n', k, ...
            norm(c.'*VPrim(:, k) - cr.'*VrPrim(:, k), 2))
        fprintf(1, '--------------------------------------------------\n');
    end

    fprintf(1, '2. QUADRATIC, RIGHT TANGENTIAL CONDITIONS.\n')
    fprintf(1, '--------------------------------------------------\n');
    for i = 1:r
        for j = 1:r
            fprintf(1, 'ABSOLUTE ERROR IN CONDITION (%d, %d): %d\n', i, j, ...
                norm(MVec.'*kron(VPrim(:, i), VPrim(:, j)) -  MrVec.'*kron(VrPrim(:, i), VrPrim(:, j)), 2))
            fprintf(1, '--------------------------------------------------\n');
        end
    end

    % Only check mixed conditions if doing mixed projection.
    if strcmp(opts.interpolation, 'mixed')
        fprintf(1, '3. MIXED, HERMITE TANGENTIAL CONDITIONS.\n')
        fprintf(1, '--------------------------------------------------\n');
        for k = 1:r
            roQuadTerm = 0;
            foQuadTerm = 0;
            for i = 1:r
                foQuadTerm = foQuadTerm + qoLeftTangents(k, i)*MVec.'*kron(-((shifts(k)*E - A)\(E*VPrim(:, k))), VPrim(:, i)) + ...
                    qoLeftTangents(i, k)*MVec.'*kron(VPrim(:, i), -((shifts(k)*E - A)\VPrim(:, k)));
    
                roQuadTerm = roQuadTerm + qoLeftTangents(k, i)*MrVec.'*kron(-((shifts(k)*Er - Ar)\VrPrim(:, k)), VrPrim(:, i)) + ...
                    qoLeftTangents(i, k)*MrVec.'*kron(VrPrim(:, i), -((shifts(k)*Er - Ar)\(Er*VrPrim(:, k))));
            end
            foLinTerm = -loLeftTangents(k)*c.'*((shifts(k)*E - A)\(E*VPrim(:, k)));
            roLinTerm = -loLeftTangents(k)*cr.'*((shifts(k)*Er - Ar)\(Er*VrPrim(:, k)));
    
            fprintf(1, 'ABSOLUTE ERROR IN CONDITION %d: %d\n', k, ...
                norm(foLinTerm + foQuadTerm - (roLinTerm + roQuadTerm), 2))
            fprintf(1, '--------------------------------------------------\n');
        end
    end
end

end