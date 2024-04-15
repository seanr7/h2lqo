function [Wprim, Vprim, Worth, Vorth, H_shifts, pW, pV] = interpolatory_solves(...
    Mso, Dso, Kso, Qfo, bso, shifts, r, opts)
% INTERPOLATORY_SOLVES Driver function for computing interpolatory
% projections of frequency-domain linear quadratic output systems.
%
%
% DESCRIPTION:
%   Function to compute interpolatory projection matrices Vorth, Worth, for
%   use in the model reduction frequency domain linear quadratic output 
%   systems. The methodologies used are outlined in Section 5 of
%   "Interpolatory model order reduction of large-scale dynamical systems
%   with root mean squared error measures"
%   2n x q primitive interpolatiton bases are computed using linear solves 
%   according to the input opts.proj, or passed; q = length(shifts), q>=r.
%   Then, the primitive bases
%   are compressed according to the selection procedure opts.compress.
%
%   The left 2n x r interpolatory projection matrix is computed according 
%   to
%
%       Vprim(:, k) = (shifts(k)*Efo - Afo)\bfo; (1)
%
%   If performing Galerkin projection, take Wprim = Vprim.
%   If performing Petrov-Galerkin projection, right 2n x r interpolatory
%   projection matrix is computed according to
%
%       Wprim(:, k) = (shifts(k)*Efo - Afo)'\(Qfo*
%                       ((shifts(k)*Efo - Afo)\bfo); (2)
%
%   where
%
%       Efo = [I  0;   0     Mso]; (3)
%       Afo = [I  0;  -Kso  -Dso]; (4)
%       bfo = [0; bso];            (5)
%
%   The linear solves are not computed directly as in (1), (2). Instead, 
%   they are computed using the specialized companion function 
%   'so_structured_solve.m', which leverages the underlying second-order
%   structure of the first-order system matrices in (3) - (5). 
%   See 'help so_structured_solve' for details.
%
%   If opts.compress == Linfty, primitive basis vectors are selected
%   iteratively via a greedy search; new points are chosen where the
%   pointwise absolute error is the largest. Once selected, these 2n x r
%   matrices are orthonormalized via an economy-size QR to give Vorth and
%   Worth.
% 
%   If opts.compress == avg, primitive basis vectors are `averaged' using a
%   pivoted QR decomposition. The leading r <= q orthonormal columns are
%   chosen as Vrorth, Wrorth.
%   It is assumed that the eigenvalues of (s*Efo-Afo) lie in the open left
%   half-plane, and that the points shifts are placed along the imaginary
%   axis.
%
% INPUTS:
%   Mso    - sparse second order mass matrix with dimensions n x n in 
%            (3)
%   Dso    - sparse second order damping matrix with dimensions n x n 
%            in (4)
%   Kso    - sparse second order stiffness matrix with dimensions n x n 
%            in (4)
%   bso    - sparse second order input matrix with dimensions n x 1 in 
%            (5)
%   Qfo    - sparse symmetric quadratic output matrices with dimensions 
%            2*n x 2*n, such that Qfo(1:n, 1:n) is nonzero, Qfo is zero 
%            elsewhere
%   shifts - q >= r complex interpolation points along iR
%   r      - order of reduction
%   opts   - structure, containing the following optional entries:
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | proj            | option to type of projection to construct Vprim,  |
%   |                 | Wprim, according to                               |
%   |                 | 'g'  - Build Vprim, Wprim as in (1) according to  |
%   |                 |        Galerkin projection.                       |
%   |                 | 'pg' - Build Vprim, Wprim as in (2) according to  |
%   |                 |        Petrov-Galerkin projection.                |
%   |                 | (default 'g')                                     |
%   +-----------------+---------------------------------------------------+
%   | compress        | option to compress primitive bases Vprim, Wprim   |
%   |                 | 'avg'    - choose r >= q columns via pivoted QR   |
%   |                 | 'Linfty' - choose r >= q columns via greedy search|     
%   |                 | (default 'avg')                                   |
%   +-----------------+---------------------------------------------------+
%   | recomp_bases    | bool, do we recompute the matrices Vr, Wr, or are |
%   |                 | they passed as opts?                              |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | recomp_evals    | bool, do we recompute the transfer function       |
%   |                 | evaluations used in the Linfty search, or are they|
%   |                 | passed as opts?                                   |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | Vprim           | Left n x r interpolatory projection matrix        | 
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%   | Wprim           | Right n x r interpolatory projection matrix       | 
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%   | Hshifts         | Full-order ransfer function evaluations at given  |
%   |                 | shifts, as q x 1 array                            | 
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Wprim    - Primitive 2n x r left interpolatory projection matrix, 
%              constructed according to opts.proj
%   Vprim    - Primitive 2nn x r right interpolatory projection matrix, 
%              constructed according to opts.proj
%   Worth    - Orthonormalized 2n x r left interpolatory projection matrix, 
%              compressed according to opts.compress
%   Vorth    - Orthonormalized 2nn x r right interpolatory projection matrix, 
%              compressed according to opts.compress
%   pV       - Columns of Vprim selected by opts.compress
%   pV       - Columns of Wprim selected by opts.compress
%   H_shifts - Transfer function evaluations of full-order model (Efo, Afo, 
%              bfo, Qfo) computed at inputted shifts.

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model order reduction of large-scale 
% dynamical systems with root mean squared error measures"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

% Virginia Tech, USA
% Last editied: 4/15/2024
%%
% Grab dimensions.
q =  max(size(shifts));
[n, ~] = size(Mso);

% Check and set inputs
if (nargin < 7) 
    opts = struct(); % Empty struct 
end

if ~isfield(opts, 'proj')
    opts.proj = 'g'; 
end
if ~isfield(opts, 'compress')
    fprintf(1, 'Setting default value for opts.compression = avg\n')
    fprintf(1, '------------------------------------------------\n')
    opts.compress = 'avg'; 
end
if ~isfield(opts, 'recomp_bases')
    fprintf(1, 'Setting default value for opts.recomp_bases = 0\n')
    fprintf(1, '------------------------------------------------\n')
    opts.recomp_bases = 0; 
end
if ~isfield(opts, 'recomp_evals')
    fprintf(1, 'Setting default value for opts.recomp_evals = 0\n')
    fprintf(1, '------------------------------------------------\n')
    opts.recomp_evals = 0;
end
if ~isfield(opts, 'Vprim')
    fprintf(1, 'Setting default value for opts.Vprim = []\n')
    fprintf(1, '------------------------------------------------\n')
    opts.Vprim = [];
end
if ~isfield(opts, 'Wprim')
    fprintf(1, 'Setting default value for opts.Wprim = []\n')
    fprintf(1, '------------------------------------------------\n')
    opts.Wprim = [];
end
if ~isfield(opts, 'H_shifts')
    fprintf(1, 'Setting default value for opts.H_shifts = []\n')
    fprintf(1, '------------------------------------------------\n')
    opts.H_shifts = [];
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

%% Methods.
% If opts.recomp_bases, pre-compute Vprim and Wprim
if opts.recomp_bases 
    % Space allocation
    Vprim = zeros(2*n, q);   Wprim = zeros(2*n, q);
    % Galerkin projection
    if strcmp(opts.proj, 'g')
        % In Galerkin projection, primitive bases are identical in (1)
        linear_solves = tic;
        fprintf(1, 'Computing model reduction bases via Galerkin projection \n')
        fprintf(1, '--------------------------------------------------------\n')
        for k = 1:q
            fprintf(1, 'Current iterate is k = %d\n', k)
            fprintf(1, '-------------------------\n')
            % Option 0 in 'so_structured_solve.m' implements (-shifts(k) * E - A)\b)
            v           = so_structured_solve(Mso, Dso, Kso, bso, -shifts(k), 0, 1);
            Vprim(:, k) = v;
            Wprim(:, k) = v;
        end
        fprintf(1, 'Primitive bases computed in %.2f s\n', toc(linear_solves))
        fprintf(1, '----------------------------------\n')
        % Fail-safe; comment out if you don't want to save your bases mid
        % iteration
        % fprintf(1, 'Saving your bases!\n')
        % fprintf(1, '------------------\n')
        % filename = 'results/prim_bases_g.mat';
        % save(filename, 'Vprim', 'Wprim');
    end
    if strcmp(opts.proj, 'pg')
        % In Petrov-Galerkin projection, Vprim, Wprim are as in (1), (2)
        linear_solves = tic;
        fprintf(1, 'Computing model reduction bases via Petrov-Galerkin projection \n')
        fprintf(1, '---------------------------------------------------------------\n')
        for k = 1:q
            fprintf(1, 'Current iterate is k = %d\n', k)
            fprintf(1, '-------------------------\n')
            % Option 0 in 'so_structured_solve.m' implements (-shifts(k) * E - A)\b)
            v           = so_structured_solve(Mso, Dso, Kso, bso, -shifts(k), 0, 1);
            Vprim(:, k) = v;
            % Option 1 in 'so_structured_solve.m' implements ((-shifts(k) * E - A)')\(Q*v))
            w           = so_structured_solve(Mso, Dso, Kso, Qfo*v, -shifts(k), 1, 1);
            Wprim(:, k) = w;
        end
        fprintf(1, 'Primitive bases computed in %.2f s\n', toc(linear_solves))
        fprintf(1, '----------------------------------\n')
        % Fail-safe; comment out if you don't want to save your bases mid
        % iteration
        % fprintf(1, 'Saving your bases!\n')
        % fprintf(1, '------------------\n')
        % filename = 'results/prim_bases_pg.mat';
        % save(filename, 'Vprim', 'Wprim');
    end
% Otherwise, pre-computed bases are passed as arguments.
else 
    fprintf(1, 'Primitive bases passed as args; not recomputing\n')
    fprintf(1, '-----------------------------------------------\n')
    Vprim = opts.Vprim; 
    Wprim = opts.Wprim;
end

% If opts.recomp_evals, compute transfer function evaluations
% Cheap, because Vprim already given
if opts.recomp_evals
    compute_tf = tic;
    fprintf(1, 'Computing transfer function values at given shifts \n')
    fprintf(1, '---------------------------------------------------\n')
    H_shifts = zeros(q, 1);
    for k = 1:q
        H_shifts(k) = Vprim(:, k)'*Qfo*Vprim(:, k);
    end
    fprintf(1, 'H_shifts computed in %.2f s\n', toc(compute_tf))
    fprintf(1, '----------------------------\n')
else 
    fprintf(1, 'Transfer function values passed as args; not recomputing \n')
    fprintf(1, '---------------------------------------------------------\n')
    H_shifts = opts.H_shifts;
end

if strcmp(opts.compress, 'avg')
    compression = tic;
    fprintf(1, 'Computing orthonormalized model reduction bases via pivoted QR\n')
    fprintf(1, '--------------------------------------------------------------\n')
    [Vorth, ~, pV] = qr(Vprim, 'vector', 'econ');   
    [Worth, ~, pW] = qr(Wprim, 'vector', 'econ');
    % Grab r leading columns orthonormal columns from pivoted QR of Vprim
    % and Wprim
    Vorth = Vorth(:, 1:r);   Worth = Worth(:, 1:r); 
    fprintf(1, 'Vorth and Worth computed in %.2f s\n', toc(compression))
    fprintf(1, '----------------------------------\n')
end
if strcmp(opts.compress, 'Linfty')
    compression = tic;
    fprintf(1, 'Computing orthonormalized model reduction bases via greedy Linfty search\n')
    fprintf(1, '------------------------------------------------------------------------\n')
    % Space for indices
    p = zeros(r, 1); 
    % Columns of Vprim, Wprim, chosen where abs transfer function error is
    % maximized (among the discrete locations in shifts)

    % p1 is just the maximal value of the full-order transfer function,
    % since the reduced-model is `zero' at start
    [~, p1] = max(abs(H_shifts));   p(1) = p1;
    % Project down
    Er = Wprim(:, p1)' * Efo * Vprim(:, p1);  Ar = Wprim(:, p1)' * Afo * Vprim(:, p1); 
    Qr = Vprim(:, p1)' * Qfo * Vprim(:, p1);  br = Wprim(:, p1)' * bfo;
    for k = 1:r-1
        iter_compression = tic;
        fprintf(1, 'Current iterate of greedy Linfty search is k = %d  \n', k)
        fprintf(1, '---------------------------------------------------\n')
        % Evaluate the previous interpolatory reduced model at the
        % locations in shifts for Linfty error estimation
        Hr_shifts = zeros(q, 1);
        for j = 1:q
            tmpr = (shifts(j)*Er - Ar)\br;
            Hr_shifts(j) = br'*(((shifts(j)*Er - Ar)')\(Qr*tmpr)); 
        end
        % Approximate Linfty error, choose next column where error is max
        Linfty_error = abs(H_shifts - Hr_shifts);  
        % Set error at previously chosen points to zero
        for i = 1:k
            Linfty_error(p(i)) = 0;
        end
        [~, pk] = max(Linfty_error);
        p(k + 1) = pk;
        % Orthonormalize columns chosen by greedy search; avoids ill-cond
        [Worth, ~] = qr(Wprim(:, p(1:k+1)), 'econ');    
        [Vorth, ~] = qr(Vprim(:, p(1:k+1)), 'econ');
        % Next Proj LQO-ROM; grab columns stored in P
        Er = Worth'*Efo*Vorth;  Ar = Worth'*Afo*Vorth; 
        Qr = Vorth'*Qfo*Vorth;  br = Worth'*bfo;
        fprintf(1, 'Search iteration done in %.2f s\n', toc(iter_compression))
        fprintf(1, '-------------------------------\n')
    end
    pW = p; pV = p; % Save indices
    fprintf(1, 'Vorth and Worth computed in %.2f s\n', toc(compression))
    fprintf(1, '----------------------------------\n')
end
fprintf('Outputting projection matrices\n')
fprintf(1, '---------------------------\n')
end

