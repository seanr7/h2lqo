function [WPrim, VPrim, WOrth, VOrth, HShifts, pW, pV] = interpolatory_solves(...
                                                                    Mso, ...
                                                                    Dso, ...
                                                                    Kso, ...
                                                                    Qfo, ...
                                                                    Bso, ...
                                                                    shifts, ...
                                                                    r, ...
                                                                    opts)
%INTERPOLATORY_SOLVES Driver function for computing interpolatory
% projections of single input, single quadratic output linear systems
%
% SYNTAX: 
%   [WPrim, VPrim, WOrth, VOrth, HShifts, pW, pV] = interpolatory_solves(Mso, ...
%                                           Dso, Kso, Qfo, Bso, shifts, r, ...
%                                           opts)
%   [WPrim, VPrim, WOrth, VOrth, HShifts, pW, pV] = interpolatory_solves(Mso, ...
%                                           Dso, Kso, Qfo, Bso, shifts, r)
%                                                                       
%
% DESCRIPTION:
%   This function computes interpolatory projection matrices Vorth, Worth, 
%   for use in the model reduction of LQO systems with a single input and
%   single quadratic output. The methodologies used are outlined in Section 
%   5 of "Interpolatory model reduction of dynamical systems with root mean 
%   squared error" [Reiter and Werner, 2024].
%   2n x q primitive interpolatiton bases are computed using linear solves 
%   according to the input opts.proj, or passed; q = length(shifts), q>=r.
%   The primitive bases are compressed according to the selection procedure 
%   opts.compress.
%
%   The left 2n x r model reduction basis is computed as
%
%       VPrim(:, k) = (shifts(k)*Efo - Afo)\bfo; (1)
%
%   If performing Galerkin projection, take WPrim = VPrim.
%   If performing Petrov-Galerkin projection, right 2n x r model reduction
%   basis is computed as
%
%       WPrim(:, k) = (shifts(k)*Efo - Afo)'\(Qfo*
%                       ((shifts(k)*Efo - Afo)\bfo); (2)
%
%   where
%
%       Efo = [I  0;   0     Mso]; (3)
%       Afo = [I  0;  -Kso  -Dso]; (4)
%       Bfo = [0; Bso];            (5)
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
%   matrices are orthonormalized via an economy-size QR to give VOrth and
%   WOrth.
% 
%   If opts.compress == avg, primitive basis vectors are `averaged' using a
%   pivoted QR decomposition. The leading r <= q orthonormal columns are
%   chosen as VOrth, WOrth.
%
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
%   Bso    - sparse second order input matrix with dimensions n x 1 in 
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
%   | recompBases     | bool, do we recompute the matrices Vr, Wr, or are |
%   |                 | they passed as opts?                              |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | recompEvals     | bool, do we recompute the transfer function       |
%   |                 | evaluations used in the Linfty search, or are they|
%   |                 | passed as opts?                                   |
%   |                 | (default 1)                                       |
%   +-----------------+---------------------------------------------------+
%   | VPrim           | Left n x r interpolatory projection matrix        | 
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%   | WPrim           | Right n x r interpolatory projection matrix       | 
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%   | HShifts         | Full-order ransfer function evaluations at given  |
%   |                 | shifts, as q x 1 array                            | 
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   WPrim   - Primitive 2n x r left interpolatory model reduction basis,
%             constructed according to opts.proj
%   VPrim   - Primitive 2nn x r right interpolatory model reduction basis,
%             constructed according to opts.proj
%   WOrth   - Orthonormalized 2n x r left interpolatory model reduction 
%             basis, compressed according to opts.compress
%   VOrth   - Orthonormalized 2nn x r right interpolatory model reduction 
%             basis, compressed according to opts.compress
%   pV      - Columns of VPrim selected by opts.compress
%   pV      - Columns of WPrim selected by opts.compress
%   HShifts - Transfer function evaluations of full-order model (Efo, Afo, 
%             Bfo, Qfo) computed at inputted shifts.

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model reduction of dynamical systems 
% with root mean squared error"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

% Virginia Tech, USA
% Last editied: 9/12/2024
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of solves.
q      =  max(size(shifts));
[n, ~] = size(Mso);

% Set opts.
if (nargin < 7) 
    opts = struct(); 
end
if ~isfield(opts, 'proj')
    opts.proj = 'g'; 
end
if ~isfield(opts, 'compress')
    fprintf(1, 'DEFAULT VALUE FOR opts.compress = avg\n')
    fprintf(1, '------------------------------------------------\n')
    opts.compress = 'avg'; 
end
if ~isfield(opts, 'recompBases')
    fprintf(1, 'DEFAULT VALUE FOR opts.recomp_bases = 0\n')
    fprintf(1, '------------------------------------------------\n')
    opts.recompBases = 0; 
end
if ~isfield(opts, 'recompEvals')
    fprintf(1, 'DEFAULT VALUE FOR opts.recomp_evals = 0\n')
    fprintf(1, '------------------------------------------------\n')
    opts.recompEvals = 0;
end
if ~isfield(opts, 'VPrim')
    fprintf(1, 'NO LEFT PRIMITIVE BASIS SET opts.Vprim = []\n')
    fprintf(1, '------------------------------------------------\n')
    opts.VPrim = [];
end
if ~isfield(opts, 'WPrim')
    fprintf(1, 'NO LEFT PRIMITIVE BASIS SET opts.Wprim = []\n')
    fprintf(1, '------------------------------------------------\n')
    opts.WPrim = [];
end
if ~isfield(opts, 'HShifts')
    fprintf(1, 'DEFAULT VALUE FOR opts.HShifts = []\n')
    fprintf(1, '------------------------------------------------\n')
    opts.HShifts = [];
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
% METHODS.                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If not passed, primitive bases are computed.
if opts.recompBases 
    % Space allocation.
    VPrim = zeros(2*n, q);   WPrim = zeros(2*n, q);

    % Option 1. Galerkin projection, VPrim = WPrim.
    if strcmp(opts.proj, 'g')
        linearSolves = tic;
        % Time all solves.
        fprintf(1, 'COMPUTING BASES VIA GALERKIN PROJECTION. \n')
        fprintf(1, '----------------------------------------------\n')

        for k = 1:q
            fprintf(1, 'CURRENT ITERATE IS k = %d\n', k)
            fprintf(1, '----------------------------------------------\n')
            % Solve v = (shifts(k) * E - A)\b.
            v           = so_structured_solve(Mso, Dso, Kso, Bso, shifts(k), 0, 1);
            VPrim(:, k) = v;
            WPrim(:, k) = v;
        end

        fprintf(1, 'PRIMITIVE BASES COMPUTED IN %.2f s\n', toc(linearSolves))
        fprintf(1, '----------------------------------------------\n')
        % `Fail safe' to save bases mid iteration.

        % fprintf(1, 'SAVING YOUR BASES!\n')
        % fprintf(1, '------------------\n')
        % filename = 'results/prim_bases_g.mat';
        % save(filename, 'Vprim', 'Wprim');
    end

    % Option 2. Petrov-Galerkin projection, VPrim != WPrim.
    if strcmp(opts.proj, 'pg')
        linearSolves = tic;
        % Time all solves.
        fprintf(1, 'COMPUTING BASES VIA PETROV-GALERKIN PROJECTON.\n')
        fprintf(1, '----------------------------------------------\n')

        for k = 1:q
            fprintf(1, 'Current iterate is k = %d\n', k)
            fprintf(1, '-------------------------\n')
            % Solve v = (shifts(k) * E - A)\b.
            v           = so_structured_solve(Mso, Dso, Kso, Bso, shifts(k), 0, 1);
            VPrim(:, k) = v;

            % Solve w = (shifts(k) * E - A)'\(Q*v).
            tmp         = Qfo*v;
            w           = so_structured_solve(Mso, Dso, Kso, tmp(1:n, 1), shifts(k), 1, 1);
            WPrim(:, k) = w;
        end
        fprintf(1, 'PRIMITIVE BASES COMPUTED IN %.2f s\n', toc(linearSolves))
        fprintf(1, '----------------------------------------------\n')
        % `Fail safe' to save bases mid iteration.

        % fprintf(1, 'SAVING YOUR BASES!\n')
        % fprintf(1, '------------------\n')
        % filename = 'results/prim_bases_pg.mat';
        % save(filename, 'Vprim', 'Wprim');
    end
else 
    fprintf(1, 'VPrim, WPrim, PASSED AS ARGUMENTS.\n')
    fprintf(1, '----------------------------------------------\n')
    VPrim = opts.VPrim; 
    WPrim = opts.WPrim;
end

% Transfer function evaluations for selection.
if opts.recompEvals
    computeTf = tic;
    fprintf(1, 'COMPUTING TRANSFER FUNCTION VALUES.\n')
    fprintf(1, '----------------------------------------------\n')
    HShifts = zeros(q, 1);
    for k = 1:q
        HShifts(k) = VPrim(:, k)'*Qfo*VPrim(:, k);
    end
    fprintf(1, 'VALUES COMPUTED IN %.2f s\n', toc(computeTf))
    fprintf(1, '----------------------------------------------\n')
else 
    fprintf(1, 'HShifts PASSED AS ARGUMENTS.\n')
    fprintf(1, '----------------------------------------------\n')
    HShifts = opts.HShifts;
end

% Option 1. Compress bases via pivoted QR.
if strcmp(opts.compress, 'avg')
    compression = tic;
    fprintf(1, '`avg` COMPRESSION'.\n')
    fprintf(1, '----------------------------------------------\n')
    [VOrth, ~, pV] = qr(VPrim, 'vector', 'econ');   
    [WOrth, ~, pW] = qr(WPrim, 'vector', 'econ');

    % Leading columns.
    VOrth = VOrth(:, 1:r);   WOrth = WOrth(:, 1:r); 
    fprintf(1, 'VOrth, WOrth COMPUTED IN%.2f s\n', toc(compression))
    fprintf(1, '----------------------------------------------\n')
end
if strcmp(opts.compress, 'Linfty')
    compression = tic;
    fprintf(1, '`Linfty` COMPRESSION.\n')
    fprintf(1, '----------------------------------------------\n')
    % Greedily choose next column based on the shift at where the error is
    % highest.
    p = zeros(r, 1); 

    % Note: p1 is just the maximal value of the full-order transfer 
    % function, since the reduced-model is `zero' at start.
    [~, p1] = max(abs(HShifts));   p(1) = p1;
    
    % Project.
    Er = WPrim(:, p1)'*Efo*VPrim(:, p1);  Ar = WPrim(:, p1)'*Afo*VPrim(:, p1); 
    Qr = VPrim(:, p1)'*Qfo*VPrim(:, p1);  Br = WPrim(:, p1)'*Bfo;
    for k = 1:r-1
        timeCompression = tic;
        fprintf(1, 'CURRENT ITERATE OF Linfty SEARCH IS k = %d\n', k)
        fprintf(1, '----------------------------------------------\n')
        
        % Error evaluation.
        HrShifts = zeros(q, 1);
        for j = 1:q
            tmpr = (shifts(j)*Er - Ar)\Br;
            HrShifts(j) = Br'*(((shifts(j)*Er - Ar)')\(Qr*tmpr)); 
        end
        % Approximate Linfty error, choose next column where error is max
        Linfty_error = abs(HShifts - HrShifts);  
        % Set error at previously chosen points to zero
        for i = 1:k
            Linfty_error(p(i)) = 0;
        end
        [~, pk] = max(Linfty_error);
        p(k + 1) = pk;

        % Orthogonalize to avoid ill-conditioning in intermediate
        % projections.
        [WOrth, ~] = qr(WPrim(:, p(1:k+1)), 'econ');    
        [VOrth, ~] = qr(VPrim(:, p(1:k+1)), 'econ');

        Er = WOrth'*Efo*VOrth;  Ar = WOrth'*Afo*VOrth; 
        Qr = VOrth'*Qfo*VOrth;  Br = WOrth'*Bfo;
        fprintf(1, 'SEARCH ITERATION FINISHED IN %.2f s\n', toc(timeCompression))
        fprintf(1, '----------------------------------------------\n')
    end

    % Indices picked by greedy search.
    pW = p; pV = p; 
    fprintf(1, 'VOrth, WOrth COMPUTED IN %.2f s\n', toc(compression))
    fprintf(1, '----------------------------------------------\n')
end
fprintf('OUTPUTTING BASES.\n')
fprintf(1, '----------------------------------------------\n')
end