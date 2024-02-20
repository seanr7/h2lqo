function [Wprim, Vprim, Worth, Vorth, H_shifts, pW, pV, opts] = interpolatory_solves(E, A, B, Q, shifts, r, opts)
% Author: Sean Reiter (seanr7@vt.edu)

% Function to compute linear solves for interpolatory model-order reduction
% of linear quadratic-output (LQO) dynamical systems. `Primitive' bases are 
% computed via linear solves at the inputted shifts; projection matrices 
% are compressed according to input ops.

% Inputs:
% @param E, A, B, C, Q: Parameters of full-order LQO model
% @param shifts:        Complex frequencies at which to interpolate the FOM
% @param opts:          Struct containing input opts

% Outputs:
% @param Wprim:    `Primitive' left model reduction basis; contruction
%                  according to opts.proj
% @param Vprim:    `Primitive' right model reduction basis; contruction
%                  according to opts.proj
% @param Worth:    `Compressed' left model reduction basis; contruction
%                  according to opts.compression
% @param Vorth:    `Compressed' left model reduction basis; contruction
%                  according to opts.compression
% @param H_shifts: The quadratic-output transfer function evaluated at the
%                  inputted shifts
% @param: pV, pW:  Indices of shifts included in projection subspaces;
%                  outputted for the sake of checking interpolatory cons

% Opts:
% @param compression:   How to compute Wroth, Vrorth, from Wrprim, Vrprim
%                       Options
%                           'Linfty': Build orth bases by selectiong r <= q 
%                                     primitive basis vectors via a greedy
%                                     selection; subsequent vectors are
%                                     selected at shifts whre the transfer
%                                     function error is maximized (recycles
%                                     linear solves)
%                           'avg':    Compute orth bases via a pivoted QR
%                                     applied to prim bases; keep r <=
%                                     q leading columns
% @param proj:          What kind of projection do we use?
%                           'g':      Primitive bases taken as Wrprim =
%                                     Vrprim
%                           'pg':     Primitive bases are computed
%                                     differently; Wrprim encodes
%                                     information of QO matrix Q
% @param recomp_bases:  Bool; do we need to recompute primitive bases, or
%                       are they passed as an arg?
% @param recomp_tf:     Bool; do we need to recompute tf values, or
%                       are they passed as an arg?
% @param Vprim, Wprim:  Precommputed primitive bases
% @param H_shifts:      Precommputed tf values
%%
% Default opts
if ~isfield(opts, 'proj')
    fprintf('Setting default value for opts.proj = g\n')
    opts.proj = 'g'; % Galerkin projection
end
if ~isfield(opts, 'compression')
    fprintf('Setting default value for opts.compression = avg\n')
    opts.compression = 'avg'; % Orth via Pivoted QR
end
if ~isfield(opts, 'recomp_bases')
    fprintf('Setting default value for opts.recomp_bases = 0\n')
    opts.recomp_bases = 0; % False; recomp primitive bases
end
if ~isfield(opts, 'recomp_tf')
    fprintf('Setting default value for opts.recomp_tf = 0\n')
    opts.recomp_tf = 0; % False; recomp primitive bases
end
if ~isfield(opts, 'Vprim')
    fprintf('Setting default value for opts.Vprim = []\n')
    opts.Vprim = []; % No pre-computed bases
end
if ~isfield(opts, 'Wprim')
    fprintf('Setting default value for opts.Wprim = []\n')
    opts.Wprim = [];
end
if ~isfield(opts, 'H_shifts')
    fprintf('Setting default value for opts.H_shifts = []\n')
    opts.H_shifts = [];
end

q =  max(size(shifts));
[n, ~] = size(A);
if opts.recomp_bases % If we do compute the primitive bases
    Vprim = zeros(n, q);   Wprim = zeros(n, q);
    % Compute primitive bases
    if strcmp(opts.proj, 'g')
        % Primitive bases are identical as in Galerkin projection
        linear_solves = tic;
        fprintf('Computing model reduction bases via Galerkin projection \n')
        for k = 1:q
            tmp = (shifts(k) * E - A)\B;
            Vprim(:, k) = tmp;  Wprim(:, k) = tmp;
        end
        fprintf('Vr and Wr computed in %.2f s\n', toc(linear_solves))
    end
    if strcmp(opts.proj, 'pg')
        % Primitive bases are different as in Petrov-Galerkin projection
        % Left projection matrix encodes QO matrix
        linear_solves = tic;
        fprintf('Computing model reduction bases via Petrov-Galerkin projection \n')
        for k = 1:q
            tmp = (shifts(k) * E - A)\B;
            Vprim(:, k) = tmp;
            Wprim(:, k) = ((shifts(k) * E - A)' \ (Q * tmp));
        end
        fprintf('Vr and Wr computed in %.2f s\n', toc(linear_solves))
    end
    % Just to be safe ...
    fprintf('Saving your bases!\n')
    filename = 'prim_bases.mat';
    save(filename, 'Vprim', 'Wprim');
else % Bases are given in opts
    fprintf('Primitive bases passed as args; not recomputing \n')
    Vprim = opts.Vprim; Wprim = opts.Wprim;
end

if opts.recomp_tf
    % Next, compute H at shifts; use precomputed solves
    compute_tf = tic;
    fprintf('Computing transfer function values at given shifts \n')
    H_shifts = zeros(q, 1);
    for k = 1:q
        H_shifts(k) = Vprim(:, k)' * Q * Vprim(:, k);
    end
    fprintf('H_shifts computed in %.2f s\n', toc(compute_tf))
else % Tf values given in opts
    fprintf('Tf values passed as args; not recomputing \n')
    H_shifts = opts.H_shifts;
end

if strcmp(opts.compression, 'avg')
    compression = tic;
    fprintf('Compressing projection bases via pivoted-QR \n')
    [Vorth, ~, pV] = qr(Vprim, 'vector', 'econ');   [Worth, ~, pW] = qr(Wprim, 'vector', 'econ');
    % Grab r `leading' columns of primitive bases according to pivoted QR
    Vorth = Vorth(:, 1:r);   Worth = Worth(:, 1:r); % Double check this
    fprintf('Vr and Wr orthonormalized in %.2f s\n', toc(compression))
end
if strcmp(opts.compression, 'Linfty')
    overall_compression = tic;
    fprintf('Compressing projection bases via greedy Linfty search \n')
    % Choose r columns greedily where error is maximized
    p = zeros(r, 1); % Space for indices chosen via greedy search
    % First index is just max magnitude tf value, since `LQO-ROM' is zero
    [~, p1] = max(abs(H_shifts));    p(1) = p1;
    % Project down
    Er = Wprim(:, p1)' * E * Vprim(:, p1);  Ar = Wprim(:, p1)' * A * Vprim(:, p1); 
    Qr = Vprim(:, p1)' * Q * Vprim(:, p1);  Br = Wprim(:, p1)' * B;
    for k = 1:r-1
        iter_compression = tic;
        fprintf('Starting k = %d iteration of greedy Linfty search \n', k)
        % Eval. LQO-ROM and compute error (order is k at each iter)
        Hr_shifts = zeros(q, 1); % Eval at q shifts
        for j = 1:q
            tmpr = (shifts(j) * Er - Ar)\Br;
            Hr_shifts(j) = Br'* (((shifts(j) * Er - Ar)') \ (Qr * tmpr)); 
        end
        Linfty_error = abs(H_shifts - Hr_shifts); % Compute Linfty_error
        % Choose next index based on where error is biggest 
        % Hacky, but set previously chosen indices to 0
        for i = 1:k
            Linfty_error(p(i)) = 0;
        end
        [~, pk] = max(Linfty_error);
        p(k+1) = pk;
        % QR compressed basis to avoid ill conditioning
        [Worth, ~] = qr(Wprim(:, p(1:k+1)), 'econ');    [Vorth, ~] = qr(Vprim(:, p(1:k+1)), 'econ');
        % Next Proj LQO-ROM; grab columns stored in P
        Er = Worth' * E * Vorth;  Ar = Worth' * A * Vorth; 
        Qr = Vorth' * Q * Vorth;  Br = Worth' * B;
        fprintf('Search iteration done in %.2f s\n', toc(iter_compression))
    end
    % Worth = Wprim(:, p(1:k+1));    Vorth = Vprim(:, p(1:k+1));
    pW = p; pV = p; % Save indices
    fprintf('Vr and Wr orthonormalized in %.2f s\n', toc(overall_compression))
end
fprintf('Outputting\n')
end

