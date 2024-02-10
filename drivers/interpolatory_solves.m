function [Wprim, Vprim, Worth, Vorth, H_shifts, opts] = interpolatory_solves(E, A, B, Q, shifts, r, opts)
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

%%
if isempty(fieldnames(opts))
    % Default opts
    opts.proj = 'g'; % Galerkin projection
    opts.compression = 'avg'; % Orth via Pivoted QR
end
if length(fields(opts)) > 2
    error('Too many fields in input opts; only 2 allowed')
else
    if length(fields(opts)) == 1
        if isfield(opts, 'proj')
            opts.compression = 'avg';
        else % isfield(opts, 'compression')
            if isfield(opts, 'compression')
                opts.proj = 'g';
            else
                error('Unrecognized field in input opts')
            end
        end
    end
end

q =  max(size(shifts));
Vprim = zeros(q, 1);   Wprim = zeros(q, 1);
% Compute primitive bases
if strcmp(opts.proj, 'p')
    % Primitive bases are identical as in Galerkin projection
    for k = 1:q
        tmp = (shifts(k) * E - A)\B;
        Vprim(:, k) = tmp;  Wprim(:, k) = tmp;
    end
if strcmp(opts.proj, 'pg')
    % Primitive bases are different as in Petrov-Galerkin projection
    % Left projection matrix encodes QO matrix
    for k = 1:q
        tmp = (shifts(k) * E - A)\B;
        Vprim(:, k) = tmp;
        Wprim(:, k) = ((shifts(k) * E - A)' \ (Q * tmp));
    end
else
    error('NotImplementedError')
end
end

% Next, compute H at shifts; use precomputed solves
H_shifts = zeros(q, 1);
for k = 1:q
    H_shifts(k) = Vprim(:, k)' * Q * Vprim(:, k);
end

if strcmp(opts.compression, 'avg')
    [Vorth, ~, pV] = qr(Vprim, 'vector');   [Worth, ~, pW] = qr(Wprim, 'vector');
    % Grab r `leading' columns of primitive bases according to pivoted QR
    Vorth = Vorth(:, 1:r);   Worth = Worth(:, 1:r); % Double check this
if strcmp(opts.compression, 'Linfty')
    % Choose r columns greedily where error is maximized
    p = zeros(k, 1); % Space for indices chosen via greedy search
    % First index is just max magnitude tf value, since `LQO-ROM' is zero
    [~, p1] = max(abs(H_shifts));    p(1) = p1;
    % Project down
    Er = Wprim(:, p1)' * E * Vprim(:, p1);  Ar = Wprim(:, p1)' * A * Vprim(:, p1); 
    Qr = Vprim(:, p1)' * Q * Vprim(:, p1);  Br = Wprim(:, p1)' * B;
    for k = 1:r 
        % Eval. LQO-ROM and compute error (order is k at each iter)
        Hr_shifts = zeros(q, 1); % Eval at q shifts
        for j = 1:q
            tmpr = (shifts(j) * Er - Ar)\Br;
            Hr_shifts(j) = (((shifts(j) * Er - Ar)') \ (Qr * tmpr)); 
        end
        Linfty_error = abs(Hshifts - Hr_shifts); % Compute Linfty_error
        % Choose next index based on where error is biggest 
        [pk, ~] = mamx(Linfty_error);
        p(k) = pk;
        % Next Proj LQO-ROM; grab columns stored in P
        Er = Wprim(:, p)' * E * Vprim(:, p);  Ar = Wprim(:, p)' * A * Vprim(:, p); 
        Qr = Vprim(:, p)' * Q * Vprim(:, p);  Br = Wprim(:, p)' * B;
    end
    % Now, grab compressed basis
    Worth = Wprim(:, p);    Vorth = Vprim(:, p);
else
    error('NotImplementedError')
end
end

end

