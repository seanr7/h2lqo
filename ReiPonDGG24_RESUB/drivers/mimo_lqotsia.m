function [Ar, Br, Cr, Mr, info] = mimo_lqotsia(A, B, C, M, r, opts)
%MIMO_LQOTSIA Two-sided iteration algorithm for model-order reudction of
% linear systems with multiple quadratic outputs
%
% SYNTAX:
%   [Ar, Br, Cr, Mr, info] = mimo_lqotsia(A, B, C, M, r, opts)
%
% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Ar, Br, Cr, Mr)
%   using the two-sided iteration algorithm given in 
%   "$\mathcal{H}_2 optimal model reduction of linear systems with 
%   multiple quadratic outputs"
%   At each iteration, two Sylvester equations
%
%       A*X + X*Ar' + B*Br' = 0                                   (1)
%       A'*Z + Z*Ar - 2*M1*X*M1r - ... - 2*Mp*X*Mpr - C*Cr' = 0   (2)
%
%   are solved, and a linear quadratic output reduced model is obtained via
%   projeciton W = Z and V = X.
%   It is assumed that the eigenvalues of (s*I-A) lie in the open left
%   half-plane.
%
% INPUTS:
%   A    - state matrix with dimensions n x n in (1)
%   B    - input matrix with dimensions n x m in (1)
%   C    - linear output matrix with dimensions p x n in (2)
%          if empty set to zeros(p, n)
%   M    - 3d-array of (symmetric) quadratic output matrices with 
%          dimensions p x n x n in (2)
%          if empty set to zeros(n, n)
%   r    - order of reduction
%   opts - structure, containing the following optional entries:
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | tol             | nonnegative scalar, tolerance for convergence     |
%   |                 | based on the H2error of the reduced model         |
%   |                 | (default 10e-4)                                   |
%   +-----------------+---------------------------------------------------+
%   | maxiter         | positive integer, maximum number of iteration     |
%   |                 | steps                                             |
%   |                 | (default 100)                                     |
%   +-----------------+---------------------------------------------------+
%   | Ar              | initial state matrix                              |
%   |                 | (default -diag(logspace(0,4,r))                   |
%   +-----------------+---------------------------------------------------+
%   | Br              | initial input matrix                              |
%   |                 | (default eye(n, m))                               |
%   +-----------------+---------------------------------------------------+
%   | Cr              | initial linear-output matrix                      |
%   |                 | (default eye(p, r))                               |
%   +-----------------+---------------------------------------------------+
%   | Mr              | initial quadratic-output matrices                 |
%   |                 | (default repmat(eye(n,n), 1, 1, p)                |
%   +-----------------+---------------------------------------------------+
%   | storeBases      | option to store final model reduction bases       |
%   |                 | (default false)                                   |
%   +-----------------+---------------------------------------------------+
%   | convMonitoring  | string, for monitoring convergence using exact or |
%   |                 | approximate calculation of the squared H2 error   |
%   |                 | (default 'approximate')                           |
%   +-----------------+---------------------------------------------------+
%   | fomNorm         | optional argument, value of the H2 norm of the    |
%   |                 | fom; use in computing a hierarchy of ROMs         |
%   |                 | (default [])                                      |
%   +-----------------+---------------------------------------------------+
%
% OUTPUTS:
%   Ar   - reduced state matrix with dimensions r x r in (1)
%   Br   - reduced descriptor matrix with dimensions r x m in (1)
%   Cr   - reduced linear output matrix with dimensions p x r in (2)
%          If C is zero then Cr is zeros(p, r)
%   Mr   - 3d-array of reduced (symmetric) quadratic output matrices 
%          with dimensions p x r x r in (2) 
%          If M is zero then Mr is zeros(r, r)
%   info - structure, containing the following information for monitoring
%          convergence
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | errors          | history of squared relative H2 errors throughout  |
%   |                 | the iteration                                     |
%   |                 | (tails if opts.convMonitoring == 'approximate')   |
%   +-----------------+---------------------------------------------------+
%   | tails           | variable parts of the relative H2 error           |
%   +-----------------+---------------------------------------------------+ 
%   | changeInErrors  | change in squared H2 errors throughout            |
%   |                 | (changeInTails if opts.convMonitoring ==          |
%   |                 |  'approximate')                                   |
%   +-----------------+---------------------------------------------------+
%   | changeInTails   | change in tails of squared H2 errors throughout   |
%   +-----------------+---------------------------------------------------+ 
%   | Wt              | final left model reduction basis                  |
%   +-----------------+---------------------------------------------------+
%   | V               | final right model reduction basis                 |
%   +-----------------+---------------------------------------------------+ 
%

%
% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 9/23/2024
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem dimensions.
n = size(A, 1);
m = size(B, 2);
if isempty(C)
    try
        [~, ~, p] = size(M, 1);
    catch
        p = 1;
    end
else
    p = size(C, 1);
end

% Check and set input opts.
if (nargin < 6) 
    opts = struct();
end

if ~isfield(opts, 'tol')
    opts.tol = 10e-6;
end
if ~isfield(opts, 'maxIter')
    opts.maxIter = 100;
end
if ~isfield(opts, 'Ar')
    opts.Ar = -diag(logspace(0, 4, r)); 
end
if ~isfield(opts, 'Br')
    opts.Br = eye(r, m);  
end
if ~isfield(opts, 'Cr')
    if isempty(C)
        opts.Cr = zeros(p, r);
    else
        opts.Cr = eye(p, r);
    end
end
if ~isfield(opts, 'Mr')
    if isempty(M)
        opts.Mr = repmat(zeros(r, r), 1, 1, p);
    else
        opts.Mr = repmat(eye(r, r), 1, 1, p);
    end
end
if ~isfield(opts, 'storeBases')
    opts.storeBases = false;
end
if ~isfield(opts, 'convMonitoring')
    opts.convMonitoring = 'approximate';
end
if ~isfield(opts, 'fomNorm')
    opts.fomNorm = [];
end

assert(strcmp(opts.convMonitoring, 'approximate') || strcmp(opts.convMonitoring, 'exact'), ...
    'Unrecognized input for convergence monitoring!')

% Use bool to avoid repeated strcmps. 
if strcmp(opts.convMonitoring, 'exact')
    exactConv = true;
else
    exactConv = false; 
end

% Set linear output to zero if not included.
if isempty(C)
    C = zeros(p, n);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ITERATION.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overallStart = tic;
fprintf(1, 'INITIALIZING ALGORITHM.\n');
fprintf(1, '--------------------------------------------------\n');
% Initial reduced model.
Ar = opts.Ar;   Br = opts.Br;   Cr = opts.Cr;   Mr = opts.Mr;  

% Pre-compute FOM H2 norm (if desired).
if exactConv
    if ~isempty(opts.fomNorm)
        % If FOM H2 norm is passed.
        h = opts.fomNorm;
        fprintf(1, 'H2 NORM OF THE FULL-ORDER MODEL PASSED AS AN ARG.\n');
        fprintf(1, '--------------------------------------------------\n');
    else
        normStart = tic;
        fprintf(1, 'COMPUTING H2 NORM OF FULL-ORDER MODEL.\n');
        fprintf(1, '--------------------------------------------------\n');
        P   = lyap(A, B*B');
        rhs = C'*C; 
        for i = 1:p
            rhs = rhs + M(:, :, i)*P*M(:, :, i);
        end
        Q  = lyap(A', rhs);
        
        % (Squared) H2 norm of full-order model.
        h  = trace(B'*Q*B); 
        fprintf(1, 'H2 NORM COMPUTED IN %d s\n', toc(normStart));
        fprintf(1, '--------------------------------------------------\n');
    end
else
    fprintf(1, 'H2 NORM NOT COMPUTED OR PASSED; MONITORING TAILS FOR CONVERGENCE.\n');
    fprintf(1, '-----------------------------------------------------------------\n');
end

% (Squared) relative H2 error due to initial reduced model.
Pr    = lyap(Ar, Br*Br'); 
X     = mess_sylvester_sparse_dense(A, 'N', Ar, 'T', B*Br', speye(n, n), eye(r, r));
rhsQr = Cr'*Cr; 
rhsY  = -C'*Cr;
for i = 1:p
    rhsQr = rhsQr + Mr(:, :, i)*Pr*Mr(:, :, i);
    rhsY  = rhsY - M(:, :, i)*X*Mr(:, :, i);
end
Qr = lyap(Ar', rhsQr);
Y  = mess_sylvester_sparse_dense(A, 'T', Ar, 'N', rhsY, speye(n, n), eye(r, r));

% (Squared) relative H2 error from trace formula.
iterate         = 1;
tails(iterate)  = trace(Br'*Qr*Br) + 2*trace(B'*Y*Br); % `Tail' of H2 error
if exactConv
    errors(iterate) = (h + tails(iterate))/h; % Relative, squared H2 error
else
    errors(iterate) = tails(iterate);
end

% Set iteration count and tolerance to engage `while'.
changeInErrors(iterate) = opts.tol + 1;
changeInTails(iterate)  = opts.tol + 1;

% If opts.convMonitoring == 'exact', changeInErrors will be exact change in
% the error. Else, it will be changeInTails. 
while (changeInErrors(iterate) > opts.tol && iterate <= opts.maxIter)
    thisIterStart = tic; 
    fprintf(1, 'CURRENT ITERATE IS k = %d\n', iterate);
    fprintf(1, '--------------------------------------------------\n');

    % Solve for right projection matrix V = X in (1).
    % XCheck = lyap(A, Ar', B*Br');
    X = mess_sylvester_sparse_dense(A, 'N', Ar, 'T', B*Br', speye(n, n), eye(r, r));

    % Solve for left projection matrix W = Z in (2).
    rhs = -C'*Cr;
    for i = 1:p
        rhs = rhs - 2*M(:, :, i)*X*Mr(:, :, i);
    end
    % ZCheck = lyap(full(A'), Ar, rhs);
    Z = mess_sylvester_sparse_dense(A, 'T', Ar, 'N', rhs, speye(n, n), eye(r, r));
    
    % Orthonormalize projection matrice.
    [V, ~] = qr(X, "econ");     [W, ~] = qr(Z, "econ");

    % Compute reduced model via projection.
    Wt = (W.'*V)\(W.');
    Ar = Wt*A*V;   Br = Wt*B;   Cr = C*V;   
    for i = 1:p
        Mr(:, :, i) = V'*M(:, :, i)*V;    
    end
  
    iterate = iterate + 1;  
    fprintf(1, 'COMPUTING H2 ERROR AT CURRENT ITERATE.\n');
    fprintf(1, '--------------------------------------------------\n');
    Pr    = lyap(Ar, Br*Br'); 
    rhsQr = Cr'*Cr; 
    rhsY  = -C'*Cr;
    for i = 1:p
        % Note: No factor of 2 in right hand side for Y or Qr.
        rhsQr = rhsQr + Mr(:, :, i)*Pr*Mr(:, :, i); 
        % X saved above from Sylvester solve.
        rhsY  = rhsY - M(:, :, i)*X*Mr(:, :, i);  
    end
    Qr = lyap(Ar', rhsQr);
    Y  = mess_sylvester_sparse_dense(A, 'T', Ar, 'N', rhsY, speye(n, n), eye(r, r));

    % Trace formula.
    tails(iterate)  = trace(Br'*Qr*Br) + 2*trace(B'*Y*Br); % `Tail' of H2 error
    if exactConv
        errors(iterate) = abs((h + tails(iterate))/h); % Relatve, squared H2 error
        fprintf(1, 'RELATIVE, SQUARED H2 ERROR OF CURRENT MODEL ITERATE IS: e(k) = ||G - Gr(k)||_H2^2/||G||_H2^2 = %.12f\n', ...
            errors(iterate));
    else
        errors(iterate) = tails(iterate);
    end

    % Convergence tracking.
    changeInErrors(iterate) = abs(errors(iterate) - errors(iterate-1))/abs(errors(1));
    changeInTails(iterate)  = abs(tails(iterate) - tails(iterate-1))/abs(tails(1));
    if exactConv
        fprintf('RELATIVE CHANGE IN SQUARED H2 ERRORS IS |e(k) - e(k-1)|/e(1) = %.12f \n', changeInErrors(iterate))
    else
        fprintf('CHANGE IN VARIABLE TAILS OF THE SQUARED H2 ERROR IS |t(k) - t(k-1)|/t(1) = %.12f \n', changeInErrors(iterate))
    end
    fprintf(1, '--------------------------------------------------\n');
    
    % Timing.
    fprintf(1, 'CURRENT ITERATE FINISHED IN %.2f s.\n', toc(thisIterStart));
    fprintf(1, 'END OF CURRENT ITERATE k = %d.\n', iterate);
    fprintf(1, '--------------------------------------------------\n');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TERMINATION.                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
info                = struct();
info.errors         = errors;    
info.tails          = tails;
info.changeInErrors = changeInErrors;
info.changeInTails  = changeInTails;
if opts.storeBases
    info.Wt  = Wt;
    info.V   = V;
end

if iterate == (opts.maxIter + 1)
    fprintf('ALGORITHM HAS TERMINATED DUE TO REACHING MAX NO. OF ITERATIONS; TOTAL TIME ELAPSED IS %.2f s\n', toc(overallStart))
     fprintf(1, '--------------------------------------------------\n');
else
    fprintf('ALGORITHM HAS CONVERGED IN %d ITERATIONS.\n', iterate)
    fprintf('TOTAL TIME ELAPSED IS %.2f s\n', toc(overallStart))
    fprintf(1, '--------------------------------------------------\n');
end
end