function [Ar, Br, Cr, Mr, info] = mimolqo_bt(A, B, C, M, r)
% MIMOLQO_BT Balanced truncation algorithm for model-order reudction of
% linear systems with multiple quadratic outputs
%
% DESCRIPTION:
%   Computes a linear quadratic output reduced model (Ar, Br, Cr, M1r, ...,
%   Mpr) using the balanced truncation algorithm given in:
%   "Gramians, Energy Functionals, and Balanced Truncation for Linear 
%    Dynamical Systems with Quadratic Outputs" by Peter Benner, 
%   Pawan Goyal, and Igor Pontes Duff. [BenGP22]
% 
%   Square root factors of the quadratic output system Gramians P = U*U'
%   and Q = L*L', which satisfy
%
%       A*P + P*A' + B*B' = 0                               (1)
%       A'*Q + Q*A + C'*C + M1*P*M1 + ... + 2*Mp*P*Mp = 0   (2)
%
%   are solved for. A singular value decomposition (svd) of U'*L is
%   computed
%
%       U'*L = [Z1, Z2] * diag(S1, S2) * [Y1, Y2]'          (3)
%
%   Where the matrices in the svd are partitioned according to the
%   reduction order 1 < r <= n. A reduced order model is obtained using the
%   projection matrices
%
%       Vr = U*Z1*S1^(-1/2) and Wr = L*Y1*S1^(-1/2)         (4)
%
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
% 
% OUTPUTS:
%   Ar    - reduced state matrix with dimensions r x r in (1)
%   Br    - reduced descriptor matrix with dimensions r x m in (1)
%   Cr    - reduced linear output matrix with dimensions p x r in (2)
%           If C is zero then Cr is zeros(p, r)
%   Mr    - 3d-array of reduced (symmetric) quadratic output matrices with 
%           dimensions p x r x r in (2) 
%           If M is zero then Mr is zeros(r, r)
%   info  - output info, containing the following 
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | svs             | singular values of the matrix U'*L                |
%   +-----------------+---------------------------------------------------+
%   | leftproj        | left projection matrix in the order reduction     |
%   +-----------------+---------------------------------------------------+
%   | rightproj       | right projection matrix in the order reduction    |
%   +-----------------+---------------------------------------------------+
%

%
% Copyright (c) 2024 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 3/27/2024
%

%%
% Grab state, input, output dimensions
% Check input matrices.
n = size(A, 1);
if isempty(C)
    try
        [~, ~, p] = size(M, 1);
    catch
        p = 1;
    end
    pureqo = true;
    C      = zeros(p, r);
    Cr     = zeros(p, r);
else
    p      = size(C, 1);
    pureqo = false;
end

%% Begin algorithm.
overall_start = tic;
fprintf(1, 'Beginning algorithm\n')
fprintf(1, '-------------------\n');

lyap_time = tic; 
fprintf(1, 'Solving for Gramians from Lyapunov equations\n')
fprintf(1, '--------------------------------------------\n');
% Solve for reachability Gramian P
P = lyap(A, B*B');

% Solve for quadratic output observability Gramian Q
if ~pureqo % If output is not purely quadratic
    rhs = C'*C;
else
    rhs = zeros(n, n);
end
for i = 1:p
    rhs = rhs + M(:, :, i)*P*M(:, :, i);
end
Q = lyap(A', rhs);
fprintf(1, 'Gramians solved for in %.2f s\n', toc(lyap_time))

fprintf(1, 'Computing projection matrices\n')
fprintf(1, '---------------------------- \n');

% Compute square root factors of Gramians, and take the svd of U'*L
try
    U = chol(P);   
    U = U';
catch
    U = chol(P + eye(n, n)*10e-10);   
    U = U';
end
try
    L = chol(Q);   
    L = L';
catch
    L = chol(Q + eye(n, n)*10e-10);    
    L = L';
end
[Z, S, Y] = svd(U'*L);

% Compute projection matrices
V = U*Z(:, 1:r)*S(1:r, 1:r)^(-1/2); % Right
W = L*Y(:, 1:r)*S(1:r, 1:r)^(-1/2);% Left

% Compute reduced order model via projection
Ar  = W'*A*V;   Br = W'*B;   
% Linear output
if ~pureqo
    Cr = C*V;   
end
% Quadratic output
Mr = repmat(zeros(r, r), 1, 1, p);
for i = 1:p
    Mr(:, :, i) = V'*M(:, :, i)*V;
end

% Output info
info           = struct();
info.svs       = diag(S);
info.leftproj  = W;
info.rightproj = V;

fprintf(1, 'Algorithm has completed\n')
fprintf(1, 'Total time elapsed is %.2f s\n', toc(overall_start))
fprintf(1, '---------------------------------------\n')
end