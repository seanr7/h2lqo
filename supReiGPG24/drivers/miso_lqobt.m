function [Er, Ar, Br, cr, Mr, info] = miso_lqobt(E, A, B, c, M, r)
% MISO_LQOBT Balanced truncation algorithm for model-order reudction of
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
%       A*P*E' + E*P*A' + B*B' = 0                               (1)
%       A'*Q*E + E'*Q*A + c*c' + M*P*M= 0   (2)
%
%   are solved for. A singular value decomposition (svd) of U'*E*L is
%   computed
%
%       U'*E*L = [Z1, Z2] * diag(S1, S2) * [Y1, Y2]'          (2)
%
%   Where the matrices in the svd are partitioned according to the
%   reduction order 1 < r <= n. A reduced order model is obtained using the
%   projection matrices
%
%       Vr = U*Z1*S1^(-1/2) and Wr = L*Y1*S1^(-1/2)         (3)
%
%   It is assumed that the eigenvalues of (s*E - A) lie in the open left
%   half-plane, and that E is nonsingular.
%
% INPUTS:
%   E    - descriptor matrix with dimensions n x n in (1)
%   A    - state matrix with dimensions n x n in (1)
%   B    - input matrix with dimensions n x m in (1)
%   c    - linear output matrix with dimensions n x 1 in (2)
%   M    - (symmetric) quadratic output matrix n x n
%   r    - order of reduction
% 
% OUTPUTS:
%   Er    - reduced descriptor matrix with dimensions r x r
%   Ar    - reduced state matrix with dimensions r x r 
%   Br    - reduced descriptor matrix with dimensions r x m 
%   cr    - reduced linear output matrix with dimensions r x 1 
%   Mr    - reduced (symmetric) quadratic output matrix with dimensions 
%           r x r
%   info  - output info, containing the following 
%   +-----------------+---------------------------------------------------+
%   |    PARAMETER    |                     MEANING                       |
%   +-----------------+---------------------------------------------------+
%   | svs             | singular values of the matrix U'*L                |
%   +-----------------+---------------------------------------------------+
%   | leftProj        | left projection matrix in the order reduction     |
%   +-----------------+---------------------------------------------------+
%   | rightProj       | right projection matrix in the order reduction    |
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(A, 1);
m = size(B, 2);
p = 1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
overallStart = tic;

lyapTime = tic; 
fprintf(1, 'SOLVING FOR GRAMIANS FROM LYAPUNOV EQUATIONS.\n')
fprintf(1, '---------------------------------------------\n');
% Solve for reachability Gramian P.
% P = lyap(A, B*B', [], E);
% U = mess_lyap(A, B, [], [], E);

% Solve for quadratic output observability Gramian Q.
% L = mess_lyap(A', [c, M*U], [], [], E');
% Q = lyap(A', (c*c' + M*P*M), [], E');
fprintf(1, 'GRAMIANS SOLVED FOR IN %.2f s\n', toc(lyapTime))

fprintf(1, 'COMPUTING PROJECTION MATRICES.\n')
fprintf(1, '---------------------------- -\n');

% Compute square root factors of Gramians, and take the svd of U'*L.
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
[Z, S, Y] = svd(U'*E*L);

% Compute projection matrices.
V = U*Z(:, 1:r)*S(1:r, 1:r)^(-1/2); % Right
W = L*Y(:, 1:r)*S(1:r, 1:r)^(-1/2);% Left

% Compute reduced order model via projection.
Er  = eye(r, r);
Ar  = W.'*A*V;   Br = W.'*B;   
% Linear output.
cr = V.'*c;
% Quadratic output.
Mr = V'*M*V;

% Output info.
info           = struct();
info.svs       = diag(S);
info.leftproj  = W;
info.rightproj = V;

fprintf(1, 'ALGORITHM HAS COMPLETED.\n')
fprintf(1, 'TOTAL ELAPSED TIME IS %.2f s\n', toc(overallStart))
fprintf(1, '---------------------------------------\n')
end