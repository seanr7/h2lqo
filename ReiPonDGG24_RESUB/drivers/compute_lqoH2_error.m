function error = compute_lqoH2_error(A, B, C, M, Ar, Br, Cr, Mr, fomH2norm)
% COMPUTE_LQOH2_ERROR function to compute H2 error between two
% multiple-input, single-output linear quadratic-output (LQO) systems.
%
% SYNTAX:
%   error = compute_lqoH2_error(A, B, C, M, Ar, Br, Cr, Mr, fomH2norm)
%
% DESCRIPTION:
%   Computes the H2 error between two LQO systems specified by:
%       (1) Sys    = (A, B, C, M)
%       (2) SysRed = (Ar, Br, Cr, Mr)
%   using the Gramian- and Sylvester equation-based formulas (40) and (42)
%   for the squared H2 error taken from the companion paper.
%
% INPUTS:
%   A         - state matrix with dimensions n x n in (1)
%   B         - input matrix with dimensions n x m in (1)
%   C         - linear output matrix with dimensions 1 x n in (1)
%   M         - (symmetric) quadratic output matrix with dimensions n x n 
%               in (1)
%   Ar        - state matrix with dimensions r x r in (2)
%   Br        - input matrix with dimensions r x m in (2)
%   Cr        - linear output matrix with dimensions 1 x r in (2)
%   Mr        - (symmetric) quadratic output matrix with dimensions r x r 
%               in (2)
%   fomH2norm - H2 norm of full-order model Sys in (1)
%
% OUTPUTS:
%   error - computed H2 (relative) error between Sys and SysRed

%
% Copyright (c) 2025 Sean Reiter
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%
% Virginia Tech, USA
% Last editied: 6/23/2025
%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUTS.                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[n, ~] = size(A);
[r, ~] = size(Ar);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ERROR.                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X  = mess_sylvester_sparse_dense(A, 'N', Ar, 'T', B*Br', speye(n, n), eye(r, r));
Pr = lyap(Ar, Br*Br'); 
Y  = mess_sylvester_sparse_dense(A, 'T', Ar, 'N', -C'*Cr - M*X*Mr, speye(n, n), eye(r, r));
Qr = lyap(Ar', Cr'*Cr + Mr*Pr*Mr);

% Trace formula.
error = sqrt(abs(fomH2norm^2 + trace(Br'*Qr*Br) + 2*trace(B'*Y*Br)))/fomH2norm;
end