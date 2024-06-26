function [v] = so_structured_solve(Mso, Dso, Kso, bso, s, struct_rhs, time)
% SO_STRUCTURED_SOLVE Function to implement a linear solve for descriptor
% and mass matrices with second order system structure.
%
%
% DESCRIPTION:
%   Function to compute a single (2n x 2n) linear system solve of the form
%
%       v = (s*E - A)\b;    (0a)
%    or w = ((s*E - A)')\b; (0b)
%
%   where the mass matrix (A), descriptor matrix (E), and right hand side
%   (b) are obtained from the first order realization of a second order 
%   dynamical system, and thus have the particular structure
%
%       E = [I  0;   0     Mso];       (1)
%       A = [I  0;  -Kso  -Dso];       (2)
%       b = [0; bso]; or b = [bso; 0]; (3)
%
%   Via the Woodbury matrix identity and the inverse formula of a 2 x 2
%   block matrix, v is instead computed in an equivalent way using only 
%   n x n linear solves.
%   Option 1: v = (s*E - A)\b for b = [0; bso], then
%       
%       z = (s*Mso + Dso)\bso;                              (4a)
%       v = [(1/s)*(z - ((s^2)*Mso + s*Dso + Kso)\(Kso*z)); (4b)
%            s*((s^2)*Mso + s*Dso + Kso)\bso];
% 
%   If w = ((s*E - A)')\b for b = [bso; 0], then
%
%       z = ((conj(s)^2)*Mso + conj(s)*Dso + Kso)\(Kso*bso); (5a)
%       v = [(1/conj(s))*(bso - Kso*z);                      (5b)
%            z];               
%
%   It is assumed that the complex shift s is not a pole of the matrix
%   pencil (s*E - A) and (s*M + D), and that s is strictly nonzero.
%
% INPUTS:
%   Mso       - sparse second order mass matrix with dimensions n x n in 
%                (1)
%   Dso       - sparse second order damping matrix with dimensions n x n 
%                in (2)
%   Kso       - sparse second order stiffness matrix with dimensions n x n 
%                in (2)
%   bso       - sparse second order input matrix with dimensions n x 1 in 
%                (3)
%   s         - complex shift in linear solve
%   solve_opt - boolean, do we solve system (0a) or (0b)?
%                  0 if v = (s*E - A)\b with b = [0;   bso];                            
%                  1 if w = ((s*E - A)')\b with b = [bso; 0];  
%   time       - optional boolean argument to print time required
%
% OUTPUTS:
%   v - sparse solution to the linear system (0a) or (0b) with dimensions 
%       2n x 1 computed accoding to (4a) and (4b) or (5a) and (5b)

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model order reduction of large-scale 
% dynamical systems with root mean squared error measures"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

% Virginia Tech, USA
% Last editied: 4/12/2024
%%
% Check and set inputs
if nargin < 7
    % Default is to not time solutions
    time = 0;
end
if nargin < 6
    error('Must specify structure of the right hand side!\n')
end
[n, ~] = size(Mso);

%%
if time
    tic
    fprintf(1, 'Initiating structured linear solve\n')
    fprintf(1, '----------------------------------\n')
end

% Structured solve; option 1 in (4a), (4b)
if struct_rhs == 0 % if v = (s*E - A)\b with b = [0;   bso]
    z  = (s*Mso + Dso)\bso; 
    v1 = (1/s).*(z - ((s^2).*Mso + s.*Dso + Kso)\(Kso*z));
    v2 = s.*(((s^2).*Mso + s.*Dso + Kso)\bso);
else % if w = ((s*E - A)')\b with b = [bso; 0]
    sconj = conj(s);    
    z     = ((sconj^2).*Mso + sconj.*Dso + Kso)\bso;
    v1    = (1/sconj).*(bso - Kso*z);
    v2    = z;
end
v             = spalloc(2*n, 1, nnz(v1) + nnz(v2));
v(1:n, :)     = v1;
v(n+1:2*n, :) = v2;

if time
    fprintf(1, 'Structured solve finished in %.2f s\n',toc)
    fprintf(1, '-----------------------------------------\n');
end
