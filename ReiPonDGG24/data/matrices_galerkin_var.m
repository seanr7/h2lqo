function [Agal,Bgal,Cgal,Ngal,Hgal] = matrices_galerkin_var(degr,sigma)
% (semi-)analytical computation of matrices in Galerkin system
% with variance as single output
%
% Original Galerkin system :
%   d/dt E*v = A*v + B*u
%          w = C*v
% Differential equation for variance:
%   d/dt Var = 2 [ v' * A' * (Ccut'Ccut) * v + u * B' * (Ccut'Ccut) * v ]
%
% Input parameters:
%  degr    degree in polynomical chaos expansion
%  sigma   magnitude of random perturbation
%

% evaluation of Galerkin system with gPC coefficients as output
[Ahat,Bhat,Chat,Ehat] = matrices_galerkin(degr,sigma);
dim = length(Ahat);

% invert matrices
Ahat = Ehat\Ahat;
Bhat = Ehat\Bhat;

% system matrix
Agal = sparse(dim+1,dim+1);
Agal(1:dim,1:dim) = Ahat;

% input matrix
Bgal = sparse(dim+1,1);
Bgal(1:dim) = Bhat;

% output matrix
Cgal = sparse(1,dim+1);
Cgal(end) = 1.;

% matrix-matrix product (used twice)
% first row of Chat cut off
Ccut = Chat(2:end,:);
M = Ccut'*Ccut;

% matrix of bilinear term
Ngal = sparse(1,dim+1);
Ngal(1,1:dim) = 2*Bhat'*M;

% matrix of quadratic term
% symmetrisation is used
A = Ahat'*M;
Hgal = sparse(dim+1,dim+1);
Hgal(1:dim,1:dim) = A + A';

%return
end