function [A,B,M] = mat_galerkin_ex(degr,sigma)

% evaluation of Galerkin system with gPC coefficients as output
[Ahat,Bhat,Chat,Ehat] = matrices_galerkin(degr,sigma);
dim = length(Ahat);

% invert matrices
A = Ehat\Ahat;
B = Ehat\Bhat;

% matrix-matrix product (used twice)
% first row of Chat cut off
Ccut = Chat(2:end,:);
M = Ccut'*Ccut;

end

