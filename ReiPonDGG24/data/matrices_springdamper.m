function [A,B,C,E] = matrices_springdamper(n,p)
% matrices of linear dynamical system for mass-spring-damper
%    E(p) * x' = A(p) * x + B(p) * u
%           y  = C * x
% The dimension of the system is 2*n.
% The system is single-input-single-output.
% Output is the position of the last mass.
%
% The equations read as
%  m_1 z_1'' = d_2 (z_2'-z_1') 
%              + c_1 (u-z_1) + c_2 (z_2-z_1) + c_{n+1} (z_n-z_1)
%  m_j z_j'' = d_j (z_{j-1}'-z_j') + d_{j+1} (z_{j+1}'-z_j') 
%              + c_j (z_{j-1}-z_j) + c_{j+1} (z_{j+1}-z_j)   for j=2,...,n-1
%  m_n z_n'' = d_n (z_{n-1}'-z_n') - d_{n+1} z_n' 
%              + c_n (z_{n-1}-z_n) + c_{n+1} (z_1-z_n) - c_{n+2} z_n.
% The second-order ODEs are transformed into a first-order system.
%
% The physical parameters are arranged in the vector p:
%  (m_1,...,m_n,c_1,...,c_{n+2},d_2,...,d_{n+1})
% The number of parameters is 3*n+2.
%
% Literature: B. Lohmann, R. Eid, Efficient order reduction of parametric
% and nonlinear models by superposition of locally reduced models, 2009.

% rearrange parameters
[mp,np] = size(p);
if ((mp~=(3*n+2))|(np~=1))
   disp('Wrong dimension of parameter vector!')
   A = []; B =[]; C = []; E = [];
   return
end
mass = p(1:n);
spring = p((n+1):(2*n+2));
damper = zeros(n+1,1);
damper(2:end) = p((2*n+3):end);

% definition of second-order system M * w'' = S * w' + T * w
T = zeros(n);
T(1,1) = - spring(1) - spring(2) - spring(n+1);
T(1,2) = spring(2);
T(1,n) = spring(n+1);
for j = 2:(n-1)
   T(j,j-1) = spring(j);
   T(j,j) = - spring(j) - spring(j+1);
   T(j,j+1) = spring(j+1);
end
T(n,1) = spring(n+1);
T(n,n-1) = spring(n);
T(n,n) = - spring(n) - spring(n+1) - spring(n+2);

S = zeros(n);
S(1,1) = - damper(2);
S(1,2) = damper(2);
for j = 2:(n-1)
   S(j,j-1) = damper(j);
   S(j,j) = - damper(j) - damper(j+1);
   S(j,j+1) = damper(j+1);
end
S(n,n-1) = damper(n);
S(n,n) = - damper(n) - damper(n+1);

M = diag(mass);

% matrix A in dynamical system
A = zeros(2*n);
A(1:n,(n+1):(2*n)) = eye(n);
A((n+1):(2*n),1:n) = T;
A((n+1):(2*n),(n+1):(2*n)) = S;

% input matrix
B = zeros(2*n,1);
B(n+1) = spring(1);

% output matrix
C = zeros(1,2*n);
C(n) = 1.;

% matrix E in dynamical system
E = eye(2*n);
E((n+1):(2*n),(n+1):(2*n)) = M;


return
end

