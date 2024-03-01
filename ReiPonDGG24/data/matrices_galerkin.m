function [Agal,Bgal,Cgal,Egal] = matrices_galerkin(degr,sigma)
% (semi-)analytical computation of matrix in Galerkin system

% internal parameters
n = 4;           % number of oscillators
ndim = 2*n;      % dimension of dynamical system
npar = 14;       % number of parameters

% Parameters from article Lohmann/Eid
mass = [1 5 25 125]';
spring = [27 9 3 1 2 3]';
damper = [0.1 0.4 1.6 1.]';
pmean = [mass; spring; damper;];

% turn off warnings
warning off

% initialisation of arrays
Acell = cell(1,npar);
Bcell = cell(1,npar);
Ccell = cell(1,npar);
Ecell = cell(1,npar);

% computation of mean values
[Amean,Bmean,Cmean,Emean] = matrices_springdamper(n,pmean);

% computation of differnce matrices
eta = 1e-8;
for k = 1:npar
   
   pvec = pmean;
   pvec(k) = pvec(k) + eta;
   
   [deltaA,deltaB,deltaC,deltaE] = matrices_springdamper(n,pvec);
   
   Acell{k} = sigma*pmean(k)*(deltaA - Amean)/eta; 
   Bcell{k} = sigma*pmean(k)*(deltaB - Bmean)/eta;
   Ccell{k} = sigma*pmean(k)*(deltaC - Cmean)/eta; 
   Ecell{k} = sigma*pmean(k)*(deltaE - Emean)/eta;
   
end

% table from Gauss-Legendre quadrature
k = 10;   % number of nodes
d = 3;    % degree of polynomials
[nodes,weights,bpoly] = gauss_legendre_1d(d,k);
Smat = zeros(d+1);
for j1 = 1:(d+1)
   for j2 = 1:(d+1)
      for l = 1:k
         Smat(j1,j2) = Smat(j1,j2) ...
	      + weights(l)*nodes(l)*bpoly(j1,l)*bpoly(j2,l);
      end
      if (abs(Smat(j1,j2))<1e-14)
         Smat(j1,j2) = 0.;
      end
   end
end

% determination of multi-indices: yields variable mindex
mindex = multiindex(npar,degr);
% determination of number of basis polynomials
nbpoly = size(mindex,1);
Agal = sparse(nbpoly*ndim,nbpoly*ndim);
Bgal = sparse(nbpoly*ndim,nbpoly);
Cgal = sparse(nbpoly,nbpoly*ndim);
Egal = sparse(nbpoly*ndim,nbpoly*ndim);

% loop over parameters
for k = 1:npar
   
   % matrix < p_k phi_i phi_j >
   Tmat = zeros(nbpoly,nbpoly);
   for i = 1:nbpoly
      for j = 1:nbpoly
	     ind = 0.;
         if (i~=j)
            for l = 1:npar
               if (l~=k)
                  if (mindex(i,2+l)~=mindex(j,2+l))
                     ind = ind + 1;
                  end
               end
            end
         end
	     if (ind==0)
            Tmat(i,j) = Smat(mindex(i,2+k)+1,mindex(j,2+k)+1);
         end
      end
   end
   Tmat = sparse(Tmat);
   
   % add Kronecker product
   deltaA = Acell{k};
   Agal = Agal + sparse(kron(Tmat,sparse(deltaA)));
   deltaB = Bcell{k};
   Bgal = Bgal + sparse(kron(Tmat,sparse(deltaB)));
   deltaC = Ccell{k};
   Cgal = Cgal + sparse(kron(Tmat,sparse(deltaC)));
   deltaE = Ecell{k};
   Egal = Egal + sparse(kron(Tmat,sparse(deltaE)));
   
end

% finalize matrices - add block diagonal
Agal = kron(speye(nbpoly),sparse(Amean)) + Agal;
Bgal = kron(speye(nbpoly),sparse(Bmean)) + Bgal;
Bgal = Bgal(:,1);
Cgal = kron(speye(nbpoly),sparse(Cmean)) + Cgal;
Egal = kron(speye(nbpoly),sparse(Emean)) + Egal;


return
end