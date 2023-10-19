function [Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(A, b, c, K, ...
    lambda_prev, phi_prev, kappa_prev, itermax, eps, plot_conv)
% Author: Sean Reiter
% Date Modified: 4-12-23
% 
% Please report any issues with this code to seanr7@vt.edu
%
% Function to implement an iterative rational krylov algorithm (IRKA) for
% H2-optimal interpolatory MOR of LTI systems with linear quadratic output
% (LQO). 
%
% Inputs:
% :param (A, b, c, K): ((n,n), (n,1), (1,n), (n,n)) State-space 
%                      realization of the LQO FOM.
% :param lambda_prev : Initial interpolation points, (r,1) array
% :param phi_prev    : Initial 1st-order residues, (r,1) array
% :param kappa_prev  : Initial 2nd-order residues, (r,r) array
% :param itermax     : Max no. of iterations to run
% :param eps         : Convergence tol
% :param plot_conv   : Bool; do you want to plot convergence of the 
%                           poles? 
% 
% Outputs:
% :param (Ar, br, cr, Kr): ((r,r), (r,1), (1,r), (r,r)) State-space 
%                           realization of the LQO-IRKA ROM.
% :param lambda          : Converged interpolation points, (r,1) array
% :param phi             : Converged 1st-order residues, (r,1) array
% :param kappa           : Converged 2nd-order residues, (r,r) array

%%
% Build initial projection bases 
r = length(lambda_prev);     [~,n] = size(A);
if nargin == 9 % Default is not to plot
    plot_conv = 0;
end
if nargin == 8 % If no tolerance set
    eps = 10e-12;
end
if nargin == 7 % If no max iterations set
    itermax = 100;
end

I = eye(n,n);
iter = 1;   err(iter) = eps + 1; % So while loop initiates

while err(iter) > eps && iter < itermax
    % In building Wr/Vr, need reduced poles AND residues
    % Fill out left/right projection matrix Wr/Vr
    Vr = zeros(n,r);     Wr = zeros(n,r);
    % TODO: This part below can be optimized so as to not repeat any linear
    % solves 
    for k = 1:r
        % Right
        Vr(:,k) = (-lambda_prev(k)*I-A)\b;     
        % Left
        Wr(:,k) = phi_prev(k)*((-lambda_prev(k)*I-A')\c');
        tmp = (-lambda_prev(k)*I-A')\K';
        for i = 1:r
            Wr(:,k) = Wr(:,k) + kappa_prev(i,k)*tmp*((-lambda_prev(i)*I-A)\b);
        end
    end
    [Vr,~] = qr(Vr, "econ");     [Wr,~] = qr(Wr, "econ");
    % Compute ROM, check for convergence
    Wrt = Wr';  Vrt = Vr';
    Ar = (Wrt*Vr)\(Wrt*A*Vr);   br = (Wrt*Vr)\(Wrt*b);
    cr = c*Vr;  Kr = Vrt*K*Vr;
    % New interpolation points + residues
    [Xr, Lr] = eig(Ar);     lambda = diag(Lr);
    % 1st-order residues
    phi = (cr*Xr)'.*(Xr\br);
    % 2nd-order residues
    kappa = Xr'*Kr*Xr;

    % Track convergence of poles 
    iter = iter + 1; 
    err(iter) = max(abs(lambda - lambda_prev));
    % Overwrite reduced-order poles/residues from previous iteration
    lambda_prev = lambda;   phi_prev = phi;    kappa_prev = kappa;
end
if plot_conv == true
    semilogy([1:iter], err, '-o', LineWidth=1.5)
    xlim([1,iter])
    xlabel('Iteration count')
    ylabel('Magnitude of change in \lambda')
end
e      nd