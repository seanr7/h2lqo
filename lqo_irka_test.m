% Test lqo-irka
clc
clear all
close all

%% example from lqo Overleaf doc

A = -diag([1,2,3,4]);    b = [1;1;2;2];  c = [2,2,3,3];
K = [20, -10, 0, 0;
    -10, 20, -10, 0;
    0, -10, 20, -10;
    0, 0, -10, 20];

n = 4;  r = 2; % 
lambda_prev = rand(r,1);  phi_prev = rand(r,1);   kappa_prev = rand(r,r);
[Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(A, b, c, K, lambda_prev, phi_prev, kappa_prev, 100, 10e-16, 1);

% Are interpolation conditions satisfied?
I = eye(n,n);   Ir = eye(r,r);
H1 = @(s) c*((s*I-A)\b); 
H2 = @(s1, s2) b'*((s1*I-A)'\K)*((s2*I-A)\b);

H1r = @(s) cr*((s*Ir-Ar)\br); 
H2r = @(s1, s2) br'*((s1*Ir-Ar)'\Kr)*((s2*Ir-Ar)\br); 

l1 = lambda(1); l2 = lambda(2);
fprintf('1st-order optimality conditions, 2r in total')
for i = 1:r
    [H1(-lambda(1)) - H1r(-lambda(1))]
end

fprintf('2nd-order optimality conditions, r^2 in total')
for i = 1:r
    for j = 1:r
        [H2(-lambda(i),-lambda(j)) - H2r(-lambda(i),-lambda(j))]
    end
    
end

% Derivatives?
H1_prime = @(s) -c*((s*I-A)\((s*I-A)\b));
H1r_prime = @(s) -cr*((s*Ir-Ar)\((s*Ir-Ar)\br));

H2_s1prime = @(s1,s2) -b'*((s1*I-A)'\((s1*I-A)'\K))*((s2*I-A)\b); % Partial deriv w.r.t s1
H2r_s1prime = @(s1,s2) -br'*((s1*Ir-Ar)'\((s1*Ir-Ar)'\Kr))*((s2*Ir-Ar)\br);   % Partial deriv w.r.t s1

H2_s2prime = @(s1,s2) -b'*((s1*I-A)'\K)*((s2*I-A)\((s2*I-A)\b)); % Partial deriv w.r.t s2
H2r_s2prime = @(s1,s2) -br'*((s1*Ir-Ar)'\Kr)*((s2*Ir-Ar)\((s2*Ir-Ar)\br)); % Partial deriv w.r.t s2

% So, strictly derivatives are not interpolated... But are the mixed
% conditions satisfied? 

fprintf('Mixed linear + quadratic optimality conditions, 2r in total')

for i = 1:r
    ROpart = 2*phi(i)*H1r_prime(-lambda(i));
    FOpart = 2*phi(i)*H1_prime(-(lambda(i)));
    for j = 1:r
        ROpart = ROpart + kappa(i,j)*H2r_s1prime(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2r_s2prime(-lambda(j),-lambda(i));
        FOpart = FOpart + kappa(i,j)*H2_s1prime(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2_s2prime(-lambda(j),-lambda(i));
    end
    [FOpart - ROpart]
end

% fprintf('Linear derivative interpolation?')
% [H1_prime(-l1)-H1r_prime(-l1), H1_prime(-l2)-H1r_prime(-l2)]
% 
% fprintf('Quadratic derivative interpolation? H2^(1,0)')
% [H2_s1prime(-l1,-l2)-H2r_s1prime(-l1,-l2)]
% [H2_s1prime(-l2,-l1)-H2r_s1prime(-l2,-l1)]
% 
% fprintf('Quadratic derivative interpolation? H2^(0,1)')
% [H2_s2prime(-l1,-l2)-H2r_s2prime(-l1,-l2)]
% [H2_s2prime(-l2,-l1)-H2r_s2prime(-l2,-l1)]

%% a random example...
n = 10;
r = 4;

% TODO: Somewhere under the hood, I must be assuming implicitly that A is a
% diagonal matrix at the start, because the code fails if that is not the
% case...

A = diag(rand(n,1));
b = rand(n,1);  c = rand(1,n);  
K = rand(n,n);  K = (K+K')/2;

lambda_prev = rand(r,1);  phi_prev = rand(r,1);   kappa_prev = rand(r,r);
[Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(A, b, c, K, lambda_prev, phi_prev, kappa_prev, 100, 10e-12, 1);

% Are interpolation conditions satisfied?
I = eye(n,n);   Ir = eye(r,r);
H1 = @(s) c*((s*I-A)\b); 
H2 = @(s1, s2) b'*((s1*I-A)'\K)*((s2*I-A)\b);

H1r = @(s) cr*((s*Ir-Ar)\br); 
H2r = @(s1, s2) br'*((s1*Ir-Ar)'\Kr)*((s2*Ir-Ar)\br); 

l1 = lambda(1); l2 = lambda(2);
fprintf('1st-order optimality conditions, 2r in total')
for i = 1:r
    [H1(-lambda(1)) - H1r(-lambda(1))]
end

fprintf('2nd-order optimality conditions, r^2 in total')
for i = 1:r
    for j = 1:r
        [H2(-lambda(i),-lambda(j)) - H2r(-lambda(i),-lambda(j))]
    end
    
end

% Derivatives?
H1_prime = @(s) -c*((s*I-A)\((s*I-A)\b));
H1r_prime = @(s) -cr*((s*Ir-Ar)\((s*Ir-Ar)\br));

H2_s1prime = @(s1,s2) -b'*((s1*I-A)'\((s1*I-A)'\K))*((s2*I-A)\b); % Partial deriv w.r.t s1
H2r_s1prime = @(s1,s2) -br'*((s1*Ir-Ar)'\((s1*Ir-Ar)'\Kr))*((s2*Ir-Ar)\br);   % Partial deriv w.r.t s1

H2_s2prime = @(s1,s2) -b'*((s1*I-A)'\K)*((s2*I-A)\((s2*I-A)\b)); % Partial deriv w.r.t s2
H2r_s2prime = @(s1,s2) -br'*((s1*Ir-Ar)'\Kr)*((s2*Ir-Ar)\((s2*Ir-Ar)\br)); % Partial deriv w.r.t s2

% So, strictly derivatives are not interpolated... But are the mixed
% conditions satisfied? 

fprintf('Mixed linear + quadratic optimality conditions, 2r in total')

for i = 1:r
    ROpart = 2*phi(i)*H1r_prime(-lambda(i));
    FOpart = 2*phi(i)*H1_prime(-(lambda(i)));
    for j = 1:r
        ROpart = ROpart + kappa(i,j)*H2r_s1prime(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2r_s2prime(-lambda(j),-lambda(i));
        FOpart = FOpart + kappa(i,j)*H2_s1prime(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2_s2prime(-lambda(j),-lambda(i));
    end
    [FOpart - ROpart]
end

%% lqo-irka

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
end