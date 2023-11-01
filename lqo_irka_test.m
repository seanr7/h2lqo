% Test lqo-irka
clc
clear all
close all

%%
addpath('/Users/seanr/Desktop/h2lqo/benchmarks/')
load('heat-cont.mat')
A = full(A);    B = full(B);    C = full(C);
[n,~] = size(A);    p = 1; % relevant dimensions
E = eye(n, n);
b = B(:, p);  c = C(p, :); % If using ISS model, make SISO
K = (diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1));

r = 10; % 
lambda_prev = -logspace(-2,4,r)';  phi_prev = rand(r,1);   tmp = rand(r,r);
kappa_prev = (tmp+tmp')/2;
% lambda_prev = -10*rand(r,1);  phi_prev = rand(r,1);   kappa_prev = rand(r,r);
[Er, Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(E, A, b, c, K, ...
    lambda_prev, phi_prev, kappa_prev, 100, 10e-8, 1);

%% H2 errors
addpath('/Users/seanr/Documents/MATLAB/Research/QuadBT/code4Sean_Oct10th_2022',...
    '/Users/seanr/Documents/MATLAB/Research/QuadBT/Quad-PRBT')

load('heat-cont.mat')
% load('iss1R.mat')
A = full(A);    B = full(B);    C = full(C);
[n,~] = size(A);    p = 1; % relevant dimensions
E = eye(n, n);
b = B(:, p);  c = C(p, :); % If using ISS model, make SISO
M = (diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1));
M = M*10e-2;

testcases = 10;
lambda_init = -10*rand(testcases,1);  phi_init = rand(testcases,1); 
tmp = rand(testcases,testcases);
kappa_init = (tmp+tmp')/2;
count = 1;
for r = 2:2:testcases
    figure 
    [Ar, br, cr, Mr, lambda, phi, kappa] = lqo_irka(A, b, c, M, lambda_init(1:r), phi_init(1:r), kappa_init(1:r,1:r), 100, 10e-8, 1);
    % Compute H2 error
    Aerr_irka = blkdiag(A, Ar);
    berr_irka = [b; br];
    cerr_irka = [c, -cr];
    Merr_irka = blkdiag(M, -Mr);
    Perr_irka = lyap(Aerr_irka, berr_irka*berr_irka');
    Qerr_irka = lyap(Aerr_irka', cerr_irka'*cerr_irka + Merr_irka*Perr_irka*Merr_irka);
    H2_errors_irka(count) = sqrt(berr_irka'*Qerr_irka*berr_irka);
% 
% 
%     P_qbt = lyap(A,b*b'); % OF P
%     U = chol(P_qbt+eye(n,n)*10e-16);
%     Q_qbt = lyap(A',c'*c + M*P_qbt*M); % OF Q
%     L = chol(Q_qbt+eye(n,n)*10e-16);
%     
%     [Z,S,Y] = svd(U'*L2);
%     hsv = diag(S);
%     
%     Z1 = Z(:,1:r);
%     Y1 = Y(:,1:r);
%     S1 = S(1:r,1:r);
%     
%     W = L*Y1/(sqrt(S1));
%     V = U*Z1/(sqrt(S1));
%     
%     Ar = W'*A*V;
%     br = W'*b;
%     cr = c*V;
%     Mr = V'*M*V;
%     Aerr_qbt = blkdiag(A, Ar);
%     berr_qbt = [b; br];
%     cerr_qbt = [c, -cr];
%     Merr_qbt = blkdiag(M, -Mr);
%     Perr_qbt = lyap(Aerr_qbt, berr_qbt*berr_qbt');
%     Qerr_qbt = lyap(Aerr_qbt', cerr_qbt'*cerr_qbt + Merr_qbt*Perr_qbt*Merr_qbt);
%     H2_errors_qbt(count) = sqrt(berr_qbt'*Qerr_qbt*berr_qbt);
%     [Ar, br, cr, Mr, interp_pts] = twosided_lqo(A, b, c, M, lambda_init(1:r), phi_init(1:r), kappa_init(1:r,1:r), 10e-8);
%     Aerr_tsia = blkdiag(A, Ar);
%     berr_tsia = [b; br];
%     cerr_tsia = [c, -cr];
%     Merr_tsia = blkdiag(M, -Mr);
%     Perr_tsia = lyap(Aerr_tsia, berr_tsia*berr_tsia');
%     Qerr_tsia = lyap(Aerr_tsia', cerr_tsia'*cerr_tsia + Merr_tsia*Perr_tsia*Merr_tsia);
%     H2_errors_tsia(count) = sqrt(berr_tsia'*Qerr_tsia*berr_tsia);
    count = count + 1;
end
figure
golden_ratio = (sqrt(5)+1)/2;
axes('position', [.125 .15 .75 golden_ratio-1])
semilogy([2:2:testcases], H2_errors_irka, '-o','markersize', 10, 'linewidth',2)
hold on
% semilogy([2:2:testcases], H2_errors_qbt, '-*','markersize', 10, 'linewidth',2)
xlabel('Order of reduction, $r$','FontSize',14,'Interpreter','latex')
ylabel('$\|\mathcal{S}-\mathcal{S}_r\|_{\mathcal{H}_2}$','FontSize',14, 'Interpreter','latex')
legend('LQO-IRKA', 'LQO-TSIA')

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

% function [Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(A, b, c, K, ...
%     lambda_prev, phi_prev, kappa_prev, itermax, eps, plot_conv)
% % Author: Sean Reiter
% % Date Modified: 4-12-23
% % 
% % Please report any issues with this code to seanr7@vt.edu
% %
% % Function to implement an iterative rational krylov algorithm (IRKA) for
% % H2-optimal interpolatory MOR of LTI systems with linear quadratic output
% % (LQO). 
% %
% % Inputs:
% % :param (A, b, c, K): ((n,n), (n,1), (1,n), (n,n)) State-space 
% %                      realization of the LQO FOM.
% % :param lambda_prev : Initial interpolation points, (r,1) array
% % :param phi_prev    : Initial 1st-order residues, (r,1) array
% % :param kappa_prev  : Initial 2nd-order residues, (r,r) array
% % :param itermax     : Max no. of iterations to run
% % :param eps         : Convergence tol
% % :param plot_conv   : Bool; do you want to plot convergence of the 
% %                           poles? 
% % 
% % Outputs:
% % :param (Ar, br, cr, Kr): ((r,r), (r,1), (1,r), (r,r)) State-space 
% %                           realization of the LQO-IRKA ROM.
% % :param lambda          : Converged interpolation points, (r,1) array
% % :param phi             : Converged 1st-order residues, (r,1) array
% % :param kappa           : Converged 2nd-order residues, (r,r) array
% 
% %%
% % Build initial projection bases 
% r = length(lambda_prev);     [~,n] = size(A);
% if nargin == 9 % Default is not to plot
%     plot_conv = 0;
% end
% if nargin == 8 % If no tolerance set
%     eps = 10e-12;
% end
% if nargin == 7 % If no max iterations set
%     itermax = 100;
% end
% 
% I = eye(n,n);
% iter = 1;   err(iter) = eps + 1; % So while loop initiates
% 
% while err(iter) > eps && iter < itermax
%     % In building Wr/Vr, need reduced poles AND residues
%     % Fill out left/right projection matrix Wr/Vr
%     Vr = zeros(n,r);     Wr = zeros(n,r);
%     % TODO: This part below can be optimized so as to not repeat any linear
%     % solves 
%     for k = 1:r
%         % Right
%         Vr(:,k) = (-lambda_prev(k)*I-A)\b;     
%         % Left
%         Wr(:,k) = phi_prev(k)*((-lambda_prev(k)*I-A')\c');
%         tmp = (-lambda_prev(k)*I-A')\K';
%         for i = 1:r
%             Wr(:,k) = Wr(:,k) + kappa_prev(i,k)*tmp*((-lambda_prev(i)*I-A)\b);
%         end
%     end
%     [Vr,~] = qr(Vr, "econ");     [Wr,~] = qr(Wr, "econ");
%     % Compute ROM, check for convergence
%     Wrt = Wr';  Vrt = Vr';
%     Ar = (Wrt*Vr)\(Wrt*A*Vr);   br = (Wrt*Vr)\(Wrt*b);
%     cr = c*Vr;  Kr = Vrt*K*Vr;
%     % New interpolation points + residues
%     [Xr, Lr] = eig(Ar);     lambda = diag(Lr);
%     % 1st-order residues
%     phi = (cr*Xr)'.*(Xr\br);
%     % 2nd-order residues
%     kappa = Xr'*Kr*Xr;
% 
%     % Track convergence of poles 
%     iter = iter + 1; 
%     err(iter) = max(abs(lambda - lambda_prev));
%     % Overwrite reduced-order poles/residues from previous iteration
%     lambda_prev = lambda;   phi_prev = phi;    kappa_prev = kappa;
% end
% if plot_conv == true
%     golden_ratio = (sqrt(5)+1)/2;
%     axes('position', [.125 .15 .75 golden_ratio-1])
%     semilogy([1:iter], err, '-o','linewidth',2)
%     xlim([1,iter])
%     xlabel('Iteration count','FontSize',14,'Interpreter','latex')
%     ylabel('Magnitude of change in $\lambda(\mathbf{A})$','FontSize',14, 'Interpreter','latex')
% end
% end

%% TSIA
function [Ar, br, cr, Kr, interp_pts] = twosided_lqo(A, b, c, K, lambda, phi, kappa, eps)
% Author: Sean Reiter
% Date Modified: 3-27-23
% 
% Please report any issues with this code to seanr7@vt.edu
%
% Function to implement two-sided iterative algorithm for MOR of LQO
% systems in: [Gosea/Antoulas, '19]: `` A two-sided iterative framework for
% model reduction of linear systems with quadratic output''.
%
% Inputs:
% :param (A, b, c, K): ((n,n), (n,1), (1,n), (n,n)) State-space realization 
%                      of LQO FOM.
% :param lambda      : Left interpolation points
% :param mu          : Right interpolation points
% :param eps         : Convergence tol

%%
% Build initial projection bases 
r = length(lambda);     [~,n] = size(A);
% if mod(r,2) ~= 0 % For now, assume even no. of interpolation points
%     error('An even number if interpolation points is assumed')
% end
% if r ~= length(mu)
%     error('Must have an equal number of left/right interpolation points')
% end

% Vr = zeros(n,r);    Wr = zeros(n,r);
I = eye(n,n);
% for k = 1:(r/2)
%     Vr(:,2*(k-1)+1) = (lambda(2*(k-1)+1)*I-A)\b;     Vr(:,2*k) = (lambda(2*k)*I-A)\b; 
%     Wr(:,2*(k-1)+1) = (mu(2*(k-1)+1)*I-A')\c';       Wr(:,2*k) = ((mu(2*k)*I-A')\K)*((lambda(2*(k-1)+1)*I-A)\b); 
% end
% % Initial projection step 
% Wrt = Wr';  Vrt = Vr';
% Ar = (Wrt*Vr)\(Wrt*A*Vr);   br = (Wrt*Vr)\(Wrt*b);
% cr = c*Vr;  Kr = Vrt*K*Vr;
% Initialise
Vr = zeros(n,r);     Wr = zeros(n,r);
% TODO: This part below can be optimized so as to not repeat any linear
% solves 
for k = 1:r
    % Right
    Vr(:,k) = (-lambda(k)*I-A)\b;     
    % Left
    Wr(:,k) = phi(k)*((-lambda(k)*I-A')\c');
    tmp = (-lambda(k)*I-A')\K';
    for i = 1:r
        Wr(:,k) = Wr(:,k) + kappa(i,k)*tmp*((-lambda(i)*I-A)\b);
    end
end
[Vr,~] = qr(Vr, "econ");     [Wr,~] = qr(Wr, "econ");
% Compute ROM, check for convergence
Wrt = Wr';  Vrt = Vr';
Ar = (Wrt*Vr)\(Wrt*A*Vr);   br = (Wrt*Vr)\(Wrt*b);
cr = c*Vr;  Kr = Vrt*K*Vr;
lambda = eig(Ar);

% Iterate until eigenvalues converge 
interp_pts = [lambda];

while max(abs(lambda-interp_pts)) > eps
    Arprev = Ar;
    Brprev = Br;
    Crprev = Cr;
    Mrprev = Mr;
    
    X = sylvester(A,Ar',-B*Br');
    Y = sylvester(A',Ar,-M'*X*Mr-C'*Cr);

    %V = orth(X);
    %W = orth(Y);
    V = X;
    W = Y;
    W = (W'*V)\W';

    Ar = W*A*V;
    Br = W*B;
    Cr = C*V;
    Mr = V'*M*V;
    
    [U,V] = eig(Ar);
    Kr = Mr(:)';


    Ar = U^(-1)*Ar*U;
    Br = U^(-1)*Br;
    Cr = Cr*U;
    Kr = Kr*kron(U,U);
 
    Db = diag(Br);
    Ar = Db^(-1)*Ar*Db;
    br = Db^(-1)*Br;
    cr = Cr*Db;
    Kr = Kr*kron(Db,Db);

    Kr = reshape(Kr,2*r,2*r);
    interp_pts = eig(Ar);
end

% while max(abs(lambda-interp_pts)) > eps
%     lambda = eig(Ar);
%     % Solve + orthonormalize projection matrices
%     X = sylvester(A, Ar, -b*br');    Y = sylvester(A', Ar, -(K'*X*Kr + c'*cr));
%     [Vr,~] = qr(X);     [Wr,~] = qr(Y);
%     % Compute ROM, check for convergence
%     Wrt = Wr';  Vrt = Vr';
%     Ar = (Wrt*Vr)\(Wrt*A*Vr);   br = (Wrt*Vr)\(Wrt*b);
%     cr = c*Vr;  Kr = Vrt*K*Vr;
%     interp_pts = eig(Ar);
% end
end