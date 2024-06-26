% Unit tests for LQO IRKA
clc
clear all
close all

%% 2nd-order -> 1st-order test probel
n1 = 10; alpha=.002; beta=alpha; v = 5;

[M, D, K]=triplechain_MSD(n1, alpha, beta, v);

M = full(M);
D = full(D);
K = full(K);

O  = zeros(size(K,1),1);
Cv = O';
Cp = ones(1,size(K,1));
B  = ones(size(K,1),1);
C  = [Cp, Cv];

[n, ~] = size(M);

E_qo = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is (sparse) mass matrix

A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(D) + n);  % A_qo = [0, I; -K, -D]
A_qo(1:n, n+1:2*n) = speye(n); 
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -D; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B;
% No scalar output in this example; only QO

C_qo = zeros(1, 2*n);
% Our QO (Quadratic-output) matrix Q
Q_qo = blkdiag(eye(n), zeros(n));

%%

% r = 4;
% interp_pts = -10 * rand(r, 1);   FO_res = rand(r, 1);    tmp = rand(r, 1);
% SO_res = (tmp + tmp')/2;
% max_iter = 100; tol = 10e-8;   plot = true;

% Run iteration + plot conv
% [E_qo_r, A_qo_r, B_qo_r, ~, Q_qo_r, conv_nodes, conv_FO_res, conv_SO_res, pole_history] = lqo_irka(E_qo, A_qo, B_qo, [], Q_qo, ...
%     interp_pts, FO_res, SO_res, max_iter, tol, plot);

addpath("ReiW24\drivers\")

r = 4;
[E_qo_r, A_qo_r, B_qo_r, ~, Q_qo_r, info] = sisolqo_irka(E_qo, A_qo, B_qo, [], Q_qo, r);
conv_nodes = info.poles;    conv_SO_res = info.sores;
% Check interpolation conditions
H2 = @(s1, s2) B_qo' * (((s1 * E_qo - A_qo)')\Q_qo) * ((s2 * E_qo - A_qo)\B_qo);
H2r = @(s1, s2) B_qo_r' * (((s1 * E_qo_r - A_qo_r)')\Q_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r);

fprintf('2nd-order optimality conditions, %d^2 in total', r)
for i = 1:r
    for j = 1:r
        H2(-(conv_nodes(i)),  -(conv_nodes(j))) - H2r(-(conv_nodes(i)), -(conv_nodes(j)))
    end
    
end
% 
% % Partial deriv wrt first argument of H2, H2r
H2_prime_s1 = @(s1, s2) -B_qo' * (((s1 * E_qo - A_qo)'\E_qo') * ((s1 * E_qo - A_qo)'\Q_qo)) * ((s2 * E_qo - A_qo)\B_qo); 
H2_prime_s2 = @(s1, s2) -B_qo' * ((s1 * E_qo - A_qo)'\Q_qo) * ((s2 * E_qo - A_qo)\E_qo) * ((s2 * E_qo - A_qo)\B_qo); 
H2r_prime_s1 = @(s1, s2) -B_qo_r' * (((s1 * E_qo_r - A_qo_r)'\E_qo_r') * ((s1 * E_qo_r - A_qo_r)'\Q_qo_r)) * ((s2 * E_qo_r- A_qo_r)\B_qo_r); 
H2r_prime_s2 = @(s1, s2) -B_qo_r' * ((s1 * E_qo_r - A_qo_r)'\Q_qo_r) * (((s2 * E_qo_r - A_qo_r)\E_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r)); 
% 
fprintf('Mixed linear + quadratic optimality conditions, %d in total', r)

for i = 1:r
    ro_side = 0;
    fo_side = 0;
    for j = 1:r
        ro_side = ro_side + 2 * conv_SO_res(j, i) * H2r_prime_s1(-conv_nodes(i), -conv_nodes(j));
        fo_side = fo_side + 2 * conv_SO_res(j, i) * H2_prime_s1(-conv_nodes(i), -conv_nodes(j));
    end
    fo_side - ro_side
end

%% 1a. Toy model
n = 8; 
% Random scalar + QO component
A = -10 * diag(rand(n, 1));
% A = [-1 + 1i, 0; 0, -1 - 1i];
% A = diag([-1 + 1i, -1 - 1i, -2 + 1i, -2 - 1i]);
b = rand(n, 1);     c = rand(1, n);
% QO matrix
M = 10*(diag(2 * ones(1, n)) + diag(-1 * ones(1, n-1), 1) + diag(-1 * ones(1, n - 1), -1));

% To find good interpolation points, look at the poles
r = 4;
interp_pts = -10 * rand(r, 1);   FO_res = rand(r, 1);    tmp = rand(r, 1);
SO_res = (tmp + tmp')/2;
max_iter = 100; tol = 10e-10;   plot = true;

% First, E = I case
E = eye(n, n);

% Run iteration + plot conv
[Er, Ar, br, cr, Mr, conv_nodes, conv_FO_res, conv_SO_res, V_r] = lqo_irka(E, A, b, c, M, ...
    interp_pts, FO_res, SO_res, max_iter, tol, plot);


% Check interpolation conditions
H1 = @(s) c * ((s * E - A)\b); 
H2 = @(s1, s2) b.' * ((s1 * E - A).'\M) * ((s2 * E - A)\b);

H1r = @(s) cr * ((s * Er - Ar)\br); 
H2r = @(s1, s2) br.' * ((s1 * Er - Ar).'\Mr) * ((s2 * Er - Ar)\br); 

fprintf('1st-order optimality conditions, %d in total', r)
for i = 1:r
    H1(-conv_nodes(i)) - H1r(-conv_nodes(i))
end

fprintf('2nd-order optimality conditions, %d^2 in total', r)
for i = 1:r
    for j = 1:r
        H2(-conv_nodes(i),  -conv_nodes(j)) - H2r(-conv_nodes(i), -conv_nodes(j))
    end
    
end

% Derivs of H1, H1r
H1_prime = @(s) -c * (((s * E - A)\E) * ((s * E - A)\b));
H1r_prime = @(s) -cr * (((s * Er - Ar)\Er) * ((s * Er - Ar)\br));

% Partial deriv wrt first argument of H2, H2r
H2_prime_s1 = @(s1, s2) -b.' * (((s1 * E - A).'\E.') * ((s1 * E - A).'\M)) * ((s2 * E - A)\b); 
H2_prime_s2 = @(s1, s2) -b.' * ((s1 * E - A).'\M) * ((s2 * E - A)\E) * ((s2 * E - A)\b); 
H2r_prime_s1 = @(s1, s2) -br.' * (((s1 * Er - Ar).'\Er.') * ((s1 * Er - Ar).'\Mr)) * ((s2 * Er - Ar)\br);   % Partial deriv w.r.t s1
H2r_prime_s2 = @(s1, s2) -br.' * ((s1 * Er - Ar).'\Mr) * (((s2 * Er - Ar)\Er) * ((s2 * Er - Ar)\br)); 


fprintf('Mixed linear + quadratic optimality conditions, %d in total', r)

for i = 1:r
    ro_side = conv_FO_res(i) * H1r_prime(-(conv_nodes(i)));
    fo_side = conv_FO_res(i) * H1_prime(-(conv_nodes(i)));
    for j = 1:r
        ro_side = ro_side + 2 * conv_SO_res(i, j) * H2r_prime_s1(-conv_nodes(i), -conv_nodes(j));
        fo_side = fo_side + 2 * conv_SO_res(i, j) * H2_prime_s1(-conv_nodes(i), -conv_nodes(j));
        % % 
        % % ro_side = ro_side + conv_SO_res(i, j) * H2r_prime_s1(-conv_nodes(i), -conv_nodes(j)) + ...
        % %             conv_SO_res(j, i) * H2r_prime_s2(-conv_nodes(j), -conv_nodes(i));
        % % fo_side = fo_side + conv_SO_res(i, j) * H2_prime_s1(-conv_nodes(i), -conv_nodes(j)) + ...
        % %             conv_SO_res(j, i) * H2_prime_s2(-conv_nodes(j), -conv_nodes(i));
    end
    fo_side - ro_side
end

% If everything stays real, the above works... - SR (1-30-24)
% conv_nodes
%% Now E ~= I, but nonsingular 
E = diag(rand(n, 1));
% x = -2 + 1i*[1, 2, 3, 4];
% A = diag([x, conj(x)]);
% Run iteration + plot conv
[Er, Ar, br, cr, Mr, conv_nodes, conv_FO_res, conv_SO_res] = lqo_irka(E, A, b, c, M, ...
    interp_pts, FO_res, SO_res, max_iter, tol, plot);

% Check interpolation conditions; redefine functions 
H1 = @(s) c * ((s * E - A)\b); 
H2 = @(s1, s2) b.' * ((s1 * E - A).'\M) * ((s2 * E - A)\b);

H1r = @(s) cr * ((s * Er - Ar)\br); 
H2r = @(s1, s2) br.' * ((s1 * Er - Ar).'\Mr) * ((s2 * Er - Ar)\br); 

fprintf('1st-order optimality conditions, %d in total', r)
for i = 1:r
    H1(-(conv_nodes(i))) - H1r(-(conv_nodes(i)))
end

fprintf('2nd-order optimality conditions, %d^2 in total', r)
for i = 1:r
    for j = 1:r
        H2(-(conv_nodes(i)),  -(conv_nodes(j))) - H2r(-(conv_nodes(i)), -(conv_nodes(j)))
    end
    
end

% Derivs of H1, H1r
H1_prime = @(s) -c * (((s * E - A)\E) * ((s * E - A)\b));
H1r_prime = @(s) -cr * (((s * Er - Ar)\Er) * ((s * Er - Ar)\br));

% Partial deriv wrt first argument of H2, H2r
H2_prime_s1 = @(s1, s2) -b' * (((s1 * E - A)'\E') * ((s1 * E - A)'\M)) * ((s2 * E - A)\b); 
H2_prime_s2 = @(s1, s2) -b' * ((s1 * E - A)'\M) * ((s2 * E - A)\E) * ((s2 * E - A)\b); 
H2r_prime_s1 = @(s1, s2) -br' * (((s1 * Er - Ar)'\Er') * ((s1 * Er - Ar)'\Mr)) * ((s2 * Er - Ar)\br);   % Partial deriv w.r.t s1
H2r_prime_s2 = @(s1, s2) -br' * ((s1 * Er - Ar)'\Mr) * (((s2 * Er - Ar)\Er) * ((s2 * Er - Ar)\br)); 



fprintf('Mixed linear + quadratic optimality conditions, %d in total', r)

for i = 1:r
    ro_side = conv_FO_res(i) * H1r_prime(-(conv_nodes(i)));
    fo_side = conv_FO_res(i) * H1_prime(-(conv_nodes(i)));
    for j = 1:r
        % Commented out code works if real-valued matrices
        % ro_side = ro_side + 2 * conv_SO_res(i, j) * H2r_prime_s1(-conv_nodes(i), -conv_nodes(j));
        % fo_side = fo_side + 2 * conv_SO_res(i, j) * H2_prime_s1(-conv_nodes(i), -conv_nodes(j));
        ro_side = ro_side + conv_SO_res(i, j) * H2r_prime_s1(-conv_nodes(i), -conv_nodes(j)) ...
            + conv_SO_res(j, i) * H2r_prime_s2(-conv_nodes(j), -conv_nodes(i));
        fo_side = fo_side + conv_SO_res(i, j) * H2_prime_s1(-conv_nodes(i), -conv_nodes(j)) ...
            + conv_SO_res(j, i) * H2_prime_s2(-conv_nodes(j), -conv_nodes(i));
    end
    fo_side - ro_side
end


% eig(Er\Ar)
%% 
addpath('/Users/seanr/Desktop/h2lqo/benchmarks/')
load('heat-cont.mat')
A = full(A);    b = full(B);    C = full(C);
[n,~] = size(A);    p = 1; % relevant dimensions
E = eye(n, n);
b = b(:, p);  c = C(p, :); % If using ISS model, make SISO
M = (diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1));

r = 10; % 
interp_pts = -logspace(-3,5,r)';  FO_res = rand(r,1);   tmp = rand(r,r);
SO_res = (tmp+tmp')/2;
[Er, Ar, br, cr, Mr, conv_nodes, conv_FO_res, conv_SO_res] = lqo_irka(E, A, b, c, M, ...
    interp_pts, FO_res, SO_res, max_iter, tol, plot);

%% H2 errors
addpath('/Users/seanr/Documents/MATLAB/Research/QuadBT/code4Sean_Oct10th_2022',...
    '/Users/seanr/Documents/MATLAB/Research/QuadBT/Quad-PRBT')

load('heat-cont.mat')
% load('iss1R.mat')
A = full(A);    b = full(B);    C = full(C);
[n,~] = size(A);    p = 1; % relevant dimensions
E = eye(n, n);
b = b(:, p);  c = C(p, :); % If using ISS model, make SISO
M = (diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1));
M = M*10e-2;

testcases = 10;
count = 1;
max_iter = 100; tol = 10e-8;    plot = 0;
for r = 2:2:testcases
    interp_pts = -logspace(-3,5,r)';  FO_res = rand(r,1);   tmp = rand(r,r);
    SO_res = (tmp+tmp')/2;
    [Er, Ar, br, cr, Mr, conv_nodes, conv_FO_res, conv_SO_res] = lqo_irka(E, A, b, c, M, ...
    interp_pts, FO_res, SO_res, max_iter, tol, plot);
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
E = eye(4,4);

n = 4;  r = 2; % 
lambda_prev = rand(r,1);  phi_prev = rand(r,1);   kappa_prev = rand(r,r);
[Er, Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(E, A, b, c, K, lambda_prev, phi_prev, kappa_prev, 100, 10e-16, 1);

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

H2_prime_s1 = @(s1,s2) -b'*((s1*I-A)'\((s1*I-A)'\K))*((s2*I-A)\b); % Partial deriv w.r.t s1
H2r_prime_s1 = @(s1,s2) -br'*((s1*Ir-Ar)'\((s1*Ir-Ar)'\Kr))*((s2*Ir-Ar)\br);   % Partial deriv w.r.t s1

H2_s2prime = @(s1,s2) -b'*((s1*I-A)'\K)*((s2*I-A)\((s2*I-A)\b)); % Partial deriv w.r.t s2
H2r_s2prime = @(s1,s2) -br'*((s1*Ir-Ar)'\Kr)*((s2*Ir-Ar)\((s2*Ir-Ar)\br)); % Partial deriv w.r.t s2

% So, strictly derivatives are not interpolated... But are the mixed
% conditions satisfied? 

fprintf('Mixed linear + quadratic optimality conditions, 2r in total')

for i = 1:r
    ROpart = 2*phi(i)*H1r_prime(-lambda(i));
    FOpart = 2*phi(i)*H1_prime(-(lambda(i)));
    for j = 1:r
        ROpart = ROpart + kappa(i,j)*H2r_prime_s1(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2r_s2prime(-lambda(j),-lambda(i));
        FOpart = FOpart + kappa(i,j)*H2_prime_s1(-lambda(i),-lambda(j)) + ...
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

E = eye(n,n);
A = diag(-rand(n,1));
b = rand(n,1);  c = rand(1,n);  
K = rand(n,n);  K = (K+K')/2;

lambda_prev = rand(r,1);  phi_prev = rand(r,1);   kappa_prev = rand(r,r);
[Er, Ar, br, cr, Kr, lambda, phi, kappa] = lqo_irka(E, A, b, c, K, lambda_prev, phi_prev, kappa_prev, 100, 10e-12, 1);

% Are interpolation conditions satisfied?
I = eye(n,n);   Ir = eye(r,r);
H1 = @(s) c*((s*E-A)\b); 
H2 = @(s1, s2) b'*((s1*E-A)'\K)*((s2*E-A)\b);

H1r = @(s) cr*((s*Er-Ar)\br); 
H2r = @(s1, s2) br'*((s1*Er-Ar)'\Kr)*((s2*Er-Ar)\br); 

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
H1_prime = @(s) -c*((s*E-A)\((s*E-A)\b));
H1r_prime = @(s) -cr*((s*Er-Ar)\((s*Er-Ar)\br));

H2_prime_s1 = @(s1,s2) -b'*((s1*E-A)'\((s1*E-A)'\K))*((s2*E-A)\b); % Partial deriv w.r.t s1
H2r_prime_s1 = @(s1,s2) -br'*((s1*Er-Ar)'\((s1*Er-Ar)'\Kr))*((s2*Er-Ar)\br);   % Partial deriv w.r.t s1

H2_s2prime = @(s1,s2) -b'*((s1*E-A)'\K)*((s2*E-A)\((s2*E-A)\b)); % Partial deriv w.r.t s2
H2r_s2prime = @(s1,s2) -br'*((s1*Er-Ar)'\Kr)*((s2*Er-Ar)\((s2*Er-Ar)\br)); % Partial deriv w.r.t s2

% So, strictly derivatives are not interpolated... But are the mixed
% conditions satisfied? 

fprintf('Mixed linear + quadratic optimality conditions, 2r in total')

for i = 1:r
    ROpart = 2*phi(i)*H1r_prime(-lambda(i));
    FOpart = 2*phi(i)*H1_prime(-(lambda(i)));
    for j = 1:r
        ROpart = ROpart + kappa(i,j)*H2r_prime_s1(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2r_s2prime(-lambda(j),-lambda(i));
        FOpart = FOpart + kappa(i,j)*H2_prime_s1(-lambda(i),-lambda(j)) + ...
            kappa(j,i)*H2_s2prime(-lambda(j),-lambda(i));
    end
    [FOpart - ROpart]
end

