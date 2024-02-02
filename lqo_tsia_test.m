%% 
% Test script for Linear Quadratic Output (LQO) Two-Sided Iteration
% Algorithm (TSIA)

% Requires: MESS (Matrix Equation Sparse Solvers) toolbox 
clc 
clear all

%% 
% (From MESS) Test problem is a mass spring damper system of three coupled 
% mass spring damper chains with proportional damping;
% System is given in second-order (SO) and converted to first-order (FO);
%   @param n0: order of SO system is n = 3 * n0 + 1
%   @param a, b; proportional damping coeffs., D = a * M + b * K
%   @param viscosity of added dampers

n0 = 10; a = .002; b = a; v = 5;
[M, D, K] = triplechain_MSD(n0, a, b, v);

n = 3 * n0 + 1;
O = zeros(n ,1);
% Velcoity & position output matrices
Cv = O;    Cp = ones(n, 1);
B = ones(n, 1); % Input matrix

% Convert SO system to FO-LQO system
% Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo = spalloc(2*n, 2*n, nnz(M) + n); 
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is (sparse) mass matrix

% State matrix; A_qo = [0, I; -K, -D]
A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(D) + n);  
A_qo(1:n, n+1:2*n) = speye(n); 
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -stiffness matrix
A_qo(n+1:2*n, n+1:2*n) = -D; % (2, 2) block is -damping matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B;

QO_only = 0; % Optional; set to true if linear component of output is present
if QO_only
    C_qo = [];
else
    C_qo = [Cp; Cv];
end

% QO matrix Q
Q_qo = blkdiag(eye(n), zeros(n));

%%
% Set-up input opts
r = 10; % Order of LQO-ROM
itermax = 100;  eps = 10e-8;    plot_conv = true;

% Initial ROM
E_qo_init = eye(r, r);
A_qo_init = -10 * diag(logspace(-1, 2, r));
B_qo_init = ones(r, 1);    
if QO_only
    C_qo_init = [];
else
    C_qo_init = ones(r, 1);
end
Q_qo_init = eye(r, r);

% Perform two-sided iteration algorithm
[E_qo_r, A_qo_r, B_qo_r, C_qo_r, Q_qo_r, pole_history, FOres_history, SOres_history]  ...
    = lqo_tsia(E_qo, A_qo, B_qo, C_qo, Q_qo, E_qo_init, A_qo_init, B_qo_init, C_qo_init, Q_qo_init, itermax, eps, plot_conv);

%%
% Now, check the expected interpolation conditions ... 
% Need to look at poles + residues from second to last iteration, since
% those were used to build the final ROM
[~, no_iters] = size(pole_history);

% Full-order transfer functions
if ~QO_only
    H1 = @(s) C_qo' * ((s * E_qo - A_qo)\B_qo); 
end
H2 = @(s1, s2) B_qo' * ((s1 * E_qo' - A_qo')\Q_qo) * ((s2 * E_qo - A_qo)\B_qo);

if ~QO_only
    H1_r = @(s) C_qo_r' * ((s * E_qo_r - A_qo_r)\B_qo_r); 
end
H2_r = @(s1, s2) B_qo_r' * ((s1 * E_qo_r' - A_qo_r')\Q_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r);

if ~QO_only
    fprintf('Checking 1st-order optimality conditions, %d in total', r)
    for i = 1:r
        H1(-pole_history(i, no_iters - 1)) - H1_r(-pole_history(i, no_iters - 1))
    end
end

fprintf('Checking 2nd-order optimality conditions, %d^2 in total', r)
for i = 1:r
    for j = 1:r
        H2(-pole_history(i, no_iters - 1),  -pole_history(j, no_iters - 1)) ...
            - H2_r(-pole_history(i, no_iters - 1), -pole_history(j, no_iters - 1))
    end
end

% Derivatives of H1(s), H1_r(s)
H1_prime = @(s) - C_qo' * (((s * E_qo - A_qo)\E_qo) * ((s * E_qo - A_qo)\B_qo));
H1_r_prime = @(s) - C_qo_r' * (((s * E_qo_r - A_qo_r)\E_qo_r) * ((s * E_qo_r - A_qo_r)\B_qo_r));

% Partial derivatives of H2(s1, s2), H2_r(s1, s2), w.r.t both arguments
H2_prime_s1 = @(s1, s2) -B_qo' * (((s1 * E_qo' - A_qo')\E_qo') * ...
    ((s1 * E_qo' - A_qo')\Q_qo)) * ((s2 * E_qo - A_qo)\B_qo); 
H2_prime_s2 = @(s1, s2) -B_qo' * ((s1 * E_qo' - A_qo')\Q_qo) * ((s2 * E_qo - A_qo)\E_qo) ...
    * ((s2 * E_qo - A_qo)\B_qo); 
H2_r_prime_s1 = @(s1, s2) -B_qo_r' * (((s1 * E_qo_r' - A_qo_r')\E_qo_r') * ...
    ((s1 * E_qo_r' - A_qo_r')\Q_qo_r)) * ((s2 * E_qo_r - A_qo_r)\B_qo_r); 
H2_r_prime_s2 = @(s1, s2) -B_qo_r' * ((s1 * E_qo_r' - A_qo_r')\Q_qo_r) * ...
    ((s2 * E_qo_r - A_qo_r)\E_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r); 



fprintf('Checking Hermite optimality conditions, %d in total', r)
for i = 1:r
    if ~QO_only
        ro_side = -FOres_history(i, no_iters - 1) * H1_r_prime(-pole_history(i, no_iters - 1));
        fo_side = -FOres_history(i, no_iters - 1) * H1_prime(-pole_history(i, no_iters - 1));
    else
        ro_side = 0;
        fo_side = 0;
    end
    for j = 1:r
        % ro_side = ro_side + 2*SOres_history(i, j, no_iters) * H2_r_prime_s2(-pole_history(i, no_iters - 1), -pole_history(j, no_iters - 1));
        % fo_side = fo_side + 2*SOres_history(i, j, no_iters) * H2_prime_s2(-pole_history(i, no_iters - 1), -pole_history(j, no_iters - 1));
        
        ro_side = ro_side + SOres_history(i, j, no_iters - 1) * H2_r_prime_s1(-pole_history(i, no_iters - 1), -pole_history(j, no_iters - 1)) + ...
                    SOres_history(j, i, no_iters - 1) * H2_r_prime_s2(-pole_history(j, no_iters - 1), -pole_history(i, no_iters - 1));
        fo_side = fo_side + SOres_history(i, j, no_iters - 1) * H2_prime_s1(-pole_history(i, no_iters - 1), -pole_history(j, no_iters - 1)) + ...
                    SOres_history(j, i, no_iters - 1) * H2_prime_s2(-pole_history(j, no_iters - 1), -pole_history(i, no_iters - 1));
    end
    fo_side - ro_side
end