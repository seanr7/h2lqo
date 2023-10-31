function [A_r, b_r, c_r, M_r, L, phi, mu] = lqo_irka(A, b, c, M, ...
    L_prev, phi_prev, mu_prev, itermax, eps, plot_conv)
% Author: Sean Reiter (seanr7@vt.edu)

% Apply the Iterative Rational Krylov Algorithm (IRKA) to Linear Quadratic
% Output (LQO) systems.
 
% Paramters
% ---------
% @param (A, b, c, M) \in (Rnn, Rn1, R1n, Rnn):
%   State-space description of an n-dimensional LQO model. Assumed to be
%   asymptotically stable.
%   M is assumed to be symmetric; M = M'. 
%   A is assumed to be diagonalizable.

% @param L_prev:
%   Initial selection of interpolation points in open left half-plane.
%   Stored as an (r x 1) array.
%   Note: These are the `reduced-order poles' at initialization; 
%   i.e., Ar = diag(L_prev).

% @param phi_prev, mu_prev:
%   Initial selection of 1st and 2nd-order (complex) residues.
%   Stored as (r x 1) and (r x r) arrays, respectively. 
%   Note: In the basis of Ar's eigenvectors, these are the linear and
%   quadratic output matrices; i.e., c_r = phi_prev, M_r = mu_prev.

% @param itermax:
%   Max no. of iterations for LQO-IRKA to run. (Default itermax = 100)

% @param eps:
%   Convergence tol; if (max(eig(A_r) - L_prev) < eps), converged. 
%   (Default eps = 10e-8)

% @param plot_conv:
%   Bool; do we plot the convergence of the reduced-order poles eig(A_r)?
%   (Default plot_conv = False)

% Outputs
% ---------
% @param (A_r, b_r, c_r, M_r) \in (Rrr, Rr1, R1r, Rrr):
%   State-space realization of the locally H2-optimal r-dimensional
%   LQO-ROM, obtained via Petrov-Galerkin (PG) projection.

% @param L:
%   Converged interpolation points, as an (r x 1) array.

% @param phi, mu:
%   Converged 1st and 2nd-order residues, as (r x 1) and (r x r) arrays.

%% LQO-IRKA
% Take input data; Compute initial LQO-ROM via PG-proj
r = length(L_prev);     [~, n] = size(A); % ROM and FOM dimensions
I = eye(n, n);

% Check iteration/convergence tolerances; set to defaults if unspecified
if nargin == 9 % (Default is to not plot convergence of the RO-poles)
    plot_conv = 0;
end
if nargin == 8 % (No convergence tolerance set)
    eps = 10e-8;
end
if nargin == 7 % (No max number of iterations set)
    itermax = 100;
end

% Counter + tolerance to enter while
iter = 1;   err(iter) = eps + 1; 
while (err(iter) > eps && iter < itermax)
    % Pre-allocate space for right, left PG-proj bases V_r, W_r
    % (Note: Construction requires RO-poles and residues @ each iter)
    V_r = zeros(n, r);     W_r = zeros(n, r);
    % TODO: Optimize to minimze no. of linear solves
    for k = 1:r % Fill out columns of V_r, W_r
        % 1. Construction of Vr enforces r + r^2 interpolation conditions:
        %   @math: H1(-L_prev(k)) = H1r(-L_prev(k)), k = 1, ..., r
        %   @math: H2(-L_prev(i), -L_prev(j)) = H2(-L_prev(i), -L_prev(j)), 
        %           i, j = 1, ..., r
        V_r(:, k) = (-L_prev(k) * I - A)\b; 
        % 2. Construction of W_r enforces r `mixed' Hermite conditions:
        %   @math: phi_prev(k) * H1^(1)(-L_prev(k)) + \sum_{j = 1}^{r} ...
        %          mu_prev(k, j) * H2^(1, 0)(-L_prev(k), -L_prev(j)) = ...
        %          phi_prev(k) * H1_r^(1)(-L_prev(k)) + \sum_{j = 1}^{r} ...
        %          mu_prev(k, j) * H2_r^(1, 0)(-L_prev(k), -L_prev(j)) = ...
        W_r(:, k) = phi_prev(k)*((-L_prev(k) * I - A')\c');
        tmp = (-L_prev(k) * I - A')\M; 
        for i = 1:r
            W_r(:, k) = W_r(:, k) + mu_prev(i, k)* tmp * ...
                            ((-L_prev(i) * I - A)\b);
        end
    end
    % Orthonormalize projection matrices
    [V_r, ~] = qr(V_r, "econ");     [W_r, ~] = qr(W_r, "econ");
    % Compute LQO-ROM via PG-proj and monitor convergence (E_r = eye(r, r))
    W_rt = W_r';  V_rt = V_r';
    A_r = (W_rt * V_r)\(W_rt * A * V_r);   b_r = (W_rt * V_r)\(W_rt * b);
    c_r = c * V_r;  M_r = V_rt * M * V_r;
    
    % Diagonalize Ar; Get RO-poles ...
    [X_r, L_r] = eig(A_r);     L = diag(L_r);
    % ... + 1st, 2nd-order residues of H1r, H2r
    phi = (c_r * X_r)'.* (X_r\b_r);    mu = X_r' * M_r * X_r;

    % Track convergence of poles; Update iter count + compute max deviation 
    % between poles across iterations
    iter = iter + 1;    err(iter) = max(abs(L - L_prev));
    % Overwrite RO-poles + residues
    L_prev = L;   phi_prev = phi;    mu_prev = mu;
end
% If plotting convergence
if plot_conv == true
    semilogy(1:iter, err, '-o', LineWidth=1.5)
    xlim([1,iter])
    xlabel('Iteration count')
    ylabel('Magnitude of change in \lambda')
end
end