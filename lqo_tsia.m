function [Er, Ar, Br, Cr, Mr, pole_history, FOres_history, SOres_history] = ...
    lqo_tsia(E, A, B, C, M, Er, Ar, Br, Cr, Mr, itermax, eps, plot_conv, compute_res)
% Author: Sean Reiter (seanr7@vt.edu)

% Apply the Two-sided Iteration Algorithm (TSIA) to Linear Quadratic
% Output (LQO) systems @cite [GosA19]
 
% Paramters
% ---------
% @param (E, A, b, c, M) \in (Rnn, Rnn, Rn1, R1n, Rnn):
%   State-space description of an n-dimensional LQO model. Assumed to be
%   asymptotically stable.
%   M is assumed to be symmetric; M = M'. 
%   A is assumed to be diagonalizable.
%   E (mass) is assumed to be nonsingular.

% @param (Er, Ar, Br, Cr, Mr) \in (Err, Rrr, Rr1, R1r, Rrr)
%   Initial LQO-ROM to initialize TSIA. Assumed to be asymptotically
%   stable.
%   Mr is assumed to be symmetric; M_init = M_init'. 
%   Ar is assumed to be diagonalizable.
%   Er (mass) is assumed to be nonsingular.

% @param itermax:
%   Max no. of iterations for algorithm to run. (Default itermax = 100)

% @param eps:
%   Convergence tol; if (max(eig(A_r) - L_prev) < eps), converged. 
%   (Default eps = 10e-8)

% @param plot_conv:
%   Bool; do we plot the convergence of the reduced-order poles eig(A_r)?
%   (Default plot_conv = False)

% @param compute_res
%   Bool; do we compute the 1st and 2nd-order residues of the LQO-ROM?
%   (Default compute_res = False)

% Outputs
% ---------
% @param (Er, Ar, br, cr, Mr) \in (Err, Rrr, Rr1, R1r, Rrr):
%   State-space realization of the locally H2-optimal r-dimensional
%   LQO-ROM, obtained via Petrov-Galerkin (PG) projection.

% @param pole_history:
%   Reduced-order poles (mirros of interpolation poitns), as an 
%   (r x iter) array.

% @param FOres_history, SOres_history
%   History of 1st and 2nd-order residues; for purpose of checking
%   interpolation conditions, only.

%%
% Check iteration/convergence tolerances; set to defaults if unspecified
if nargin == 13 % (Default is to not plot convergence of the RO-poles)
    compute_res = 0;
end
if nargin == 12 % (Default is to not plot convergence of the RO-poles)
    plot_conv = 0;
end
if nargin == 11 % (No convergence tolerance set)
    eps = 10e-8;
end
if nargin == 10 % (No max number of iterations set)
    itermax = 100;
end

[r, ~] = size(Ar); % Dimension of ROM
% If no c term inputted; output is purely quadratic
pure_QO = false;
if isempty(C)
    pure_QO = true;
    Cr = []; 
    FOres_history = [];
end

% Start the clock on the total iteration
overall_timer_start = tic;
fprintf('Beginning LQO-TSIA\n')

% Counter + tolerance to enter while
iter = 1;   err(iter) = eps + 1; 

[Xr, Lr] = eig(Ar, Er); 
pole_history(:, iter) = diag(Lr); % Save initial poles
% Compute initial poles + residues
if ~compute_res
    FOres_history = [];     SOres_history = [];
else
    % Now, compute and save residues for the sake of checking interpolation
    if ~pure_QO % Compute FO residues
        FOres = (Cr' * Xr) * diag(Xr\(Er\Br));
        FOres_history(:, iter) = FOres.';
    end
    SOres = diag(Xr\(Er\Br)).' * (Xr.' * Mr * Xr) * diag(Xr\(Er\Br));
    SOres_history(:, :, iter) = SOres;
end

while (err(iter) > eps && iter <= itermax)
    % Start the clock
    thisiter_timer_start = tic;
    fprintf('Beginning next TSIA iteration; current iterate is k = %d\n', iter)
    
    % Solve sparse-dense Sylvester equations; 
    % X \in Rnr satisfies
    %   @math: A*X*Er' + E*X*Ar' + B*Br' = 0
    % For some reason; Convergence is VERY poor if -B*Br' not present...
    [X, ~, ~, ~, ~] = mess_sylvester_sparse_dense(A, 'N', Ar, 'T', -B*Br', E, Er);
    % And, Z \in \Rnr satisfies 
    %   @math: A'*Z*Er + E'*Z*Ar - 2*M*X*Mr - C*Cr' = 0
    rhs = -2*M*X*Mr';
    if ~pure_QO
        rhs = rhs - C*Cr';
    end
    [Z, ~, ~, ~, ~] = mess_sylvester_sparse_dense(A, 'T', Ar, 'N', rhs, E, Er);
    
    % Orthonormalize projection matrices
    [Vr, ~] = qr(X, "econ");     [Wr, ~] = qr(Z, "econ");
    % Compute LQO-ROM via PG-proj 
    Er = Wr'*E*Vr;   Ar = Wr'*A*Vr;   Br = Wr'*B;
    if ~pure_QO % Compute reduced linear output term
        Cr = (C'*Vr)';  
    end
    Mr = Vr'*M*Vr;

    % End the clock
    fprintf('End of current TSIA iterate k = %d\n', iter)
    fprintf('Current TSIA iteration finished in %.2f s\n',toc(thisiter_timer_start))

    iter = iter + 1;    
    % Get poles (interpolation points) from previous iteration
    [Xr, Lr] = eig(Ar, Er); 
    % [Xr, Lr] = eig(Er\Ar); 
    pole_history(:, iter) = diag(Lr);
    err(iter) = max(abs(pole_history(:, iter) - pole_history(:, iter - 1)));
    fprintf('Largest magnitude change in interpolation points is %.2f \n', err(iter))

    % Compute residues
    if compute_res
        % Now, compute and save residues for the sake of checking interpolation
        if ~pure_QO % Compute FO residues
            FOres = (Cr' * Xr) * diag(Xr\(Er\Br));
            FOres_history(:, iter) = FOres.';
        end
        SOres = diag(Xr\(Er\Br)).' * (Xr.' * Mr * Xr) * diag(Xr\(Er\Br));
        SOres_history(:, :, iter) = SOres;
    end
end

if iter == (itermax + 1)
    fprintf('TSIA has terminated due to reaching the max no. of iterations; total time elapsed is %.2f s\n', toc(overall_timer_start))
else
    fprintf('TSIA has converged in %d iterations\n', iter)
    fprintf('Total time elapsed is %.2f s\n', toc(overall_timer_start))
end

% If plotting convergence
if plot_conv == true
    semilogy(1:iter, err, '-o', LineWidth=1.5)
    xlim([1,iter])
    xlabel('Iteration count')
    ylabel('Magnitude of change in \lambda')
end
end