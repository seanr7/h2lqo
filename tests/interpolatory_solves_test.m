clc
clear all
close all

addpath('drivers/')
addpath('tests/')
%% 2nd-order -> 1st-order test probel
n1 = 30; alpha=.002; beta=alpha; v = 5;

[M, D, K]=triplechain_MSD(n1, alpha, beta, v);

O = zeros(size(K,1),1);
Cv = O';
Cp = ones(1,size(K,1));
B = ones(size(K,1),1);
C = [Cp, Cv];

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

r = 10; q = 250;
shifts = 2*pi*1i*linspace(1,251, 251); % 1 - 251 Hz
fprintf('1. Linfty and Galerkin projection\n')
opts.compression = 'Linfty';
opts.proj = 'g';
[Wprim, Vprim, Worth, Vorth, H_shifts, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r = Worth'*E_qo*Vorth; A_qo_r = Worth'*A_qo*Vorth; 
Q_qo_r = Vorth'*Q_qo*Vorth; B_qo_r = Worth'*B_qo;
% Check interpolation conditions
H2 = @(s1, s2) B_qo' * ((s1 * E_qo - A_qo)'\Q_qo) * ((s2 * E_qo - A_qo)\B_qo);

H2r = @(s1, s2) B_qo_r' * ((s1 * E_qo_r - A_qo_r)'\Q_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r);

fprintf(', %d^2 Interpolation conditions', r)
for i = 1:r
    for j = 1:r
        fprintf('Due to V')
        H2(shifts(pV(i)),  shifts(pV(j))) - H2r(shifts(pV(i)), shifts(pV(j)))
    end
end
for i = 1:r
    for j = 1:r
        fprintf('Due to W')
        H2(shifts(pW(i)),  shifts(pW(j))) - H2r(shifts(pW(i)), shifts(pW(j)))
    end
end

%%
fprintf('2. Linfty and Petrov-Galerkin projection\n')
opts.compression = 'Linfty';
opts.proj = 'pg';
[Wprim, Vprim, Worth, Vorth, H_shifts, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r = Worth'*E_qo*Vorth; A_qo_r = Worth'*A_qo*Vorth; 
Q_qo_r = Vorth'*Q_qo*Vorth; B_qo_r = Worth'*B_qo;
% Check interpolation conditions
H2 = @(s1, s2) B_qo' * ((s1 * E_qo - A_qo)'\Q_qo) * ((s2 * E_qo - A_qo)\B_qo);

H2r = @(s1, s2) B_qo_r' * ((s1 * E_qo_r - A_qo_r)'\Q_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r);

fprintf(', %d^2 Interpolation conditions', r)
for i = 1:r
    for j = 1:r
        fprintf('Due to V')
        H2(shifts(pV(i)),  shifts(pV(j))) - H2r(shifts(pV(i)), shifts(pV(j)))
    end
end
for i = 1:r
    for j = 1:r
        fprintf('Due to W')
        H2(shifts(pW(i)),  shifts(pW(j))) - H2r(shifts(pW(i)), shifts(pW(j)))
    end
end

%%
fprintf('3. avg and Galerkin projection\n')
opts.compression = 'avg';
opts.proj = 'g';
[Wprim, Vprim, Worth, Vorth, H_shifts, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r = Worth'*E_qo*Vorth; A_qo_r = Worth'*A_qo*Vorth; 
Q_qo_r = Vorth'*Q_qo*Vorth; B_qo_r = Worth'*B_qo;
% Check interpolation conditions
H2 = @(s1, s2) B_qo' * ((s1 * E_qo - A_qo)'\Q_qo) * ((s2 * E_qo - A_qo)\B_qo);

H2r = @(s1, s2) B_qo_r' * ((s1 * E_qo_r - A_qo_r)'\Q_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r);

fprintf(', %d^2 Interpolation conditions', r)
for i = 1:r
    for j = 1:r
        fprintf('Due to V')
        H2(shifts(pV(i)),  shifts(pV(j))) - H2r(shifts(pV(i)), shifts(pV(j)))
    end
end
for i = 1:r
    for j = 1:r
        fprintf('Due to W')
        H2(shifts(pW(i)),  shifts(pW(j))) - H2r(shifts(pW(i)), shifts(pW(j)))
    end
end

%%
fprintf('4. avg and Petrov-Galerkin projection\n')
opts.compression = 'avg';
opts.proj = 'g';
[Wprim, Vprim, Worth, Vorth, H_shifts, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r = Worth'*E_qo*Vorth; A_qo_r = Worth'*A_qo*Vorth; 
Q_qo_r = Vorth'*Q_qo*Vorth; B_qo_r = Worth'*B_qo;
% Check interpolation conditions
H2 = @(s1, s2) B_qo' * ((s1 * E_qo - A_qo)'\Q_qo) * ((s2 * E_qo - A_qo)\B_qo);

H2r = @(s1, s2) B_qo_r' * ((s1 * E_qo_r - A_qo_r)'\Q_qo_r) * ((s2 * E_qo_r - A_qo_r)\B_qo_r);

fprintf(', %d^2 Interpolation conditions', r)
for i = 1:r
    for j = 1:r
        fprintf('Due to V')
        H2(shifts(pV(i)),  shifts(pV(j))) - H2r(shifts(pV(i)), shifts(pV(j)))
    end
end
for i = 1:r
    for j = 1:r
        fprintf('Due to W')
        H2(shifts(pW(i)),  shifts(pW(j))) - H2r(shifts(pW(i)), shifts(pW(j)))
    end
end