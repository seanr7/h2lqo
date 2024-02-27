%%
% Author: Sean Reiter (seanr7@vt.edu)
clear
close all
%% 
% 'rms' - evaluate root mean square z-displacement of all nodes on the
%         plate surface

% Objective: frequency-domain rms evaluations of plate model in
%            `plateTVA_n201900m1q28278.mat'. Available at
%            https://zenodo.org/record/7671686

% To compute rms `output', need to simulate a SIMO model, with 28278
% outputs. Instead, treat as an linear quadratic output system, 
% and compute rms evaluations directly

addpath('drivers/')
addpath('data/')
fprintf('Loading plateTVA model...\n')
load('data/plateTVA_n201900m1q28278_fo')
n_nodes = full(sum(sum(C)));

%% 
% Convert plate model to FO (first-order) from SO (second-order)
% Model is given in SO-form
% Necessarily, need to conver to FO to do LQO_IRKA for now
fprintf('Converting 2nd-order LTI system to a 1st-order LQO system\n')
tic
[n, ~] = size(M);

E_qo = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is (sparse) mass matrix

A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(E) + n);  % A_qo = [0, I; -K, -E]
A_qo(1:n, n+1:2*n) = speye(n); % (1, 2) block of A_qo
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B;
% No scalar output in this example; only QO

% Our `M' matrix (i.e., the quadratic output matrix) is C' * C
Q_qo = spalloc(2*n, 2*n, nnz(C' * C));
Q_qo(1:n, 1:n) = C' * C; 
fprintf('1st-order LQO realization built in %.2f s\n',toc)

%% 
% Simulate full-order model 
% 250 frequencies to sample at from 0hz - 250hz (s) are given in '*.mat' file 
% Instead, use 500 frequences from 0hz - 500hz

s = 1i*linspace(1,2*pi*501, 501);% Comment out to keep original freq range

% recompute = true;
recompute = false;
if recompute == true
    fprintf('Beginning full-order simulation. Estimated time of completion is %.2f s\n', 250*15.72)
    overall_start = tic;
    res = zeros(1,length(s));
    for ii=1:length(s)
        fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
        current_iter = tic;
        tmp = (s(ii) * E_qo - A_qo) \ B_qo;
        % res(ii) is H2(s(ii), s(ii)) = sqrt(((s(ii) * E_qo - A_qo)') \ Q_qo ((s(ii) * E_qo - A_qo) \ eye(n, n)))
        % i.e., sqrt() of quadratic-output transfer function
        res(ii) = sqrt((tmp' * Q_qo * tmp) / n_nodes); % Q: So really, we want sqrt(H_2(s(ii), s(ii))/n_nodes ? (Just to put it in my language..)
        fprintf('Iteration of full-order simulation finished in %.2f s\n',toc(current_iter))
    end
    fprintf('Full-order simulation finished; total time of completion is %.2f s\n', toc(overall_start))
    f = imag(s)/2/pi;
    mag = 10*log10(abs(res)/1e-9);
    frpintf('Saving full-order simulation data\n')
    filename = 'FOsim_data.mat';
    save(filename,'res','f','mag')
    movefile FOsim_data.mat data/FOsim_data.mat

else
    fprintf('Not re-running the full-order simulation; loading saved data from file FOSIM_data.mat\n')
    % load('FOSIM_data.mat')
end

% figure('name','Transfer function')
% plot(f,mag)
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')

%% 
% Now, compute LQO-ROM via 
%   1. Linfty + Galerkin projection
%   2. Linfty + Petrov-Galerkin projection
%   3. Avg + Galerkin projection
%   4. Avg + Petrov-Galerkin projection


shifts = 1i*linspace(1,2*pi*251, 250);% Comment out to keep original freq rang
r = 50; % Order
fprintf('1. Computing LQO-ROM via Linfty sampling and Galerkin projection\n')
opts.compression = 'Linfty';
opts.proj = 'g';
% Load bases and tf values
load('prim_bases_g') % Saved on remote server
load('H_shifts')
opts.recomp_bases = 0;
opts.recomp_tf = 0;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
opts.H_shifts = [];
[~, ~, Worth, Vorth, H_shifts, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r50_Linfty_g = Worth'*E_qo*Vorth; A_qo_r50_Linfty_g = Worth'*A_qo*Vorth; 
Q_qo_r50_Linfty_g = Vorth'*Q_qo*Vorth; B_qo_r50_Linfty_g = Worth'*B_qo;

% Save shifts 
% save('H_shifts.mat', 'H_shifts');

filename = 'plateTVAlqo_r50_Linfty_g.mat';
save(filename, 'E_qo_r50_Linfty_g', 'A_qo_r50_Linfty_g', 'B_qo_r50_Linfty_g', 'Q_qo_r50_Linfty_g', ...
    'pW', 'pV') 
movefile plateTVAlqo_r50_Linfty_g.mat data/plateTVAlqo_r50_Linfty_g.mat


%%
fprintf('2. Computing LQO-ROM via Linfty sampling and Petrov-Galerkin projection\n')
opts.compression = 'Linfty';
opts.proj = 'pg';
% Load bases and tf values
load('prim_bases_g') % Saved on remote server
opts.recomp_bases = 0; 
opts.recomp_tf = 0;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
load('H_shifts') % Just to be safe
opts.H_shifts = H_shifts;
[~, ~, Worth, Vorth, ~, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r50_Linfty_pg = Worth'*E_qo*Vorth; A_qo_r50_Linfty_pg = Worth'*A_qo*Vorth; 
Q_qo_r50_Linfty_pg = Vorth'*Q_qo*Vorth; B_qo_r50_Linfty_pg = Worth'*B_qo;


filename = 'plateTVAlqo_r50_Linfty_pg.mat';
save(filename, 'E_qo_r50_Linfty_pg', 'A_qo_r50_Linfty_pg', 'B_qo_r50_Linfty_pg', 'Q_qo_r50_Linfty_pg', ...
    'pW', 'pV') 
movefile plateTVAlqo_r50_Linfty_pg.mat data/plateTVAlqo_r50_Linfty_pg.mat

%%
fprintf('3. Computing LQO-ROM via avg (pivoted QR) sampling and Galerkin projection\n')
opts.compression = 'avg';
opts.proj = 'g';
% Grab primitive bases and tf values from previous iteration
load('prim_bases_g') % Saved on remote server
opts.recomp_bases = 0;
opts.recomp_tf = 0;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
load('H_shifts') % Just to be safe
opts.H_shifts = H_shifts;
[~, ~, Worth, Vorth, ~, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r50_avg_g = Worth'*E_qo*Vorth; A_qo_r50_avg_g = Worth'*A_qo*Vorth; 
Q_qo_r50_avg_g = Vorth'*Q_qo*Vorth; B_qo_r50_avg_g = Worth'*B_qo;


filename = 'plateTVAlqo_r50_avg_g.mat';
save(filename, 'E_qo_r50_avg_g', 'A_qo_r50_avg_g', 'B_qo_r50_avg_g', 'Q_qo_r50_avg_g',  ...
    'pW', 'pV')  
movefile plateTVAlqo_r50_avg_g.mat data/plateTVAlqo_r50_avg_g.mat

%%
fprintf('4. Computing LQO-ROM via avg (pivoted QR) sampling and Petrov-Galerkin projection\n')
opts.compression = 'avg';
opts.proj = 'pg';
% Grab primitive bases and tf values from previous iteration
load('prim_bases_pg') % Saved on remote server
opts.recomp_bases = 0;
opts.recomp_tf = 0;
opts.Vprim = Vprim;
opts.Wprim = Wprim;
load('H_shifts') % Just to be safe
opts.H_shifts = H_shifts;
[~, ~, Worth, Vorth, ~, pW, pV, opts] = interpolatory_solves(E_qo, A_qo, B_qo, Q_qo, shifts, r, opts);
% Compute LQO-ROM
E_qo_r50_avg_pg = Worth'*E_qo*Vorth; A_qo_r50_avg_pg = Worth'*A_qo*Vorth; 
Q_qo_r50_avg_pg = Vorth'*Q_qo*Vorth; B_qo_r50_avg_pg = Worth'*B_qo;


filename = 'plateTVAlqo_r50_avg_pg.mat';
save(filename, 'E_qo_r50_avg_pg', 'A_qo_r50_avg_pg', 'B_qo_r50_avg_pg', 'Q_qo_r50_avg_pg', ...
    'pW', 'pV') 
movefile plateTVAlqo_r50_avg_pg.mat data/plateTVAlqo_r50_avg_pg.mat


% %%
% % Now, do reduced-order simulations
% 
% fprintf('Beginning reduced-order simulation\n')
% overall_start = tic;
% 
% res_ro = zeros(1,length(s));
% for ii=1:length(s)
%     fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
%     current_iter = tic;
%     tmp = (s(ii) * E_qo_r - A_qo_r) \ B_qo_r;
%     res_ro(ii) = sqrt((tmp' * Q_qo_r * tmp) / n_nodes);
%     fprintf('Iteration of reduced-order simulation finished in %.2f s\n',toc(current_iter))
% end
% fprintf('Reduced-order simulation finished; time of completion is %.2f s/n', toc(overall_start))
% 
% f_ro = imag(s)/2/pi;
% mag_ro = 10*log10(abs(res_ro)/1e-9);
% 
% fprintf('Saving reduced-order simulation data')
% filename = 'ROsim_plateTVA_r50_lqoirka.mat';
% movefile ROsim_plateTVA_r50_lqoirka.mat data/ROsim_plateTVA_r50_lqoirka.mat
% save(filename,'res_ro','f_ro','mag_ro')
