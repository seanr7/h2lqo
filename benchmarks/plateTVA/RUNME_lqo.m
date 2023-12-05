%%
% Author: Sean Reiter (seanr7@vt.edu)
clear
close all
%% 
% 'rms' - evaluate root mean square z-displacement of all nodes on the
%         plate surface

% End goal: frequency-domain rms evaluations of plate model in
%           `plateTVA_n201900m1q28278.mat'. Available at
%           https://zenodo.org/record/7671686

% To compute rms `output', need to simulation a SIMO model, with 28278
% outputs. Instead, treat as an LQO model, and compute these evaluations
% directly

fprintf('Loading plateTVA model...\n')
load('plateTVA_n201900m1q28278_full')
n_nodes = full(sum(sum(C)));

%% Convert plate model to FO (first-order) from SO (second-order)
% Model is given in SO-form
% Necessarily, need to conver to FO to do LQO_IRKA for now
fprintf('Converting SO realization to a FO-LQO system\n')
tic
[n, ~] = size(M);

E_qo = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is mass matrix

A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(E) + n);  % A_qo = [0, I; -K, -E]
A_qo(1:n, 1:n) = speye(n); % (1, 1) block of A_qo
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
B_qo(n+1:2*n, :) = B;
% No scalar output in this example; only QO

% Our `M' matrix (i.e., the quadratic output matrix) is C' * C
M_qo = spalloc(2*n, 2*n, nnz(C' * C));
M_qo(1:n, 1:n) = C' * C; % Double check this...
fprintf('FO-LQO realization built in %.2f s/n',toc)

%% Sample H2(s1, s2) (QO-tf)
% frequencies to sample at (s) are given in '*.mat' file 

recompute = true;
fprintf('Beginning full-order simulation. Estimated time of completion is %.2f s/n', 250*15.72)
% recompute = false;
overall_start = tic;
if recompute == true
    res = zeros(1,length(s));
    for ii=1:length(s)
        fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
        current_iter = tic;
        tmp = (s(ii) * E_qo - A_qo) \ B_qo;
        res(ii) = sqrt(tmp' * M_qo * tmp) / n_nodes; % Q: So really, we want sqrt(H_2(s(ii), s(ii))/n_nodes ? (Just to put it in my language..)
        fprintf('Iteration of FO-sim finished in %.2f s\n',toc(current_iter))
    end
end

fprintf('Full-order sim finished; time of completion is %.2f s/n', toc(overall_start))

f = imag(s)/2/pi;
mag = 10*log10(abs(res)/1e-9);

fprintf('Saving FO simulation data')
filename = 'FOSIM_data.mat';
save(filename,'res','f','mag')

% figure('name','Transfer function')
% plot(f,mag)
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [dB]')

%% Do via LQO-IRKA
addpath('~/h2lqo')
% addpath('~/Desktop/h2lqo')

fprintf('Beginning construction of the LQO-ROM via LQO-IRKA')

r = 10; % r = in Steffen's; but lets see if we can run this for now..
poles_prev = -logspace(-2, 4, r)'; % Spread 
tmp = rand(r, r);
SO_res_prev = (tmp+tmp')/2; % Try this since M_qo is I?
[E_qo_r, A_qo_r, B_qo_r, ~, M_qo_r, poles, ~, SO_res] = lqo_irka(E_qo, A_qo, B_qo, ...
    [], M_qo, poles_prev, [], SO_res_prev, 100, 10e-8, 1);

filename = 'plateTVAlqo_r250.m';
save(filename, 'E_qo_r', 'A_qo_r', 'B_qo_r', 'M_qo_r', 'poles', 'SO_res') 

% Compute H2 error
% A_qo_err = spalloc(2*n + r, 2*n + r, nnz(A_qo) + nnz(A_qo_r));
% A_qo_err(1:2*n, 1:2*n) = A_qo;  
% A_qo_err(2*n + 1:2*n + r, 2*n + 1:2*n + r) = A_qo_r;
% % blkdiag(A_qo, A_qo_r);
% B_qo_err = spalloc(2*n + r, 1, nnz(B_qo) + nnz(B_qo_r));
% B_qo_err(1:2*n, 1) = B_qo;  
% B_qo_err(2*n + 1:2*n + r, 1) = B_qo_r;
% % [B_qo; B_qo_r];
% M_qo_err = spalloc(2*n + r, 2*n + r, nnz(M_qo) + nnz(M_qo_r));
% M_qo_err(1:2*n, 1:2*n) = M_qo;  
% M_qo_err(2*n + 1:2*n + r, 2*n + 1:2*n + r) = M_qo_r;
% % blkdiag(M_qo, -M_qo_r);
% fprintf('Sanity; is M_qo_r symm?')
% norm(M_qo_r - M_qo_r', 2)
% 
% P_qo_err = lyap(A_qo_err, B_qo_err * B_qo_err');
% Q_qo_err = lyap(A_qo_err', M_qo_err * P_qo_err * M_qo_err);
% fprintf('H2 error')
% sqrt(B_qo_err' * Q_qo_err * B_qo_err)

fprintf('Beginning reduced-order simulation.')
% recompute = false;
overall_start = tic;

res_ro = zeros(1,length(s));
for ii=1:length(s)
    fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    tmp = (s(ii) * speye(r, r) - A_qo_r) \ B_qo_r;
    res_ro(ii) = sqrt(tmp' * M_qo_r * tmp) / n_nodes; % Q: So really, we want sqrt(H_2(s(ii), s(ii))/n_nodes ? (Just to put it in my language..)
    fprintf('Iteration of RO-sim finished in %.2f s\n',toc(current_iter))
end
fprintf('Reduced-order sim finished; time of completion is %.2f s/n', toc(overall_start))

% Save output filefprintf('Full-order sim finished; time of completion is %.2f s/n', toc(overall_start))

f_ro = imag(s)/2/pi;
mag_ro = 10*log10(abs(res_ro)/1e-9);

fprintf('Saving RO simulation data')
filename = 'ROsim_data.mat';
save(filename,'res_ro','f_ro','mag_ro')
