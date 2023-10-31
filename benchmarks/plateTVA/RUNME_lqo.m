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

load('plateTVA_n201900m1q28278_full')
n_nodes = full(sum(sum(C)));

%% Convert plate model to FO (first-order) from SO (second-order)
% Model is given in SO-form
% Necessarily, need to conver to FO to do LQO_IRKA for now
[n, ~] = size(M);

E_qo = spalloc(2*n, 2*n, nnz(M) + n); % Descriptor matrix; E_qo = [I, 0: 0, M]
E_qo(1:n, 1:n) = speye(n); % (1, 1) block
E_qo(n+1:2*n, n+1:2*n) = M; % (2, 2) block is mass matrix

A_qo = spalloc(2*n, 2*n, nnz(K) + nnz(E) + n);  % A_qo = [0, I; -K, -E]
A_qo(1:n, 1:n) = speye(n); % (1, 1) block of A_qo
A_qo(n+1:2*n, 1:n) = -K;  % (2, 1) block is -damping matrix
A_qo(n+1:2*n, n+1:2*n) = -E; % (2, 2) block is -stiffness matrix

B_qo = spalloc(2*n, 1, nnz(B)); % B_qo = [0; B];
% No scalar output in this example; only QO

% Our `M' matrix (i.e., the quadratic output matrix) is just the (n_nodes x n_nodes) identity
M_qo = speye(n_nodes);

%% Sample H2(s1, s2) (QO-tf)
% frequencies to sample at (s) are given in '*.mat' file 

recompute = true;
if recompute == true
    res = zeros(1,length(s));
    for ii=1:length(s)
        fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
        tic
        tmp = (s(ii) * E_qo - A_qo) \ B_qo;
        print(size(tmp))
        res(ii) = sqrt(tmp' * M_qo * tmp) / n_nodes; % Q: So really, we want sqrt(H_2(s(ii), s(ii))/n_nodes ? (Just to put it in my language..)
        fprintf('finished in %.2f s\n',toc)
    end
end

f = imag(s)/2/pi;
mag = 10*log10(abs(res)/1e-9);

save(rms_lqo,"res",'-mat')

figure('name','Transfer function')
plot(f,mag)
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

%% Do via LQO-IRKA