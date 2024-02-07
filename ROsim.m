clc
clear all

%%
% Load data
load('ROsim_data_r_10.mat')
load('plateTVAlqo_r_10.mat')
% load('FOSIM_data.mat')
load('plateTVA_n201900m1q28278_full.mat')

%%

f = imag(s)/2/pi;
mag = 10*log10(abs(res)/1e-9);

n_nodes = full(sum(sum(C)));
fprintf('Beginning reduced-order simulation\n')
% recompute = false;
overall_start = tic;

res_ro = zeros(1,length(s));
for ii=1:length(s)
    fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    tmp = (s(ii) * E_qo_r - A_qo_r) \ B_qo_r;
    res_ro(ii) = sqrt((tmp' * Q_qo_r * tmp)/ n_nodes) ; % Q: So really, we want sqrt(H_2(s(ii), s(ii))/n_nodes ? (Just to put it in my language..)
    fprintf('Iteration of RO-sim finished in %.2f s\n',toc(current_iter))
end
fprintf('Reduced-order sim finished; time of completion is %.2f s/n', toc(overall_start))

% Save output filefprintf('Full-order sim finished; time of completion is %.2f s/n', toc(overall_start))

f_ro = imag(s)/2/pi;
mag_ro = 10*log10(abs(res_ro)/1e-9);


plot(f_ro,mag_ro)
hold on
plot(f, mag)