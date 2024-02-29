clc
clear all
% addpath('path/to/data/')
addpath('data/')

%%
% r = 25 examples
load('data/plateTVA_r25_redux_lqoirka.mat')
load('data/plateTVAlqo_r25_Linfty_g.mat')
load('data/plateTVAlqo_r25_Linfty_pg.mat')
load('data/plateTVAlqo_r25_avg_g.mat')
load('data/plateTVAlqo_r25_avg_pg.mat')

% And, FOsim data
load('data/simdata/FOSIM_data.mat')
% Also the C matrix to get n_nodes
load('data/plateTVA_n201900m1q28278_fo.mat','C')
s = 1i*linspace(0,2*pi*250, 500); % Frequencies used in FOsim
s_hz = imag(s)/2/pi; % In Hz, if you so desire
mag = 10*log10(abs(res)/1e-9); % Magnitude of FO transfer function

n_nodes = full(sum(sum(C)));


%%
% r = 25 simulation
n_nodes = full(sum(sum(C)));
fprintf('Beginning reduced-order simulation\n')
overall_start = tic;

res_r25_irka = zeros(1,length(s));
res_r25_Linfty_g = zeros(1, length(s));
res_r25_Linfty_pg = zeros(1, length(s));
res_r25_avg_g = zeros(1, length(s));
res_r25_avg_pg = zeros(1, length(s));
for ii=1:length(s)
    fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    tmp_r25_irka = (s(ii) * E_qo_r - A_qo_r) \ B_qo_r;
    tmp_r25_Linfty_g = (s(ii) * E_qo_r25_Linfty_g - A_qo_r25_Linfty_g) \ B_qo_r25_Linfty_g;
    tmp_r25_Linfty_pg = (s(ii) * E_qo_r25_Linfty_pg - A_qo_r25_Linfty_pg) \ B_qo_r25_Linfty_pg;
    tmp_r25_avg_g = (s(ii) * E_qo_r25_avg_g - A_qo_r25_avg_g) \ B_qo_r25_avg_g;
    tmp_r25_avg_pg = (s(ii) * E_qo_r25_avg_pg - A_qo_r25_avg_pg) \ B_qo_r25_avg_pg;
    % Transfer function evals
    res_r25_irka(ii) = sqrt((tmp_r25_irka' * Q_qo_r * tmp_r25_irka)/ n_nodes);
    res_r25_Linfty_g(ii) = sqrt((tmp_r25_Linfty_g' * Q_qo_r25_Linfty_g * tmp_r25_Linfty_g)/ n_nodes);
    res_r25_Linfty_pg(ii) = sqrt((tmp_r25_Linfty_pg' * Q_qo_r25_Linfty_pg * tmp_r25_Linfty_pg)/ n_nodes);
    res_r25_avg_g(ii) = sqrt((tmp_r25_avg_g' * Q_qo_r25_avg_g * tmp_r25_avg_g)/ n_nodes);
    res_r25_avg_pg(ii) = sqrt((tmp_r25_avg_pg' * Q_qo_r25_avg_pg * tmp_r25_avg_pg)/ n_nodes);
    fprintf('Iteration of RO-sim finished in %.2f s\n',toc(current_iter))
end
fprintf('Reduced-order sim finished; time of completion is %.2f s\n', toc(overall_start))

mag_r25_irka = 10*log10(abs(res_r25_irka)/1e-9);
mag_r25_Linfty_g = 10*log10(abs(res_r25_Linfty_g)/1e-9);
mag_r25_Linfty_pg = 10*log10(abs(res_r25_Linfty_pg)/1e-9);
mag_r25_avg_g = 10*log10(abs(res_r25_avg_g)/1e-9);
mag_r25_avg_pg = 10*log10(abs(res_r25_avg_pg)/1e-9);
 
write = 1;
if write
    % Store data
    magmatrix = [s_hz', mag', mag_r25_irka', mag_r25_Linfty_g', mag_r25_Linfty_pg', mag_r25_avg_g', mag_r25_avg_pg'];
    filename = 'outputs_for_paper/r25_mag.dat';
    dlmwrite(filename, magmatrix, ...
        'delimiter', '\t', 'precision', 8);
    errmatrix = [s_hz', abs(mag-mag_r25_irka)./abs(mag)', abs(mag-mag_r25_Linfty_g)./abs(mag)', abs(mag-mag_r25_Linfty_pg)./abs(mag)', abs(mag-mag_r25_avg_g)./abs(mag)', abs(mag-mag_r25_avg_pg)./abs(mag)'];
    filename = 'outputs_for_paper/r25_err.dat';
    dlmwrite(filename, errmatrix, ...
        'delimiter', '\t', 'precision', 8);
end

% plot(s_hz,mag_ro10)
% hold on
% plot(f, mag)

ColMat = zeros(6,3);

%ColMat(1,:) = [1 0.6 0.4];
ColMat(1,:) = [ 0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [  0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];
ColMat(6,:) = [1 0.4 0.6];

% f=figure;
% f.Position = [476 445 700 280];
% Make aspect ration `golden'
%%
figure
golden_ratio = (sqrt(5)+1)/2;

plot(s_hz, mag,'-o','color',ColMat(1,:),'markersize',10,LineWidth=1);hold on;
plot(s_hz, mag_r25_irka, '-.','color',ColMat(2,:),LineWidth=1);
plot(s_hz, mag_r25_Linfty_g,'--','color', ColMat(3,:),LineWidth=1);
plot(s_hz, mag_r25_Linfty_pg,'--.','color', ColMat(4,:),LineWidth=1);
plot(s_hz, mag_r25_avg_g,'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, mag_r25_avg_pg,'--.','color', ColMat(6,:),LineWidth=1);

legend('Full-order', 'r=25, IRKA', 'r=25 Linfty-g', 'r=25 Linfty-pg', 'r=25 avg-g', 'r=25 avg-pg')
%%

% Plot erro
% f=figure;
% f.Position = [476 445 700 280];
% Make aspect ration `golden'
figure

semilogy(s_hz, abs(mag-mag_r25_irka)./abs(mag), '-.','color',ColMat(2,:),LineWidth=1);hold on
semilogy(s_hz, abs(mag-mag_r25_Linfty_g)./abs(mag),'--','color', ColMat(3,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_Linfty_pg)./abs(mag),'--.','color', ColMat(4,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_avg_g)./abs(mag),'--.','color', ColMat(5,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_avg_pg)./abs(mag),'--.','color', ColMat(6,:),LineWidth=1);

legend('r=25, IRKA', 'r=25 Linfty-g', 'r=25 Linfty-pg', 'r=25 avg-g', 'r=25 avg-pg')

%%
% r = 50 examples
% load('data/plateTVA_r50_redux_lqoirka.mat')
load('data/plateTVAlqo_r50_Linfty_g.mat')
load('data/plateTVAlqo_r50_Linfty_pg.mat')
load('data/plateTVAlqo_r50_avg_g.mat')
load('data/plateTVAlqo_r50_avg_pg.mat')

% And, FOsim data
load('data/simdata/FOSIM_data.mat')
% Also the C matrix to get n_nodes
load('data/plateTVA_n201900m1q28278_fo.mat','C')
s = 1i*linspace(0,2*pi*250, 500); % Frequencies used in FOsim
s_hz = imag(s)/2/pi; % In Hz, if you so desire
mag = 10*log10(abs(res)/1e-9); % Magnitude of FO transfer function

n_nodes = full(sum(sum(C)));


%%
% r = 25 simulation
n_nodes = full(sum(sum(C)));
fprintf('Beginning reduced-order simulation\n')
overall_start = tic;

% res_r50_irka = zeros(1,length(s));
res_r50_Linfty_g = zeros(1, length(s));
res_r50_Linfty_pg = zeros(1, length(s));
res_r50_avg_g = zeros(1, length(s));
res_r50_avg_pg = zeros(1, length(s));
for ii=1:length(s)
    fprintf('Frequency step %d, f=%.2f Hz ... ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    % tmp_r50_irka = (s(ii) * E_qo_r - A_qo_r) \ B_qo_r;
    tmp_r50_Linfty_g = (s(ii) * E_qo_r50_Linfty_g - A_qo_r50_Linfty_g) \ B_qo_r50_Linfty_g;
    tmp_r50_Linfty_pg = (s(ii) * E_qo_r50_Linfty_pg - A_qo_r50_Linfty_pg) \ B_qo_r50_Linfty_pg;
    tmp_r50_avg_g = (s(ii) * E_qo_r50_avg_g - A_qo_r50_avg_g) \ B_qo_r50_avg_g;
    tmp_r50_avg_pg = (s(ii) * E_qo_r50_avg_pg - A_qo_r50_avg_pg) \ B_qo_r50_avg_pg;
    % Transfer function evals
    % res_r50_irka(ii) = sqrt((tmp_r50_irka' * Q_qo_r * tmp_r50_irka)/ n_nodes);
    res_r50_Linfty_g(ii) = sqrt((tmp_r50_Linfty_g' * Q_qo_r50_Linfty_g * tmp_r50_Linfty_g)/ n_nodes);
    res_r50_Linfty_pg(ii) = sqrt((tmp_r50_Linfty_pg' * Q_qo_r50_Linfty_pg * tmp_r50_Linfty_pg)/ n_nodes);
    res_r50_avg_g(ii) = sqrt((tmp_r50_avg_g' * Q_qo_r50_avg_g * tmp_r50_avg_g)/ n_nodes);
    res_r50_avg_pg(ii) = sqrt((tmp_r50_avg_pg' * Q_qo_r50_avg_pg * tmp_r50_avg_pg)/ n_nodes);
    fprintf('Iteration of RO-sim finished in %.2f s\n',toc(current_iter))
end
fprintf('Reduced-order sim finished; time of completion is %.2f s\n', toc(overall_start))

% mag_r50_irka = 10*log10(abs(res_r50_irka)/1e-9);
mag_r50_Linfty_g = 10*log10(abs(res_r50_Linfty_g)/1e-9);
mag_r50_Linfty_pg = 10*log10(abs(res_r50_Linfty_pg)/1e-9);
mag_r50_avg_g = 10*log10(abs(res_r50_avg_g)/1e-9);
mag_r50_avg_pg = 10*log10(abs(res_r50_avg_pg)/1e-9);
 
% plot(s_hz,mag_ro10)
% hold on
% plot(f, mag)

ColMat = zeros(6,3);

%ColMat(1,:) = [1 0.6 0.4];
ColMat(1,:) = [ 0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [  0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];
ColMat(6,:) = [1 0.4 0.6];

% f=figure;
% f.Position = [476 445 700 280];
% Make aspect ration `golden'
%%
figure
golden_ratio = (sqrt(5)+1)/2;

plot(s_hz, mag,'-o','color',ColMat(1,:),'markersize',10,LineWidth=1);hold on;
% plot(s_hz, mag_r50_irka, '-.','color',ColMat(2,:),LineWidth=1);
plot(s_hz, mag_r50_Linfty_g,'--','color', ColMat(3,:),LineWidth=1);
plot(s_hz, mag_r50_Linfty_pg,'--.','color', ColMat(4,:),LineWidth=1);
plot(s_hz, mag_r50_avg_g,'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, mag_r50_avg_pg,'--.','color', ColMat(6,:),LineWidth=1);

legend('Full-order', 'r=50 Linfty-g', 'r=50 Linfty-pg', 'r=50 avg-g', 'r=50 avg-pg')
% legend('Full-order', 'r=25, IRKA', 'r=25 Linfty-g', 'r=25 Linfty-pg', 'r=25 avg-g', 'r=25 avg-pg')
%%

figure
golden_ratio = (sqrt(5)+1)/2;

% plot(s_hz, abs(mag-mag_r50_irka), '-.','color',ColMat(2,:),LineWidth=1);hold on
plot(s_hz, abs(mag-mag_r50_Linfty_g),'--','color', ColMat(3,:),LineWidth=1);hold on
plot(s_hz, abs(mag-mag_r50_Linfty_pg),'--.','color', ColMat(4,:),LineWidth=1);
plot(s_hz, abs(mag-mag_r50_avg_g),'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, abs(mag-mag_r50_avg_pg),'--.','color', ColMat(6,:),LineWidth=1);

% legend('r=25, IRKA', 'r=25 Linfty-g', 'r=25 Linfty-pg', 'r=25 avg-g', 'r=25 avg-pg')
legend('r=50 Linfty-g', 'r=50 Linfty-pg', 'r=50 avg-g', 'r=50 avg-pg')