%% RUNME
% Script file to run all experiments on the vibro-acoustic plate (plateTVA) 
% data set.

%
% This file is part of the archive Code and Results for Numerical 
% Experiments in "Interpolatory model order reduction of large-scale 
% dynamical systems with root mean squared error measures"
% Copyright (c) 2024 Sean Reiter, Steffen W. R. Werner
% All rights reserved.
% License: BSD 2-Clause license (see COPYING)
%

clc;
clear all;
close all;

% Get and set all paths
[rootpath, filename, ~] = fileparts(mfilename('fullpath'));
loadname            = [rootpath filesep() ...
    'data' filesep() filename];
savename            = [rootpath filesep() ...
    'results' filesep() filename];

% Add paths to drivers and data
addpath([rootpath, '/drivers'])
addpath([rootpath, '/data'])

% Write .log file, put in `out' folder
if exist([savename '.log'], 'file') == 2
    delete([savename '.log']);
end
outname = [savename '.log']';

diary(outname)
diary on; 

fprintf(1, ['SCRIPT: ' upper(filename) '\n']);
fprintf(1, ['========' repmat('=', 1, length(filename)) '\n']);
fprintf(1, '\n');

%% Load base data for r = 25.
fprintf(1, 'Load problem data for order r = 25 reduced models\n');
fprintf(1, '-------------------------------------------------\n');

% Order r = 25 reduced models
load('results/plateTVAlqo_r25_lqoirka.mat')
load('results/plateTVAlqo_r25_Linfty_g.mat')
load('results/plateTVAlqo_r25_Linfty_pg.mat')
load('results/plateTVAlqo_r25_avg_g.mat')
load('results/plateTVAlqo_r25_avg_pg.mat')

% Simulation data for full-order plateTVA model
load('results/fosim_data.mat')
% Load C matrix from original model to get n_nodes used in calculating rms
% response
load('data/plateTVA_n201900m1q28278.mat','C')
n_nodes = full(sum(sum(C)));

% Compute frequencies used in the simulation
s    = 1i*linspace(0,2*pi*250, 500); 
s_hz = imag(s)/2/pi; 

fprintf(1, '\n');

%% Simulate rms response r = 25.
fprintf(1, 'Simulate rms response for order r = 25 reduced models\n');
fprintf(1, '-----------------------------------------------------\n');

% Allocate space for results
res_r25_irka      = zeros(1,length(s));
res_r25_Linfty_g  = zeros(1, length(s));
res_r25_Linfty_pg = zeros(1, length(s));
res_r25_avg_g     = zeros(1, length(s));
res_r25_avg_pg    = zeros(1, length(s));

fprintf(1, 'Beginning reduced order simulations\n')
fprintf(1, '-----------------------------------\n');
overall_start = tic;

for ii=1:length(s)
    fprintf(1, 'Frequency step %d, f=%.2f Hz ...\n ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    tmp_r25_irka      = (s(ii)*Efo_r - Afo_r)\Bfo_r;
    tmp_r25_Linfty_g  = (s(ii)*Efo_r25_Linfty_g - Afo_r25_Linfty_g)\Bfo_r25_Linfty_g;
    tmp_r25_Linfty_pg = (s(ii)*Efo_r25_Linfty_pg - Afo_r25_Linfty_pg)\Bfo_r25_Linfty_pg;
    tmp_r25_avg_g     = (s(ii)*Efo_r25_avg_g - Afo_r25_avg_g)\Bfo_r25_avg_g;
    tmp_r25_avg_pg    = (s(ii)*Efo_r25_avg_pg - Afo_r25_avg_pg)\Bfo_r25_avg_pg;
    % Evaluate rms response
    res_r25_irka(ii)      = sqrt((tmp_r25_irka'*Qfo_r*tmp_r25_irka)/n_nodes);
    res_r25_Linfty_g(ii)  = sqrt((tmp_r25_Linfty_g'*Qfo_r25_Linfty_g*tmp_r25_Linfty_g)/n_nodes);
    res_r25_Linfty_pg(ii) = sqrt((tmp_r25_Linfty_pg'*Qfo_r25_Linfty_pg*tmp_r25_Linfty_pg)/n_nodes);
    res_r25_avg_g(ii)     = sqrt((tmp_r25_avg_g' * Qfo_r25_avg_g * tmp_r25_avg_g)/n_nodes);
    res_r25_avg_pg(ii)    = sqrt((tmp_r25_avg_pg'*Qfo_r25_avg_pg*tmp_r25_avg_pg)/n_nodes);
    fprintf(1, 'Current iteration of reduced order simulation finished in %.2f s\n',toc(current_iter))
    fprintf(1, '----------------------------------------------------------------------\n');
end
fprintf(1, 'Reduced order simulations finished in %.2f s\n', toc(overall_start))
fprintf(1, '--------------------------------------------------\n');

% Compute magnitudes to plot
% To plot magnitudes in absolute value, uncomment the code below
% mag               = abs(res);
% mag_r25_irka      = abs(res_r25_irka);
% mag_r25_Linfty_g  = abs(res_r25_Linfty_g);
% mag_r25_Linfty_pg = abs(res_r25_Linfty_pg);
% mag_r25_avg_g     = abs(res_r25_avg_g);
% mag_r25_avg_pg    = abs(res_r25_avg_pg);

% To plot magnitudes in dB, uncomment the code below
mag               = 10*log10(abs(res)/1e-9); 
mag_r25_irka      = 10*log10(abs(res_r25_irka)/1e-9);
mag_r25_Linfty_g  = 10*log10(abs(res_r25_Linfty_g)/1e-9);
mag_r25_Linfty_pg = 10*log10(abs(res_r25_Linfty_pg)/1e-9);
mag_r25_avg_g     = 10*log10(abs(res_r25_avg_g)/1e-9);
mag_r25_avg_pg    = 10*log10(abs(res_r25_avg_pg)/1e-9);
 
write = 1;
if write
    % Store data
    magmatrix = [s_hz', mag', mag_r25_irka', mag_r25_Linfty_g', mag_r25_Linfty_pg', ...
        mag_r25_avg_g', mag_r25_avg_pg'];
    dlmwrite('results/r25_mag.dat', magmatrix, ...
        'delimiter', '\t', 'precision', 8);
    errmatrix = [s_hz', (abs(mag-mag_r25_irka)./abs(mag))', (abs(mag-mag_r25_Linfty_g)./abs(mag))', ...
        (abs(mag-mag_r25_Linfty_pg)./abs(mag))', (abs(mag-mag_r25_avg_g)./abs(mag))', ...
        (abs(mag-mag_r25_avg_pg)./abs(mag))'];
    dlmwrite('results/r25_err.dat', errmatrix, ...
        'delimiter', '\t', 'precision', 8);
end

%% Plot transfer function errors and magnitude for r = 25 reduced models.
fprintf(1, 'Plotting response magnitude and error for order r = 25 reduced models\n');
fprintf(1, '---------------------------------------------------------------------\n');

% Plot colors
ColMat = zeros(6,3);
ColMat(1,:) = [ 0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [  0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];
ColMat(6,:) = [1 0.4 0.6];

figure(1)
% Pointwise L_infty error
subplot(2,1,1)
semilogy(s_hz, abs(mag-mag_r25_irka)./abs(mag), '-.','color',ColMat(2,:),LineWidth=1);hold on
semilogy(s_hz, abs(mag-mag_r25_Linfty_g)./abs(mag),'--','color', ColMat(3,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_Linfty_pg)./abs(mag),'--','color', ColMat(4,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_avg_g)./abs(mag),'--.','color', ColMat(5,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_avg_pg)./abs(mag),'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

% Magnitude
subplot(2,1,2)
plot(s_hz, mag,'--o','color',ColMat(1,:),'markersize',4,LineWidth=1);hold on;
plot(s_hz, mag_r25_irka, '-.','color',ColMat(2,:),LineWidth=1);
plot(s_hz, mag_r25_Linfty_g,'--','color', ColMat(3,:),LineWidth=1);
plot(s_hz, mag_r25_Linfty_pg,'--','color', ColMat(4,:),LineWidth=1);
plot(s_hz, mag_r25_avg_g,'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, mag_r25_avg_pg,'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Full-order', 'r=25, IRKA', 'r=25 Linfty-g', 'r=25 Linfty-pg', ...
    'r=25 avg-g', 'r=25 avg-pg', 'Orientation', 'horizontal', 'Location', ...
    'southoutside', 'NumColumns', 3)

% Overwrite figure
saveas(figure(1), 'results/r25_plots.png')

%% Load base data for r = 50.
fprintf(1, 'Load problem data for order r = 50 reduced models\n');
fprintf(1, '-------------------------------------------------\n');

% order r = 50 reduced models
load('results/plateTVAlqo_r50_lqoirka.mat')
load('results/plateTVAlqo_r50_Linfty_g.mat')
load('results/plateTVAlqo_r50_Linfty_pg.mat')
load('results/plateTVAlqo_r50_avg_g.mat')
load('results/plateTVAlqo_r50_avg_pg.mat')

% All other data already loaded
fprintf(1, '\n');

%% Simulate rms response r = 50.
fprintf(1, 'Simulate rms response for order r = 50 reduced models\n');
fprintf(1, '-----------------------------------------------------\n');

% Allocate space for results
res_r50_irka      = zeros(1,length(s));
res_r50_Linfty_g  = zeros(1, length(s));
res_r50_Linfty_pg = zeros(1, length(s));
res_r50_avg_g     = zeros(1, length(s));
res_r50_avg_pg    = zeros(1, length(s));

fprintf(1, 'Beginning reduced order simulations\n')
fprintf(1, '-------------------------------\n');
overall_start = tic;

for ii=1:length(s)
    fprintf(1, 'Frequency step %d, f=%.2f Hz ...\n ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    tmp_r50_irka      = (s(ii)*Efo_r - Afo_r)\Bfo_r;
    tmp_r50_Linfty_g  = (s(ii)*Efo_r50_Linfty_g - Afo_r50_Linfty_g)\Bfo_r50_Linfty_g;
    tmp_r50_Linfty_pg = (s(ii)*Efo_r50_Linfty_pg - Afo_r50_Linfty_pg)\Bfo_r50_Linfty_pg;
    tmp_r50_avg_g     = (s(ii)*Efo_r50_avg_g - Afo_r50_avg_g)\Bfo_r50_avg_g;
    tmp_r50_avg_pg    = (s(ii)*Efo_r50_avg_pg - Afo_r50_avg_pg)\Bfo_r50_avg_pg;
    % Evaluate rms response
    res_r50_irka(ii)      = sqrt((tmp_r50_irka'*Qfo_r*tmp_r50_irka)/n_nodes);
    res_r50_Linfty_g(ii)  = sqrt((tmp_r50_Linfty_g'*Qfo_r50_Linfty_g*tmp_r50_Linfty_g)/n_nodes);
    res_r50_Linfty_pg(ii) = sqrt((tmp_r50_Linfty_pg'*Qfo_r50_Linfty_pg*tmp_r50_Linfty_pg)/n_nodes);
    res_r50_avg_g(ii)     = sqrt((tmp_r50_avg_g' * Qfo_r50_avg_g * tmp_r50_avg_g)/n_nodes);
    res_r50_avg_pg(ii)    = sqrt((tmp_r50_avg_pg'*Qfo_r50_avg_pg*tmp_r50_avg_pg)/n_nodes);
    fprintf(1, 'Current iteration of reduced order simulation finished in %.2f s\n',toc(current_iter))
    fprintf(1, '----------------------------------------------------------------------\n');
end
fprintf(1, 'Reduced order simulations finished in %.2f s\n', toc(overall_start))
fprintf(1, '--------------------------------------------------\n');

% Compute magnitudes to plot
% To plot magnitudes in absolute value, uncomment the code below
% mag_r50_irka      = abs(res_r50_irka);
% mag_r50_Linfty_g  = abs(res_r50_Linfty_g);
% mag_r50_Linfty_pg = abs(res_r50_Linfty_pg);
% mag_r50_avg_g     = abs(res_r50_avg_g);
% mag_r50_avg_pg    = abs(res_r50_avg_pg);

% To plot magnitudes in dB, uncomment the code below
mag_r50_irka      = 10*log10(abs(res_r50_irka)/1e-9);
mag_r50_Linfty_g  = 10*log10(abs(res_r50_Linfty_g)/1e-9);
mag_r50_Linfty_pg = 10*log10(abs(res_r50_Linfty_pg)/1e-9);
mag_r50_avg_g     = 10*log10(abs(res_r50_avg_g)/1e-9);
mag_r50_avg_pg    = 10*log10(abs(res_r50_avg_pg)/1e-9);
 
write = 1;
if write
    % Store data
    magmatrix = [s_hz', mag', mag_r50_irka', mag_r50_Linfty_g', mag_r50_Linfty_pg', ...
        mag_r50_avg_g', mag_r50_avg_pg'];
    dlmwrite('results/r50_mag.dat', magmatrix, ...
        'delimiter', '\t', 'precision', 8);
    errmatrix = [s_hz', (abs(mag-mag_r50_irka)./abs(mag))', (abs(mag-mag_r50_Linfty_g)./abs(mag))', ...
        (abs(mag-mag_r50_Linfty_pg)./abs(mag))', (abs(mag-mag_r50_avg_g)./abs(mag))', ...
        (abs(mag-mag_r50_avg_pg)./abs(mag))'];
    dlmwrite('results/r50_err.dat', errmatrix, ...
        'delimiter', '\t', 'precision', 8);
end

%% Plot transfer function errors and magnitude for r = 50 reduced models.
fprintf(1, 'Plotting response magnitude and error for order r = 50 reduced models\n');
fprintf(1, '---------------------------------------------------------------------\n');
figure(2)
% Pointwise L_infty error
subplot(2,1,1)
semilogy(s_hz, abs(mag-mag_r50_irka)./abs(mag), '-.','color',ColMat(2,:),LineWidth=1);hold on
semilogy(s_hz, abs(mag-mag_r50_Linfty_g)./abs(mag),'--','color', ColMat(3,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r50_Linfty_pg)./abs(mag),'--','color', ColMat(4,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r50_avg_g)./abs(mag),'--.','color', ColMat(5,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r50_avg_pg)./abs(mag),'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

% Magnitude
subplot(2,1,2)
plot(s_hz, mag,'--o','color',ColMat(1,:),'markersize',4,LineWidth=1);hold on;
plot(s_hz, mag_r50_irka, '-.','color',ColMat(2,:),LineWidth=1);
plot(s_hz, mag_r50_Linfty_g,'--','color', ColMat(3,:),LineWidth=1);
plot(s_hz, mag_r50_Linfty_pg,'--','color', ColMat(4,:),LineWidth=1);
plot(s_hz, mag_r50_avg_g,'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, mag_r50_avg_pg,'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Full-order', 'r=50, IRKA', 'r=50 Linfty-g', 'r=50 Linfty-pg', ...
    'r=50 avg-g', 'r=50 avg-pg', 'Orientation', 'horizontal', 'Location', ...
    'southoutside', 'NumColumns', 3)

% Overwrite figure
saveas(figure(2), 'results/r50_plots.png')

%% Load base data for r = 75.
fprintf(1, 'Load problem data for order r = 75 reduced models\n');
fprintf(1, '-------------------------------------------------\n');

% order r = 75 reduced models
load('results/plateTVAlqo_r75_lqoirka_restart2.mat')
load('results/plateTVAlqo_r75_Linfty_g.mat')
load('results/plateTVAlqo_r75_Linfty_pg.mat')
load('results/plateTVAlqo_r75_avg_g.mat')
load('results/plateTVAlqo_r75_avg_pg.mat')

% All other data already loaded
fprintf(1, '\n');

%% Simulate rms response r = 75.
fprintf(1, 'Simulate rms response for order r = 75 reduced models\n');
fprintf(1, '-----------------------------------------------------\n');

% Allocate space for results
res_r75_irka      = zeros(1,length(s));
res_r75_Linfty_g  = zeros(1, length(s));
res_r75_Linfty_pg = zeros(1, length(s));
res_r75_avg_g     = zeros(1, length(s));
res_r75_avg_pg    = zeros(1, length(s));

fprintf(1, 'Beginning reduced order simulations\n')
fprintf(1, '-------------------------------\n');
overall_start = tic;

for ii=1:length(s)
    fprintf(1, 'Frequency step %d, f=%.2f Hz ...\n ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    tmp_r75_irka      = (s(ii)*Efo_r - Afo_r)\Bfo_r;
    tmp_r75_Linfty_g  = (s(ii)*Efo_r75_Linfty_g - Afo_r75_Linfty_g)\Bfo_r75_Linfty_g;
    tmp_r75_Linfty_pg = (s(ii)*Efo_r75_Linfty_pg - Afo_r75_Linfty_pg)\Bfo_r75_Linfty_pg;
    tmp_r75_avg_g     = (s(ii)*Efo_r75_avg_g - Afo_r75_avg_g)\Bfo_r75_avg_g;
    tmp_r75_avg_pg    = (s(ii)*Efo_r75_avg_pg - Afo_r75_avg_pg)\Bfo_r75_avg_pg;
    % Evaluate rms response
    res_r75_irka(ii)      = sqrt((tmp_r75_irka'*Qfo_r*tmp_r75_irka)/n_nodes);
    res_r75_Linfty_g(ii)  = sqrt((tmp_r75_Linfty_g'*Qfo_r75_Linfty_g*tmp_r75_Linfty_g)/n_nodes);
    res_r75_Linfty_pg(ii) = sqrt((tmp_r75_Linfty_pg'*Qfo_r75_Linfty_pg*tmp_r75_Linfty_pg)/n_nodes);
    res_r75_avg_g(ii)     = sqrt((tmp_r75_avg_g' * Qfo_r75_avg_g * tmp_r75_avg_g)/n_nodes);
    res_r75_avg_pg(ii)    = sqrt((tmp_r75_avg_pg'*Qfo_r75_avg_pg*tmp_r75_avg_pg)/n_nodes);
    fprintf(1, 'Current iteration of reduced order simulation finished in %.2f s\n',toc(current_iter))
    fprintf(1, '----------------------------------------------------------------------\n');
end
fprintf(1, 'Reduced order simulations finished in %.2f s\n', toc(overall_start))
fprintf(1, '--------------------------------------------------\n');

% Compute magnitudes to plot
% To plot magnitudes in absolute value, uncomment the code below
% mag_r75_irka      = abs(res_r75_irka);
% mag_r75_Linfty_g  = abs(res_r75_Linfty_g);
% mag_r75_Linfty_pg = abs(res_r75_Linfty_pg);
% mag_r75_avg_g     = abs(res_r75_avg_g);
% mag_r75_avg_pg    = abs(res_r75_avg_pg);

% To plot magnitudes in dB, uncomment the code below
mag_r75_irka      = 10*log10(abs(res_r75_irka)/1e-9);
mag_r75_Linfty_g  = 10*log10(abs(res_r75_Linfty_g)/1e-9);
mag_r75_Linfty_pg = 10*log10(abs(res_r75_Linfty_pg)/1e-9);
mag_r75_avg_g     = 10*log10(abs(res_r75_avg_g)/1e-9);
mag_r75_avg_pg    = 10*log10(abs(res_r75_avg_pg)/1e-9);
 
write = 1;
if write
    % Store data
    magmatrix = [s_hz', mag', mag_r75_irka', mag_r75_Linfty_g', mag_r75_Linfty_pg', ...
        mag_r75_avg_g', mag_r75_avg_pg'];
    dlmwrite('results/r75_mag.dat', magmatrix, ...
        'delimiter', '\t', 'precision', 8);
    errmatrix = [s_hz', (abs(mag-mag_r75_irka)./abs(mag))', (abs(mag-mag_r75_Linfty_g)./abs(mag))', ...
        (abs(mag-mag_r75_Linfty_pg)./abs(mag))', (abs(mag-mag_r75_avg_g)./abs(mag))', ...
        (abs(mag-mag_r75_avg_pg)./abs(mag))'];
    dlmwrite('results/r75_err.dat', errmatrix, ...
        'delimiter', '\t', 'precision', 8);
end

%% Plot transfer function errors and magnitude for r = 75 reduced models.
fprintf(1, 'Plotting response magnitude and error for order r = 75 reduced models\n');
fprintf(1, '---------------------------------------------------------------------\n');
figure(3)
% Pointwise L_infty error
subplot(2,1,1)
semilogy(s_hz, abs(mag-mag_r75_irka)./abs(mag), '-.','color',ColMat(2,:),LineWidth=1);hold on
semilogy(s_hz, abs(mag-mag_r75_Linfty_g)./abs(mag),'--','color', ColMat(3,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r75_Linfty_pg)./abs(mag),'--','color', ColMat(4,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r75_avg_g)./abs(mag),'--.','color', ColMat(5,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r75_avg_pg)./abs(mag),'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

% Magnitude
subplot(2,1,2)
plot(s_hz, mag,'--o','color',ColMat(1,:),'markersize',4,LineWidth=1);hold on;
plot(s_hz, mag_r75_irka, '-.','color',ColMat(2,:),LineWidth=1);
plot(s_hz, mag_r75_Linfty_g,'--','color', ColMat(3,:),LineWidth=1);
plot(s_hz, mag_r75_Linfty_pg,'--','color', ColMat(4,:),LineWidth=1);
plot(s_hz, mag_r75_avg_g,'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, mag_r75_avg_pg,'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Full-order', 'r=75, IRKA', 'r=75 Linfty-g', 'r=75 Linfty-pg', ...
    'r=75 avg-g', 'r=75 avg-pg', 'Orientation', 'horizontal', 'Location', ...
    'southoutside', 'NumColumns', 3)

% Overwrite figure
saveas(figure(3), 'results/r75_plots.png')

%% Load base data for r = 100.
fprintf(1, 'Load problem data for order r = 100 reduced models\n');
fprintf(1, '--------------------------------------------------\n');

% order r = 100 reduced models
load('results/plateTVAlqo_r100_lqoirka_restart2.mat')
load('results/plateTVAlqo_r100_Linfty_g.mat')
load('results/plateTVAlqo_r100_Linfty_pg.mat')
load('results/plateTVAlqo_r100_avg_g.mat')
load('results/plateTVAlqo_r100_avg_pg.mat')

% All other data already loaded
fprintf(1, '\n');

%% Simulate rms response r = 100.
fprintf(1, 'Simulate rms response for order r = 100 reduced models\n');
fprintf(1, '-----------------------------------------------------\n');

% Allocate space for results
res_r100_irka      = zeros(1,length(s));
res_r100_Linfty_g  = zeros(1, length(s));
res_r100_Linfty_pg = zeros(1, length(s));
res_r100_avg_g     = zeros(1, length(s));
res_r100_avg_pg    = zeros(1, length(s));

fprintf(1, 'Beginning reduced order simulations\n')
fprintf(1, '-------------------------------\n');
overall_start = tic;

for ii=1:length(s)
    fprintf(1, 'Frequency step %d, f=%.2f Hz ...\n ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    tmp_r100_irka      = (s(ii)*Efo_r - Afo_r)\Bfo_r;
    tmp_r100_Linfty_g  = (s(ii)*Efo_r100_Linfty_g - Afo_r100_Linfty_g)\Bfo_r100_Linfty_g;
    tmp_r100_Linfty_pg = (s(ii)*Efo_r100_Linfty_pg - Afo_r100_Linfty_pg)\Bfo_r100_Linfty_pg;
    tmp_r100_avg_g     = (s(ii)*Efo_r100_avg_g - Afo_r100_avg_g)\Bfo_r100_avg_g;
    tmp_r100_avg_pg    = (s(ii)*Efo_r100_avg_pg - Afo_r100_avg_pg)\Bfo_r100_avg_pg;
    % Evaluate rms response
    res_r100_irka(ii)      = sqrt((tmp_r100_irka'*Qfo_r*tmp_r100_irka)/n_nodes);
    res_r100_Linfty_g(ii)  = sqrt((tmp_r100_Linfty_g'*Qfo_r100_Linfty_g*tmp_r100_Linfty_g)/n_nodes);
    res_r100_Linfty_pg(ii) = sqrt((tmp_r100_Linfty_pg'*Qfo_r100_Linfty_pg*tmp_r100_Linfty_pg)/n_nodes);
    res_r100_avg_g(ii)     = sqrt((tmp_r100_avg_g' * Qfo_r100_avg_g * tmp_r100_avg_g)/n_nodes);
    res_r100_avg_pg(ii)    = sqrt((tmp_r100_avg_pg'*Qfo_r100_avg_pg*tmp_r100_avg_pg)/n_nodes);
    fprintf(1, 'Current iteration of reduced order simulation finished in %.2f s\n',toc(current_iter))
    fprintf(1, '----------------------------------------------------------------------\n');
end
fprintf(1, 'Reduced order simulations finished in %.2f s\n', toc(overall_start))
fprintf(1, '--------------------------------------------------\n');

% Compute magnitudes to plot
% To plot magnitudes in absolute value, uncomment the code below
% mag_r100_irka      = abs(res_r100_irka);
% mag_r100_Linfty_g  = abs(res_r100_Linfty_g);
% mag_r100_Linfty_pg = abs(res_r100_Linfty_pg);
% mag_r100_avg_g     = abs(res_r100_avg_g);
% mag_r100_avg_pg    = abs(res_r100_avg_pg);

% To plot magnitudes in dB, uncomment the code below
mag_r100_irka      = 10*log10(abs(res_r100_irka)/1e-9);
mag_r100_Linfty_g  = 10*log10(abs(res_r100_Linfty_g)/1e-9);
mag_r100_Linfty_pg = 10*log10(abs(res_r100_Linfty_pg)/1e-9);
mag_r100_avg_g     = 10*log10(abs(res_r100_avg_g)/1e-9);
mag_r100_avg_pg    = 10*log10(abs(res_r100_avg_pg)/1e-9);

write = 1;
if write
    % Store data
    magmatrix = [s_hz', mag', mag_r100_irka', mag_r100_Linfty_g', mag_r100_Linfty_pg', ...
        mag_r100_avg_g', mag_r100_avg_pg'];
    dlmwrite('results/r100_mag.dat', magmatrix, ...
        'delimiter', '\t', 'precision', 8);
    errmatrix = [s_hz', (abs(mag-mag_r100_irka)./abs(mag))', (abs(mag-mag_r100_Linfty_g)./abs(mag))', ...
        (abs(mag-mag_r100_Linfty_pg)./abs(mag))', (abs(mag-mag_r100_avg_g)./abs(mag))', ...
        (abs(mag-mag_r100_avg_pg)./abs(mag))'];
    dlmwrite('results/r100_err.dat', errmatrix, ...
        'delimiter', '\t', 'precision', 8);
end

%% Plot transfer function errors and magnitude for r = 100 reduced models.
fprintf(1, 'Plotting response magnitude and error for order r = 100 reduced models\n');
fprintf(1, '----------------------------------------------------------------------\n');
figure(4)
% Pointwise L_infty error
subplot(2,1,1)
semilogy(s_hz, abs(mag-mag_r100_irka)./abs(mag), '-.','color',ColMat(2,:),LineWidth=1);hold on
semilogy(s_hz, abs(mag-mag_r100_Linfty_g)./abs(mag),'--','color', ColMat(3,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r100_Linfty_pg)./abs(mag),'--','color', ColMat(4,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r100_avg_g)./abs(mag),'--.','color', ColMat(5,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r100_avg_pg)./abs(mag),'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')

% Magnitude
subplot(2,1,2)
plot(s_hz, mag,'--o','color',ColMat(1,:),'markersize',4,LineWidth=1);hold on;
plot(s_hz, mag_r100_irka, '-.','color',ColMat(2,:),LineWidth=1);
plot(s_hz, mag_r100_Linfty_g,'--','color', ColMat(3,:),LineWidth=1);
plot(s_hz, mag_r100_Linfty_pg,'--','color', ColMat(4,:),LineWidth=1);
plot(s_hz, mag_r100_avg_g,'--.','color', ColMat(5,:),LineWidth=1);
plot(s_hz, mag_r100_avg_pg,'--.','color', ColMat(6,:),LineWidth=1);
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
legend('Full-order', 'r=100, IRKA', 'r=100 Linfty-g', 'r=100 Linfty-pg', ...
    'r=100 avg-g', 'r=100 avg-pg', 'Orientation', 'horizontal', 'Location', ...
    'southoutside', 'NumColumns', 3)

% Overwrite figure
saveas(figure(4), 'results/r100_plots.png')

%% Compute single error values
fprintf(1, 'Computing relative H2 errors\n')
fprintf(1, '-----------------------------\n')

fprintf(1, 'Order r = %d\n', 25)
fprintf(1, '--------------\n')
fprintf(1, 'Relative H2 error due to IRKA     : %.10f \n', sum((abs(mag-mag_r25_irka)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-g : %.10f \n', sum((abs(mag-mag_r25_Linfty_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-pg: %.10f \n', sum((abs(mag-mag_r25_Linfty_pg)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-g    : %.10f \n', sum((abs(mag-mag_r25_avg_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-pg   : %.10f \n', sum((abs(mag-mag_r25_avg_pg)))/sum(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, 'Order r = %d\n', 50)
fprintf(1, '--------------\n')
fprintf(1, 'Relative H2 error due to IRKA     : %.10f \n', sum((abs(mag-mag_r50_irka)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-g : %.10f \n', sum((abs(mag-mag_r50_Linfty_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-pg: %.10f \n', sum((abs(mag-mag_r50_Linfty_pg)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-g    : %.10f \n', sum((abs(mag-mag_r50_avg_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-pg   : %.10f \n', sum((abs(mag-mag_r50_avg_pg)))/sum(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, 'Order r = %d\n', 75)
fprintf(1, '--------------\n')
fprintf(1, 'Relative H2 error due to IRKA     : %.10f \n', sum((abs(mag-mag_r75_irka)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-g : %.10f \n', sum((abs(mag-mag_r75_Linfty_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-pg: %.10f \n', sum((abs(mag-mag_r75_Linfty_pg)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-g    : %.10f \n', sum((abs(mag-mag_r75_avg_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-pg   : %.10f \n', sum((abs(mag-mag_r75_avg_pg)))/sum(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, 'Order r = %d\n', 100)
fprintf(1, '--------------\n')
fprintf(1, 'Relative H2 error due to IRKA     : %.10f \n', sum((abs(mag-mag_r100_irka)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-g : %.10f \n', sum((abs(mag-mag_r100_Linfty_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to Linfty-pg: %.10f \n', sum((abs(mag-mag_r100_Linfty_pg)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-g    : %.10f \n', sum((abs(mag-mag_r100_avg_g)))/sum(abs(mag)))
fprintf(1, 'Relative H2 error due to avg-pg   : %.10f \n', sum((abs(mag-mag_r100_avg_pg)))/sum(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')


%%
fprintf(1, 'Computing relative Hinfty errors\n')
fprintf(1, '-----------------------------\n')

fprintf(1, 'Order r = %d\n', 25)
fprintf(1, '--------------\n')
fprintf(1, 'Relative Hinfty error due to IRKA     : %.10f \n', max((abs(mag-mag_r25_irka)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-g : %.10f \n', max((abs(mag-mag_r25_Linfty_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-pg: %.10f \n', max((abs(mag-mag_r25_Linfty_pg)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-g    : %.10f \n', max((abs(mag-mag_r25_avg_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-pg   : %.10f \n', max((abs(mag-mag_r25_avg_pg)))/max(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, 'Order r = %d\n', 50)
fprintf(1, '--------------\n')
fprintf(1, 'Relative Hinfty error due to IRKA     : %.10f \n', max((abs(mag-mag_r50_irka)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-g : %.10f \n', max((abs(mag-mag_r50_Linfty_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-pg: %.10f \n', max((abs(mag-mag_r50_Linfty_pg)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-g    : %.10f \n', max((abs(mag-mag_r50_avg_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-pg   : %.10f \n', max((abs(mag-mag_r50_avg_pg)))/max(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, 'Order r = %d\n', 75)
fprintf(1, '--------------\n')
fprintf(1, 'Relative Hinfty error due to IRKA     : %.10f \n', max((abs(mag-mag_r75_irka)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-g : %.10f \n', max((abs(mag-mag_r75_Linfty_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-pg: %.10f \n', max((abs(mag-mag_r75_Linfty_pg)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-g    : %.10f \n', max((abs(mag-mag_r75_avg_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-pg   : %.10f \n', max((abs(mag-mag_r75_avg_pg)))/max(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

fprintf(1, 'Order r = %d\n', 100)
fprintf(1, '--------------\n')
fprintf(1, 'Relative Hinfty error due to IRKA     : %.10f \n', max((abs(mag-mag_r100_irka)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-g : %.10f \n', max((abs(mag-mag_r100_Linfty_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to Linfty-pg: %.10f \n', max((abs(mag-mag_r100_Linfty_pg)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-g    : %.10f \n', max((abs(mag-mag_r100_avg_g)))/max(abs(mag)))
fprintf(1, 'Relative Hinfty error due to avg-pg   : %.10f \n', max((abs(mag-mag_r100_avg_pg)))/max(abs(mag)))
fprintf(1, '------------------------------------------------------------\n')

%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off
