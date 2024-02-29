%% RUNME
% Script file to run all experiments on the vibro-acoustic plate (plateTVA) 
% data set.

% This file is part of the archive Code, Data and Results for Numerical 
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
fullpath = matlab.desktop.editor.getActiveFilename;
[rootpath, filename, ~] = fileparts(fullpath(3:end));
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

% order r = 25 reduced models
load('results/plateTVAlqo_r25_lqoirka.mat')
load('results/plateTVAlqo_r25_Linfty_g.mat')
load('results/plateTVAlqo_r25_Linfty_pg.mat')
load('results/plateTVAlqo_r25_avg_g.mat')
load('results/plateTVAlqo_r25_avg_pg.mat')

% Simulation data for full-order plateTVA model
load('data/fosim_data.mat')
% Load C matrix from original model to get n_nodes used in calculating rms
% response
load('data/plateTVA_n201900m1q28278_fo.mat','C')
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
fprintf(1, '-------------------------------\n');
overall_start = tic;

for ii=1:length(s)
    fprintf(1, 'Frequency step %d, f=%.2f Hz ...\n ',ii,imag(s(ii))/2/pi)
    current_iter = tic;
    % Tmp solves
    tmp_r25_irka = (s(ii)*E_qo_r - A_qo_r)\B_qo_r;
    tmp_r25_Linfty_g = (s(ii)*E_qo_r25_Linfty_g - A_qo_r25_Linfty_g)\B_qo_r25_Linfty_g;
    tmp_r25_Linfty_pg = (s(ii)*E_qo_r25_Linfty_pg - A_qo_r25_Linfty_pg)\B_qo_r25_Linfty_pg;
    tmp_r25_avg_g = (s(ii)*E_qo_r25_avg_g - A_qo_r25_avg_g)\B_qo_r25_avg_g;
    tmp_r25_avg_pg = (s(ii)*E_qo_r25_avg_pg - A_qo_r25_avg_pg)\B_qo_r25_avg_pg;
    % Evaluate rms response
    res_r25_irka(ii) = sqrt((tmp_r25_irka'*Q_qo_r*tmp_r25_irka)/n_nodes);
    res_r25_Linfty_g(ii) = sqrt((tmp_r25_Linfty_g'*Q_qo_r25_Linfty_g*tmp_r25_Linfty_g)/n_nodes);
    res_r25_Linfty_pg(ii) = sqrt((tmp_r25_Linfty_pg'*Q_qo_r25_Linfty_pg*tmp_r25_Linfty_pg)/n_nodes);
    res_r25_avg_g(ii) = sqrt((tmp_r25_avg_g' * Q_qo_r25_avg_g * tmp_r25_avg_g)/n_nodes);
    res_r25_avg_pg(ii) = sqrt((tmp_r25_avg_pg'*Q_qo_r25_avg_pg*tmp_r25_avg_pg)/n_nodes);
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
    filename = 'results/r25_mag.dat';
    dlmwrite(filename, magmatrix, ...
        'delimiter', '\t', 'precision', 8);
    errmatrix = [s_hz', (abs(mag-mag_r25_irka)./abs(mag))', (abs(mag-mag_r25_Linfty_g)./abs(mag))', ...
        (abs(mag-mag_r25_Linfty_pg)./abs(mag))', (abs(mag-mag_r25_avg_g)./abs(mag))', ...
        (abs(mag-mag_r25_avg_pg)./abs(mag))'];
    filename = 'results/r25_err.dat';
    dlmwrite(filename, errmatrix, ...
        'delimiter', '\t', 'precision', 8);
end

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
axes('position', [.125 .15 .75 golden_ratio-1])

plot(s_hz, mag,'--o','color',ColMat(1,:),'markersize',4,LineWidth=1);hold on;
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
axes('position', [.125 .15 .75 golden_ratio-1])

semilogy(s_hz, abs(mag-mag_r25_irka)./abs(mag), '-.','color',ColMat(2,:),LineWidth=1);hold on
semilogy(s_hz, abs(mag-mag_r25_Linfty_g)./abs(mag),'--','color', ColMat(3,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_Linfty_pg)./abs(mag),'--.','color', ColMat(4,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_avg_g)./abs(mag),'--.','color', ColMat(5,:),LineWidth=1);
semilogy(s_hz, abs(mag-mag_r25_avg_pg)./abs(mag),'--.','color', ColMat(6,:),LineWidth=1);

legend('r=25, IRKA', 'r=25 Linfty-g', 'r=25 Linfty-pg', 'r=25 avg-g', 'r=25 avg-pg')
diary off