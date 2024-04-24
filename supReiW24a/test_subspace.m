%%
% Script to test second order structured solves
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

%% Test difference in subspaces.
% Load bases computed with first-order solves.
load('results/prim_bases_g_exactsolves.mat', 'Vprim')
load('results/prim_bases_pg_exactsolves.mat', 'Wprim')
Vprim_exact = Vprim;
Wprim_exact = Wprim;

% Load bases computed with second-order solves.
load('results/prim_bases_g.mat', 'Vprim_g')
load('results/prim_bases_pg.mat', 'Wprim_pg')
Vprim = Vprim_g;
Wprim = Wprim_pg;

r = 25;
% Orthogonalize via pivoted QR and compare. 
compression = tic;
fprintf(1, 'Computing orthonormalized model reduction bases via pivoted QR\n')
fprintf(1, '--------------------------------------------------------------\n')
[Vorth_exact, ~, pV_exact] = qr(Vprim_exact, 'vector', 'econ');   
[Worth_exact, ~, pW_exact] = qr(Wprim_exact, 'vector', 'econ');
% Grab r leading columns orthonormal columns from pivoted QR of Vprim
% and Wprim
Vorth_exact = Vorth_exact(:, 1:r);   Worth_exact = Worth_exact(:, 1:r); 
fprintf(1, 'Vorth_exact and Worth_exact computed in %.2f s\n', toc(compression))
fprintf(1, '----------------------------------\n')

% Orthogonalize via pivoted QR and compare. 
compression = tic;
fprintf(1, 'Computing orthonormalized model reduction bases via pivoted QR\n')
fprintf(1, '--------------------------------------------------------------\n')
[Vorth, ~, pV] = qr(Vprim, 'vector', 'econ');   
[Worth, ~, pW] = qr(Wprim, 'vector', 'econ');
% Grab r leading columns orthonormal columns from pivoted QR of Vprim
% and Wprim
Vorth = Vorth(:, 1:r);   Worth = Worth(:, 1:r); 
fprintf(1, 'Vorth_exact and Worth_exact computed in %.2f s\n', toc(compression))
fprintf(1, '----------------------------------\n')

fprintf(1, 'Relative error in Vorth: %.16f \n', norm(Vorth - Vorth_exact,2)/norm(Vorth_exact,2))
fprintf(1, 'Relative error in Vprim: %.16f \n', norm(Vprim - Vprim_exact,2)/norm(Vprim_exact,2))
fprintf(1, 'Relative error in Worth: %.16f \n', norm(Worth - Worth_exact,2)/norm(Worth_exact,2))
fprintf(1, 'Relative error in Wprim: %.16f \n', norm(Wprim - Wprim_exact,2)/norm(Wprim_exact,2))

for i = 1:250
fprintf(1, 'Relative error in column i of Vorth: %.16f \n', norm(Vorth(:, i) - Vorth_exact(:, i),2)/norm(Vorth_exact,2))
fprintf(1, 'Relative error in column i of Vprim: %.16f \n', norm(Vprim(:, i) - Vprim_exact(:, i),2)/norm(Vprim_exact,2))
fprintf(1, 'Relative error in column i of Worth: %.16f \n', norm(Worth(:, i) - Worth_exact(:, i),2)/norm(Worth_exact,2))
fprintf(1, 'Relative error in column i of Wprim: %.16f \n', norm(Wprim(:, i) - Wprim_exact(:, i),2)/norm(Wprim_exact,2))
end
%% Finished script.
fprintf(1, 'FINISHED SCRIPT.\n');
fprintf(1, '================\n');
fprintf(1, '\n');

diary off