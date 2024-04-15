%%
% Script to test lqoirka on a toy problem
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

%% Load base data.
% Convert second order mass spring damper system to first order linear
% quadratic output system.

fprintf(1, 'Creating toy model...\n')
fprintf(1, '-------------------------\n');

n1              = 10; 
alpha           = .002; 
beta            = alpha; 
v               = 5;
[Mso, Dso, Kso] = triplechain_MSD(n1, alpha, beta, v);

[n, ~] = size(Mso);
cso    = ones(1, n);
bso    = ones(n, 1);

%% Convert plate model to first-order from second-order.
fprintf(1, 'Converting second-order realization to first-order linear quadratic output system\n')
tic

Efo                   = spalloc(2*n, 2*n, nnz(Mso) + n); % Descriptor matrix; Efo = [I, 0: 0, Mso]
Efo(1:n, 1:n)         = speye(n);                        % (1, 1) block
Efo(n+1:2*n, n+1:2*n) = Mso;                               % (2, 2) block is (sparse) mass matrix

Afo                   = spalloc(2*n, 2*n, nnz(Kso) + nnz(Dso) + n); % Afo = [0, I; -Kso, -Dso]
Afo(1:n, n+1:2*n)     = speye(n);                                   % (1, 2) block of Afo
Afo(n+1:2*n, 1:n)     = -Kso;                                       % (2, 1) block is -Kso
Afo(n+1:2*n, n+1:2*n) = -Dso;                                       % (2, 2) block is -Dso 

bfo             = spalloc(2*n, 1, nnz(bso)); % bfo = [0; bso];
bfo(n+1:2*n, :) = bso;       

% Our quadratic output matrix is C' * C
Qfo           = spalloc(2*n, 2*n, nnz(cso' * cso));
Qfo(1:n, 1:n) = cso' * cso; 
fprintf(1, 'First-order realization built in %.2f s\n',toc)
fprintf(1, '--------------------------------------------\n');
% Will use for comparison, and checking interpolation conditions

%% Compute reduced-order model and check interpolation conditions.
r = 8;
% Opts
opts          = struct();
opts.plotconv = 1;
opts.tol      = 10e-7;
[Efo_r, Afo_r, bfo_r, Qfo_r, info]    = lqoirka(Mso, Dso, Kso, bso, Qfo, r, opts);
pole_hist                             = info.pole_hist;    
conv_nodes                            = pole_hist(:, end-1);
convres                               = info.sores;

% Check interpolation conditions
H2  = @(s1, s2) bfo' * (((s1 * Efo - Afo)')\Qfo) * ((s2 * Efo - Afo)\bfo);
H2r = @(s1, s2) bfo_r' * (((s1 * Efo_r - Afo_r)')\Qfo_r) * ((s2 * Efo_r - Afo_r)\bfo_r);

fprintf('2nd-order optimality conditions, %d^2 in total\n', r)
fprintf('----------------------------------------------\n')
for i = 1:r
    for j = 1:r
        H2(-(conv_nodes(i)),  -(conv_nodes(j))) - H2r(-(conv_nodes(i)), -(conv_nodes(j)))
    end
    
end
 
% Partial deriv wrt first argument of H2, H2r
H2_prime_s1  = @(s1, s2) -bfo' * (((s1 * Efo - Afo)'\Efo') * ((s1 * Efo - Afo)'\Qfo)) * ((s2 * Efo - Afo)\bfo); 
H2_prime_s2  = @(s1, s2) -bfo' * ((s1 * Efo - Afo)'\Qfo) * ((s2 * Efo - Afo)\Efo) * ((s2 * Efo - Afo)\bfo); 
H2r_prime_s1 = @(s1, s2) -bfo_r' * (((s1 * Efo_r - Afo_r)'\Efo_r') * ((s1 * Efo_r - Afo_r)'\Qfo_r)) * ((s2 * Efo_r- Afo_r)\bfo_r); 
H2r_prime_s2 = @(s1, s2) -bfo_r' * ((s1 * Efo_r - Afo_r)'\Qfo_r) * (((s2 * Efo_r - Afo_r)\Efo_r) * ((s2 * Efo_r - Afo_r)\bfo_r)); 
% 
fprintf('Hermite quadratic optimality conditions, %d in total\n', r)
fprintf('----------------------------------------------------\n')

for i = 1:r
    ro_side = 0;
    fo_side = 0;
    for j = 1:r
        ro_side = ro_side + 2 * convres(i, j) * H2r_prime_s1(-conv_nodes(i), -conv_nodes(j));
        fo_side = fo_side + 2 * convres(i, j) * H2_prime_s1(-conv_nodes(i), -conv_nodes(j));
    end
    fo_side - ro_side
end
