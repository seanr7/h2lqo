%% RUNME_ADVECDIFF
% Script file to run all experiments on advection diffusion model in
% "$\mathcal{H}_2%-optimal model reduction of linear systems with multiple
% quadratic outputs"
% Copyright (c) 2024 Sean Reiter
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
addpath([rootpath, '/QBDynamicsQBOutput_MATLAB'])

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
fprintf(1, 'Load problem data\n');
fprintf(1, '-----------------\n');

% TODO: Check best practice for using someone else's code in code you want
% to publish ... 
nx   = 300;  % No. of spatial grid points
diff = 1e-2; % Diffusion parameter
adv  = 1;    % Advection parameter

% Advection-diffusion problem taken from [Diaz et al., 2024]
[E, A, B, ~] = AdvDiff_Probgen_1D(nx, adv, diff);  
E = full(E);     % E = I in this example
A = full(A);    
B = full(B);   
Clin = eye(nx, nx); % Linear output matrix
Mquad = eye(nx, nx); % Quadratic output matrix

fprintf(1, '\n');

%% Run ode15 and simulate output.
fprintf(1, 'Solve AdvecDiff problem via ode15\n');
fprintf(1, '---------------------------------\n');

Tfin = 2;               % final time for simulation
nt   = round(Tfin*128); % number of time steps

% Set system inputs 
u0 = @(t) 0.5*(cos(pi*t)+1); % Dirichlet input on left boundary
u1 = @(t) zeros(size(t));    % Neumann input on right boundary

ode_rtol = 1e-6; 
tsteps   = linspace(0, Tfin, nt+1); % Time-steps
v0       = zeros(nx, 1);            % Initial condition 
options  = odeset('AbsTol',1.e-2/nx^2, 'RelTol',ode_rtol, ...
                  'Mass', E, 'Jacobian', A, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiff = @(t,y)(A*y + B*[u0(t);u1(t)]);
[t, v]   = ode15s(fAdvDiff, tsteps, v0, options);
% Note, v is nt \times nx
% Compute quadratic cost by 1/(2*n) * ||x(t) - \bf{1}||_2^2
h = 1/nx;
output = (h/2) * sum((Clin*v' - ones(nx,1)*ones(1,length(t)) ).^2, 1);

% plot state
xplot = linspace(0, 1, nx+1); 
vfull = [u0(t), v];  % Adds back in Dirichlet BC data for full state
% figure(1)
% mesh(xplot, t, vfull); 
% xlabel('$x$'); 
% ylabel('$t$')
% zlabel('$v(x, t)$')
% title('State (original)')

% fprintf(1, '\n');
% figure(2)
% plot(t, output, '-'); hold on
% xlabel('$t$'); 
% ylabel('$z(t)$')
% hold on 

fprintf(1, 'Simulate full output\n');
fprintf(1, '--------------------\n');
% Output has mixed linear and quadratic terms
Clin = -h*ones(1, nx);
y    = Clin*v';
for tt = 1:length(t)
    y(:,tt) = y(:,tt) + (h/2)*v(tt, :)*Mquad*v(tt, :)';
end
y = y + (1/2) * ones(1, nt+1);

figure(1)
% Define colormat
ColMat(1,:) = [ 0.8500    0.3250    0.0980];
ColMat(2,:) = [0.3010    0.7450    0.9330];
ColMat(3,:) = [  0.9290    0.6940    0.1250];
ColMat(4,:) = [0.4660    0.6740    0.1880];
ColMat(5,:) = [0.4940    0.1840    0.5560];


% Make aspect ration `golden'
figure
% golden_ratio = (sqrt(5)+1)/2;
% axes('position', [.125 .15 .75 golden_ratio-1])
subplot(2,1,1)
plot(t, y, '-','color',ColMat(1,:), LineWidth=1.5)
hold on
grid on
% xlabel('$t$','interpreter','latex'); 
% ylabel('$y(t)$','interpreter','latex')

%% Run algorithm.
fprintf(1, 'Computed reduced-order models via LQO-TSIA\n');
fprintf(1, '------------------------------------------\n');

r = 25;
[Er, Ar, Br, Clinr, Mquadr, poles] = mimolqo_tsia(E, A, B, Clin, Mquad, r);
%%
pole_hist = poles;
maxiter = max(size(poles));
for k = 2:maxiter
    err(k) = abs(max(poles(:, k)-poles(:,k-1)));
end

figure
golden_ratio = (sqrt(5)+1)/2;
axes('position', [.125 .15 .75 golden_ratio-1])
semilogy([1:maxiter], err, '-o', 'color',ColMat(3,:), LineWidth=1.5)
xlim([2,maxiter])
xlabel('$k$', 'interpreter', 'latex')
%% Run ode15 and simulate reduced output.
fprintf(1, 'Solve reduced AdvecDiff problem via ode15\n');
fprintf(1, '-----------------------------------------\n');

ode_rtol = 1e-6; 
nxr      = r;
vr0      = zeros(nxr, 1); % Initial condition 
options  = odeset('AbsTol',1.e-2/nxr^2, 'RelTol',ode_rtol, ...
                  'Mass', Er, 'Jacobian', Ar, ...
                  'MStateDependence', 'none', 'Stats','on');
fAdvDiffred = @(tr,yr)(Ar*yr + Br*[u0(tr);u1(tr)]);
[tr, vr]    = ode15s(fAdvDiffred, tsteps, vr0, options);
% Note, vr is nt \times nxr

fprintf(1, 'Simulate reduced output\n');
fprintf(1, '-----------------------\n');
yr = Clinr*vr';
for tt = 1:length(tr)
    yr(:,tt) = yr(:,tt) + (h/2)*vr(tt, :)*Mquadr*vr(tt, :)';
end
yr = yr + (1/2) * ones(1, nt+1);

plot(t, yr, '-.', 'color',ColMat(3,:), LineWidth=1.5); hold on
lgd = legend('$y(t)$', '$\hat{y}(t)$','interpreter','latex','FontName','Arial',...
    'location', 'northeast');
fontsize(lgd,10,'points')

subplot(2,1,2)
plot(t, abs(y - yr), '-','color',ColMat(3,:),LineWidth=1.5); hold on
% xlabel('$t$'); 
% ylabel('$y(t)-\hat{y}(t)$','interpreter','latex')
% figure(3)
% [~, totaliter] = size(poles);   totaliter = totaliter - 1;
% err = max(abs(poles - poles))
% plot()
%%

% %% Hierarchy
% addpath('supReiW24/drivers/')
% % load('data/heat-cont.mat')
% % % A = full(A);    B = full(B);    Clin = full(C);    E = eye(200,200);
% % Clin = C;
% % [n,~] = size(A);    p = 1; % relevant dimensions
% % B = B(:, p);  Clin = Clin(p, :); % If using ISS model, make SISO
% % Mquad = (diag(2*ones(1,n)) + diag(-1*ones(1,n-1),1) + diag(-1*ones(1,n-1),-1));
% % Mquad = Mquad*10e-2;
% % Elin = speye(200, 200);
% % nx = 200;
% B = B(:,1);
% for r = 2:2:20  
%     fprintf(1, 'r=%d', r)
%     [Er, Ar, Br, Clinr, Mquadr, poles] = sisolqo_irka(E, A, B, Clin, Mquad, r);
%     % [Er, Ar, Br, Clinr, Mquadr, poles] = mimolqo_tsia(E, A, B, Clin, Mquad, r);
%     Aerr = [A, zeros(nx, r); zeros(r, nx), Ar]; Berr = [B; Br];
%     Cerr = [Clin, - Clinr]; Merr = [Mquad, zeros(nx, r); zeros(r, nx), -Mquadr];
%     Perr = lyap(Aerr, Berr*Berr');
%     Qerr = lyap(Aerr', Cerr'*Cerr + Merr*Perr*Merr);
%     k=r/2;  H2errs(k) = sqrt(trace(Berr'*Qerr*Berr));
%     H2errs(k)
% end
% %%
% plot([2, 4, 6, 8, 10, 12, 14, 16, 18, 20],abs(H2errs))