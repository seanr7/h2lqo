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

load('plateTVA_n201900m1q28278')
n_nodes = full(sum(sum(C)));
