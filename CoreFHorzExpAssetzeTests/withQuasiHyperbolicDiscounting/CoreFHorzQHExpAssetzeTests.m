% Tests of the core VFI Toolkit FHorz ExpAssetze commands under QUASI-HYPERBOLIC discounting.
% Runs the CoreFHorzExpAssetzeTests models (Naive + Sophisticated), plus exponential cross-tests:
%   (i)   Naive: the continuation value equals the exponential value function (any beta0)
%   (ii)  Naive, beta0=1: both the main value fn AND the continuation value equal exponential
%   (iii) Sophisticated, beta0=1: both the main value fn AND the continuation value equal exponential
%
% experienceassetze: aprime depends on (d2,a2,z,e), so z and e are ALWAYS present -> every subcode
% is z_e (there is no z_noe). QH+ExpAssetze supports the BASELINE only for a single a1 (there are no
% DC1/GI1 QH-ExpAssetze raws; DC/GI require two standard assets -- see the _with2A1 subcodes).
% No QH+ExpAssetze+semiz family exists.
%
% This suite lives in the parent CoreFHorzExpAssetzeTests/ folder and reuses its setup and ReturnFns.

addpath('../CoreFHorzExpAssetzeTests_Setup/')
addpath('../CoreFHorzExpAssetze_ReturnFns/')
addpath('./CoreFHorzQHExpAssetzeTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzQHExpAssetzeTests_subcodes/With2A1_subcodes/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetze_setup

Params.beta0=0.9; % additional today-tomorrow (present-bias) discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};

%% withA1 (single standard asset): baseline + lowmemory, Naive & Sophisticated
%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=1;
output=CoreFHorzQHExpAssetze_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=2;
output=CoreFHorzQHExpAssetze_d1_z_e_nosemiz_withA1(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1 (two standard assets): base + DC2A + GI2A + DC2A_GI2A, Naive & Sophisticated.
% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1). The second
% standard endogenous state a1_2 is BINARY (a capped high-return asset). a1main is kept modest.
% THESE EXERCISE THE 12 (currently non-compliant) QH DC2A/GI2A/DC2A_GI2A ExpAssetze raws.

n_a_2A1=[51,13]; % [a1_1, a2]; the binary a1_2 is added inside the subcode -> [51,2,13]
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3; a2_grid];
n_a_2A1_notsobig=[151,13];
a_grid_2A1_notsobig=[5*linspace(0,1,n_a_2A1_notsobig(1))'.^3; a2_grid];

%% with2A1, without d1, with z, with e
figure_c=3;
output=CoreFHorzQHExpAssetze_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e
figure_c=4;
output=CoreFHorzQHExpAssetze_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
