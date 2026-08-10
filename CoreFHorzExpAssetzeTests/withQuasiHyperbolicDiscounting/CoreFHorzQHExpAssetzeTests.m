% Tests of the core VFI Toolkit FHorz ExpAssetze commands under QUASI-HYPERBOLIC discounting.
% Runs the CoreFHorzExpAssetzeTests models (Naive + Sophisticated), plus exponential cross-tests:
%   (i)   Naive: the continuation value equals the exponential value function (any beta0)
%   (ii)  Naive, beta0=1: both the main value fn AND the continuation value equal exponential
%   (iii) Sophisticated, beta0=1: both the main value fn AND the continuation value equal exponential
%
% experienceassetze: aprime depends on (d2,a2,z,e), so z and e are ALWAYS present -> every subcode
% is z_e (there is no z_noe). Nosemiz: single a1 supports the baseline only (DC/GI require two
% standard assets -- see the _with2A1 subcodes). The QH+ExpAssetze+semiz family (figs 3,4) runs
% single a1, base/DC1/GI1/DC1_GI1 x {d1,nod1} x {Naive,Sophisticated}, with the full lowmemory
% ladder. The QH+ExpAssetze+semiz+with2A1 family (figs 7,8) runs base/DC2A/GI2A/DC2A_GI2A
% x {d1,nod1} x {Naive,Sophisticated}, also with the full lowmemory ladder.
%
% This suite lives in the parent CoreFHorzExpAssetzeTests/ folder and reuses its setup and ReturnFns.

addpath('../CoreFHorzExpAssetzeTests_Setup/')
addpath('../CoreFHorzExpAssetze_ReturnFns/')
addpath('./CoreFHorzQHExpAssetzeTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzQHExpAssetzeTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzQHExpAssetzeTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzQHExpAssetzeTests_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAssetze_ReturnFns/Semiz_ReturnFns/')

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

%% Semiz variant (STAGE 2a): QH+semiz+experienceassetze, d1 BASE method only
%% without d1, with z, with e, with semiz
n_a_notsobig=[151,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=3;
output=CoreFHorzQHExpAssetze_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, with semiz
figure_c=4;
output=CoreFHorzQHExpAssetze_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);



%% with2A1 (two standard assets): base + DC2A + GI2A + DC2A_GI2A, Naive & Sophisticated.
% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1). The second
% standard endogenous state a1_2 is BINARY (a capped high-return asset). a1main is kept modest.
% THESE EXERCISE THE 12 (currently non-compliant) QH DC2A/GI2A/DC2A_GI2A ExpAssetze raws.

n_a_2A1=[51,13]; % [a1_1, a2]; the binary a1_2 is added inside the subcode -> [51,2,13]
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3; a2_grid];
n_a_2A1_notsobig=[151,13];
a_grid_2A1_notsobig=[5*linspace(0,1,n_a_2A1_notsobig(1))'.^3; a2_grid];

%% with2A1, without d1, with z, with e
figure_c=5;
output=CoreFHorzQHExpAssetze_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e
figure_c=6;
output=CoreFHorzQHExpAssetze_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);



%% with2A1 + semiz: QH+semiz+experienceassetze with TWO standard endogenous assets
% Exercises the QH ExpAssetzeSemiExo DC2A / GI2A / DC2A_GI2A code paths (length(n_a1)>1).
% Reuses the with2A1 grids (n_a_2A1 etc, the binary a1_2 is added inside the subcode).

%% with2A1, without d1, with z, with e, with semiz
figure_c=7;
output=CoreFHorzQHExpAssetze_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e, with semiz
figure_c=8;
output=CoreFHorzQHExpAssetze_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
