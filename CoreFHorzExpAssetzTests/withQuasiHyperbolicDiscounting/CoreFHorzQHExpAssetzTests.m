% Tests of the core VFI Toolkit FHorz ExpAssetz commands under QUASI-HYPERBOLIC discounting.
% Runs the same models as CoreFHorzExpAssetzTests, but with quasi-hyperbolic discounting -- each
% subcode is run once for Naive and once for Sophisticated, with a continuation-value (Valt) check
% beside every value-fn (V) check. Each subcode also ends with cross-tests versus exponential
% discounting:
%   (i)   Naive: the continuation value equals the exponential value function (any beta0)
%   (ii)  Naive, beta0=1: both the main value fn AND the continuation value equal exponential
%   (iii) Sophisticated, beta0=1: both the main value fn AND the continuation value equal exponential
%
% NOTE ON COVERAGE: QH+experienceassetz supports the BASELINE for a single standard asset a1, and
% the DC2A/GI2A/DC2A_GI2A paths only with TWO standard endogenous assets (length(n_a1)>1). There is
% no QH+experienceassetz+semiz family. So this suite has:
%  - withA1 subcodes (single a1): baseline + lowmemory, Naive & Sophisticated;
%  - with2A1 subcodes (two a1, a binary a1_2 spliced in): base + DC2A + GI2A + DC2A_GI2A + lowmemory.
% No semiz variants.

% This suite lives in the parent CoreFHorzExpAssetzTests/ folder and reuses its setup and
% ReturnFns (same models -- QH only changes discounting), so those are one level up ('../').
addpath('../CoreFHorzExpAssetzTests_Setup/')
addpath('../CoreFHorzExpAssetz_ReturnFns/')
% Cross-tests reference experienceassete ReturnFns (in the sibling CoreFHorzExpAsseteTests/, two levels up)
addpath('../../CoreFHorzExpAsseteTests/CoreFHorzExpAssete_ReturnFns/')

addpath('./CoreFHorzQHExpAssetzTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzQHExpAssetzTests_subcodes/With2A1_subcodes/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetz_setup

Params.beta0=0.9; % additional today-tomorrow (present-bias) discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};

%% withA1 (single standard asset): baseline + lowmemory, Naive & Sophisticated
%% without d1, with z, without e
figure_c=1;
output=CoreFHorzQHExpAssetz_nod1_z_noe_nosemiz_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e
figure_c=2;
output=CoreFHorzQHExpAssetz_d1_z_noe_nosemiz_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e
figure_c=3;
output=CoreFHorzQHExpAssetz_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
figure_c=4;
output=CoreFHorzQHExpAssetz_d1_z_e_nosemiz_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);


%% with2A1 (two standard endogenous assets): base + DC2A + GI2A + DC2A_GI2A + lowmemory, Naive & Sophisticated
% Splice a binary second standard asset a1_2 in each subcode -> n_a=[a1_1,2,a2] -> length(n_a1)>1 -> DC2A/GI2A path
n_a_2A1=[51,n_a_justexpasset]; % [a1_1, a2]; the binary a1_2 is added inside the subcode -> [51,2,13]
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3; a2_grid];
n_a_2A1_notsobig=[151,n_a_justexpasset];
a_grid_2A1_notsobig=[5*linspace(0,1,n_a_2A1_notsobig(1))'.^3; a2_grid];

%% with2A1, without d1, with z, without e
figure_c=5;
output=CoreFHorzQHExpAssetz_nod1_z_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, without e
figure_c=6;
output=CoreFHorzQHExpAssetz_d1_z_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, without d1, with z, with e
figure_c=7;
output=CoreFHorzQHExpAssetz_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e
figure_c=8;
output=CoreFHorzQHExpAssetz_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
