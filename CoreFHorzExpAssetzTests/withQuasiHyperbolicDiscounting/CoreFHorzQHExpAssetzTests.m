% Tests of the core VFI Toolkit FHorz ExpAssetz commands under QUASI-HYPERBOLIC discounting.
% Runs the same models as CoreFHorzExpAssetzTests, but with quasi-hyperbolic discounting -- each
% subcode is run once for Naive and once for Sophisticated, with a continuation-value (Valt) check
% beside every value-fn (V) check. Each subcode also ends with cross-tests versus exponential
% discounting:
%   (i)   Naive: the continuation value equals the exponential value function (any beta0)
%   (ii)  Naive, beta0=1: both the main value fn AND the continuation value equal exponential
%   (iii) Sophisticated, beta0=1: both the main value fn AND the continuation value equal exponential
%
% COVERAGE: this suite now mirrors the baseline CoreFHorzExpAssetzTests one-for-one, figs 1-24 in
% the same order:
%   figs 1-4:   noa1, nosemiz  (experience asset a2 is the only endogenous state)
%   figs 5-8:   noa1, semiz
%   figs 9-12:  withA1, nosemiz (single standard asset a1)
%   figs 13-16: withA1, semiz
%   figs 17-20: with2A1, nosemiz (two standard endogenous assets; a binary a1_2 spliced in)
%   figs 21-24: with2A1, semiz
% Within each group the order is: nod1_z_noe, d1_z_noe, nod1_z_e, d1_z_e.
%
% Which solution methods each subcode runs (mirroring its baseline counterpart):
%   noa1 (figs 1-8):     base only (no a1 for divide-and-conquer / grid interpolation to act on)
%   withA1 nosemiz (9-12): base only (there are no DC1/GI1 QH-ExpAssetz raws)
%   withA1 semiz (13-16):  base + DC1 + GI1 + DC1_GI1
%   with2A1 (17-24):       base + DC2A + GI2A + DC2A_GI2A
%
% STATUS (2026-08-19): COMPLETE. This bank was written test-first; the toolkit has now caught up
% on all of it. ExperienceAssetz carries 104 QH raws (40 nosemiz + 64 semiz), the
% QuasiHyperbolicExpAssetzSemiExo dispatcher and its {DC,GI,DC_GI} sub-dispatchers all exist, and
% no 'not yet implemented' error survives anywhere in the family.
%
% The 2026-08-19 run: 1268 checks, 24 of 24 figures, every one zero, no errors. That covers all
% four solver tiers x nod1/with-d1 x no-e/with-e x Naive/Sophisticated, including the with-e 2A
% paths at lowmemory 3 and the beta0=1 degeneracy checks that collapse QH onto the exponential.
% (Toolkit commit 2f291d33.)
%
% The subcodes here take figure_c but draw no figures (alternative preferences change nothing beyond
% Policy, so there are no StationaryDist/AllStats/LifeCycleProfiles/SimPanel blocks), so only the
% diary is saved.

%% Diary of the command window output (written to the parent bank's TestOutput folder)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzQHExpAssetzTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzQHExpAssetzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzQHExpAssetzTestsdiary.txt

% This suite lives in the parent CoreFHorzExpAssetzTests/ folder and reuses its setup and
% ReturnFns (same models -- QH only changes discounting), so those are one level up ('../').
addpath('../CoreFHorzExpAssetzTests_Setup/')
addpath('../CoreFHorzExpAssetz_ReturnFns/')
addpath('../CoreFHorzExpAssetz_ReturnFns/Noa1_ReturnFns/')
addpath('../CoreFHorzExpAssetz_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzExpAssetz_ReturnFns/Semiz_ReturnFns/')
% Cross-tests reference experienceassete ReturnFns (in the sibling CoreFHorzExpAsseteTests/, two levels up)
addpath('../../CoreFHorzExpAsseteTests/CoreFHorzExpAssete_ReturnFns/')

addpath('./CoreFHorzQHExpAssetzTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzQHExpAssetzTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzQHExpAssetzTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzQHExpAssetzTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzQHExpAssetzTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzQHExpAssetzTests_subcodes/With2A1_subcodes/Semiz_subcodes/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetz_setup

Params.beta0=0.9; % additional today-tomorrow (present-bias) discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};

%% ================= noa1 (figs 1-8): experience asset a2 is the ONLY endogenous state =================
% Base method only (no DC/GI/DC+GI -- irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.
% Previously test-first; the toolkit now supports QH+experienceassetz noa1, and these PASS.

%% without d1, with z, without e, noa1, nosemiz
figure_c=1;
output=CoreFHorzQHExpAssetz_nod1_z_noe_nosemiz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, noa1, nosemiz
figure_c=2;
output=CoreFHorzQHExpAssetz_d1_z_noe_nosemiz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, noa1, nosemiz
figure_c=3;
output=CoreFHorzQHExpAssetz_nod1_z_e_nosemiz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, noa1, nosemiz
figure_c=4;
output=CoreFHorzQHExpAssetz_d1_z_e_nosemiz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, without e, noa1, semiz
figure_c=5;
output=CoreFHorzQHExpAssetz_nod1_z_noe_semiz_noa1(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, noa1, semiz
figure_c=6;
output=CoreFHorzQHExpAssetz_d1_z_noe_semiz_noa1(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, noa1, semiz
figure_c=7;
output=CoreFHorzQHExpAssetz_nod1_z_e_semiz_noa1(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, noa1, semiz
figure_c=8;
output=CoreFHorzQHExpAssetz_d1_z_e_semiz_noa1(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);


%% ================= withA1 (figs 9-16): single standard endogenous asset a1 =================
% nosemiz (figs 9-12): baseline method + lowmemory, Naive & Sophisticated.

%% without d1, with z, without e
figure_c=9;
output=CoreFHorzQHExpAssetz_nod1_z_noe_nosemiz_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e
figure_c=10;
output=CoreFHorzQHExpAssetz_d1_z_noe_nosemiz_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e
figure_c=11;
output=CoreFHorzQHExpAssetz_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
figure_c=12;
output=CoreFHorzQHExpAssetz_d1_z_e_nosemiz_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% withA1 + semiz (figs 13-16): base + DC1 + GI1 + DC1_GI1 + lowmemory, Naive & Sophisticated
% Previously test-first; the toolkit now supports QH+experienceassetz+semiz, and these PASS.

%% without d1, with z, without e, with semiz
figure_c=13;
output=CoreFHorzQHExpAssetz_nod1_z_noe_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, with semiz
figure_c=14;
output=CoreFHorzQHExpAssetz_d1_z_noe_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, with semiz
figure_c=15;
output=CoreFHorzQHExpAssetz_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, with semiz
figure_c=16;
output=CoreFHorzQHExpAssetz_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);


%% ================= with2A1 (figs 17-24): two standard endogenous assets =================
% Splice a binary second standard asset a1_2 in each subcode -> n_a=[a1_1,2,a2] -> length(n_a1)>1
% -> DC2A/GI2A/DC2A_GI2A path. base + DC2A + GI2A + DC2A_GI2A + lowmemory, Naive & Sophisticated.
n_a_2A1=[51,n_a_justexpasset]; % [a1_1, a2]; the binary a1_2 is added inside the subcode -> [51,2,13]
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3; a2_grid];
n_a_2A1_notsobig=[151,n_a_justexpasset];
a_grid_2A1_notsobig=[5*linspace(0,1,n_a_2A1_notsobig(1))'.^3; a2_grid];

%% with2A1, without d1, with z, without e
figure_c=17;
output=CoreFHorzQHExpAssetz_nod1_z_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, without e
figure_c=18;
output=CoreFHorzQHExpAssetz_d1_z_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, without d1, with z, with e
figure_c=19;
output=CoreFHorzQHExpAssetz_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e
figure_c=20;
output=CoreFHorzQHExpAssetz_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1 + semiz (figs 21-24)
% Reuses the 2A1 grids (n_a_2A1 etc.) from the nosemiz with2A1 section above, with the *semiz d-grids.
% Implemented and GPU-green as of 2026-08-19 (toolkit commit 2f291d33): the 24 semiz 2A raws
% (DC2A / GI2A / DC2A_GI2A x nod1/with-d1 x no-e/with-e x Naive/Sophisticated) plus the 2A branch
% of each of the three semiz sub-dispatchers.

%% with2A1, without d1, with z, without e, with semiz
figure_c=21;
output=CoreFHorzQHExpAssetz_nod1_z_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, without e, with semiz
figure_c=22;
output=CoreFHorzQHExpAssetz_d1_z_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, without d1, with z, with e, with semiz
figure_c=23;
output=CoreFHorzQHExpAssetz_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e, with semiz
figure_c=24;
output=CoreFHorzQHExpAssetz_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

diary off
