% Tests of the core VFI Toolkit FHorz ExpAssete commands under QUASI-HYPERBOLIC discounting.
% Runs the CoreFHorzExpAsseteTests models (Naive + Sophisticated), plus exponential cross-tests:
%   (i)   Naive: the continuation value equals the exponential value function (any beta0)
%   (ii)  Naive, beta0=1: both the main value fn AND the continuation value equal exponential
%   (iii) Sophisticated, beta0=1: both the main value fn AND the continuation value equal exponential
%
% experienceassete: aprime depends on (d2,a2,e), so e is ALWAYS present -> every subcode is _e.
% Same 24 combinations, in the same figure order, as the baseline CoreFHorzExpAsseteTests:
%   figs 1-4:   noa1 (the experience asset a2 is the only endogenous state), base method only
%   figs 5-8:   noa1 + semiz, base method only
%   figs 9-12:  withA1 (single standard asset a1), base / DC1 / GI1 / DC1_GI1
%   figs 13-16: withA1 + semiz, base / DC1 / GI1 / DC1_GI1
%   figs 17-20: with2A1 (binary second standard asset), base / DC2A / GI2A / DC2A_GI2A
%   figs 21-24: with2A1 + semiz, base / DC2A / GI2A / DC2A_GI2A
% Each subcode runs its methods at every valid lowmemory level (set by which shocks are
% present: {e}->{0,1}, {z,e} or {semiz,e}->{0,1,2}, {semiz,z,e}->{0,1,2,3}).
%
% STATUS (2026-08-19): the toolkit NOW HAS quasi-hyperbolic support for experienceassete. The
% 128 QH raws live in ValueFnIter/FHorz/ExperienceAssete/QuasiHyperbolic/, and
% ValueFnIter_Case1_FHorz routes to ValueFnIter_FHorz_QuasiHyperbolicExpAssete (nosemiz) and
% ValueFnIter_FHorz_QuasiHyperbolicExpAsseteSemiExo (semiz). Committed in toolkit b6bc5acf.
%
% This bank was written test-first, but the solvers landed before it was last run: the
% 2026-08-18 run is FULLY GREEN -- 1624 checks, every one zero, no errors. It covers Naive and
% Sophisticated on V / Valt / Policy / Policyalt across each variant's lowmemory ladder, plus 24
% ValueFnFromPolicy oracle checks (including the Valt reconstruction for Naive).
%
% One known toolkit limit, which does NOT affect any figure here:
% ValueFnIter_FHorz_QuasiHyperbolicExpAsseteSemiExo_{DC,GI,DC_GI} still error on N_a1==0
% ('Have not implemented experience assets with semi-exogenous shocks, without also having a
% standard asset'). The noa1+semiz block (figs 5-8) is base-method only, so it never reaches
% those dispatchers.
%
% This suite lives in the parent CoreFHorzExpAsseteTests/ folder and reuses its setup and ReturnFns.
%
% The subcodes here take figure_c but draw no figures, so only the diary is saved.

%% Diary of the command window output (written to the parent bank's TestOutput folder)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzQHExpAsseteTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzQHExpAsseteTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzQHExpAsseteTestsdiary.txt

addpath('../CoreFHorzExpAsseteTests_Setup/')
addpath('../CoreFHorzExpAssete_ReturnFns/')
addpath('../CoreFHorzExpAssete_ReturnFns/Noa1_ReturnFns/')
addpath('../CoreFHorzExpAssete_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzExpAssete_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzExpAssete_ReturnFns/With2A1_ReturnFns/')
addpath('../CoreFHorzExpAssete_ReturnFns/With2A1_ReturnFns/Semiz_ReturnFns/')
% The baseline bank also adds the experienceassetz ReturnFns (for its cross-tests); kept here
% so this suite sees exactly the same path as the baseline bank.
addpath('../../CoreFHorzExpAssetzTests/CoreFHorzExpAssetz_ReturnFns/')
addpath('./CoreFHorzQHExpAsseteTests_subcodes/')
addpath('./CoreFHorzQHExpAsseteTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzQHExpAsseteTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzQHExpAsseteTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzQHExpAsseteTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzQHExpAsseteTests_subcodes/With2A1_subcodes/Semiz_subcodes/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssete_setup

Params.beta0=0.9; % additional today-tomorrow (present-bias) discount factor
vfoptionsbaseline.QHadditionaldiscount='beta0';

%% ================= WITHOUT a1 (figs 1-8): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1), so these run the base method at every valid
% lowmemory level, for Naive then Sophisticated, plus the exponential cross-tests.
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

%% noa1 nosemiz (4 variants)

%% without d1, without z, with e, noa1
figure_c=1;
output=CoreFHorzQHExpAssete_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, noa1
figure_c=2;
output=CoreFHorzQHExpAssete_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, noa1
figure_c=3;
output=CoreFHorzQHExpAssete_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, noa1
figure_c=4;
output=CoreFHorzQHExpAssete_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% noa1 semiz (4 variants)

%% without d1, without z, with e, noa1, semiz
figure_c=5;
output=CoreFHorzQHExpAssete_nod1_noz_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, noa1, semiz
figure_c=6;
output=CoreFHorzQHExpAssete_d1_noz_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, noa1, semiz
figure_c=7;
output=CoreFHorzQHExpAssete_nod1_z_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, noa1, semiz
figure_c=8;
output=CoreFHorzQHExpAssete_d1_z_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);


%% ================= WITH a1 (figs 9-16) =================
% base / DC1 / GI1 / DC1_GI1, Naive and Sophisticated, at every valid lowmemory level.

%% without d1, without z, with e
figure_c=9;
output=CoreFHorzQHExpAssete_nod1_noz_e(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=10;
output=CoreFHorzQHExpAssete_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=11;
output=CoreFHorzQHExpAssete_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=12;
output=CoreFHorzQHExpAssete_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% Semiz variants (figs 13-16)

%% without d1, without z, with e, with semiz
figure_c=13;
output=CoreFHorzQHExpAssete_nod1_noz_e_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=14;
output=CoreFHorzQHExpAssete_d1_noz_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, with semiz
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=15;
output=CoreFHorzQHExpAssete_nod1_z_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, with semiz
n_a_notsobig=[151,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=16;
output=CoreFHorzQHExpAssete_d1_z_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);


%% ================= With TWO standard endogenous assets: a = [a1_1, a1_2 (binary), a2 (experienceassete)] (figs 17-20) =================
% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1: the
% first standard endogenous state is divide-conquered, the rest are folded/brute-forced).
% The second standard endogenous state a1_2 is BINARY (a capped high-return asset).
% a1main is kept modest here because the binary second asset doubles the a-grid.

% n_a_2A1=[a1, a1_2, a2] and a_grid_2A1 come from the setup; a1_2 is a genuine multi-point
% second standard asset, so the subcodes take n_a/a_grid as given and build nothing.
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

%% with2A1, without d1, without z, with e
figure_c=17;
output=CoreFHorzQHExpAssete_nod1_noz_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, without z, with e
figure_c=18;
output=CoreFHorzQHExpAssete_d1_noz_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, without d1, with z, with e
figure_c=19;
output=CoreFHorzQHExpAssete_nod1_z_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e
figure_c=20;
output=CoreFHorzQHExpAssete_d1_z_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% ================= with2A1 + semiz tier (figs 21-24) =================
% Same 2A1 grids as figs 17-20 (n_a_2A1 etc. defined above), but now with the semiz
% d-grids (d3 search effort drives the semi-exogenous employment state).

%% with2A1+semiz, without d1, without z, with e
figure_c=21;
output=CoreFHorzQHExpAssete_nod1_noz_e_with2A1_semiz(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1+semiz, with d1, without z, with e
figure_c=22;
output=CoreFHorzQHExpAssete_d1_noz_e_with2A1_semiz(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1+semiz, without d1, with z, with e
figure_c=23;
output=CoreFHorzQHExpAssete_nod1_z_e_with2A1_semiz(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1+semiz, with d1, with z, with e
figure_c=24;
output=CoreFHorzQHExpAssete_d1_z_e_with2A1_semiz(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

diary off
