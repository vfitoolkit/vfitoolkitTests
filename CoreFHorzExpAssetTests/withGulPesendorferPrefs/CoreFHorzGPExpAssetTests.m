% GUL-PESENDORFER preference tests of the core VFI Toolkit FHorz ExpAsset commands.
% Mirrors the nosemiz 1A1 tiers of CoreFHorzExpAssetTests.m (and its QH twin
% CoreFHorzQHExpAssetTests.m, the donor bank): noa1 nosemiz figs 1-8, withA1 nosemiz
% figs 9-16, over with/without d1, with/without z, with/without e, and every solution
% method the baseline runs (withA1: base / DC1 / GI1 / DC1_GI1; noa1: base only).
%
% Gul-Pesendorfer temptation and self-control:
%   V_j = max_{d,a1'} [ u + v + beta*E V_{j+1} ] - max_{d,a1'} v
% so Policy maximizes the tempted objective and V nets off the self-control cost. The
% temptation fn v is declared as vfoptions.temptationFn (same input signature convention
% as the ReturnFn, own parameters; here deliberately different: temptation is over
% consumption only, no leisure term). Baseline: v = lambdaGP*u_c(c) + shiftGP with
% lambdaGP=0.1, shiftGP=0.
%
% Each subcode runs the methods with their lowmemory ladders, a ValueFnFromPolicy oracle
% on every method, and V_Jplus1 sections. The cross-test files carry the lambda=0 anchors
% (GP with lambda=0 == standard ExpAsset, at every tier), the hand-rolled brute force
% (the ExpAsset a2' expectation assembled by hand from aprimeFn), constant-v/shift
% invariance, age-dependent lambda_j, age-dependent pi_z_J (the KLM2021 configuration),
% V_Jplus1 at all tiers, and the V_GP<=V_standard sign check.
%
% STATE (2026-09-14, see GulPesendorfer_ExpAsset_proposal.md in the toolkit repo): GPU-green
% against toolkit 6de29e94 (the GulPesendorferExpAsset family, 40 raws). 508 checks: 482
% at/below the ULP floor, the other 24 at 5.6e-9/7.5e-9 = 3-4 ULP of the poor-corner
% V~1.1e7 in op-order-different comparisons (FromPolicy oracle, constant-v, shift
% invariance in the z figs) — accepted as the ULP floor scaled to this family's V. Run 1
% found the empty-feasible-set degeneracy (a2=0 kills earnings, so some states have no
% feasible choice; GP's V=max(-Inf)-max(-Inf)=NaN corrupted continuations), fixed by the
% MostTempting==-Inf -> 0 guard in every raw.
%
% Memory note: GP holds the return matrix AND its temptation twin simultaneously (~2x the
% core solver's memory), so every big/notsobig grid below is tightened relative to the QH
% twin's (which are themselves at the GPU's edge in places).
%
% The subcodes draw no figures, so only the diary is saved.
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzGPExpAssetTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzGPExpAssetTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzGPExpAssetTestsdiary.txt

addpath('./CoreFHorzGPExpAssetTests_subcodes/') % the GPTemptationFn_* wrappers
addpath('./CoreFHorzGPExpAssetTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzGPExpAssetTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzGPExpAssetTests_subcodes/CrossTests/')
addpath('../CoreFHorzExpAssetTests_Setup/')
addpath('../CoreFHorzExpAsset_ReturnFns/')
addpath('../CoreFHorzExpAsset_ReturnFns/Noa1_ReturnFns/')

%% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAsset_setup

%% Gul-Pesendorfer: the temptation parameters
% The temptation fns live in CoreFHorzGPExpAssetTests_subcodes/ (GPTemptationFn_*), one
% per case signature; each subcode declares its own vfoptions.temptationFn (mirroring how
% the ReturnFn is declared per-subcode)
Params.lambdaGP=0.1; % temptation strength: v = lambdaGP*u_c(c) + shiftGP
Params.shiftGP=0; % constant shift of the temptation utility; should never change anything

% vfoptions.exoticpreferences='GulPesendorfer' is set inside each subcode

%% ================= WITHOUT a1 (figs 1-8): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

%% without d1, without z, without e, noa1, nosemiz
figure_c=1;
output=CoreFHorzGPExpAsset_nod1_noz_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, without e, noa1, nosemiz
figure_c=2;
output=CoreFHorzGPExpAsset_d1_noz_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, without e, noa1, nosemiz
figure_c=3;
output=CoreFHorzGPExpAsset_nod1_z_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, noa1, nosemiz
figure_c=4;
output=CoreFHorzGPExpAsset_d1_z_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, without z, with e, noa1, nosemiz
figure_c=5;
output=CoreFHorzGPExpAsset_nod1_noz_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, noa1, nosemiz
figure_c=6;
output=CoreFHorzGPExpAsset_d1_noz_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, noa1, nosemiz
figure_c=7;
output=CoreFHorzGPExpAsset_nod1_z_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, noa1, nosemiz
figure_c=8;
output=CoreFHorzGPExpAsset_d1_z_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% ================= WITH a1 (figs 9-16) =================
% Reset the setup (the without-a1 half above left the workspace alone, but re-run for safety/independence)
CoreFHorzExpAsset_setup
Params.lambdaGP=0.1;
Params.shiftGP=0;

% GP big-grid rule: the temptation twin ~doubles solver memory, so each case's GI grid is
% tightened relative to its QH twin ([1001,13] -> [501,13]; [301,13] -> [201,13];
% [201,13] -> [151,13])
n_a_GPbig=[501,n_a_justexpasset];
a1_grid_GPbig=5*linspace(0,1,n_a_GPbig(1))'.^3;
a_grid_GPbig=[a1_grid_GPbig;a2_grid];

%% without d1, without z, without e, without semiz
figure_c=9;
output=CoreFHorzGPExpAsset_nod1_noz_noe_nosemiz(n_d_withoutd1,n_a,n_a_GPbig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_GPbig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, without e, without semiz
figure_c=10;
output=CoreFHorzGPExpAsset_d1_noz_noe_nosemiz(n_d_withd1,n_a,n_a_GPbig,n_z,N_j,d_grid_withd1,a_grid,a_grid_GPbig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, without e, without semiz
figure_c=11;
output=CoreFHorzGPExpAsset_nod1_z_noe_nosemiz(n_d_withoutd1,n_a,n_a_GPbig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_GPbig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, without semiz
n_a_notsobig=[201,n_a_justexpasset]; % QH twin uses [301,13]; GP rule tightens it
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=12;
output=CoreFHorzGPExpAsset_d1_z_noe_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, without z, with e, without semiz
figure_c=13;
output=CoreFHorzGPExpAsset_nod1_noz_e_nosemiz(n_d_withoutd1,n_a,n_a_GPbig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_GPbig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, without semiz
n_a_notsobig=[201,n_a_justexpasset]; % QH twin uses [301,13]; GP rule tightens it
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=14;
output=CoreFHorzGPExpAsset_d1_noz_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, without semiz
n_a_notsobig=[201,n_a_justexpasset]; % QH twin uses [301,13]; GP rule tightens it
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=15;
output=CoreFHorzGPExpAsset_nod1_z_e_nosemiz(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, without semiz
n_a_notsobig=[151,n_a_justexpasset]; % QH twin uses [201,13]; GP rule tightens it
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=16;
output=CoreFHorzGPExpAsset_d1_z_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% The cross tests (see the comments at the top of the cross-test subcodes for what they cover)
output=GPExpAsset_CrossTests_nod1(n_d_withoutd1,n_a,n_a_GPbig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_GPbig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=GPExpAsset_CrossTests_d1(n_d_withd1,n_a,n_a_GPbig,n_z,N_j,d_grid_withd1,a_grid,a_grid_GPbig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Done
diary off
