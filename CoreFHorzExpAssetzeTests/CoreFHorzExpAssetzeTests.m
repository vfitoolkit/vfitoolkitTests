% Implement tests of the core VFI Toolkit FHorz with ExpAssetze commands
% experienceassetze: aprime depends on (d2,a2,z,e), so z and e are always present
% with/without d1
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% This is all done with a1 (standard endogenous state alongside the experienceassetze)
% Semiz tests are skipped (experienceassetze does not yet support semiz)


addpath('./CoreFHorzExpAssetzeTests_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_Setup/')
addpath('./CoreFHorzExpAssetze_ReturnFns/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/CrossTests/')


%% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetze_setup

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=1;
output=CoreFHorzExpAssetze_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=2;
output=CoreFHorzExpAssetze_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% Cross-tests
% CrossTest 1: 'fake' experienceassetze that ignores e vs actual experienceassetz (should match)
output=CoreFHorzExpAssetze_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassetze that ignores z vs actual experienceassete (should match)
output=CoreFHorzExpAssetze_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests2_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 3: 'fake' experienceassetze that ignores both z and e vs plain experienceasset (should match)
output=CoreFHorzExpAssetze_CrossTests3_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests3_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 4: experienceassetze with iid-markov z + e vs experienceassete with 2-dim e (should match)
output=CoreFHorzExpAssetze_CrossTests4_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests4_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% THINGS NOT CHECKED
% Tests with semiz (experienceassetze does not yet support semiz)
% Tests with the experienceassetze as the only asset (no a1)
