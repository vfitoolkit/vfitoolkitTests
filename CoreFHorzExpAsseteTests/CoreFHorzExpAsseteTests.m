% Implement tests of the core VFI Toolkit FHorz with ExpAssete commands
% experienceassete: aprime depends on (d2,a2,e), so e is always present
% with/without d1
% with/without z
% with/without divide-and-conquer
% with/without grid interpolation (NOT YET IMPLEMENTED in toolkit -- tests will throw)
% with/without low memory (where appropriate)
%
% This is all done with a1 (standard endogenous state alongside the experienceassete)
% Semiz tests are skipped (experienceassete does not yet support semiz)


addpath('./CoreFHorzExpAsseteTests_subcodes/')
addpath('./CoreFHorzExpAsseteTests_Setup/')
addpath('./CoreFHorzExpAssete_ReturnFns/')
addpath('./CoreFHorzExpAsseteTests_subcodes/CrossTests/')
% Cross-tests compare against experienceassetz and experienceasset, so need their ReturnFns
addpath('../CoreFHorzExpAssetzTests/CoreFHorzExpAssetz_ReturnFns/')


%% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssete_setup

%% without d1, without z, with e
figure_c=1;
output=CoreFHorzExpAssete_nod1_noz_e(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=2;
output=CoreFHorzExpAssete_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=3;
output=CoreFHorzExpAssete_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=4;
output=CoreFHorzExpAssete_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% Cross-tests
% CrossTest 1: experienceassete with iid e vs experienceassetz with iid-markov z (should match)
output=CoreFHorzExpAssete_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassete that ignores e vs plain experienceasset (should match)
output=CoreFHorzExpAssete_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests2_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% THINGS NOT CHECKED
% Tests with semiz (experienceassete does not yet support semiz)
% Tests with the experienceassete as the only asset (no a1)
