% Implement tests of the core VFI Toolkit FHorz with ExpAssetze commands
% experienceassetze: aprime depends on (d2,a2,z,e), so z and e are always present
% with/without d1
% with/without semiz
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% NOTE ON NAMING: Every subcode here is run WITH a1 (a standard endogenous state
% alongside the experienceassetze), and so carries the '_withA1' suffix, following the
% naming convention of the neighbouring CoreFHorzRiskyAssetTests.
% Unlike RiskyAsset, there is NO 'without a1' counterpart: experienceassetze without a1
% has DELIBERATELY NOT been implemented in VFI Toolkit (an experience asset is only
% meaningful alongside a standard endogenous asset). Hence there is nothing to
% distinguish against, but the '_withA1' suffix is kept to make this explicit.


addpath('./CoreFHorzExpAssetzeTests_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/WithA1_subcodes/') % all subcodes are 'with a1' (see note above)
addpath('./CoreFHorzExpAssetzeTests_Setup/')
addpath('./CoreFHorzExpAssetze_ReturnFns/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetze_setup

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=1;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=2;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_withA1(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);


%% Cross-tests
% CrossTest 1: 'fake' experienceassetze that ignores e vs actual experienceassetz (should match)
output=CoreFHorzExpAssetze_CrossTests_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassetze that ignores z vs actual experienceassete (should match)
output=CoreFHorzExpAssetze_CrossTests2_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests2_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 3: 'fake' experienceassetze that ignores both z and e vs plain experienceasset (should match)
output=CoreFHorzExpAssetze_CrossTests3_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests3_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 4: experienceassetze with iid-markov z + e vs experienceassete with 2-dim e (should match)
output=CoreFHorzExpAssetze_CrossTests4_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests4_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);





%% Semiz variants
addpath('./CoreFHorzExpAssetzeTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetze_ReturnFns/Semiz_ReturnFns/')

%% without d1, with z, with e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=3;
output=CoreFHorzExpAssetze_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, with semiz
n_a_notsobig=[151,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=4;
output=CoreFHorzExpAssetze_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% Semiz cross-tests

% CrossTest1+semiz: 'fake' experienceassetze+semiz that ignores e vs experienceassetz+semiz
output=CoreFHorzExpAssetze_CrossTests_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest2+semiz: 'fake' experienceassetze+semiz that ignores z vs experienceassete+semiz
output=CoreFHorzExpAssetze_CrossTests2_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests2_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest3+semiz: 'fake' experienceassetze+semiz that ignores both z and e vs plain experienceasset+semiz
output=CoreFHorzExpAssetze_CrossTests3_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests3_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest4+semiz: experienceassetze+semiz with iid-markov z + e vs experienceassete+semiz with 2-dim e
output=CoreFHorzExpAssetze_CrossTests4_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests4_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);




%% With TWO standard endogenous assets: a = [a1_1, a1_2 (binary), a2 (experienceassetze)]
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A1_subcodes/Semiz_subcodes/')

% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1: the
% first standard endogenous state is divide-conquered, the rest are folded/brute-forced).
% The second standard endogenous state a1_2 is BINARY (a capped high-return asset).
% a1main is kept modest here because the binary second asset doubles the a-grid.

n_a_2A1=[51,13]; % [a1_1, a2]; the binary a1_2 is added inside the subcode -> [51,2,13]
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3; a2_grid];
n_a_2A1_notsobig=[151,13];
a_grid_2A1_notsobig=[5*linspace(0,1,n_a_2A1_notsobig(1))'.^3; a2_grid];

%% with2A1, without d1, with z, with e
figure_c=5;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1, with d1, with z, with e
figure_c=6;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with2A1 + semiz --- COMMENTED OUT AS NOT YET BUILD

% figure_c=7;
% output=CoreFHorzExpAssetze_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% figure_c=8;
% output=CoreFHorzExpAssetze_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

