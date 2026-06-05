% Implement lots of tests of the core VFI Toolkit FHorz commands
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz


addpath('../CoreFHorzTests/CoreFHorzTests_Setup/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/')

addpath('./QuasiHyperbolicDiscountingFHorzTests_subcodes/')
addpath('./QuasiHyperbolicDiscountingFHorzTests_subcodes/CrossTests/')

%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};

% vfoptions.exoticpreferences='QuasiHyperbolic';
% vfoptions.quasi_hyperbolic='Naive';
% vfoptions.quasi_hyperbolic='Sophisticated';

%% The only functions worth testing are the value fn ones, as after you have Policy everything else is anyway ignoring the quasi-hyperbolic discounting

%% without d, without z, without e, without semiz
figure_c=1;
output=QHDFHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

% Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
% all just 1e-3, so irrelevant. Is because interpolation creates tiny amount of variance where the is none.
% (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

%% with d, without z, without e, without semiz
figure_c=2;
output=QHDFHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

% Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
% all just 1e-3, so irrelevant.  Is because interpolation creates tiny amount of variance where the is none.
% (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

%% without d, with z, without e, without semiz
figure_c=3;
output=QHDFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, with z, without e, without semiz
figure_c=4;
output=QHDFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d, without z, with e, without semiz
figure_c=5;
output=QHDFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, without z, with e, without semiz
figure_c=6;
output=QHDFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d, with z, with e, without semiz
figure_c=7;
output=QHDFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, with z, with e, without semiz
figure_c=8;
output=QHDFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
% Note: These are likely redundant, as already compared them all to exponential discounting and did these cross-tests for exponentinal
% discounting. But it is going to test the very unlikely case that they still somehow manage to differ when beta0 is involved

output=QHDFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=QHDFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% all looking good :)














%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

% d1 is a decision variable that is not in the SemiExoStateFn

addpath('../CoreFHorzTests/CoreFHorzTests_Setup/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/')

addpath('./QuasiHyperbolicDiscountingFHorzTests_subcodes/')
addpath('./QuasiHyperbolicDiscountingFHorzTests_subcodes/CrossTests/')

addpath('./QuasiHyperbolicDiscountingFHorzTests_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorz_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount={'beta0'};


%% without d1, without z, without e, with semiz
figure_c=9;
output=QHDFHorz_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, without e, with semiz
figure_c=10;
output=QHDFHorz_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, without e, with semiz
figure_c=11;
output=QHDFHorz_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good :)

%% with d1, with z, without e, with semiz
figure_c=12;
output=QHDFHorz_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, without z, with e, with semiz
figure_c=13;
output=QHDFHorz_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, with e, with semiz
figure_c=14;
output=QHDFHorz_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, with e, with semiz
figure_c=15;
output=QHDFHorz_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, with e, with semiz
figure_c=16;
output=QHDFHorz_d1_z_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
% Note: These are likely redundant, as already compared them all to exponential discounting and did these cross-tests for exponentinal
% discounting. But it is going to test the very unlikely case that they still somehow manage to differ when beta0 is involved

output=QHDFHorz_CrossTests_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=QHDFHorz_CrossTests_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% All looks good!

%% Now some further cross-tests, using a semi-exo that is really just a markov

output=QHDFHorz_CrossTests2_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=QHDFHorz_CrossTests2_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% All looks good!

