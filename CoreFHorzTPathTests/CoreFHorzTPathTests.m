% Implement lots of tests of the core VFI Toolkit FHorz commands
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz


%%
addpath('./CoreFHorzTPathTests_subcodes/')
addpath('./CoreFHorzTPathTests_Setup/')
addpath('./CoreFHorzTPath_ReturnFns/')
addpath('./CoreFHorzTPathTests_subcodes/CrossTests/')

%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPath_setup

%%
z_grid=z_grid.*ones(1,N_j,'gpuArray');
pi_z=pi_z.*ones(1,1,N_j,'gpuArray');
vfoptionsbaseline.pi_e=vfoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');
vfoptionsbaseline.e_grid=vfoptionsbaseline.e_grid.*ones(1,N_j,'gpuArray');
simoptionsbaseline.pi_e=simoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');
simoptionsbaseline.pi_e=simoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');

%% without d, without z, without e, without semiz
figure_c=1;
output=CoreFHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, without z, without e, without semiz
figure_c=2;
output=CoreFHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d, with z, without e, without semiz
figure_c=3;
output=CoreFHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, with z, without e, without semiz
figure_c=4;
output=CoreFHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d, without z, with e, without semiz
figure_c=5;
output=CoreFHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, without z, with e, without semiz
figure_c=6;
output=CoreFHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d, with z, with e, without semiz
figure_c=7;
output=CoreFHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d, with z, with e, without semiz
figure_c=8;
output=CoreFHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
output=CoreFHorzTPath_CrossTests_nod_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

output=CoreFHorzTPath_CrossTests_d_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

% all looking good :)


%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

% %% That is all the without semiz, now with semiz
% % From here on, it is the eight with semiz
% % From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)
% 
% % d1 is a decision variable that is not in the SemiExoStateFn
% 
% addpath('./CoreFHorzTPathTests_subcodes/')
% addpath('./CoreFHorzTPathTests_Setup/')
% addpath('./CoreFHorzTPath_ReturnFns/')
% addpath('./CoreFHorzTPathTests_subcodes/CrossTests/')
% 
% addpath('./CoreFHorzTPathTests_subcodes/Semiz_subcodes/')
% addpath('./CoreFHorzTPath_ReturnFns/Semiz_ReturnFns/')
% % Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
% CoreFHorzTPath_setup
% 
% % For models without d1, use:
% % n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% % For models with d1, use:
% % n_d_semiz and d_grid_semiz (as n_d and d_grid)
% 
% %% without d1, without z, without e, with semiz
% figure_c=9;
% output=CoreFHorzTPath_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% with d1, without z, without e, with semiz
% figure_c=10;
% output=CoreFHorzTPath_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% without d1, with z, without e, with semiz
% figure_c=11;
% output=CoreFHorzTPath_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good :)
% 
% %% with d1, with z, without e, with semiz
% figure_c=12;
% output=CoreFHorzTPath_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% without d1, without z, with e, with semiz
% figure_c=13;
% output=CoreFHorzTPath_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% with d1, without z, with e, with semiz
% figure_c=14;
% output=CoreFHorzTPath_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% without d1, with z, with e, with semiz
% figure_c=15;
% output=CoreFHorzTPath_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% with d1, with z, with e, with semiz
% figure_c=16;
% n_a_notsobig=751; % To avoid out-of-memory errors
% a_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
% 
% output=CoreFHorzTPath_d1_z_e_semiz(n_d_semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% % looks good
% 
% %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
% output=CoreFHorzTPath_CrossTests_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% 
% output=CoreFHorzTPath_CrossTests_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% 
% % All looks good!
% 
% %% Now some further cross-tests, using a semi-exo that is really just a markov
% 
% output=CoreFHorzTPath_CrossTests2_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% 
% output=CoreFHorzTPath_CrossTests2_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% 
% % All looks good!

%% Done! Damn that was a lot of tests. Glad that is over.




%% THINGS NOT CHECKED
% Check using two decision variables in the semiz codes (both for d1 and for d2, and without d1)


