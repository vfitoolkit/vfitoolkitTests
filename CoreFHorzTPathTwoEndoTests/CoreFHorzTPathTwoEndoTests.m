% Implement core TPath tests of the VFI Toolkit FHorz commands with two endogenous states.
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz [NOT YET IMPLEMENTED]


%%
addpath('./CoreFHorzTPathTwoEndoTests_subcodes/')
addpath('./CoreFHorzTPathTwoEndoTests_Setup/')
addpath('./CoreFHorzTPathTwoEndo_ReturnFns/')

%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPathTwoEndo_setup

%% without d, without z, without e, without semiz
figure_c=1;
output=CoreFHorzTPathTwoEndo_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, without z, without e, without semiz
figure_c=2;
output=CoreFHorzTPathTwoEndo_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d, with z, without e, without semiz
figure_c=3;
output=CoreFHorzTPathTwoEndo_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, with z, without e, without semiz
figure_c=4;
output=CoreFHorzTPathTwoEndo_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d, without z, with e, without semiz
figure_c=5;
output=CoreFHorzTPathTwoEndo_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, without z, with e, without semiz
figure_c=6;
output=CoreFHorzTPathTwoEndo_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d, with z, with e, without semiz
figure_c=7;
output=CoreFHorzTPathTwoEndo_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, with z, with e, without semiz
figure_c=8;
output=CoreFHorzTPathTwoEndo_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
