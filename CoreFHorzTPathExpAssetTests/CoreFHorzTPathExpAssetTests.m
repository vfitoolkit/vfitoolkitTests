% Implement core TPath tests of the VFI Toolkit FHorz commands with ExperienceAsset.
% with/without d1
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz [NOT YET IMPLEMENTED]
% noa1 [NOT YET IMPLEMENTED]


%%
addpath('./CoreFHorzTPathExpAssetTests_subcodes/')
addpath('./CoreFHorzTPathExpAssetTests_Setup/')
addpath('./CoreFHorzTPathExpAsset_ReturnFns/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPathExpAsset_setup

%% without d1, without z, without e, without semiz
figure_c=1;
output=CoreFHorzTPathExpAsset_nod1_noz_noe_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, without e, without semiz
figure_c=2;
output=CoreFHorzTPathExpAsset_d1_noz_noe_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, without e, without semiz
figure_c=3;
output=CoreFHorzTPathExpAsset_nod1_z_noe_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, without e, without semiz
figure_c=4;
output=CoreFHorzTPathExpAsset_d1_z_noe_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, without z, with e, without semiz
figure_c=5;
output=CoreFHorzTPathExpAsset_nod1_noz_e_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, without z, with e, without semiz
figure_c=6;
output=CoreFHorzTPathExpAsset_d1_noz_e_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d1, with z, with e, without semiz
figure_c=7;
output=CoreFHorzTPathExpAsset_nod1_z_e_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d1, with z, with e, without semiz
figure_c=8;
output=CoreFHorzTPathExpAsset_d1_z_e_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
