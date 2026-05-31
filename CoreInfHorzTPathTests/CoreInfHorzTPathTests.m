% Implement core TPath tests of the VFI Toolkit InfHorz commands.
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% No semiz (skipped intentionally for this build).


%%
addpath('./CoreInfHorzTPathTests_subcodes/')
addpath('./CoreInfHorzTPathTests_Setup/')
addpath('./CoreInfHorzTPath_ReturnFns/')

%% Setup so that use the same d,a,z,e in all the models that use them
CoreInfHorzTPath_setup

%% without d, without z, without e
figure_c=1;
output=CoreInfHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, without z, without e
figure_c=2;
output=CoreInfHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d, with z, without e
figure_c=3;
output=CoreInfHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, with z, without e
figure_c=4;
output=CoreInfHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d, without z, with e
figure_c=5;
output=CoreInfHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, without z, with e
figure_c=6;
output=CoreInfHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% without d, with z, with e
figure_c=7;
output=CoreInfHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);

%% with d, with z, with e
figure_c=8;
output=CoreInfHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
