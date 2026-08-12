% Implement core TPath tests of the VFI Toolkit FHorz commands with ExperienceAssetz.
% experienceassetz: aprime depends on (d2,a2,z), so z is always present
% with/without d1
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz [NOT YET IMPLEMENTED]
% noa1 [NOT YET IMPLEMENTED]
%
% This is all done with a1 (standard endogenous state alongside the experienceassetz)
% experienceassetz without a1 is not supported by VFI Toolkit


%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathExpAssetzTestsdiary.txt

%%
addpath('./CoreFHorzTPathExpAssetzTests_subcodes/')
addpath('./CoreFHorzTPathExpAssetzTests_Setup/')
addpath('./CoreFHorzTPathExpAssetz_ReturnFns/')
addpath('./CoreFHorzTPathExpAssetzTests_subcodes/CrossTests/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPathExpAssetz_setup

%% without d1, with z, without e, without semiz
figure_c=1;
output=CoreFHorzTPathExpAssetz_nod1_z_noe_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% look good

%% with d1, with z, without e, without semiz
figure_c=2;
% Following the TPath ExpAsset tests, use a smaller 'big' grid to avoid out-of-memory
n_a_notsobig=[501,n_a_justexpasset]; % to test Grid Interpolation
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];
output=CoreFHorzTPathExpAssetz_d1_z_noe_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% Some of the lowmemory are not quite right, seems to be just the Policy (V
% seems fine); at first I thought L2flag, but appears to impact the DC
% (without GI) so that is not the reason.

%% without d1, with z, with e, without semiz
figure_c=3;
% experienceassetz machinery carries z-dependence, so use the smaller 'big' grid to avoid out-of-memory
n_a_notsobig=[301,n_a_justexpasset]; % to test Grid Interpolation
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];
output=CoreFHorzTPathExpAssetz_nod1_z_e_nosemiz(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% look good

%% with d1, with z, with e, without semiz
figure_c=4;
n_a_notsobig=[301,n_a_justexpasset]; % to test Grid Interpolation
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];
output=CoreFHorzTPathExpAssetz_d1_z_e_nosemiz(T,PricePath,ParamPath,n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% Some of the lowmemory are not quite right, seems to be just the Policy (V
% seems fine); at first I thought L2flag, but appears to impact the DC
% (without GI) so that is not the reason.

%% Cross-tests
% CrossTest 2: 'fake' experienceassetz that ignores z vs plain experienceasset (should match)
output=CoreFHorzTPathExpAssetz_CrossTest2_nod1(T,PricePath,ParamPath,n_d_withoutd1,n_a,n_z,N_j,d_grid_withoutd1,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
output=CoreFHorzTPathExpAssetz_CrossTest2_d1(T,PricePath,ParamPath,n_d_withd1,n_a,n_z,N_j,d_grid_withd1,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);

diary off
