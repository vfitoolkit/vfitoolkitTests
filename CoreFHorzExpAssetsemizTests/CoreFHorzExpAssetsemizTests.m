% Implement tests of the core VFI Toolkit FHorz with ExpAssetsemiz commands
% experienceassetsemiz: aprime depends on (d2,a2,semiz), so semiz is always present
% (semiz drives the experience asset). semiz is the semi-exogenous state, so the
% semiz decision d3 is also always present.
% with/without d1
% with/without ordinary z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% This is all done with a1 (standard endogenous state alongside the experienceassetsemiz)
% experienceassetsemiz without a1 is not supported by VFI Toolkit

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt

addpath('./CoreFHorzExpAssetsemizTests_subcodes/')
addpath('./CoreFHorzExpAssetsemizTests_Setup/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/')
addpath('./CoreFHorzExpAssetsemizTests_subcodes/CrossTests/')


%% Setup so that use the same d,a,semiz,z,e in all the models that use them
CoreFHorzExpAssetsemiz_setup

%% without d1, without z, without e
figure_c=1;
output=CoreFHorzExpAssetsemiz_nod1_noz_noe(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, without e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=2;
output=CoreFHorzExpAssetsemiz_d1_noz_noe(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, without z, with e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=3;
output=CoreFHorzExpAssetsemiz_nod1_noz_e(n_d_withoutd1,n_a,n_a_notsobig,0,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e
n_a_notsobig=[201,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=4;
output=CoreFHorzExpAssetsemiz_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, without e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=5;
output=CoreFHorzExpAssetsemiz_nod1_z_noe(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=6;
output=CoreFHorzExpAssetsemiz_d1_z_noe(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=7;
output=CoreFHorzExpAssetsemiz_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e
n_a_notsobig=[201,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=8;
output=CoreFHorzExpAssetsemiz_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)


%% Cross-tests
% CrossTest 1: experienceassetsemiz with a degenerate (d3-invariant) semiz vs experienceassetz (should match)
output=CoreFHorzExpAssetsemiz_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassetsemiz that ignores semiz vs plain experienceasset+semiz (should match)
output=CoreFHorzExpAssetsemiz_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests2_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 3: experienceassetsemiz+z vs experienceassetz with combined z=[semiz,z] (pins bothz ordering)
output=CoreFHorzExpAssetsemiz_CrossTests3_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests3_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

diary off
