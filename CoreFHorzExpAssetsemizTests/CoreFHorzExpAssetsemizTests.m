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
% Tier ordering follows CoreFHorzExpAssetTests and CoreFHorzExpAssetzTests:
%   noa1 (figs 1-8):    the experienceassetsemiz a2 is the only endogenous state
%   withA1 (figs 9-16): a1 (standard endogenous state) alongside the experienceassetsemiz
%   with2A1 (figs 17-24): two standard endogenous assets, a = [a1_1, a1_2 (binary), a2]

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


%% ================= noa1 (figs 1-8): the experience asset a2 is the ONLY endogenous state =================
% TEST-FIRST: these error when run, as the ExpAssetsemiz noa1 raws do not exist in the toolkit yet
% (ValueFnIter_FHorz_ExpAssetsemiz errors on N_a1==0).
addpath('./CoreFHorzExpAssetsemizTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/Noa1_ReturnFns/')

%% without d1, without z, without e, noa1
figure_c=1;
output=CoreFHorzExpAssetsemiz_nod1_noz_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, without e, noa1
figure_c=2;
output=CoreFHorzExpAssetsemiz_d1_noz_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, without z, with e, noa1
figure_c=3;
output=CoreFHorzExpAssetsemiz_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e, noa1
figure_c=4;
output=CoreFHorzExpAssetsemiz_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, without e, noa1
figure_c=5;
output=CoreFHorzExpAssetsemiz_nod1_z_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e, noa1
figure_c=6;
output=CoreFHorzExpAssetsemiz_d1_z_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1
figure_c=7;
output=CoreFHorzExpAssetsemiz_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1
figure_c=8;
output=CoreFHorzExpAssetsemiz_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 cross-tests (TEST-FIRST: error when run, same reason as above)
% CrossTest 4: noa1 vs withA1 model where a1 (n_a1=1) is ignored (should match bit-exact)
output=CoreFHorzExpAssetsemiz_CrossTests4_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests4_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests4_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests4_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% ================= withA1 (figs 9-16): a1 alongside the experienceassetsemiz =================

%% without d1, without z, without e
figure_c=9;
output=CoreFHorzExpAssetsemiz_nod1_noz_noe(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, without e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=10;
output=CoreFHorzExpAssetsemiz_d1_noz_noe(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, without z, with e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=11;
output=CoreFHorzExpAssetsemiz_nod1_noz_e(n_d_withoutd1,n_a,n_a_notsobig,0,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e
n_a_notsobig=[201,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=12;
output=CoreFHorzExpAssetsemiz_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, without e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=13;
output=CoreFHorzExpAssetsemiz_nod1_z_noe(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=14;
output=CoreFHorzExpAssetsemiz_d1_z_noe(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=15;
output=CoreFHorzExpAssetsemiz_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e
n_a_notsobig=[201,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=16;
output=CoreFHorzExpAssetsemiz_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)


%% withA1 cross-tests
% CrossTest 1: experienceassetsemiz with a degenerate (d3-invariant) semiz vs experienceassetz (should match)
output=CoreFHorzExpAssetsemiz_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassetsemiz that ignores semiz vs plain experienceasset+semiz (should match)
output=CoreFHorzExpAssetsemiz_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests2_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 3: experienceassetsemiz+z vs experienceassetz with combined z=[semiz,z] (pins bothz ordering)
output=CoreFHorzExpAssetsemiz_CrossTests3_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests3_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% ================= With TWO standard endogenous assets (figs 17-24) =================
% a = [a1_1, a1_2 (binary), a2 (experienceassetsemiz)]
% TEST-FIRST: these error when run, as the ExpAssetsemiz DC2A/GI2A/DC2A_GI2A raws do not
% exist in the toolkit yet.
addpath('./CoreFHorzExpAssetsemizTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/With2A1_ReturnFns/')

% a1main is kept modest here because the binary second asset doubles the a-grid.
% n_a_2A1=[a1, a1_2, a2] and a_grid_2A1 come from CoreFHorzExpAssetsemiz_setup; a1_2 is a genuine
% multi-point second standard asset, so the subcodes take n_a/a_grid as given and build nothing.
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

%% with2A1, without d1, without z, without e
figure_c=17;
output=CoreFHorzExpAssetsemiz_nod1_noz_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, without z, without e
figure_c=18;
output=CoreFHorzExpAssetsemiz_d1_noz_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, without d1, without z, with e
figure_c=19;
output=CoreFHorzExpAssetsemiz_nod1_noz_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, without z, with e
figure_c=20;
output=CoreFHorzExpAssetsemiz_d1_noz_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, without d1, with z, without e
figure_c=21;
output=CoreFHorzExpAssetsemiz_nod1_z_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, without e
figure_c=22;
output=CoreFHorzExpAssetsemiz_d1_z_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, without d1, with z, with e
figure_c=23;
output=CoreFHorzExpAssetsemiz_nod1_z_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, with e
figure_c=24;
output=CoreFHorzExpAssetsemiz_d1_z_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1 cross-tests (TEST-FIRST: error when run, same reason as above)
% CrossTest 5: a degenerate second standard asset a1_2 (single point {0}) reduces the
% two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
output=CoreFHorzExpAssetsemiz_CrossTests5_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetsemiz_CrossTests5_d1_with2A1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

diary off
