% Implement lots of tests of the core VFI Toolkit FHorz with ExpAsset commands
% with/without d1
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz
%
% This is all done WITHOUT a1 first (experience asset is the only endogenous state;
% divide-and-conquer and grid interpolation layer are not relevant there), then rerun
% all 16 WITH a1 (where divide-and-conquer and grid interpolation layer also apply).

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetTestsdiary.txt

addpath('./CoreFHorzExpAssetTests_subcodes/')
addpath('./CoreFHorzExpAssetTests_Setup/')
addpath('./CoreFHorzExpAsset_ReturnFns/')
addpath('./CoreFHorzExpAssetTests_subcodes/CrossTests/')


%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzExpAsset_setup


%% ================= WITHOUT a1 (figs 1-16): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

addpath('./CoreFHorzExpAssetTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAsset_ReturnFns/Noa1_ReturnFns/')
addpath('./CoreFHorzExpAsset_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')

%% noa1 nosemiz (8 variants)

%% without d1, without z, without e, noa1, nosemiz
figure_c=1;
output=CoreFHorzExpAsset_nod1_noz_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, without e, noa1, nosemiz
figure_c=2;
output=CoreFHorzExpAsset_d1_noz_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, without e, noa1, nosemiz
figure_c=3;
output=CoreFHorzExpAsset_nod1_z_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e, noa1, nosemiz
figure_c=4;
output=CoreFHorzExpAsset_d1_z_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, without z, with e, noa1, nosemiz
figure_c=5;
output=CoreFHorzExpAsset_nod1_noz_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e, noa1, nosemiz
figure_c=6;
output=CoreFHorzExpAsset_d1_noz_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1, nosemiz
figure_c=7;
output=CoreFHorzExpAsset_nod1_z_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1, nosemiz
figure_c=8;
output=CoreFHorzExpAsset_d1_z_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 nosemiz cross-tests
% Markov-as-iid equivalence cross-test
output=CoreFHorzExpAsset_CrossTests_nod1_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests_d1_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% ExpAsset noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo model with d2
output=CoreFHorzExpAsset_CrossTests3_nod1_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests3_d1_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% no a1 vs model with a1 but where it is ignored
output=CoreFHorzExpAsset_CrossTests4_nod1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests4_d1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% noa1 semiz (8 variants)

%% without d1, without z, without e, noa1, semiz
figure_c=9;
output=CoreFHorzExpAsset_nod1_noz_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, without e, noa1, semiz
figure_c=10;
output=CoreFHorzExpAsset_d1_noz_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, without e, noa1, semiz
figure_c=11;
output=CoreFHorzExpAsset_nod1_z_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e, noa1, semiz
figure_c=12;
output=CoreFHorzExpAsset_d1_z_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, without z, with e, noa1, semiz
figure_c=13;
output=CoreFHorzExpAsset_nod1_noz_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e, noa1, semiz
figure_c=14;
output=CoreFHorzExpAsset_d1_noz_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1, semiz
figure_c=15;
output=CoreFHorzExpAsset_nod1_z_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1, semiz
figure_c=16;
output=CoreFHorzExpAsset_d1_z_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 semiz cross-tests
% Markov-as-iid equivalence cross-test
output=CoreFHorzExpAsset_CrossTests_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% semiz state as a plain Markov should give same answer
output=CoreFHorzExpAsset_CrossTests2_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests2_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% ExpAsset noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo model with d2
output=CoreFHorzExpAsset_CrossTests3_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests3_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% no a1 vs model with a1 but where it is ignored
output=CoreFHorzExpAsset_CrossTests4_nod1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests4_d1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);




%% ================= WITH a1 (figs 17-32) =================
% Reset the setup (the without-a1 half above left the workspace alone, but re-run for safety/independence)
CoreFHorzExpAsset_setup

%% a1 nosemiz (8 variants)

%% without d1, without z, without e, without semiz
figure_c=17;
output=CoreFHorzExpAsset_nod1_noz_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d1, without z, without e, without semiz
figure_c=18;
output=CoreFHorzExpAsset_d1_noz_noe_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% without d1, with z, without e, without semiz
figure_c=19;
output=CoreFHorzExpAsset_nod1_z_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d1, with z, without e, without semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=20;
output=CoreFHorzExpAsset_d1_z_noe_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% without d1, without z, with e, without semiz
figure_c=21;
output=CoreFHorzExpAsset_nod1_noz_e_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d1, without z, with e, without semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=22;
output=CoreFHorzExpAsset_d1_noz_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% without d1, with z, with e, without semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=23;
output=CoreFHorzExpAsset_nod1_z_e_nosemiz(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d1, with z, with e, without semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=24;
output=CoreFHorzExpAsset_d1_z_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
output=CoreFHorzExpAsset_CrossTests_nod1_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=CoreFHorzExpAsset_CrossTests_d1_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% all looking good :)

%% Do a test with a 'fake experience asset' and compare to a standard endogneous asset
output=CoreFHorzExpAsset_CrossTests3_nod1_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=CoreFHorzExpAsset_CrossTests3_d1_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% all looking good :)


%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

addpath('./CoreFHorzExpAssetTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAsset_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorzExpAsset_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% a1 semiz (8 variants)

%% without d1, without z, without e, with semiz
figure_c=25;
output=CoreFHorzExpAsset_nod1_noz_noe_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% looks good

%% with d1, without z, without e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=26;
output=CoreFHorzExpAsset_d1_noz_noe_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good (as good as it can be expected to given the n_a_notsobig)

%% without d1, with z, without e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=27;
output=CoreFHorzExpAsset_nod1_z_noe_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good (as good as it can be expected to given the n_a_notsobig)

%% with d1, with z, without e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=28;
output=CoreFHorzExpAsset_d1_z_noe_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good (as good as it can be expected to given the n_a_notsobig)

%% without d1, without z, with e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=29;
output=CoreFHorzExpAsset_nod1_noz_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good (as good as it can be expected to given the n_a_notsobig)

%% with d1, without z, with e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=30;
output=CoreFHorzExpAsset_d1_noz_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good (as good as it can be expected to given the n_a_notsobig)

%% without d1, with z, with e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=31;
output=CoreFHorzExpAsset_nod1_z_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good

%% with d1, with z, with e, with semiz
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=32;
output=CoreFHorzExpAsset_d1_z_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% looks good (as good as it can be expected to given the n_a_notsobig)

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
output=CoreFHorzExpAsset_CrossTests_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=CoreFHorzExpAsset_CrossTests_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% all looking good :)

%% Now some further cross-tests, using a semi-exo that is really just a markov
output=CoreFHorzExpAsset_CrossTests2_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=CoreFHorzExpAsset_CrossTests2_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% all looking good :)

%% Do a test with a 'fake experience asset' and compare to a standard endogneous asset
output=CoreFHorzExpAsset_CrossTests3_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=CoreFHorzExpAsset_CrossTests3_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% all looking good :)


%% ================= WITH 2a1 (figs 33-48): two standard endogenous assets (triggers DC2A/GI2A/DC2A_GI2A) =================
% A genuine multi-point second standard asset a1_2 (return r2) is spliced between the liquid asset
% a1 and the experience asset a2: n_a_2A1=[n_a1, n_a1_2, n_a2] (built in CoreFHorzExpAsset_setup).
% Two standard assets -> length(n_a1)>1 -> DC2A / GI2A / DC2A_GI2A.
% The 8 semiz variants (figs 41-48) exercise the ExpAssetSemiExo DC2A/GI2A/DC2A_GI2A family
% (written test-first; the toolkit raws now exist).
CoreFHorzExpAsset_setup
addpath('./CoreFHorzExpAssetTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAsset_ReturnFns/With2A1_ReturnFns/')
addpath('./CoreFHorzExpAsset_ReturnFns/With2A1_ReturnFns/Semiz_ReturnFns/')

%% 2a1 nosemiz (8 variants)
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=33;
output=CoreFHorzExpAsset_nod1_noz_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=34;
output=CoreFHorzExpAsset_d1_noz_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=35;
output=CoreFHorzExpAsset_nod1_z_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=36;
output=CoreFHorzExpAsset_d1_z_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=37;
output=CoreFHorzExpAsset_nod1_noz_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=38;
output=CoreFHorzExpAsset_d1_noz_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=39;
output=CoreFHorzExpAsset_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=40;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzExpAsset_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS

%% 2a1 nosemiz cross-tests: a degenerate second asset a1_2 (single point {0}) reduces to the with-a1 model
output=CoreFHorzExpAsset_CrossTests_nod1_nosemiz_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAsset_CrossTests_d1_nosemiz_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% 2a1 semiz (8 variants)
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=41;
output=CoreFHorzExpAsset_nod1_noz_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=42;
output=CoreFHorzExpAsset_d1_noz_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=43;
output=CoreFHorzExpAsset_nod1_z_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

figure_c=44;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzExpAsset_d1_z_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
figure_c=45;
output=CoreFHorzExpAsset_nod1_noz_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=46;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzExpAsset_d1_noz_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
figure_c=47;
output=CoreFHorzExpAsset_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
figure_c=48;
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
output=CoreFHorzExpAsset_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS


diary off

%% THINGS NOT CHECKED
% Check using two decision variables in any of d1 or d3 (the decision variables that are not in experience asset)
