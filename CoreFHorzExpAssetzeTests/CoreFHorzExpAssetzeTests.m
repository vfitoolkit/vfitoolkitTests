% Implement tests of the core VFI Toolkit FHorz with ExpAssetze commands
% experienceassetze: aprime depends on (d2,a2,z,e), so z and e are always present
% with/without d1
% with/without semiz
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% NOTE ON NAMING: Every subcode carries an explicit tag saying how many standard
% endogenous states sit alongside the experienceassetze: '_noa1' (none: the
% experience asset a2 is the only endogenous state), '_withA1' (one), or
% '_with2A1' (two). This follows the naming convention of the neighbouring
% CoreFHorzExpAssetTests/CoreFHorzRiskyAssetTests.
% (An earlier version of this suite deliberately did not implement the noa1 tier,
% on the view that an experience asset is only meaningful alongside a standard
% endogenous asset; that decision has been reversed.)
% Layout: noa1 (figs 1-4) -> withA1 (figs 5-8) -> with2A1 (figs 9-12).

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt

addpath('./CoreFHorzExpAssetzeTests_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_Setup/')
addpath('./CoreFHorzExpAssetze_ReturnFns/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetze_setup

%% ================= NOA1 (figs 1-4): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.
addpath('./CoreFHorzExpAssetzeTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetze_ReturnFns/Noa1_ReturnFns/')
addpath('./CoreFHorzExpAssetze_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')

%% noa1, without d1, with z, with e
figure_c=1;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1, with d1, with z, with e
figure_c=2;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 nosemiz cross-tests
% CrossTest 5: noa1 vs withA1 model with n_a1=1 where a1 is ignored (should match bit-exact)
output=CoreFHorzExpAssetze_CrossTests5_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests5_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% noa1 semiz
% PENDING TOOLKIT SUPPORT: these error at the ValueFnIter call until the
% ExpAssetzeSemiExo noa1 raws exist (test-first: written ahead of the toolkit code)

%% noa1, without d1, with z, with e, with semiz
figure_c=3;
output=CoreFHorzExpAssetze_nod1_z_e_semiz_noa1(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1, with d1, with z, with e, with semiz
figure_c=4;
output=CoreFHorzExpAssetze_d1_z_e_semiz_noa1(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 semiz cross-tests
% CrossTest 5 + semiz: noa1 vs withA1 model with n_a1=1 where a1 is ignored (should match bit-exact)
output=CoreFHorzExpAssetze_CrossTests5_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests5_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);




%% ================= WITH a1 (figs 5-8) =================

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=5;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=6;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_withA1(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)


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

figure_c=7;
output=CoreFHorzExpAssetze_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, with semiz
n_a_notsobig=[151,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=8;
output=CoreFHorzExpAssetze_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

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




%% ================= With TWO standard endogenous assets (figs 9-12): a = [a1_1, a1_2 (binary), a2 (experienceassetze)] =================
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A1_subcodes/Semiz_subcodes/')

% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1: the
% first standard endogenous state is divide-conquered, the rest are folded/brute-forced).
% The second standard endogenous state a1_2 is BINARY (a capped high-return asset).
% a1main is kept modest here because the binary second asset doubles the a-grid.

% n_a_2A1=[a1, a1_2, a2] and a_grid_2A1 come from the setup; a1_2 is a genuine multi-point
% second standard asset, so the subcodes take n_a/a_grid as given and build nothing.
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

%% with2A1, without d1, with z, with e
figure_c=9;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, with e
figure_c=10;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1 nosemiz cross-tests
% CrossTest 6: a degenerate second standard asset a1_2 (single point {0}) reduces the
% two-standard-asset model back to the with-a1 model
output=CoreFHorzExpAssetze_CrossTests6_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests6_d1_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% with2A1 + semiz

%% with2A1, without d1, with z, with e, with semiz
figure_c=11;
output=CoreFHorzExpAssetze_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, with e, with semiz

% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
figure_c=12;
output=CoreFHorzExpAssetze_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
% I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS

diary off
