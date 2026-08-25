% Implement tests of the core VFI Toolkit FHorz with ExpAssetz commands
% experienceassetz: aprime depends on (d2,a2,z), so z is always present
% with/without d1
% with/without e
% with/without semiz
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% NOTE ON NAMING: Subcodes carry an explicit a1-tier suffix. The '_noa1' subcodes
% (figs 1-8) run the experience asset a2 as the ONLY endogenous state (no
% divide-and-conquer/grid-interpolation blocks -- those act on a1). The '_withA1'
% subcodes (figs 9-16) run WITH a standard endogenous state a1 alongside the
% experienceassetz (where DC/GI/DC+GI also apply). Following CoreFHorzExpAssetzeTests,
% the '_with2A1' subcodes (figs 17-24) use TWO standard endogenous assets (a binary
% second asset a1_2 spliced in), which triggers the DC2A / GI2A / DC2A_GI2A code
% paths (used whenever length(n_a1)>1: the first standard endogenous state is
% divide-conquered and the rest are folded). The '_with2A1' tier is figs 17-24:
% nosemiz (figs 17-20) and semiz (figs 21-24), each covering {nod1,d1}x{noe,e}.
% TEST-FIRST: the noa1+semiz figs (5-8) and the noa1+semiz cross-tests error when
% run, as the ExpAssetzSemiExo noa1 raws do not exist in the toolkit yet. Similarly
% the with2A1+semiz figs (21-24) error when run, as the ExpAssetzSemiExo
% DC2A/GI2A/DC2A_GI2A raws are being built in parallel and do not exist yet.

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetzTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetzTestsdiary.txt

addpath('./CoreFHorzExpAssetzTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzExpAssetzTests_Setup/')
addpath('./CoreFHorzExpAssetz_ReturnFns/')
addpath('./CoreFHorzExpAssetzTests_subcodes/CrossTests/')
% Cross-tests compare against experienceassete and experienceasset, so need their ReturnFns
addpath('../CoreFHorzExpAsseteTests/CoreFHorzExpAssete_ReturnFns/')
addpath('../CoreFHorzExpAsseteTests/CoreFHorzExpAssete_ReturnFns/Noa1_ReturnFns/')


% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssetz_setup

%% ================= noa1 (figs 1-8): experience asset a2 is the ONLY endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

addpath('./CoreFHorzExpAssetzTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetzTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetz_ReturnFns/Noa1_ReturnFns/')
addpath('./CoreFHorzExpAssetz_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')

%% without d1, with z, without e, noa1, nosemiz
figure_c=1;
output=CoreFHorzExpAssetz_nod1_z_noe_nosemiz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e, noa1, nosemiz
figure_c=2;
output=CoreFHorzExpAssetz_d1_z_noe_nosemiz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1, nosemiz
figure_c=3;
output=CoreFHorzExpAssetz_nod1_z_e_nosemiz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1, nosemiz
figure_c=4;
output=CoreFHorzExpAssetz_d1_z_e_nosemiz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 nosemiz cross-tests
% CrossTest 3: noa1 vs withA1 model where a1 (n_a1=1) is ignored (should match bit-exact)
output=CoreFHorzExpAssetz_CrossTests3_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests3_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% and the same, but for the shock-leanest variant (z, without e)
output=CoreFHorzExpAssetz_CrossTests3_nod1_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests3_d1_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 1: noa1 experienceassetz with iid-markov z vs noa1 experienceassete with iid e (should match)
output=CoreFHorzExpAssetz_CrossTests_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: noa1 'fake' experienceassetz that ignores z vs noa1 plain experienceasset (should match)
output=CoreFHorzExpAssetz_CrossTests2_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests2_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% noa1 semiz variants

%% without d1, with z, without e, noa1, semiz
figure_c=5;
output=CoreFHorzExpAssetz_nod1_z_noe_semiz_noa1(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e, noa1, semiz
figure_c=6;
output=CoreFHorzExpAssetz_d1_z_noe_semiz_noa1(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1, semiz
figure_c=7;
output=CoreFHorzExpAssetz_nod1_z_e_semiz_noa1(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1, semiz
figure_c=8;
output=CoreFHorzExpAssetz_d1_z_e_semiz_noa1(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 semiz cross-tests
% CrossTest 3 + semiz: noa1 vs withA1 model where a1 (n_a1=1) is ignored (should match bit-exact)
output=CoreFHorzExpAssetz_CrossTests3_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests3_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% and the same, but for the shock-leanest variant (z, without e)
output=CoreFHorzExpAssetz_CrossTests3_nod1_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests3_d1_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 1 + semiz: noa1 experienceassetz+semiz iid-markov-z vs noa1 experienceassete+semiz iid-e (should match)
output=CoreFHorzExpAssetz_CrossTests_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2 + semiz: noa1 'fake' experienceassetz+semiz that ignores z vs noa1 plain experienceasset+semiz (should match)
output=CoreFHorzExpAssetz_CrossTests2_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests2_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% ================= withA1 (figs 9-16) =================

%% without d1, with z, without e
figure_c=9;
output=CoreFHorzExpAssetz_nod1_z_noe_nosemiz_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=10;
output=CoreFHorzExpAssetz_d1_z_noe_nosemiz_withA1(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=11;
output=CoreFHorzExpAssetz_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=12;
output=CoreFHorzExpAssetz_d1_z_e_nosemiz_withA1(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)


%% Cross-tests
% CrossTest 1: experienceassetz with iid-markov z vs experienceassete with iid e (should match)
output=CoreFHorzExpAssetz_CrossTests_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassetz that ignores z vs plain experienceasset (should match)
output=CoreFHorzExpAssetz_CrossTests2_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests2_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% Semiz variants
addpath('./CoreFHorzExpAssetzTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetz_ReturnFns/Semiz_ReturnFns/')

%% without d1, with z, without e, with semiz
figure_c=13;
output=CoreFHorzExpAssetz_nod1_z_noe_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, without e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=14;
output=CoreFHorzExpAssetz_d1_z_noe_semiz_withA1(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, with semiz
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=15;
output=CoreFHorzExpAssetz_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, with semiz
n_a_notsobig=[201,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=16;
output=CoreFHorzExpAssetz_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)


%% Semiz cross-tests
% CrossTest1+semiz: experienceassetz+semiz iid-markov-z vs experienceassete+semiz iid-e (should match)
output=CoreFHorzExpAssetz_CrossTests_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest2+semiz: 'fake' experienceassetz+semiz that ignores z vs plain experienceasset+semiz (should match)
output=CoreFHorzExpAssetz_CrossTests2_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests2_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% With TWO standard endogenous assets: a = [a1_1, a1_2 (binary), a2 (experienceassetz)] (figs 17-24)
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A1_subcodes/')

% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1: the
% first standard endogenous state is divide-conquered, the rest are folded/brute-forced).
% The second standard endogenous state a1_2 is BINARY (a capped high-return asset).
% a1main is kept modest here because the binary second asset doubles the a-grid.

% n_a_2A1=[a1, a1_2, a2] and a_grid_2A1 come from the setup; a1_2 is a genuine multi-point
% second standard asset, so the subcodes take n_a/a_grid as given and build nothing.
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

%% with2A1, without d1, with z, without e
figure_c=17;
output=CoreFHorzExpAssetz_nod1_z_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, without e
figure_c=18;
output=CoreFHorzExpAssetz_d1_z_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, without d1, with z, with e
figure_c=19;
output=CoreFHorzExpAssetz_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, with e
figure_c=20;
output=CoreFHorzExpAssetz_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1 nosemiz cross-tests
% CrossTest 4: a degenerate second standard asset a1_2 (single point {0}) reduces the
% two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
output=CoreFHorzExpAssetz_CrossTests4_nod1_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests4_d1_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% with2A1 + semiz (figs 21-24)
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
% (the semiz ReturnFns dir ./CoreFHorzExpAssetz_ReturnFns/Semiz_ReturnFns/ was already addpath-ed in the withA1 semiz section)
% Reuses the 2A1 grids (n_a_2A1 etc.) from the nosemiz with2A1 section above, with the *semiz d-grids.

%% with2A1, without d1, with z, without e, with semiz
figure_c=21;
output=CoreFHorzExpAssetz_nod1_z_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, without e, with semiz
figure_c=22;
output=CoreFHorzExpAssetz_d1_z_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, without d1, with z, with e, with semiz
figure_c=23;
output=CoreFHorzExpAssetz_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, with e, with semiz
figure_c=24;
output=CoreFHorzExpAssetz_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1 + semiz cross-tests
% CrossTest 4 + semiz: a degenerate second standard asset a1_2 (single point {0}) reduces the
% two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
output=CoreFHorzExpAssetz_CrossTests4_nod1_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetz_CrossTests4_d1_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

diary off
