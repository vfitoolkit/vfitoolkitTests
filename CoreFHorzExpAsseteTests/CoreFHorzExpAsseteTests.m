% Implement tests of the core VFI Toolkit FHorz with ExpAssete commands
% experienceassete: aprime depends on (d2,a2,e), so e is always present
% with/without d1
% with/without z
% with/without semiz
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% This is all done WITHOUT a1 first (figs 1-8: the experience asset a2 is the only
% endogenous state; divide-and-conquer and grid interpolation layer are not relevant
% there), then WITH a1 (figs 9-16: standard endogenous state alongside the
% experienceassete, where divide-and-conquer and grid interpolation layer also apply),
% then WITH TWO standard endogenous assets (figs 17-24: with2A1, a binary second
% standard asset a1_2 folded alongside the divide-conquered a1_1; triggers the
% DC2A/GI2A/DC2A_GI2A code paths; figs 17-20 nosemiz, figs 21-24 with semiz).

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAsseteTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAsseteTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAsseteTestsdiary.txt

addpath('./CoreFHorzExpAsseteTests_subcodes/')
addpath('./CoreFHorzExpAsseteTests_Setup/')
addpath('./CoreFHorzExpAssete_ReturnFns/')
addpath('./CoreFHorzExpAsseteTests_subcodes/CrossTests/')
% Cross-tests compare against experienceassetz and experienceasset, so need their ReturnFns
addpath('../CoreFHorzExpAssetzTests/CoreFHorzExpAssetz_ReturnFns/')
addpath('../CoreFHorzExpAssetzTests/CoreFHorzExpAssetz_ReturnFns/Noa1_ReturnFns/')


%% Setup so that use the same d,a,z,e in all the models that use them
CoreFHorzExpAssete_setup


%% ================= WITHOUT a1 (figs 1-8): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

addpath('./CoreFHorzExpAsseteTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssete_ReturnFns/Noa1_ReturnFns/')

%% noa1 nosemiz (4 variants)

%% without d1, without z, with e, noa1
figure_c=1;
output=CoreFHorzExpAssete_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e, noa1
figure_c=2;
output=CoreFHorzExpAssete_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1
figure_c=3;
output=CoreFHorzExpAssete_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1
figure_c=4;
output=CoreFHorzExpAssete_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 nosemiz cross-tests
% CrossTest3: noa1 vs model with a1 but where it is ignored (a1=1 degenerate)
output=CoreFHorzExpAssete_CrossTests3_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests3_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 3 at the leaner shock case (e only, no ordinary z)
output=CoreFHorzExpAssete_CrossTests3_nod1_noz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests3_d1_noz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 1 at noa1: experienceassete iid-e vs experienceassetz iid-markov-z
output=CoreFHorzExpAssete_CrossTests_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2 at noa1: fake-e-ignored experienceassete vs plain experienceasset
output=CoreFHorzExpAssete_CrossTests2_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests2_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% noa1 semiz (4 variants)
% PENDING TOOLKIT SUPPORT: ExpAssete+SemiExo+noa1 raws do not exist yet (being built
% test-first), so these currently error at the ValueFnIter call.
addpath('./CoreFHorzExpAsseteTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssete_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')

%% without d1, without z, with e, noa1, semiz
figure_c=5;
output=CoreFHorzExpAssete_nod1_noz_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e, noa1, semiz
figure_c=6;
output=CoreFHorzExpAssete_d1_noz_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, noa1, semiz
figure_c=7;
output=CoreFHorzExpAssete_nod1_z_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, noa1, semiz
figure_c=8;
output=CoreFHorzExpAssete_d1_z_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% noa1 semiz cross-tests
% CrossTest3+semiz: noa1 vs model with a1 but where it is ignored (a1=1 degenerate)
output=CoreFHorzExpAssete_CrossTests3_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests3_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 3 at the leaner shock case (e only, no ordinary z)
output=CoreFHorzExpAssete_CrossTests3_nod1_noz_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests3_d1_noz_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 1 at noa1: experienceassete iid-e vs experienceassetz iid-markov-z
output=CoreFHorzExpAssete_CrossTests_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2 at noa1: fake-e-ignored experienceassete vs plain experienceasset
output=CoreFHorzExpAssete_CrossTests2_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests2_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% ================= WITH a1 (figs 9-16) =================

%% without d1, without z, with e
figure_c=9;
output=CoreFHorzExpAssete_nod1_noz_e(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=10;
output=CoreFHorzExpAssete_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=11;
output=CoreFHorzExpAssete_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e
n_a_notsobig=[201,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=12;
output=CoreFHorzExpAssete_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)


%% Cross-tests
% CrossTest 1: experienceassete with iid e vs experienceassetz with iid-markov z (should match)
output=CoreFHorzExpAssete_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest 2: 'fake' experienceassete that ignores e vs plain experienceasset (should match)
output=CoreFHorzExpAssete_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests2_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);















%% Semiz variants
addpath('./CoreFHorzExpAsseteTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssete_ReturnFns/Semiz_ReturnFns/')


%% without d1, without z, with e, with semiz
figure_c=13;
output=CoreFHorzExpAssete_nod1_noz_e_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, without z, with e, with semiz
n_a_notsobig=[301,13]; % To avoid out-of-memory errors
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=14;
output=CoreFHorzExpAssete_d1_noz_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d1, with z, with e, with semiz
n_a_notsobig=[301,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=15;
output=CoreFHorzExpAssete_nod1_z_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d1, with z, with e, with semiz
n_a_notsobig=[151,13];
a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
a_grid_notsobig=[a1_grid_notsobig;a2_grid];

figure_c=16;
output=CoreFHorzExpAssete_d1_z_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% Semiz cross-tests

% CrossTest1+semiz: experienceassete+semiz iid-e vs experienceassetz+semiz iid-markov-z
output=CoreFHorzExpAssete_CrossTests_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% CrossTest2+semiz: 'fake' experienceassete+semiz that ignores e vs plain experienceasset+semiz
output=CoreFHorzExpAssete_CrossTests2_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests2_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);


%% ================= With TWO standard endogenous assets: a = [a1_1, a1_2 (binary), a2 (experienceassete)] (figs 17-20) =================
% PENDING TOOLKIT SUPPORT: the ExpAssete 2A raws are being built in parallel
% (test-first), so figs 17-20 error at the ValueFnIter call until they land.
addpath('./CoreFHorzExpAsseteTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssete_ReturnFns/With2A1_ReturnFns/')

% Triggers the DC2A / GI2A / DC2A_GI2A code paths (used whenever length(n_a1)>1: the
% first standard endogenous state is divide-conquered, the rest are folded/brute-forced).
% The second standard endogenous state a1_2 is BINARY (a capped high-return asset).
% a1main is kept modest here because the binary second asset doubles the a-grid.

n_a_2A1=[51,n_a_justexpasset]; % [a1_1, a2]; the binary a1_2 is added inside the subcode -> [51,2,13]
a_grid_2A1=[5*linspace(0,1,n_a_2A1(1))'.^3; a2_grid];
n_a_2A1_notsobig=[151,n_a_justexpasset];
a_grid_2A1_notsobig=[5*linspace(0,1,n_a_2A1_notsobig(1))'.^3; a2_grid];

%% with2A1, without d1, without z, with e
figure_c=17;
output=CoreFHorzExpAssete_nod1_noz_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, without z, with e
figure_c=18;
output=CoreFHorzExpAssete_d1_noz_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, without d1, with z, with e
figure_c=19;
output=CoreFHorzExpAssete_nod1_z_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1, with d1, with z, with e
figure_c=20;
output=CoreFHorzExpAssete_d1_z_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1 cross-tests
% CrossTest 4: with2A1 with degenerate a1_2 (single grid point {0}) vs plain withA1 (should match)
output=CoreFHorzExpAssete_CrossTests4_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssete_CrossTests4_d1_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% ================= with2A1 + semiz tier (figs 21-24) =================
% PENDING TOOLKIT SUPPORT: the ExpAsseteSemiExo 2A raws are being built in parallel
% (test-first), so figs 21-24 error at the ValueFnIter call until they land.
addpath('./CoreFHorzExpAsseteTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssete_ReturnFns/With2A1_ReturnFns/Semiz_ReturnFns/')

% Same 2A1 grids as figs 17-20 (n_a_2A1 etc. defined above), but now with the semiz
% d-grids (d3 search effort drives the semi-exogenous employment state).

%% with2A1+semiz, without d1, without z, with e
figure_c=21;
output=CoreFHorzExpAssete_nod1_noz_e_with2A1_semiz(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1+semiz, with d1, without z, with e
figure_c=22;
output=CoreFHorzExpAssete_d1_noz_e_with2A1_semiz(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1+semiz, without d1, with z, with e
figure_c=23;
output=CoreFHorzExpAssete_nod1_z_e_with2A1_semiz(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with2A1+semiz, with d1, with z, with e
figure_c=24;
output=CoreFHorzExpAssete_d1_z_e_with2A1_semiz(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAsseteTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% No semiz-2A1 cross-tests (matches the ExpAsset template's semiz-2A1 tier).

diary off
