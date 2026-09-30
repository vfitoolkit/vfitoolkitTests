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

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1):  noa1 nosemiz (figs 1-8)
% doPart(2):  noa1 nosemiz cross-tests
% doPart(3):  noa1 semiz (figs 9-16)
% doPart(4):  noa1 semiz cross-tests
% doPart(5):  a1 nosemiz (figs 17-24)
% doPart(6):  a1 nosemiz cross-tests
% doPart(7):  a1 nosemiz fake-experience-asset cross-tests
% doPart(8):  a1 semiz (figs 25-32)
% doPart(9):  a1 semiz cross-tests
% doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov)
% doPart(11): a1 semiz fake-experience-asset cross-tests
% doPart(12): 2a1 nosemiz (figs 33-40)
% doPart(13): 2a1 nosemiz cross-tests
% doPart(14): 2a1 semiz (figs 41-48)
% doPart(15): 2A2 cross-tests
% doPart(16): 2A2 noa1 nosemiz (figs 49-56)
% doPart(17): 2A2 noa1 semiz (figs 57-64)
% doPart(18): 2A2 a1 nosemiz (figs 65-72)   - DIES HERE, see below
% doPart(19): 2A2 a1 semiz (figs 73-80)     - never reached
%
% THE RUN OF 2026-09-21 STOPPED IN PART 18, at the first fig-65 solve:
%   Arrays have incompatible sizes of 202 and 303 in dimension 1
%   ValueFnIter_FHorz_ExpAsset_DC1_nod1_noz_raw line 172
% That is the expected state of this tier, not a regression: the 2A2 with-a1 tiers are written
% test-first (see the WITH 2A2 banner below), and the toolkit port has only reached the noa1
% tiers - parts 16 and 17 ran clean. Because the bank dies there, it never reaches CoreSummary and
% the run produces no verdict line. To get one for everything that IS implemented, set parts 18
% and 19 to zero; parts 1-17 then run to the end and the summary reports on them.
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so Fig 42
% is the same test whatever doPart says, and a png from a previous run is never overwritten by a
% different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
%
% NOTE: figs 40, 44, 46 and 48 remain commented out INSIDE parts 12 and 14 - they out-of-memory on
% this machine. doPart is for whole tiers; an individual call that cannot run stays commented.
doPart=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
% The 2A2 tiers (parts 15-19) have been checked mechanically for this: every identifier they read
% resolves to the setup file, to code outside all the if-blocks, or to their own block -- so they
% can be run on their own, in any subset. Checked for 17-19 vs 1-16, for 18/19 vs 17, and for 19
% vs 1-18. Anything added later must keep that true.

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetTestsdiary.txt
fprintf('CoreFHorzExpAssetTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
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

%% ===== doPart(1): noa1 nosemiz (figs 1-8) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): noa1 nosemiz (figs 1-8) =====\n')
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
end % doPart(1): noa1 nosemiz (figs 1-8)

%% ===== doPart(2): noa1 nosemiz cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): noa1 nosemiz cross-tests =====\n')
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
end % doPart(2): noa1 nosemiz cross-tests


%% ===== doPart(3): noa1 semiz (figs 9-16) =====
if doPart(3)==1
    fprintf('\n===== doPart(3): noa1 semiz (figs 9-16) =====\n')
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
end % doPart(3): noa1 semiz (figs 9-16)

%% ===== doPart(4): noa1 semiz cross-tests =====
if doPart(4)==1
    fprintf('\n===== doPart(4): noa1 semiz cross-tests =====\n')
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
end % doPart(4): noa1 semiz cross-tests




%% ================= WITH a1 (figs 17-32) =================
% Reset the setup (the without-a1 half above left the workspace alone, but re-run for safety/independence)
CoreFHorzExpAsset_setup

%% ===== doPart(5): a1 nosemiz (figs 17-24) =====
if doPart(5)==1
    fprintf('\n===== doPart(5): a1 nosemiz (figs 17-24) =====\n')
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
end % doPart(5): a1 nosemiz (figs 17-24)

%% ===== doPart(6): a1 nosemiz cross-tests =====
if doPart(6)==1
    fprintf('\n===== doPart(6): a1 nosemiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorzExpAsset_CrossTests_nod1_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAsset_CrossTests_d1_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(6): a1 nosemiz cross-tests

%% ===== doPart(7): a1 nosemiz fake-experience-asset cross-tests =====
if doPart(7)==1
    fprintf('\n===== doPart(7): a1 nosemiz fake-experience-asset cross-tests =====\n')
    %% Do a test with a 'fake experience asset' and compare to a standard endogneous asset
    output=CoreFHorzExpAsset_CrossTests3_nod1_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAsset_CrossTests3_d1_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(7): a1 nosemiz fake-experience-asset cross-tests


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

%% ===== doPart(8): a1 semiz (figs 25-32) =====
if doPart(8)==1
    fprintf('\n===== doPart(8): a1 semiz (figs 25-32) =====\n')
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
end % doPart(8): a1 semiz (figs 25-32)

%% ===== doPart(9): a1 semiz cross-tests =====
if doPart(9)==1
    fprintf('\n===== doPart(9): a1 semiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorzExpAsset_CrossTests_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAsset_CrossTests_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(9): a1 semiz cross-tests

%% ===== doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov) =====
if doPart(10)==1
    fprintf('\n===== doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov) =====\n')
    %% Now some further cross-tests, using a semi-exo that is really just a markov
    output=CoreFHorzExpAsset_CrossTests2_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAsset_CrossTests2_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(10): a1 semiz cross-tests 2 (semi-exo that is really a markov)

%% ===== doPart(11): a1 semiz fake-experience-asset cross-tests =====
if doPart(11)==1
    fprintf('\n===== doPart(11): a1 semiz fake-experience-asset cross-tests =====\n')
    %% Do a test with a 'fake experience asset' and compare to a standard endogneous asset
    output=CoreFHorzExpAsset_CrossTests3_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzExpAsset_CrossTests3_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % all looking good :)
end % doPart(11): a1 semiz fake-experience-asset cross-tests


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

%% ===== doPart(12): 2a1 nosemiz (figs 33-40) =====
if doPart(12)==1
    fprintf('\n===== doPart(12): 2a1 nosemiz (figs 33-40) =====\n')
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
    % COMMENTED OUT: out-of-memories on this machine, and MATLAB halts the whole script at
    % the first error, which would stop the run before it ever reaches figs 49-80.
    % output=CoreFHorzExpAsset_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(12): 2a1 nosemiz (figs 33-40)

%% ===== doPart(13): 2a1 nosemiz cross-tests =====
if doPart(13)==1
    fprintf('\n===== doPart(13): 2a1 nosemiz cross-tests =====\n')
    %% 2a1 nosemiz cross-tests: a degenerate second asset a1_2 (single point {0}) reduces to the with-a1 model
    output=CoreFHorzExpAsset_CrossTests_nod1_nosemiz_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAsset_CrossTests_d1_nosemiz_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(13): 2a1 nosemiz cross-tests

%% ===== doPart(14): 2a1 semiz (figs 41-48) =====
if doPart(14)==1
    fprintf('\n===== doPart(14): 2a1 semiz (figs 41-48) =====\n')
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
    % COMMENTED OUT: out-of-memories on this machine, and MATLAB halts the whole script at
    % the first error, which would stop the run before it ever reaches figs 49-80.
    % output=CoreFHorzExpAsset_d1_z_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=45;
    output=CoreFHorzExpAsset_nod1_noz_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=46;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % COMMENTED OUT: out-of-memories on this machine, and MATLAB halts the whole script at
    % the first error, which would stop the run before it ever reaches figs 49-80.
    % output=CoreFHorzExpAsset_d1_noz_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=47;
    output=CoreFHorzExpAsset_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=48;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % COMMENTED OUT: out-of-memories on this machine, and MATLAB halts the whole script at
    % the first error, which would stop the run before it ever reaches figs 49-80.
    % output=CoreFHorzExpAsset_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(14): 2a1 semiz (figs 41-48)




%% ================= WITH 2A2 (figs 49-80): TWO experience assets =================
% vfoptions.experienceasset=2. n_a gains a second experience-asset dimension a2_2, so the grid
% layout is [a1, a2_1, a2_2] (or [a2_1, a2_2] in the noa1 tier), and the aprimeFn takes the extra
% integer 'whicha' selector between its a2 inputs and its parameters.
% Run order mirrors the rest of the bank: noa1 first (figs 49-64, where the two experience assets
% are the only endogenous states and DC/GI are irrelevant), then with a1 (figs 65-80, three
% endogenous states, with DC1/GI1/DC1+GI1 operating on a1).
%
% AS OF WRITING most of these error out: of the 384 raws in ValueFnIter/FHorz/ExperienceAsset only
% ValueFnIter_FHorz_ExpAsset_{,nod1_}noz_raw carry an l_a2==2 arm (plus 4 QuasiHyperbolic
% siblings). The bank is the specification; the toolkit port follows it.
% The cross-tests below are the machine-precision checks; run them FIRST, since CrossTests5/6/7
% (no z, no e) sit on the two raws that already work and so are the smoke test for the whole tier.
CoreFHorzExpAsset_setup
addpath('./CoreFHorzExpAssetTests_subcodes/With2A2_subcodes/')
addpath('./CoreFHorzExpAssetTests_subcodes/With2A2_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetTests_subcodes/With2A2_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetTests_subcodes/With2A2_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAsset_ReturnFns/Noa1_ReturnFns/') % only needed while figs 1-48 are commented out; that region adds it too
addpath('./CoreFHorzExpAsset_ReturnFns/With2A2_ReturnFns/')
addpath('./CoreFHorzExpAsset_ReturnFns/With2A2_ReturnFns/Semiz_ReturnFns/')
addpath('./CoreFHorzExpAsset_ReturnFns/With2A2_ReturnFns/Noa1_ReturnFns/')
addpath('./CoreFHorzExpAsset_ReturnFns/With2A2_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(15): 2A2 cross-tests =====
if doPart(15)==1
    fprintf('\n===== doPart(15): 2A2 cross-tests =====\n')
    %% 2A2 cross-tests (machine precision; these are the real checks of the tier)
    % An inert second experience asset must reduce to the one-experience-asset model...
    output=CoreFHorzExpAsset_CrossTests5_nod1_nosemiz_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAsset_CrossTests5_d1_nosemiz_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % ...and so must an inert FIRST experience asset (the two together pin down the dim ordering)
    output=CoreFHorzExpAsset_CrossTests6_nod1_nosemiz_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAsset_CrossTests6_d1_nosemiz_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % Two live, interchangeable experience assets: V/Policy/dist must be invariant to swapping them
    output=CoreFHorzExpAsset_CrossTests7_nod1_nosemiz_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAsset_CrossTests7_d1_nosemiz_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    % Same inert-second-asset reduction, but with no standard asset at all (N_a1==0 code path)
    output=CoreFHorzExpAsset_CrossTests5_nod1_noa1_nosemiz_with2A2(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAsset_CrossTests5_d1_noa1_nosemiz_with2A2(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(15): 2A2 cross-tests

%% ===== doPart(16): 2A2 noa1 nosemiz (figs 49-56) =====
if doPart(16)==1
    fprintf('\n===== doPart(16): 2A2 noa1 nosemiz (figs 49-56) =====\n')
    %% 2A2 noa1 nosemiz (figs 49-56): the two experience assets are the only endogenous states
    % n_a_2A2_justexpasset as n_a; the n_a_big slot is unused in the noa1 tier.
    figure_c=49;
    output=CoreFHorzExpAsset_nod1_noz_noe_noa1_nosemiz_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=50;
    output=CoreFHorzExpAsset_d1_noz_noe_noa1_nosemiz_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=51;
    output=CoreFHorzExpAsset_nod1_z_noe_noa1_nosemiz_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=52;
    output=CoreFHorzExpAsset_d1_z_noe_noa1_nosemiz_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=53;
    output=CoreFHorzExpAsset_nod1_noz_e_noa1_nosemiz_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=54;
    output=CoreFHorzExpAsset_d1_noz_e_noa1_nosemiz_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=55;
    output=CoreFHorzExpAsset_nod1_z_e_noa1_nosemiz_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=56;
    output=CoreFHorzExpAsset_d1_z_e_noa1_nosemiz_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(16): 2A2 noa1 nosemiz (figs 49-56)

%% ===== doPart(17): 2A2 noa1 semiz (figs 57-64) =====
if doPart(17)==1
    fprintf('\n===== doPart(17): 2A2 noa1 semiz (figs 57-64) =====\n')
    %% 2A2 noa1 semiz (figs 57-64)
    figure_c=57;
    output=CoreFHorzExpAsset_nod1_noz_noe_noa1_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=58;
    output=CoreFHorzExpAsset_d1_noz_noe_noa1_semiz_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=59;
    output=CoreFHorzExpAsset_nod1_z_noe_noa1_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=60;
    output=CoreFHorzExpAsset_d1_z_noe_noa1_semiz_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=61;
    output=CoreFHorzExpAsset_nod1_noz_e_noa1_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=62;
    output=CoreFHorzExpAsset_d1_noz_e_noa1_semiz_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=63;
    output=CoreFHorzExpAsset_nod1_z_e_noa1_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=64;
    output=CoreFHorzExpAsset_d1_z_e_noa1_semiz_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(17): 2A2 noa1 semiz (figs 57-64)

%% ===== doPart(18): 2A2 a1 nosemiz (figs 65-72) =====
if doPart(18)==1
    fprintf('\n===== doPart(18): 2A2 a1 nosemiz (figs 65-72) =====\n')
    %% 2A2 with a1 (figs 65-72, nosemiz): three endogenous states, DC1/GI1/DC1+GI1 operate on a1
    % The 'big' a1 grid used for the with/without-grid-interpolation moment comparison shrinks as
    % shocks are added: the grid-interpolation solve builds an a1 x a1prime(fine) x a2_1 x a2_2 x
    % shocks array. Cut these further if a given GPU runs out of memory.
    n_a_2A2_notsobig=[301,n_a2_1,n_a2_2];
    a1_grid_2A2_notsobig=5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;
    a_grid_2A2_notsobig=[a1_grid_2A2_notsobig;a2_1_grid;a2_2_grid];
    figure_c=65;
    output=CoreFHorzExpAsset_nod1_noz_noe_nosemiz_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=66;
    output=CoreFHorzExpAsset_d1_noz_noe_nosemiz_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    n_a_2A2_notsobig=[201,n_a2_1,n_a2_2];
    a1_grid_2A2_notsobig=5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;
    a_grid_2A2_notsobig=[a1_grid_2A2_notsobig;a2_1_grid;a2_2_grid];
    figure_c=67;
    output=CoreFHorzExpAsset_nod1_z_noe_nosemiz_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=68;
    output=CoreFHorzExpAsset_d1_z_noe_nosemiz_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=69;
    output=CoreFHorzExpAsset_nod1_noz_e_nosemiz_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=70;
    output=CoreFHorzExpAsset_d1_noz_e_nosemiz_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    n_a_2A2_notsobig=[101,n_a2_1,n_a2_2];
    a1_grid_2A2_notsobig=5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;
    a_grid_2A2_notsobig=[a1_grid_2A2_notsobig;a2_1_grid;a2_2_grid];
    figure_c=71;
    output=CoreFHorzExpAsset_nod1_z_e_nosemiz_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=72;
    output=CoreFHorzExpAsset_d1_z_e_nosemiz_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(18): 2A2 a1 nosemiz (figs 65-72)

%% ===== doPart(19): 2A2 a1 semiz (figs 73-80) =====
if doPart(19)==1
    fprintf('\n===== doPart(19): 2A2 a1 semiz (figs 73-80) =====\n')
    %% 2A2 with a1 (figs 73-80, semiz)
    n_a_2A2_notsobig=[201,n_a2_1,n_a2_2];
    a1_grid_2A2_notsobig=5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;
    a_grid_2A2_notsobig=[a1_grid_2A2_notsobig;a2_1_grid;a2_2_grid];
    figure_c=73;
    output=CoreFHorzExpAsset_nod1_noz_noe_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=74;
    output=CoreFHorzExpAsset_d1_noz_noe_semiz_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    n_a_2A2_notsobig=[151,n_a2_1,n_a2_2];
    a1_grid_2A2_notsobig=5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;
    a_grid_2A2_notsobig=[a1_grid_2A2_notsobig;a2_1_grid;a2_2_grid];
    figure_c=75;
    output=CoreFHorzExpAsset_nod1_z_noe_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=76;
    output=CoreFHorzExpAsset_d1_z_noe_semiz_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=77;
    output=CoreFHorzExpAsset_nod1_noz_e_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=78;
    output=CoreFHorzExpAsset_d1_noz_e_semiz_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    n_a_2A2_notsobig=[101,n_a2_1,n_a2_2];
    a1_grid_2A2_notsobig=5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;
    a_grid_2A2_notsobig=[a1_grid_2A2_notsobig;a2_1_grid;a2_2_grid];
    figure_c=79;
    output=CoreFHorzExpAsset_nod1_z_e_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=80;
    output=CoreFHorzExpAsset_d1_z_e_semiz_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(19): 2A2 a1 semiz (figs 73-80)


%% One verdict for the whole run
% The bank prints over a thousand checks; this reads the diary back and says plainly whether the
% run passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzExpAssetTestsdiary.txt')

diary off

%% THINGS NOT CHECKED
% Check using two decision variables in any of d1 or d3 (the decision variables that are not in experience asset)
