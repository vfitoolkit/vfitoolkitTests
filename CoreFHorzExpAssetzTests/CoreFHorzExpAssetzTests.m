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

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1):  noa1 figs 1-4
% doPart(2):  noa1 cross-tests
% doPart(3):  noa1 figs 5-8
% doPart(4):  noa1 cross-tests 2
% doPart(5):  withA1 figs 9-12
% doPart(6):  withA1 cross-tests
% doPart(7):  withA1 figs 13-16
% doPart(8):  withA1 cross-tests 2
% doPart(9):  withA1 figs 17-20
% doPart(10): withA1 cross-tests 3
% doPart(11): withA1 figs 21-24
% doPart(12): withA1 cross-tests 4
% doPart(13): with2A2 noa1 figs 25-28       (two EXPERIENCE assets, nosemiz)
% doPart(14): with2A2 noa1 figs 29-32       (two EXPERIENCE assets, semiz)
% doPart(15): with2A2 WITH a1 figs 33-36    (a1 + two experience assets, nosemiz)
% doPart(16): with2A2 WITH a1 figs 37-40    (a1 + two experience assets, semiz)
% doPart(17): with2A2 cross-tests           (the machine-precision checks of the tier)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
%
% NOTE: fig 24 stays commented out INSIDE part 11 - it out-of-memories on this machine. doPart is for
% whole tiers; an individual call that cannot run stays commented.
% This run: all 17 parts, as the confirming full-bank run before committing. Parts 13-17 are
% green (580 checks, 0 failed, 2026-10-09). Parts 1-12 were green on an earlier build, but two
% toolkit edits have landed since (the GI per-e accumulation guard, and relocating 160 l_a2==1
% setup lines inside their guard in the SemiExo raws). Both are no-ops at l_a2==1 by inspection,
% but 15 of the touched files are exercised by parts 3,4,7,8,11,12, so measure rather than argue.
doPart=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetzTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetzTestsdiary.txt
fprintf('CoreFHorzExpAssetzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
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

%% ===== doPart(1): noa1 figs 1-4 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): noa1 figs 1-4 =====\n')
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
end % doPart(1): noa1 figs 1-4

%% ===== doPart(2): noa1 cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): noa1 cross-tests =====\n')
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
end % doPart(2): noa1 cross-tests

%% noa1 semiz variants

%% ===== doPart(3): noa1 figs 5-8 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): noa1 figs 5-8 =====\n')
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
end % doPart(3): noa1 figs 5-8

%% ===== doPart(4): noa1 cross-tests 2 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): noa1 cross-tests 2 =====\n')
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
end % doPart(4): noa1 cross-tests 2


%% ================= withA1 (figs 9-16) =================

%% ===== doPart(5): withA1 figs 9-12 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): withA1 figs 9-12 =====\n')
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
end % doPart(5): withA1 figs 9-12


%% ===== doPart(6): withA1 cross-tests =====
if doPart(6)==1
    fprintf('\n===== doPart(6): withA1 cross-tests =====\n')
    %% Cross-tests
    % CrossTest 1: experienceassetz with iid-markov z vs experienceassete with iid e (should match)
    output=CoreFHorzExpAssetz_CrossTests_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % CrossTest 2: 'fake' experienceassetz that ignores z vs plain experienceasset (should match)
    output=CoreFHorzExpAssetz_CrossTests2_nod1_withA1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests2_d1_withA1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(6): withA1 cross-tests


%% Semiz variants
addpath('./CoreFHorzExpAssetzTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetz_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(7): withA1 figs 13-16 =====
if doPart(7)==1
    fprintf('\n===== doPart(7): withA1 figs 13-16 =====\n')
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
end % doPart(7): withA1 figs 13-16


%% ===== doPart(8): withA1 cross-tests 2 =====
if doPart(8)==1
    fprintf('\n===== doPart(8): withA1 cross-tests 2 =====\n')
    %% Semiz cross-tests
    % CrossTest1+semiz: experienceassetz+semiz iid-markov-z vs experienceassete+semiz iid-e (should match)
    output=CoreFHorzExpAssetz_CrossTests_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % CrossTest2+semiz: 'fake' experienceassetz+semiz that ignores z vs plain experienceasset+semiz (should match)
    output=CoreFHorzExpAssetz_CrossTests2_nod1_semiz_withA1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests2_d1_semiz_withA1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(8): withA1 cross-tests 2


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

%% ===== doPart(9): withA1 figs 17-20 =====
if doPart(9)==1
    fprintf('\n===== doPart(9): withA1 figs 17-20 =====\n')
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
end % doPart(9): withA1 figs 17-20

%% ===== doPart(10): withA1 cross-tests 3 =====
if doPart(10)==1
    fprintf('\n===== doPart(10): withA1 cross-tests 3 =====\n')
    %% with2A1 nosemiz cross-tests
    % CrossTest 4: a degenerate second standard asset a1_2 (single point {0}) reduces the
    % two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
    output=CoreFHorzExpAssetz_CrossTests4_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests4_d1_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(10): withA1 cross-tests 3

%% with2A1 + semiz (figs 21-24)
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
% (the semiz ReturnFns dir ./CoreFHorzExpAssetz_ReturnFns/Semiz_ReturnFns/ was already addpath-ed in the withA1 semiz section)
% Reuses the 2A1 grids (n_a_2A1 etc.) from the nosemiz with2A1 section above, with the *semiz d-grids.

%% ===== doPart(11): withA1 figs 21-24 =====
if doPart(11)==1
    fprintf('\n===== doPart(11): withA1 figs 21-24 =====\n')
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
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % output=CoreFHorzExpAssetz_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(11): withA1 figs 21-24

%% ===== doPart(12): withA1 cross-tests 4 =====
if doPart(12)==1
    fprintf('\n===== doPart(12): withA1 cross-tests 4 =====\n')
    %% with2A1 + semiz cross-tests
    % CrossTest 4 + semiz: a degenerate second standard asset a1_2 (single point {0}) reduces the
    % two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
    output=CoreFHorzExpAssetz_CrossTests4_nod1_semiz_with2A1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests4_d1_semiz_with2A1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(12): withA1 cross-tests 4


%% ===== with2A2: two EXPERIENCE assets =====
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A2_subcodes/')
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A2_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A2_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A2_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzExpAssetzTests_subcodes/With2A2_subcodes/WithA1_subcodes/Semiz_subcodes/')

% n_a_2A2=[a1, a2_1, a2_2] and a_grid_2A2 come from the setup. The grid below is the FINE one: it
% feeds the subcodes' "use a really big a_grid, then the moments should be essentially the same
% with/without grid interpolation" block, so it must have MORE a1 points than n_a_2A2 (=51), not
% the same number -- an equal-sized grid there makes that comparison vacuous. 151 mirrors
% n_a_2A1_notsobig; the semiz variants get 101 because they also carry n_semiz and d3.
n_a_2A2_notsobig=[151,n_a2_1,n_a2_2];
a_grid_2A2_notsobig=[5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;a2_1_grid;a2_2_grid];
n_a_2A2_notsobig_semiz=[101,n_a2_1,n_a2_2];
a_grid_2A2_notsobig_semiz=[5*linspace(0,1,n_a_2A2_notsobig_semiz(1))'.^3;a2_1_grid;a2_2_grid];

%% ===== doPart(13): with2A2 noa1 figs 25-28 =====
if doPart(13)==1
    fprintf('\n===== doPart(13): with2A2 noa1 figs 25-28 =====\n')
    figure_c=25;
    output=CoreFHorzExpAssetz_nod1_z_noe_nosemiz_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=26;
    output=CoreFHorzExpAssetz_d1_z_noe_nosemiz_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=27;
    output=CoreFHorzExpAssetz_nod1_z_e_nosemiz_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=28;
    output=CoreFHorzExpAssetz_d1_z_e_nosemiz_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(13): with2A2 noa1 figs 25-28

%% ===== doPart(14): with2A2 noa1 figs 29-32 =====
if doPart(14)==1
    fprintf('\n===== doPart(14): with2A2 noa1 figs 29-32 =====\n')
    figure_c=29;
    output=CoreFHorzExpAssetz_nod1_z_noe_semiz_noa1_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=30;
    output=CoreFHorzExpAssetz_d1_z_noe_semiz_noa1_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=31;
    output=CoreFHorzExpAssetz_nod1_z_e_semiz_noa1_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=32;
    output=CoreFHorzExpAssetz_d1_z_e_semiz_noa1_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(14): with2A2 noa1 figs 29-32

%% ===== doPart(15): with2A2 WITH a1 figs 33-36 =====
if doPart(15)==1
    fprintf('\n===== doPart(15): with2A2 WITH a1 figs 33-36 =====\n')
    figure_c=33;
    output=CoreFHorzExpAssetz_nod1_z_noe_nosemiz_withA1_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=34;
    output=CoreFHorzExpAssetz_d1_z_noe_nosemiz_withA1_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=35;
    output=CoreFHorzExpAssetz_nod1_z_e_nosemiz_withA1_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=36;
    output=CoreFHorzExpAssetz_d1_z_e_nosemiz_withA1_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(15): with2A2 WITH a1 figs 33-36

%% ===== doPart(16): with2A2 WITH a1 figs 37-40 =====
if doPart(16)==1
    fprintf('\n===== doPart(16): with2A2 WITH a1 figs 37-40 =====\n')
    figure_c=37;
    output=CoreFHorzExpAssetz_nod1_z_noe_semiz_withA1_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig_semiz,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig_semiz,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=38;
    output=CoreFHorzExpAssetz_d1_z_noe_semiz_withA1_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig_semiz,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig_semiz,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=39;
    output=CoreFHorzExpAssetz_nod1_z_e_semiz_withA1_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig_semiz,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig_semiz,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    figure_c=40;
    output=CoreFHorzExpAssetz_d1_z_e_semiz_withA1_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig_semiz,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig_semiz,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(16): with2A2 WITH a1 figs 37-40

%% ===== doPart(17): with2A2 cross-tests =====
if doPart(17)==1
    fprintf('\n===== doPart(17): with2A2 cross-tests =====\n')
    % These build their own 2A2 grids internally, so they take the ordinary 1A n_a/a_grid.
    %% CrossTest 5: an INERT SECOND experience asset reduces to the one-experience-asset model
    output=CoreFHorzExpAssetz_CrossTests5_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests5_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    %% CrossTest 6: an INERT FIRST experience asset does too (together with 5 this pins the ordering)
    output=CoreFHorzExpAssetz_CrossTests6_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests6_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    %% CrossTest 7: SWAP SYMMETRY, which keeps both assets alive so all four corners carry weight
    output=CoreFHorzExpAssetz_CrossTests7_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests7_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    %% CrossTest 5, noa1 tier
    output=CoreFHorzExpAssetz_CrossTests5_nod1_noa1_with2A2(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetz_CrossTests5_d1_noa1_with2A2(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(17): with2A2 cross-tests

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzExpAssetzTestsdiary.txt')

diary off
