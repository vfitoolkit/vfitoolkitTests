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

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1):  NOA1 figs 1-2
% doPart(2):  NOA1 cross-tests
% doPart(3):  NOA1 figs 3-4
% doPart(4):  NOA1 cross-tests 2
% doPart(5):  WITH a1 figs 5-6
% doPart(6):  WITH a1 cross-tests
% doPart(7):  WITH a1 figs 7-8
% doPart(8):  WITH a1 cross-tests 2
% doPart(9):  With TWO standard endogenous figs 9-10
% doPart(10): With TWO standard endogenous cross-tests
% doPart(11): With TWO standard endogenous figs 11-12
% doPart(12): with2A2 NOA1 figs 13-14        (two EXPERIENCE assets, nosemiz)
% doPart(13): with2A2 NOA1 figs 15-16        (two EXPERIENCE assets, semiz)
% doPart(14): with2A2 WITH a1 figs 17-18     (a1 + two experience assets, nosemiz)
% doPart(15): with2A2 WITH a1 figs 19-20     (a1 + two experience assets, semiz)
% doPart(16): with2A2 cross-tests            (the machine-precision checks of the tier)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
% This run: parts 12-16 only, the with2A2 tier. Parts 1-10 were all green on the 2026-10-02 run
% (which is the l_a2==1 regression guard for the rewritten peel / aprimeFn arity / corner dispatch),
% so there is no reason to pay for them again. Part 11 stays off: fig 12 OOMs on this GPU in
% CreateReturnFnMatrix_ExpAsset_Disc_e. Parts are independent (no part reads another's output and
% nothing is passed through a .mat), so running this subset is safe.
doPart=[0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt
fprintf('CoreFHorzExpAssetzeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
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

%% ===== doPart(1): NOA1 figs 1-2 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): NOA1 figs 1-2 =====\n')
    %% noa1, without d1, with z, with e
    figure_c=1;
    output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% noa1, with d1, with z, with e
    figure_c=2;
    output=CoreFHorzExpAssetze_d1_z_e_nosemiz_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): NOA1 figs 1-2

%% ===== doPart(2): NOA1 cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): NOA1 cross-tests =====\n')
    %% noa1 nosemiz cross-tests
    % CrossTest 5: noa1 vs withA1 model with n_a1=1 where a1 is ignored (should match bit-exact)
    output=CoreFHorzExpAssetze_CrossTests5_nod1_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetze_CrossTests5_d1_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): NOA1 cross-tests

%% noa1 semiz
% PENDING TOOLKIT SUPPORT: these error at the ValueFnIter call until the
% ExpAssetzeSemiExo noa1 raws exist (test-first: written ahead of the toolkit code)

%% ===== doPart(3): NOA1 figs 3-4 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): NOA1 figs 3-4 =====\n')
    %% noa1, without d1, with z, with e, with semiz
    figure_c=3;
    output=CoreFHorzExpAssetze_nod1_z_e_semiz_noa1(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% noa1, with d1, with z, with e, with semiz
    figure_c=4;
    output=CoreFHorzExpAssetze_d1_z_e_semiz_noa1(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): NOA1 figs 3-4

%% ===== doPart(4): NOA1 cross-tests 2 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): NOA1 cross-tests 2 =====\n')
    %% noa1 semiz cross-tests
    % CrossTest 5 + semiz: noa1 vs withA1 model with n_a1=1 where a1 is ignored (should match bit-exact)
    output=CoreFHorzExpAssetze_CrossTests5_nod1_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetze_CrossTests5_d1_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): NOA1 cross-tests 2




%% ================= WITH a1 (figs 5-8) =================

%% ===== doPart(5): WITH a1 figs 5-6 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): WITH a1 figs 5-6 =====\n')
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
end % doPart(5): WITH a1 figs 5-6


%% ===== doPart(6): WITH a1 cross-tests =====
if doPart(6)==1
    fprintf('\n===== doPart(6): WITH a1 cross-tests =====\n')
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
end % doPart(6): WITH a1 cross-tests





%% Semiz variants
addpath('./CoreFHorzExpAssetzeTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetze_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(7): WITH a1 figs 7-8 =====
if doPart(7)==1
    fprintf('\n===== doPart(7): WITH a1 figs 7-8 =====\n')
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
end % doPart(7): WITH a1 figs 7-8

%% ===== doPart(8): WITH a1 cross-tests 2 =====
if doPart(8)==1
    fprintf('\n===== doPart(8): WITH a1 cross-tests 2 =====\n')
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
end % doPart(8): WITH a1 cross-tests 2




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

%% ===== doPart(9): With TWO standard endogenous figs 9-10 =====
if doPart(9)==1
    fprintf('\n===== doPart(9): With TWO standard endogenous figs 9-10 =====\n')
    %% with2A1, without d1, with z, with e
    figure_c=9;
    output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, with z, with e
    figure_c=10;
    output=CoreFHorzExpAssetze_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(9): With TWO standard endogenous figs 9-10

%% ===== doPart(10): With TWO standard endogenous cross-tests =====
if doPart(10)==1
    fprintf('\n===== doPart(10): With TWO standard endogenous cross-tests =====\n')
    %% with2A1 nosemiz cross-tests
    % CrossTest 6: a degenerate second standard asset a1_2 (single point {0}) reduces the
    % two-standard-asset model back to the with-a1 model
    output=CoreFHorzExpAssetze_CrossTests6_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetze_CrossTests6_d1_with2A1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(10): With TWO standard endogenous cross-tests

%% with2A1 + semiz

%% ===== doPart(11): With TWO standard endogenous figs 11-12 =====
if doPart(11)==1
    fprintf('\n===== doPart(11): With TWO standard endogenous figs 11-12 =====\n')
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
end % doPart(11): With TWO standard endogenous figs 11-12

%% ================= with2A2: TWO EXPERIENCE assets (figs 13-20) =================
% vfoptions.experienceassetze is now the integer COUNT of a2 dimensions (2), not a flag.
% Grid layout a = [a1, a2_1, a2_2], and the noa1 tier is a = [a2_1, a2_2]. n_a_2A2 /
% a_grid_2A2 / n_a_2A2_justexpasset / a_grid_2A2_justexpasset all come from the setup.
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A2_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A2_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A2_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzExpAssetzeTests_subcodes/With2A2_subcodes/WithA1_subcodes/Semiz_subcodes/')


% The 'notsobig' a1 grid for the grid-interpolation comparison has to shrink relative to the
% with2A1 tier: the arrays scale with n_a1*n_a1prime(fine)*n_a2_1*n_a2_2*n_z*n_e, and this
% family always carries BOTH z and e. Fig 12 (with2A1+semiz+d1) already out-of-memories on
% this GPU, so the semiz variants below get a smaller a1 grid again.
% This grid is the FINE one: it feeds the "use a really big a_grid, then the moments should be
% essentially the same with/without grid interpolation" block, so it must have MORE a1 points than
% n_a_2A2 (=[51,n_a2_1,n_a2_2]), not the same number. 151 mirrors n_a_2A1_notsobig; the semiz
% variants get 101 because they also carry n_semiz and d3 (fig 12, the with2A1+semiz+d1 case at
% 151, already out-of-memories on this GPU).
n_a_2A2_notsobig=[151,n_a2_1,n_a2_2];
a_grid_2A2_notsobig=[5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;a2_1_grid;a2_2_grid];
n_a_2A2_notsobig_semiz=[101,n_a2_1,n_a2_2];
a_grid_2A2_notsobig_semiz=[5*linspace(0,1,n_a_2A2_notsobig_semiz(1))'.^3;a2_1_grid;a2_2_grid];

% ===== doPart(12): with2A2 NOA1 figs 13-14 =====
if doPart(12)==1
fprintf('\n===== doPart(12): with2A2 NOA1 figs 13-14 =====\n')
% with2A2 noa1, without d1
figure_c=13;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% with2A2 noa1, with d1
figure_c=14;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(12): with2A2 NOA1 figs 13-14

% ===== doPart(13): with2A2 NOA1 figs 15-16 =====
if doPart(13)==1
fprintf('\n===== doPart(13): with2A2 NOA1 figs 15-16 =====\n')
% with2A2 noa1 + semiz, without d1
figure_c=15;
output=CoreFHorzExpAssetze_nod1_z_e_semiz_noa1_with2A2(n_d_withoutd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% with2A2 noa1 + semiz, with d1
figure_c=16;
output=CoreFHorzExpAssetze_d1_z_e_semiz_noa1_with2A2(n_d_withd1semiz,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(13): with2A2 NOA1 figs 15-16

% ===== doPart(14): with2A2 WITH a1 figs 17-18 =====
if doPart(14)==1
fprintf('\n===== doPart(14): with2A2 WITH a1 figs 17-18 =====\n')
% with2A2 with a1, without d1
figure_c=17;
output=CoreFHorzExpAssetze_nod1_z_e_nosemiz_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% with2A2 with a1, with d1
figure_c=18;
output=CoreFHorzExpAssetze_d1_z_e_nosemiz_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(14): with2A2 WITH a1 figs 17-18

% ===== doPart(15): with2A2 WITH a1 figs 19-20 =====
if doPart(15)==1
fprintf('\n===== doPart(15): with2A2 WITH a1 figs 19-20 =====\n')
% with2A2 with a1 + semiz, without d1
figure_c=19;
output=CoreFHorzExpAssetze_nod1_z_e_semiz_with2A2(n_d_withoutd1semiz,n_a_2A2,n_a_2A2_notsobig_semiz,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A2,a_grid_2A2_notsobig_semiz,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

% with2A2 with a1 + semiz, with d1
figure_c=20;
output=CoreFHorzExpAssetze_d1_z_e_semiz_with2A2(n_d_withd1semiz,n_a_2A2,n_a_2A2_notsobig_semiz,n_z,N_j,d_grid_withd1semiz,a_grid_2A2,a_grid_2A2_notsobig_semiz,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetzeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(15): with2A2 WITH a1 figs 19-20

%% ===== doPart(16): with2A2 cross-tests =====
if doPart(16)==1
fprintf('\n===== doPart(16): with2A2 cross-tests =====\n')
% These are the machine-precision checks of the with2A2 tier; the figures above are consistency
% checks only. CrossTests7 makes the SECOND experience asset inert, CrossTests8 the FIRST; the two
% together pin down the ordering of the a2 dimensions, which either alone would not. CrossTests9
% keeps BOTH assets alive and checks swap symmetry, which is the only one of the three that
% exercises the nested interpolation where all four corners carry weight.
% Run on the one-experience-asset grids: each test builds its own second a2 dimension.
output=CoreFHorzExpAssetze_CrossTests7_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests7_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests8_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests8_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests9_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests9_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
% noa1 twins of CrossTests7: the noa1 raws are a distinct code path (the aprime index IS the a2
% Kron index, with no a1prime offset), and on the plain-ExpAsset side this is the test that caught
% the first real defect.
output=CoreFHorzExpAssetze_CrossTests7_nod1_noa1_with2A2(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
output=CoreFHorzExpAssetze_CrossTests7_d1_noa1_with2A2(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(16): with2A2 cross-tests

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzExpAssetzeTestsdiary.txt')

diary off
