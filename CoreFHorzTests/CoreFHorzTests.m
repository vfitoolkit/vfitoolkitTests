% Implement lots of tests of the core VFI Toolkit FHorz commands
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1):  nosemiz (figs 1-8)
% doPart(2):  nosemiz cross-tests
% doPart(3):  semiz (figs 9-16)
% doPart(4):  semiz cross-tests
% doPart(5):  semiz cross-tests 2 (semi-exo that is really a markov)
% doPart(6):  with2A nosemiz (figs 17-24)
% doPart(7):  with2A nosemiz cross-tests
% doPart(8):  with2A semiz (figs 25-32)
% doPart(9):  with2A semiz cross-tests
% doPart(10): TestFnsToEvaluate
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so Fig 20
% is the same test whatever doPart says, and a png from a previous run is never overwritten by a
% different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTestsdiary.txt
fprintf('CoreFHorzTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTests_subcodes/')
addpath('./CoreFHorzTests_Setup/')
addpath('./CoreFHorz_ReturnFns/')
addpath('./CoreFHorzTests_subcodes/CrossTests/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

%% ===== doPart(1): nosemiz (figs 1-8) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): nosemiz (figs 1-8) =====\n')
    %% without d, without z, without e, without semiz
    figure_c=1;
    output=CoreFHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    % Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
    % all just 1e-3, so irrelevant. Is because interpolation creates tiny amount of variance where the is none.
    % (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

    %% with d, without z, without e, without semiz
    figure_c=2;
    output=CoreFHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    % Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
    % all just 1e-3, so irrelevant.  Is because interpolation creates tiny amount of variance where the is none.
    % (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

    %% without d, with z, without e, without semiz
    figure_c=3;
    output=CoreFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, with z, without e, without semiz
    figure_c=4;
    output=CoreFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, without z, with e, without semiz
    figure_c=5;
    output=CoreFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, without z, with e, without semiz
    figure_c=6;
    output=CoreFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, with z, with e, without semiz
    figure_c=7;
    output=CoreFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, with z, with e, without semiz
    figure_c=8;
    output=CoreFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good
end % doPart(1): nosemiz (figs 1-8)

%% ===== doPart(2): nosemiz cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): nosemiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % all looking good :)
end % doPart(2): nosemiz cross-tests

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

% d1 is a decision variable that is not in the SemiExoStateFn

addpath('./CoreFHorzTests_subcodes/')
addpath('./CoreFHorzTests_Setup/')
addpath('./CoreFHorz_ReturnFns/')
addpath('./CoreFHorzTests_subcodes/CrossTests/')

addpath('./CoreFHorzTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorz_ReturnFns/Semiz_ReturnFns/')
addpath('./CoreFHorzTests_subcodes/With2A_subcodes/')
addpath('./CoreFHorzTests_subcodes/With2A_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzTests_subcodes/With2A_subcodes/CrossTests/')
addpath('./CoreFHorz_ReturnFns/With2A_ReturnFns/')
addpath('./CoreFHorz_ReturnFns/With2A_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorz_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% ===== doPart(3): semiz (figs 9-16) =====
if doPart(3)==1
    fprintf('\n===== doPart(3): semiz (figs 9-16) =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=9;
    output=CoreFHorz_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, without z, without e, with semiz
    figure_c=10;
    output=CoreFHorz_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, with z, without e, with semiz
    figure_c=11;
    output=CoreFHorz_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good :)

    %% with d1, with z, without e, with semiz
    figure_c=12;
    output=CoreFHorz_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, without z, with e, with semiz
    figure_c=13;
    output=CoreFHorz_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, without z, with e, with semiz
    figure_c=14;
    output=CoreFHorz_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, with z, with e, with semiz
    figure_c=15;
    output=CoreFHorz_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, with z, with e, with semiz
    figure_c=16;
    output=CoreFHorz_d1_z_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good
end % doPart(3): semiz (figs 9-16)

%% ===== doPart(4): semiz cross-tests =====
if doPart(4)==1
    fprintf('\n===== doPart(4): semiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorz_CrossTests_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorz_CrossTests_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % All looks good!
end % doPart(4): semiz cross-tests

%% ===== doPart(5): semiz cross-tests 2 (semi-exo that is really a markov) =====
if doPart(5)==1
    fprintf('\n===== doPart(5): semiz cross-tests 2 (semi-exo that is really a markov) =====\n')
    %% Now some further cross-tests, using a semi-exo that is really just a markov

    output=CoreFHorz_CrossTests2_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorz_CrossTests2_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % All looks good!
end % doPart(5): semiz cross-tests 2 (semi-exo that is really a markov)














%% with2A: TWO standard endogenous states (triggers the DC2A/GI2A/DC2A_GI2A code paths)
% A genuine SECOND standard endogenous asset a2 (4 grid points) with preference params
% phi1,phi2. These tests were formerly the standalone CoreFHorzTwoEndoTests; they run the
% same models as above but with two standard endogenous states, sweeping DC2A / GI2A /
% DC2A_GI2A and lowmemory 0,1,2 (nosemiz) / 0,1,2,3 (semiz), plus cross-tests.
% Redefine the asset grid to two endogenous states for this section (the single-endo
% tests above have already run; TestFnsToEvaluate below still uses the single-asset n_a).
n_a_2A=[n_a,4];
n_a_2A_big=[n_a_big,4];
a2_grid_2A=[0;1;2;3];
a_grid_2A=[a_grid; a2_grid_2A];
a_grid_2A_big=[a_grid_big; a2_grid_2A];
Params.phi1=3; % second endo-state preference params
Params.phi2=0.1;
%% ===== doPart(6): with2A nosemiz (figs 17-24) =====
if doPart(6)==1
    fprintf('\n===== doPart(6): with2A nosemiz (figs 17-24) =====\n')
    %% without d, without z, without e, without semiz
    figure_c=17;
    output=CoreFHorz_nod_noz_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    % Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
    % all just 1e-3, so irrelevant. Is because interpolation creates tiny amount of variance where the is none.
    % (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

    %% with d, without z, without e, without semiz
    figure_c=18;
    output=CoreFHorz_d_noz_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    % Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
    % all just 1e-3, so irrelevant.  Is because interpolation creates tiny amount of variance where the is none.
    % (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

    %% without d, with z, without e, without semiz
    figure_c=19;
    output=CoreFHorz_nod_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, with z, without e, without semiz
    figure_c=20;
    n_a_notsobig=[501,4]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_d_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, without z, with e, without semiz
    figure_c=21;
    output=CoreFHorz_nod_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, without z, with e, without semiz
    figure_c=22;
    n_a_notsobig=[501,4]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_d_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, with z, with e, without semiz
    figure_c=23;
    output=CoreFHorz_nod_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, with z, with e, without semiz
    figure_c=24;
    n_a_notsobig=[301,4]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_d_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good
end % doPart(6): with2A nosemiz (figs 17-24)

%% ===== doPart(7): with2A nosemiz cross-tests =====
if doPart(7)==1
    fprintf('\n===== doPart(7): with2A nosemiz cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorz_CrossTests_nod_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorz_CrossTests_d_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % all looking good :)
end % doPart(7): with2A nosemiz cross-tests


%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

% d1 is a decision variable that is not in the SemiExoStateFn

% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% ===== doPart(8): with2A semiz (figs 25-32) =====
if doPart(8)==1
    fprintf('\n===== doPart(8): with2A semiz (figs 25-32) =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=25;
    output=CoreFHorz_nod1_noz_noe_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e, with semiz
    figure_c=26;
    output=CoreFHorz_d1_noz_noe_semiz_with2A(n_d_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e, with semiz
    figure_c=27;
    output=CoreFHorz_nod1_z_noe_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e, with semiz
    figure_c=28;
    n_a_notsobig=[501,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_d1_z_noe_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e, with semiz
    figure_c=29;
    output=CoreFHorz_nod1_noz_e_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e, with semiz
    figure_c=30;
    n_a_notsobig=[501,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_d1_noz_e_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e, with semiz
    figure_c=31;
    n_a_notsobig=[501,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_nod1_z_e_semiz_with2A(n_d2_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e, with semiz
    figure_c=32;
    n_a_notsobig=[251,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=CoreFHorz_d1_z_e_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(8): with2A semiz (figs 25-32)

%% ===== doPart(9): with2A semiz cross-tests =====
if doPart(9)==1
    fprintf('\n===== doPart(9): with2A semiz cross-tests =====\n')
    %% Cross-tests for semiz
    output=CoreFHorz_CrossTests_nod1_semiz_with2A(n_d_semiz,n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorz_CrossTests_d1_semiz_with2A(n_d_semiz,n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(9): with2A semiz cross-tests














%% ===== doPart(10): TestFnsToEvaluate =====
if doPart(10)==1
    fprintf('\n===== doPart(10): TestFnsToEvaluate =====\n')
    %% FnsToEvaluate-focused tests (covers every FHorz FnsToEvaluate consumer + cross-validations + analytical-truth tests)
    TestFnsToEvaluate
end % doPart(10): TestFnsToEvaluate

%% Done! Damn that was a lot of tests. Glad that is over.
%% One verdict for the whole run
% The bank prints over a thousand checks; this reads the diary back and says plainly whether the
% run passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTestsdiary.txt')

diary off




%% THINGS NOT CHECKED
% Check using two decision variables in the semiz codes (both for d1 and for d2, and without d1)


