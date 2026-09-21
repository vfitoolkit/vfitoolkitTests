% Implement lots of tests of the core VFI Toolkit FHorz commands, under QUASI-HYPERBOLIC discounting
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz
% with/without two standard endogenous assets (with2A)
%
% This is the QH mirror of CoreFHorzTests.m. Each subcode runs the model twice (Naive top,
% Sophisticated bottom) with a Valt (continuation-value) check beside every V/Policy check.
%
% TEST-FIRST STATE: nothing outstanding — everything this file tests is implemented and passing.
%
% Each subcode draws the Naive figure as figure_c and the Sophisticated figure as 100+figure_c.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): figs 1-8
% doPart(2): cross-tests
% doPart(3): figs 9-16
% doPart(4): cross-tests 2
% doPart(5): cross-tests 3
% doPart(6): figs 17-24
% doPart(7): cross-tests 4
% doPart(8): figs 25-32
% doPart(9): cross-tests 5
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzQHTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzQHTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzQHTestsdiary.txt
fprintf('CoreFHorzQHTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreFHorzTests_Setup/')
addpath('../CoreFHorz_ReturnFns/')

addpath('./withQuasiHyperbolicDiscounting_subcodes/')
addpath('./withQuasiHyperbolicDiscounting_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount='beta0';

% vfoptions.exoticpreferences='QuasiHyperbolic';
% vfoptions.quasi_hyperbolic='Naive';
% vfoptions.quasi_hyperbolic='Sophisticated';

%% The only functions worth testing are the value fn ones, as after you have Policy everything else is anyway ignoring the quasi-hyperbolic discounting

%% ===== doPart(1): figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): figs 1-8 =====\n')
    %% without d, without z, without e, without semiz
    figure_c=1;
    output=QHDFHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    % Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
    % all just 1e-3, so irrelevant. Is because interpolation creates tiny amount of variance where the is none.
    % (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

    %% with d, without z, without e, without semiz
    figure_c=2;
    output=QHDFHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    % Figure can appear to have an issue with std dev of assets, but if you look at the y-axis it is
    % all just 1e-3, so irrelevant.  Is because interpolation creates tiny amount of variance where the is none.
    % (Explanation: http://discourse.vfitoolkit.com/t/grid-interpolation-layer/394/12 )

    %% without d, with z, without e, without semiz
    figure_c=3;
    output=QHDFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% with d, with z, without e, without semiz
    figure_c=4;
    output=QHDFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% without d, without z, with e, without semiz
    figure_c=5;
    output=QHDFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% with d, without z, with e, without semiz
    figure_c=6;
    output=QHDFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% without d, with z, with e, without semiz
    figure_c=7;
    output=QHDFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% with d, with z, with e, without semiz
    figure_c=8;
    output=QHDFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good
end % doPart(1): figs 1-8

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    % Note: These are likely redundant, as already compared them all to exponential discounting and did these cross-tests for exponentinal
    % discounting. But it is going to test the very unlikely case that they still somehow manage to differ when beta0 is involved

    output=QHDFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=QHDFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): cross-tests

% all looking good :)














%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

% d1 is a decision variable that is not in the SemiExoStateFn

addpath('../CoreFHorzTests_Setup/')
addpath('../CoreFHorz_ReturnFns/')

addpath('./withQuasiHyperbolicDiscounting_subcodes/')
addpath('./withQuasiHyperbolicDiscounting_subcodes/CrossTests/')

addpath('./withQuasiHyperbolicDiscounting_subcodes/Semiz_subcodes/')
addpath('../CoreFHorz_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorz_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount='beta0';


%% ===== doPart(3): figs 9-16 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): figs 9-16 =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=9;
    output=QHDFHorz_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% with d1, without z, without e, with semiz
    figure_c=10;
    output=QHDFHorz_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% without d1, with z, without e, with semiz
    figure_c=11;
    output=QHDFHorz_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good :)

    %% with d1, with z, without e, with semiz
    figure_c=12;
    output=QHDFHorz_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% without d1, without z, with e, with semiz
    figure_c=13;
    output=QHDFHorz_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% with d1, without z, with e, with semiz
    figure_c=14;
    output=QHDFHorz_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% without d1, with z, with e, with semiz
    figure_c=15;
    output=QHDFHorz_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good

    %% with d1, with z, with e, with semiz
    figure_c=16;
    output=QHDFHorz_d1_z_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
    % looks good
end % doPart(3): figs 9-16

%% ===== doPart(4): cross-tests 2 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): cross-tests 2 =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    % Note: These are likely redundant, as already compared them all to exponential discounting and did these cross-tests for exponentinal
    % discounting. But it is going to test the very unlikely case that they still somehow manage to differ when beta0 is involved

    output=QHDFHorz_CrossTests_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=QHDFHorz_CrossTests_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): cross-tests 2

% All looks good!

%% ===== doPart(5): cross-tests 3 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): cross-tests 3 =====\n')
    %% Now some further cross-tests, using a semi-exo that is really just a markov

    output=QHDFHorz_CrossTests2_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=QHDFHorz_CrossTests2_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(5): cross-tests 3

% All looks good!







%% with2A: TWO standard endogenous states (triggers the DC2A/GI2A/DC2A_GI2A code paths), under QH
% QH mirror of CoreFHorzTests.m figs 17-32. Reuses the same With2A ReturnFns as the exponential suite.
% TEST-FIRST: the base method should pass, but the DC2A/GI2A/DC2A_GI2A solves error until plain-QH 2A
% raws are written (the QH DC/GI dispatchers currently error "only supports scalar n_a").

addpath('./withQuasiHyperbolicDiscounting_subcodes/With2A_subcodes/')
addpath('./withQuasiHyperbolicDiscounting_subcodes/With2A_subcodes/Semiz_subcodes/')
addpath('./withQuasiHyperbolicDiscounting_subcodes/With2A_subcodes/CrossTests/')
addpath('../CoreFHorz_ReturnFns/With2A_ReturnFns/')
addpath('../CoreFHorz_ReturnFns/With2A_ReturnFns/Semiz_ReturnFns/')

Params.beta0=0.9; % additional today-tomorrow discount factor
vfoptionsbaseline.QHadditionaldiscount='beta0';

% Redefine the asset grid to two endogenous states for this section
n_a_2A=[n_a,4];
n_a_2A_big=[n_a_big,4];
a2_grid_2A=[0;1;2;3];
a_grid_2A=[a_grid; a2_grid_2A];
a_grid_2A_big=[a_grid_big; a2_grid_2A];
Params.phi1=3; % second endo-state preference params
Params.phi2=0.1;

%% ===== doPart(6): figs 17-24 =====
if doPart(6)==1
    fprintf('\n===== doPart(6): figs 17-24 =====\n')
    %% without d, without z, without e, without semiz
    figure_c=17;
    output=QHDFHorz_nod_noz_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d, without z, without e, without semiz
    figure_c=18;
    output=QHDFHorz_d_noz_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d, with z, without e, without semiz
    figure_c=19;
    output=QHDFHorz_nod_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d, with z, without e, without semiz
    figure_c=20;
    n_a_notsobig=[501,4]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_d_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d, without z, with e, without semiz
    figure_c=21;
    output=QHDFHorz_nod_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d, without z, with e, without semiz
    figure_c=22;
    n_a_notsobig=[501,4]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_d_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d, with z, with e, without semiz
    figure_c=23;
    output=QHDFHorz_nod_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d, with z, with e, without semiz
    figure_c=24;
    n_a_notsobig=[301,4]; % to test Grid Interpolation
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_d_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
end % doPart(6): figs 17-24

%% ===== doPart(7): cross-tests 4 =====
if doPart(7)==1
    fprintf('\n===== doPart(7): cross-tests 4 =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=QHDFHorz_CrossTests_nod_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=QHDFHorz_CrossTests_d_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(7): cross-tests 4


%% That is all the without semiz, now with semiz (with2A)

%% ===== doPart(8): figs 25-32 =====
if doPart(8)==1
    fprintf('\n===== doPart(8): figs 25-32 =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=25;
    output=QHDFHorz_nod1_noz_noe_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d1, without z, without e, with semiz
    figure_c=26;
    output=QHDFHorz_d1_noz_noe_semiz_with2A(n_d_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d1, with z, without e, with semiz
    figure_c=27;
    output=QHDFHorz_nod1_z_noe_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d1, with z, without e, with semiz
    figure_c=28;
    n_a_notsobig=[401,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_d1_z_noe_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d1, without z, with e, with semiz
    figure_c=29;
    output=QHDFHorz_nod1_noz_e_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d1, without z, with e, with semiz
    figure_c=30;
    n_a_notsobig=[401,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_d1_noz_e_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% without d1, with z, with e, with semiz
    figure_c=31;
    n_a_notsobig=[401,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_nod1_z_e_semiz_with2A(n_d2_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)

    %% with d1, with z, with e, with semiz
    figure_c=32;
    n_a_notsobig=[251,4];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig; a2_grid_2A];
    output=QHDFHorz_d1_z_e_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Naive.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzQHTests_Fig',num2str(figure_c),'_Sophisticated.png'],'Resolution',150)
end % doPart(8): figs 25-32

%% ===== doPart(9): cross-tests 5 =====
if doPart(9)==1
    fprintf('\n===== doPart(9): cross-tests 5 =====\n')
    %% Cross-tests for semiz (with2A)
    output=QHDFHorz_CrossTests_nod1_semiz_with2A(n_d_semiz,n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=QHDFHorz_CrossTests_d1_semiz_with2A(n_d_semiz,n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(9): cross-tests 5

% Done! The QH mirror of CoreFHorzTests, including with2A.
%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreFHorzQHTestsdiary.txt')

diary off

