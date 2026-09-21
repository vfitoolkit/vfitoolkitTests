% QUASI-HYPERBOLIC discounting tests of the core VFI Toolkit FHorz ExpAsset commands.
% Mirrors CoreFHorzExpAssetTests.m one-for-one: the same 48 combinations in the same figure
% order -- noa1 nosemiz 1-8, noa1 semiz 9-16, withA1 nosemiz 17-24, withA1 semiz 25-32,
% with2A1 nosemiz 33-40, with2A1 semiz 41-48 -- over
%   with/without d1, with/without z, with/without e, with/without semiz,
%   noa1 / withA1 / with2A1, and every solution method the baseline runs.
%
% Each subcode runs Naive AND Sophisticated, the full lowmemory ladder on each method, a
% ValueFnFromPolicy oracle, and the exponential cross-tests: (i) the Naive continuation value
% equals the exponential value fn at the actual beta0; (ii)/(iii) with beta0=1 both the main
% value fn and the continuation value equal exponential, for Naive and Sophisticated.
%
% TEST-FIRST: the toolkit currently has NO quasi-hyperbolic support for experienceasset
% (ValueFnIter_Case1_FHorz has no QH branch for it), so every figure errors at its first
% ValueFnIter call. These were written ahead of the toolkit code, on purpose.
%
% The subcodes draw no figures, so only the diary is saved.
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzQHExpAssetTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzQHExpAssetTestsdiary.txt') % otherwise diary just appends to the previous run
end
%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): WITHOUT a1 figs 1-8
% doPart(2): WITHOUT a1 figs 9-16
% doPart(3): WITH a1 figs 17-24
% doPart(4): WITH a1 figs 25-32
% doPart(5): WITH 2a1 figs 33-40
% doPart(6): WITH 2a1 figs 41-48
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1];

diary ../TestOutput/CoreFHorzQHExpAssetTestsdiary.txt
fprintf('CoreFHorzQHExpAssetTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzQHExpAssetTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzQHExpAssetTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAssetTests_Setup/')
addpath('../CoreFHorzExpAsset_ReturnFns/')


%% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzExpAsset_setup

Params.beta0=0.9; % additional today-tomorrow (present-bias) discount factor
vfoptionsbaseline.QHadditionaldiscount='beta0';


%% ================= WITHOUT a1 (figs 1-16): experience asset is the only endogenous state =================
% No DC/GI/DC+GI blocks (irrelevant without a1).
% Pass n_a_justexpasset as n_a, a_grid_justexpasset as a_grid. n_a_big/a_grid_big slots unused.

addpath('./CoreFHorzQHExpAssetTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzQHExpAssetTests_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAsset_ReturnFns/Noa1_ReturnFns/')
addpath('../CoreFHorzExpAsset_ReturnFns/Noa1_ReturnFns/Semiz_ReturnFns/')

%% noa1 nosemiz (8 variants)

%% ===== doPart(1): WITHOUT a1 figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): WITHOUT a1 figs 1-8 =====\n')
    %% without d1, without z, without e, noa1, nosemiz
    figure_c=1;
    output=CoreFHorzQHExpAsset_nod1_noz_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, without z, without e, noa1, nosemiz
    figure_c=2;
    output=CoreFHorzQHExpAsset_d1_noz_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% without d1, with z, without e, noa1, nosemiz
    figure_c=3;
    output=CoreFHorzQHExpAsset_nod1_z_noe_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, with z, without e, noa1, nosemiz
    figure_c=4;
    output=CoreFHorzQHExpAsset_d1_z_noe_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% without d1, without z, with e, noa1, nosemiz
    figure_c=5;
    output=CoreFHorzQHExpAsset_nod1_noz_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, without z, with e, noa1, nosemiz
    figure_c=6;
    output=CoreFHorzQHExpAsset_d1_noz_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% without d1, with z, with e, noa1, nosemiz
    figure_c=7;
    output=CoreFHorzQHExpAsset_nod1_z_e_noa1_nosemiz(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, with z, with e, noa1, nosemiz
    figure_c=8;
    output=CoreFHorzQHExpAsset_d1_z_e_noa1_nosemiz(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
end % doPart(1): WITHOUT a1 figs 1-8

%% noa1 nosemiz cross-tests
% Markov-as-iid equivalence cross-test
% ExpAsset noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo model with d2
% no a1 vs model with a1 but where it is ignored


%% noa1 semiz (8 variants)

%% ===== doPart(2): WITHOUT a1 figs 9-16 =====
if doPart(2)==1
    fprintf('\n===== doPart(2): WITHOUT a1 figs 9-16 =====\n')
    %% without d1, without z, without e, noa1, semiz
    figure_c=9;
    output=CoreFHorzQHExpAsset_nod1_noz_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, without z, without e, noa1, semiz
    figure_c=10;
    output=CoreFHorzQHExpAsset_d1_noz_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% without d1, with z, without e, noa1, semiz
    figure_c=11;
    output=CoreFHorzQHExpAsset_nod1_z_noe_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, with z, without e, noa1, semiz
    figure_c=12;
    output=CoreFHorzQHExpAsset_d1_z_noe_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% without d1, without z, with e, noa1, semiz
    figure_c=13;
    output=CoreFHorzQHExpAsset_nod1_noz_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, without z, with e, noa1, semiz
    figure_c=14;
    output=CoreFHorzQHExpAsset_d1_noz_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% without d1, with z, with e, noa1, semiz
    figure_c=15;
    output=CoreFHorzQHExpAsset_nod1_z_e_noa1_semiz(n_d_withoutd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% with d1, with z, with e, noa1, semiz
    figure_c=16;
    output=CoreFHorzQHExpAsset_d1_z_e_noa1_semiz(n_d_withd1semiz,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1semiz,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
end % doPart(2): WITHOUT a1 figs 9-16

%% noa1 semiz cross-tests
% Markov-as-iid equivalence cross-test
% semiz state as a plain Markov should give same answer
% ExpAsset noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo model with d2
% no a1 vs model with a1 but where it is ignored




%% ================= WITH a1 (figs 17-32) =================
% Reset the setup (the without-a1 half above left the workspace alone, but re-run for safety/independence)
CoreFHorzExpAsset_setup

%% a1 nosemiz (8 variants)

%% ===== doPart(3): WITH a1 figs 17-24 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): WITH a1 figs 17-24 =====\n')
    %% without d1, without z, without e, without semiz
    figure_c=17;
    output=CoreFHorzQHExpAsset_nod1_noz_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% with d1, without z, without e, without semiz
    figure_c=18;
    output=CoreFHorzQHExpAsset_d1_noz_noe_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% without d1, with z, without e, without semiz
    figure_c=19;
    output=CoreFHorzQHExpAsset_nod1_z_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% with d1, with z, without e, without semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=20;
    output=CoreFHorzQHExpAsset_d1_z_noe_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% without d1, without z, with e, without semiz
    figure_c=21;
    output=CoreFHorzQHExpAsset_nod1_noz_e_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% with d1, without z, with e, without semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=22;
    output=CoreFHorzQHExpAsset_d1_noz_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% without d1, with z, with e, without semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=23;
    output=CoreFHorzQHExpAsset_nod1_z_e_nosemiz(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% with d1, with z, with e, without semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=24;
    output=CoreFHorzQHExpAsset_d1_z_e_nosemiz(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good
end % doPart(3): WITH a1 figs 17-24

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid

% all looking good :)

%% Do a test with a 'fake experience asset' and compare to a standard endogneous asset

% all looking good :)


%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)


addpath('../CoreFHorzExpAsset_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorzExpAsset_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% a1 semiz (8 variants)

%% ===== doPart(4): WITH a1 figs 25-32 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): WITH a1 figs 25-32 =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=25;
    output=CoreFHorzQHExpAsset_nod1_noz_noe_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % looks good

    %% with d1, without z, without e, with semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=26;
    output=CoreFHorzQHExpAsset_d1_noz_noe_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, with z, without e, with semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=27;
    output=CoreFHorzQHExpAsset_nod1_z_noe_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, with z, without e, with semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=28;
    output=CoreFHorzQHExpAsset_d1_z_noe_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, without z, with e, with semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=29;
    output=CoreFHorzQHExpAsset_nod1_noz_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% with d1, without z, with e, with semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=30;
    output=CoreFHorzQHExpAsset_d1_noz_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    % looks good (as good as it can be expected to given the n_a_notsobig)

    %% without d1, with z, with e, with semiz
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=31;
    output=CoreFHorzQHExpAsset_nod1_z_e_semiz(n_d_withoutd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    % looks good

    %% with d1, with z, with e, with semiz
    n_a_notsobig=[201,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=32;
    output=CoreFHorzQHExpAsset_d1_z_e_semiz(n_d_withd1semiz,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
end % doPart(4): WITH a1 figs 25-32

% looks good (as good as it can be expected to given the n_a_notsobig)

%% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid

% all looking good :)

%% Now some further cross-tests, using a semi-exo that is really just a markov

% all looking good :)

%% Do a test with a 'fake experience asset' and compare to a standard endogneous asset

% all looking good :)


%% ================= WITH 2a1 (figs 33-48): two standard endogenous assets (triggers DC2A/GI2A/DC2A_GI2A) =================
% A genuine multi-point second standard asset a1_2 (return r2) is spliced between the liquid asset
% a1 and the experience asset a2: n_a_2A1=[n_a1, n_a1_2, n_a2] (built in CoreFHorzExpAsset_setup).
% Two standard assets -> length(n_a1)>1 -> DC2A / GI2A / DC2A_GI2A.
% The 8 semiz variants (figs 41-48) exercise the ExpAssetSemiExo DC2A/GI2A/DC2A_GI2A family
% (written test-first; the toolkit raws now exist).
CoreFHorzExpAsset_setup
addpath('./CoreFHorzQHExpAssetTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzQHExpAssetTests_subcodes/With2A1_subcodes/Semiz_subcodes/')
addpath('../CoreFHorzExpAsset_ReturnFns/With2A1_ReturnFns/')
addpath('../CoreFHorzExpAsset_ReturnFns/With2A1_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(5): WITH 2a1 figs 33-40 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): WITH 2a1 figs 33-40 =====\n')
    %% 2a1 nosemiz (8 variants)
    n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=33;
    output=CoreFHorzQHExpAsset_nod1_noz_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=34;
    output=CoreFHorzQHExpAsset_d1_noz_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=35;
    output=CoreFHorzQHExpAsset_nod1_z_noe_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=36;
    output=CoreFHorzQHExpAsset_d1_z_noe_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=37;
    output=CoreFHorzQHExpAsset_nod1_noz_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=38;
    output=CoreFHorzQHExpAsset_d1_noz_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=39;
    output=CoreFHorzQHExpAsset_nod1_z_e_nosemiz_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=40;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    output=CoreFHorzQHExpAsset_d1_z_e_nosemiz_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(5): WITH 2a1 figs 33-40

%% 2a1 nosemiz cross-tests: a degenerate second asset a1_2 (single point {0}) reduces to the with-a1 model

%% ===== doPart(6): WITH 2a1 figs 41-48 =====
if doPart(6)==1
    fprintf('\n===== doPart(6): WITH 2a1 figs 41-48 =====\n')
    %% 2a1 semiz (8 variants)
    n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=41;
    output=CoreFHorzQHExpAsset_nod1_noz_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=42;
    output=CoreFHorzQHExpAsset_d1_noz_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    n_a_2A1_notsobig=[101,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=43;
    output=CoreFHorzQHExpAsset_nod1_z_noe_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    n_a_2A1_notsobig=[75,n_a1_2,n_a_justexpasset];
    a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
    a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

    figure_c=44;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    output=CoreFHorzQHExpAsset_d1_z_noe_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=45;
    output=CoreFHorzQHExpAsset_nod1_noz_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=46;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    output=CoreFHorzQHExpAsset_d1_noz_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    figure_c=47;
    output=CoreFHorzQHExpAsset_nod1_z_e_semiz_with2A1(n_d_withoutd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    figure_c=48;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    output=CoreFHorzQHExpAsset_d1_z_e_semiz_with2A1(n_d_withd1semiz,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1semiz,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(6): WITH 2a1 figs 41-48


%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreFHorzQHExpAssetTestsdiary.txt')

diary off

%% THINGS NOT CHECKED
% Check using two decision variables in any of d1 or d3 (the decision variables that are not in experience asset)
