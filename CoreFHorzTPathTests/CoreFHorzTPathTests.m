% Implement lots of tests of the core VFI Toolkit FHorz transition path commands, with one and two endogenous states
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
% doPart(1): figs 1-8 (each also checks AllStats and LifeCycleProfiles along the path)
% doPart(2): cross-tests
% doPart(3): figs 9-16 (each also checks AllStats and LifeCycleProfiles along the path)
% doPart(4): cross-tests 2
% doPart(5): cross-tests 3
% doPart(6): _tminus1 and _tplus1 (no figures)
% doPart(7): population demographics (no figures)
% doPart(8): z and e grids that vary over the path (no figures)
% doPart(9): localsearch, the local-search value fn step, both layouts, z/noz/e/noz_e, DC seeds (no figures, no general eqm)
% doPart(10): localsearch inside a general eqm transition path: shooting and Anderson, both layouts, DC, e, noz (no figures)
% doPart(11): shocks that vary by age, and four raw fixes (no figures)
% doPart(12): two endogenous states, without e (figs 17-20; was CoreFHorzTPathTwoEndoTests)
% doPart(13): two endogenous states, with e (figs 21-24; was CoreFHorzTPathTwoEndoTests)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathTestsdiary.txt
fprintf('CoreFHorzTPathTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTPathTests_subcodes/')
addpath('./CoreFHorzTPathTests_Setup/')
addpath('./CoreFHorzTPath_ReturnFns/')
addpath('./CoreFHorzTPathTests_subcodes/CrossTests/')
% The two endogenous state parts, doPart(12)-(13). Their ReturnFns are ReturnFn2A_*: the same model names as
% the one-asset ReturnFns, with two-asset signatures, so they cannot share a name on the path.
addpath('./CoreFHorzTPathTests_subcodes/TwoEndo_subcodes/')
addpath('./CoreFHorzTPath_ReturnFns/TwoEndo_ReturnFns/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPath_setup
% and the two endogenous state grids and parameters, under _2A names (touches no one-asset variable)
CoreFHorzTPathTwoEndo_setup

%% Uncomment to test with age-dependent shocks
% z_grid=z_grid.*ones(1,N_j,'gpuArray');
% pi_z=pi_z.*ones(1,1,N_j,'gpuArray');
% vfoptionsbaseline.pi_e=vfoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');
% vfoptionsbaseline.e_grid=vfoptionsbaseline.e_grid.*ones(1,N_j,'gpuArray');
% simoptionsbaseline.pi_e=simoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');
% simoptionsbaseline.pi_e=simoptionsbaseline.pi_e.*ones(1,N_j,'gpuArray');

%% ===== doPart(1): figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): figs 1-8 =====\n')
    %% without d, without z, without e, without semiz
    figure_c=1;
    output=CoreFHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, without z, without e, without semiz
    figure_c=2;
    output=CoreFHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, with z, without e, without semiz
    figure_c=3;
    output=CoreFHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, with z, without e, without semiz
    figure_c=4;
    output=CoreFHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, without z, with e, without semiz
    figure_c=5;
    output=CoreFHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, without z, with e, without semiz
    figure_c=6;
    output=CoreFHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d, with z, with e, without semiz
    figure_c=7;
    output=CoreFHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d, with z, with e, without semiz
    figure_c=8;
    output=CoreFHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good
end % doPart(1): figs 1-8

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorzTPath_CrossTests_nod_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzTPath_CrossTests_d_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): cross-tests

% all looking good :)


%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

% d1 is a decision variable that is not in the SemiExoStateFn

% NOTE: semiz on a TPath is handled by ValueFnOnTransPath_FHorz_SemiExo, which
% ValueFnOnTransPath_Case1_FHorz dispatches to whenever prod(vfoptions.n_semiz)>0

addpath('./CoreFHorzTPathTests_subcodes/')
addpath('./CoreFHorzTPathTests_Setup/')
addpath('./CoreFHorzTPath_ReturnFns/')
addpath('./CoreFHorzTPathTests_subcodes/CrossTests/')

addpath('./CoreFHorzTPathTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzTPath_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorzTPath_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

%% ===== doPart(3): figs 9-16 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): figs 9-16 =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=9;
    output=CoreFHorzTPath_nod1_noz_noe_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, without z, without e, with semiz
    figure_c=10;
    output=CoreFHorzTPath_d1_noz_noe_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, with z, without e, with semiz
    figure_c=11;
    output=CoreFHorzTPath_nod1_z_noe_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, with z, without e, with semiz
    figure_c=12;
    output=CoreFHorzTPath_d1_z_noe_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, without z, with e, with semiz
    figure_c=13;
    output=CoreFHorzTPath_nod1_noz_e_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, without z, with e, with semiz
    figure_c=14;
    output=CoreFHorzTPath_d1_noz_e_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% without d1, with z, with e, with semiz
    figure_c=15;
    output=CoreFHorzTPath_nod1_z_e_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good

    %% with d1, with z, with e, with semiz
    figure_c=16;
    output=CoreFHorzTPath_d1_z_e_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % looks good
end % doPart(3): figs 9-16

%% ===== doPart(4): cross-tests 2 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): cross-tests 2 =====\n')
    % Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    output=CoreFHorzTPath_CrossTests_nod1_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzTPath_CrossTests_d1_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): cross-tests 2

%% ===== doPart(5): cross-tests 3 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): cross-tests 3 =====\n')
    % Now some further cross-tests, using a semi-exo that is really just a markov
    output=CoreFHorzTPath_CrossTests2_nod1_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

    output=CoreFHorzTPath_CrossTests2_d1_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(5): cross-tests 3




%% ===== doPart(6): _tminus1 and _tplus1 =====
if doPart(6)==1
    fprintf('\n===== doPart(6): _tminus1 and _tplus1 =====\n')
    % _tminus1 and _tplus1 in GeneralEqmEqns, FnsToEvaluate and intermediateEqns, each against the same shift done by hand
    % Also a GE path whose condition uses lags, and EvalFnOnTransPath_AggVars with lags. See the header of the subcode.
    output=CoreFHorzTPath_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end % doPart(6): _tminus1 and _tplus1

%% ===== doPart(7): population demographics =====
if doPart(7)==1
    fprintf('\n===== doPart(7): population demographics =====\n')
    addpath('./CoreFHorzTPath_ReturnFns/Semiz_ReturnFns/') % the semiz model in (10); also added inside doPart(3), but parts are independent
    % Age weights (mewj) and survival (sj) that change over the path, and newborns (jequaloneDist) as a path. See the header of the subcode.
    output=CoreFHorzTPath_demographics(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,transpathoptionsbaseline);
    % (10) the age-weight and newborn checks again, on the models whose reshaping branches differ: no z, z and e, semiz
    output=CoreFHorzTPath_demographics_shapes(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,n_d2_semiz,d2_grid_semiz,Params,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(7): population demographics

%% ===== doPart(8): z and e grids that vary over the path =====
if doPart(8)==1
    fprintf('\n===== doPart(8): z and e grids that vary over the path =====\n')
    % z_grid and simoptions.e_grid with a T dimension, through ValueFnOnTransPath, AgentDistOnTransPath (fastOLG=0 and 1), one GE iteration, and the
    % AllStats, LifeCycleProfiles and AggVars path commands. See the header of the subcode.
    addpath('./CoreFHorzTPath_ReturnFns/Semiz_ReturnFns/') % the semiz models; also added inside doPart(3), but parts are independent
    output=CoreFHorzTPath_timevaryinggrids(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,n_d2_semiz,d2_grid_semiz,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(8): z and e grids that vary over the path

%% ===== doPart(9): localsearch, the local-search value fn step =====
if doPart(9)==1
    fprintf('\n===== doPart(9): localsearch, the local-search value fn step =====\n')
    %% vfoptions.localsearch=1: aprime restricted to a window around a reference, slowOLG
    % Direct calls to the single-step dispatcher on all four corners (without/with d, without/with GI),
    % every age at once. The strong check is that a window of three points centred on the exact seed
    % reference reproduces the unrestricted V exactly; the discriminator is that the same window around
    % the wrong reference does not. See the header of the subcode.
    output=CoreFHorzTPath_localsearch(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params);
end % doPart(9): localsearch, the local-search value fn step

%% ===== doPart(10): localsearch inside a general eqm transition path, by shooting =====
if doPart(10)==1
    fprintf('\n===== doPart(10): localsearch inside a general eqm transition path, by shooting =====\n')
    %% The CoreFHorzTPathAlgoTests firm model, a null reform from a bumped guess, plain shooting against
    % shooting with local search on all four corners; then ValueFnOnTransPath ignoring localsearch=1, and
    % the refusals as a set. See the header of the subcode.
    output=CoreFHorzTPath_localsearch_shooting(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end % doPart(10): localsearch inside a general eqm transition path, by shooting

%% ===== doPart(11): shocks that vary by age, and four raw fixes =====
if doPart(11)==1
    fprintf('\n===== doPart(11): shocks that vary by age, and four raw fixes =====\n')
    %% Every other part uses e and z grids and pi_e that do not vary with age, under which two of these bugs
    % gave the right numbers by accident. See the header of the subcode.
    output=CoreFHorzTPath_agedependent_rawfixes(T,PricePath,ParamPath,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames);
end % doPart(11): shocks that vary by age, and four raw fixes

%% ===== doPart(12): two endogenous states, without e (figs 17-20) =====
if doPart(12)==1
    fprintf('\n===== doPart(12): two endogenous states, without e (figs 17-20) =====\n')
    %% Two endogenous states (a1 continuous, a2 on four points), so the DC2A and GI2A code paths. These were
    % CoreFHorzTPathTwoEndoTests until 2026-10-28; their grids and parameters are the _2A ones from
    % CoreFHorzTPathTwoEndo_setup. No semiz yet.
    %% without d, without z, without e, without semiz
    figure_c=17;
    output=CoreFHorzTPathTwoEndo_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, without e, without semiz
    figure_c=18;
    output=CoreFHorzTPathTwoEndo_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, without e, without semiz
    figure_c=19;
    output=CoreFHorzTPathTwoEndo_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e, without semiz
    figure_c=20;
    output=CoreFHorzTPathTwoEndo_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(12): two endogenous states, without e (figs 17-20)

%% ===== doPart(13): two endogenous states, with e (figs 21-24) =====
if doPart(13)==1
    fprintf('\n===== doPart(13): two endogenous states, with e (figs 21-24) =====\n')
    %% As doPart(12), with an e variable. A part of its own because the with d, z and e model is the one
    % most likely to run out of memory (it uses the middle-sized a1 grid for that reason), so it can be
    % rerun on its own.
    %% without d, without z, with e, without semiz
    figure_c=21;
    output=CoreFHorzTPathTwoEndo_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, with e, without semiz
    figure_c=22;
    output=CoreFHorzTPathTwoEndo_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, with e, without semiz
    figure_c=23;
    output=CoreFHorzTPathTwoEndo_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_big_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_big_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, with e, without semiz
    figure_c=24;
    output=CoreFHorzTPathTwoEndo_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a_2A,n_a_notsobig_2A,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig_2A,z_grid,pi_z,Params_2A,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(13): two endogenous states, with e (figs 21-24)

%% Done! Damn that was a lot of tests. Glad that is over.




%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTPathTestsdiary.txt')

diary off

%% THINGS NOT CHECKED
% Check using two decision variables in the semiz codes (both for d1 and for d2, and without d1)
