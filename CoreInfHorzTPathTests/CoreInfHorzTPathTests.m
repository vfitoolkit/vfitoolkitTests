% Implement core TPath tests of the VFI Toolkit InfHorz commands.
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
% with/without semiz

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): one endogenous state (figs 1-8 live)
% doPart(2): semiz (figs 9-16), the eight semiz models (without/with d1, z, e)
% doPart(3): e cross-tests, the same path two ways: (D) e-as-z, (E) trivial e, (G) z and e as one 2-dim z (no figures)
% doPart(4): semiz cross-tests, a semiz that is really a markov vs the same process as z: alone, with z, with e, with z and e (no figures)
% doPart(5): _tminus1 and _tplus1 (no figures)
% doPart(6): localsearch, the local-search value fn step, without and with a decision variable (no figures, no general eqm)
% doPart(7): localsearch inside a real general eqm transition path, by shooting (no figures)
% doPart(8): ExogShockFn and EiidShockFn on the path, vs explicit grids (no figures)
% doPart(9): two endogenous states (figs 17-18)
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
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreInfHorzTPathTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzTPathTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreInfHorzTPathTestsdiary.txt
fprintf('CoreInfHorzTPathTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreInfHorzTPathTests_subcodes/')
addpath('./CoreInfHorzTPathTests_Setup/')
addpath('./CoreInfHorzTPath_ReturnFns/')

% Setup so that use the same d,a,z,e in all the models that use them
CoreInfHorzTPath_setup

%% ===== doPart(1): one endogenous state (figs 1-8 live) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): one endogenous state (figs 1-8 live) =====\n')
    %% without d, without z, without e
    figure_c=1;
    output=CoreInfHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, without e
    figure_c=2;
    output=CoreInfHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, without e
    figure_c=3;
    output=CoreInfHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e
    figure_c=4;
    output=CoreInfHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, with e
    figure_c=5;
    output=CoreInfHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, with e
    figure_c=6;
    output=CoreInfHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, with e
    figure_c=7;
    output=CoreInfHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, with e
    figure_c=8;
    output=CoreInfHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): one endogenous state (figs 1-8 live)

%% ===== doPart(2): semiz (figs 9-16) =====
if doPart(2)==1
    fprintf('\n===== doPart(2): semiz (figs 9-16) =====\n')
    %% The eight one-endogenous-state models again, with a semi-exogenous state semiz
    % Written 2026-10-27 test-first. semiz is employment, whose transitions depend on a binary search effort d2,
    % the last decision variable (see the semiz block at the end of the setup). Each subcode is built from the
    % doPart(1) subcode with the same z/e structure (nod1 from nod, d1 from d). The semiz fields are on
    % vfoptionsbaseline_semiz/simoptionsbaseline_semiz. For models without d1, use n_d2_semiz and d2_grid_semiz
    % (as n_d and d_grid, also in the GE block); for models with d1, use n_d_semiz and d_grid_semiz, and
    % n_d_semiz_GE and d_grid_semiz_GE in the GE block.

    %% without d1, without z, without e, with semiz
    figure_c=9;
    output=CoreInfHorzTPath_nod1_noz_noe_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d2_semiz,n_a_GE,d2_grid_semiz,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e, with semiz
    figure_c=10;
    output=CoreInfHorzTPath_d1_noz_noe_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d_semiz_GE,n_a_GE,d_grid_semiz_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e, with semiz
    figure_c=11;
    output=CoreInfHorzTPath_nod1_z_noe_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d2_semiz,n_a_GE,d2_grid_semiz,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e, with semiz
    figure_c=12;
    output=CoreInfHorzTPath_d1_z_noe_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d_semiz_GE,n_a_GE,d_grid_semiz_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e, with semiz
    figure_c=13;
    output=CoreInfHorzTPath_nod1_noz_e_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d2_semiz,n_a_GE,d2_grid_semiz,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e, with semiz
    figure_c=14;
    output=CoreInfHorzTPath_d1_noz_e_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d_semiz_GE,n_a_GE,d_grid_semiz_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e, with semiz
    figure_c=15;
    output=CoreInfHorzTPath_nod1_z_e_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_a_big,n_z,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d2_semiz,n_a_GE,d2_grid_semiz,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e, with semiz
    figure_c=16;
    output=CoreInfHorzTPath_d1_z_e_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_a_big,n_z,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,n_d_semiz_GE,n_a_GE,d_grid_semiz_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(2): semiz (figs 9-16)

%% ===== doPart(3): e cross-tests =====
if doPart(3)==1
    fprintf('\n===== doPart(3): e cross-tests =====\n')
    %% The same transition path two ways that must agree, for each way an iid e can be disguised
    % The transition path counterpart of CoreInfHorzTests' CoreInfHorz_CrossTests_{nod,d}_nosemiz_e, split
    % nod/d the same way. No figures. (D) an iid e vs the same process as a markov z with identical rows;
    % (E) a trivial e (n_e=1, e_grid=1, pi_e=1) vs no e, without and with z; (G) a z and an e vs the same
    % two shocks as one two-dimensional markov z. Each compares VPath, PolicyPath (PolicyVals under GI),
    % AgentDistPath, the AggVars along the path, and a general eqm path (Anderson), without and with GI.
    output=CoreInfHorzTPath_CrossTests_nod_nosemiz_e(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzTPath_CrossTests_d_nosemiz_e(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(3): e cross-tests

%% ===== doPart(4): semiz cross-tests =====
if doPart(4)==1
    fprintf('\n===== doPart(4): semiz cross-tests =====\n')
    %% A semiz that is really just a markov (pi_semiz ignores d2, equals pi_z, semiz_grid=z_grid) vs that process as a markov z
    % The transition path counterpart of CoreInfHorzTests' CoreInfHorz_CrossTests2_{nod1,d1}_semiz. No figures.
    % (a) alone, (b) with an additional z (vs a two-dim z), (c) with an iid e, (d) with z and e. Each compares
    % VPath, PolicyPath (PolicyVals under GI), AgentDistPath, the AggVars along the path, and a general eqm
    % path (Anderson, maxiter=5), without and with GI.
    output=CoreInfHorzTPath_CrossTests_nod1_semiz(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
    output=CoreInfHorzTPath_CrossTests_d1_semiz(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): semiz cross-tests

%% ===== doPart(5): _tminus1 and _tplus1 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): _tminus1 and _tplus1 =====\n')
    %% _tminus1 and _tplus1 in GeneralEqmEqns, FnsToEvaluate and intermediateEqns, each against the same shift done by hand
    % The InfHorz counterpart of CoreFHorzTPathTests doPart(6). No figures. Also the only place in the
    % InfHorz banks where intermediateEqns is checked for CORRECTNESS rather than merely run, and the
    % only coverage of the lag machinery in EvalFnOnTransPath_AggVars_InfHorz and
    % SimPanelValues_TransPath_InfHorz (both of which had to be given transpathoptions for it).
    output=CoreInfHorzTPath_tminus1tplus1(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(5): _tminus1 and _tplus1

%% ===== doPart(6): localsearch, the local-search value fn step, without and with a decision variable =====
if doPart(6)==1
    fprintf('\n===== doPart(6): localsearch, the local-search value fn step, without and with a decision variable =====\n')
    %% vfoptions.localsearch=1: aprime restricted to a window around an aprimeReferencePolicy
    % No figures and NO GENERAL EQM: the feature is tested here purely as a value function solver,
    % through ValueFnOnTransPath_InfHorz and one direct call to the single-step dispatcher. Working it
    % into the general eqm transition path is separate and comes later.
    % Runs on n_a as well as n_a_big: the saving the method promises is N_a/(2*nlocalsearch+1), and at
    % n_a=101 the step is kernel-launch bound, so it can only be measured on the big grid.
    % Sections (10)-(14) are the same thing with the grid interpolation layer on, where the window is
    % the same 2*nlocalsearch+1 coarse points with ngridinterp between each consecutive pair.
    % Section (15) times it against divide-and-conquer, which is the benchmark that matters: DC is
    % already in the toolkit and already exact under its own assumption, so if DC banks the saving
    % then local search has nothing to win.
    output=CoreInfHorzTPath_localsearch(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);

    %% vfoptions.localsearch=1 on a model WITH d, where the reference is conditional on d
    % The counterpart of the no-d case above. The reference is [N_d,N_a,N_z], one
    % aprime index per (d,a,z), because d enters the budget and so the optimal aprime moves with it.
    % d itself is searched in full, since an unordered d makes a window in the d index meaningless.
    % Also the only coverage of the raws' third output, the aprime chosen for EACH d, which no other
    % raw in the toolkit produces and which is what a warm start from a previous pass would use.
    output=CoreInfHorzTPath_localsearch_d(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(6): localsearch, the local-search value fn step, without and with a decision variable

%% ===== doPart(7): localsearch inside a general eqm transition path, by shooting =====
if doPart(7)==1
    fprintf('\n===== doPart(7): localsearch inside a general eqm transition path, by shooting =====\n')
    %% The first part where the restricted search sits inside a general eqm solve and so could change the answer
    % Iteration 1 is a standard sweep and seeds the reference; later iterations restrict the search
    % around the previous answer; on convergence one standard sweep re-evaluates the same path, and
    % either certifies it or supplies a fresh reference and local search resumes. Shooting only --
    % Anderson and quasi-Newton are refused by name, which is checked here, because neither can carry
    % the reference across iterations and both would silently return a wrong path if allowed to try.
    output=CoreInfHorzTPath_localsearch_shooting(T,n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(7): localsearch inside a general eqm transition path, by shooting

%% ===== doPart(8): ExogShockFn and EiidShockFn on the path =====
if doPart(8)==1
    fprintf('\n===== doPart(8): ExogShockFn and EiidShockFn on the path =====\n')
    %% Shocks given by a function, constant and with a parameter on the ParamPath, each against the same shocks as explicit grids
    % ExogShockSetup_InfHorz_TPath read an undefined ParamPathStruct in the ParamPath branches until 2026-10-27,
    % and returned an empty e grid from an EiidShockFn with constant parameters.
    output=CoreInfHorzTPath_ShockFns(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline);
end % doPart(8): ExogShockFn and EiidShockFn on the path

%% ===== doPart(9): two endogenous states (figs 17-18) =====
if doPart(9)==1
    fprintf('\n===== doPart(9): two endogenous states (figs 17-18) =====\n')
    %% with2A: TWO endogenous states (triggers the DC2A/GI2A code paths)
    % Simplified Kitao (2008) entrepreneur model: endogenous states are asset and occupation.
    % Uses its own grids and its own price/param path (both built in the setup), because the model
    % needs a calibration of its own; T is shared with the 1A subcodes.

    %% without d, with z, without e
    figure_c=17;
    output=CoreInfHorzTPath_nod_z_noe_nosemiz_with2A(T,PricePath_2A,ParamPath_2A,n_d_2A,n_a_2A,n_a_2A_big,n_z_2A,d_grid_2A,a_grid_2A,a_grid_2A_big,z_grid_2A,pi_z_2A,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_2A,n_a_2A_GE,d_grid_2A,a_grid_2A_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e
    figure_c=18;
    output=CoreInfHorzTPath_d_z_noe_nosemiz_with2A(T,PricePath_2A,ParamPath_2A,n_d_2A,n_a_2A,n_a_2A_big,n_z_2A,d_grid_2A,a_grid_2A,a_grid_2A_big,z_grid_2A,pi_z_2A,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_2A,n_a_2A_GE,d_grid_2A,a_grid_2A_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(9): two endogenous states (figs 17-18)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzTPathTestsdiary.txt')

diary off
