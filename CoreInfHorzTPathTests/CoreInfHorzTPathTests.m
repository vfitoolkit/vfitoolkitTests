% Implement core TPath tests of the VFI Toolkit InfHorz commands.
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% No semiz (skipped intentionally for this build).


%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): one endogenous state (figs 3-4 live)
% doPart(2): two endogenous states (figs 9-10)
% doPart(3): _tminus1 and _tplus1 (no figures)
% doPart(4): localsearch, the local-search value fn step (no figures, no general eqm)
% doPart(5): localsearch with a decision variable (no figures, no general eqm)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[0,0,0,1,1];

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

%% ===== doPart(1): one endogenous state (figs 3-4 live) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): one endogenous state (figs 3-4 live) =====\n')
    %% without d, without z, without e
    % figure_c=1;
    % output=CoreInfHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    % 
    % %% with d, without z, without e
    % figure_c=2;
    % output=CoreInfHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);

    %% without d, with z, without e
    figure_c=3;
    output=CoreInfHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e
    figure_c=4;
    output=CoreInfHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, with e
    % figure_c=5;
    % output=CoreInfHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    % 
    % %% with d, without z, with e
    % figure_c=6;
    % output=CoreInfHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    % 
    % %% without d, with z, with e
    % figure_c=7;
    % output=CoreInfHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
    % 
    % %% with d, with z, with e
    % figure_c=8;
    % output=CoreInfHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
end % doPart(1): one endogenous state (figs 3-4 live)

%% ===== doPart(2): two endogenous states (figs 9-10) =====
if doPart(2)==1
    fprintf('\n===== doPart(2): two endogenous states (figs 9-10) =====\n')
    %% with2A: TWO endogenous states (triggers the DC2A/GI2A code paths)
    % Simplified Kitao (2008) entrepreneur model: endogenous states are asset and occupation.
    % Uses its own grids and its own price/param path (both built in the setup), because the model
    % needs a calibration of its own; T is shared with the 1A subcodes.

    %% without d, with z, without e
    figure_c=9;
    output=CoreInfHorzTPath_nod_z_noe_nosemiz_with2A(T,PricePath_2A,ParamPath_2A,n_d_2A,n_a_2A,n_a_2A_big,n_z_2A,d_grid_2A,a_grid_2A,a_grid_2A_big,z_grid_2A,pi_z_2A,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_2A,n_a_2A_GE,d_grid_2A,a_grid_2A_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e
    figure_c=10;
    output=CoreInfHorzTPath_d_z_noe_nosemiz_with2A(T,PricePath_2A,ParamPath_2A,n_d_2A,n_a_2A,n_a_2A_big,n_z_2A,d_grid_2A,a_grid_2A,a_grid_2A_big,z_grid_2A,pi_z_2A,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_2A,n_a_2A_GE,d_grid_2A,a_grid_2A_GE,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzTPathTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(2): two endogenous states (figs 9-10)

%% ===== doPart(3): _tminus1 and _tplus1 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): _tminus1 and _tplus1 =====\n')
    %% _tminus1 and _tplus1 in GeneralEqmEqns, FnsToEvaluate and intermediateEqns, each against the same shift done by hand
    % The InfHorz counterpart of CoreFHorzTPathTests doPart(6). No figures. Also the only place in the
    % InfHorz banks where intermediateEqns is checked for CORRECTNESS rather than merely run, and the
    % only coverage of the lag machinery in EvalFnOnTransPath_AggVars_InfHorz and
    % SimPanelValues_TransPath_InfHorz (both of which had to be given transpathoptions for it).
    output=CoreInfHorzTPath_tminus1tplus1(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(3): _tminus1 and _tplus1

%% ===== doPart(4): localsearch, the local-search value fn step =====
if doPart(4)==1
    fprintf('\n===== doPart(4): localsearch, the local-search value fn step =====\n')
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
end % doPart(4): localsearch, the local-search value fn step

%% ===== doPart(5): localsearch with a decision variable =====
if doPart(5)==1
    fprintf('\n===== doPart(5): localsearch with a decision variable =====\n')
    %% vfoptions.localsearch=1 on a model WITH d, where the reference is conditional on d
    % The counterpart of doPart(4), which is the no-d case. The reference is [N_d,N_a,N_z], one
    % aprime index per (d,a,z), because d enters the budget and so the optimal aprime moves with it.
    % d itself is searched in full, since an unordered d makes a window in the d index meaningless.
    % Also the only coverage of the raws' third output, the aprime chosen for EACH d, which no other
    % raw in the toolkit produces and which is what a warm start from a previous pass would use.
    output=CoreInfHorzTPath_localsearch_d(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(5): localsearch with a decision variable

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzTPathTestsdiary.txt')

diary off
