% Implement lots of tests of the VFI Toolkit FHorz transition-path commands
% under EPSTEIN-ZIN preferences
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
% with/without fastOLG
%
% with/without semiz [NOT SUPPORTED with Epstein-Zin on a transition path: refused, see doPart(5)]
% experienceasset    [NOT SUPPORTED with Epstein-Zin on a transition path: refused, see doPart(5)]
%
% The transition-path mirror of CoreFHorzEZTests, built the way CoreFHorzTPathQHTests mirrors CoreFHorzQHTests.
% Each subcode runs its model three times, the three Epstein-Zin cases:
%   Case 1: consumption units (traditional Epstein-Zin, vfoptions.EZutils=0)
%   Case 2: utility units with positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility units with negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% The sharpest check in each is the do-nothing path: with prices and parameters held at their values in
% Params, every period of the path must equal the stationary Epstein-Zin solve exactly, which ties every
% transition-path raw to its stationary twin (and so to CoreFHorzEZTests' brute-force oracle). See the
% header of any EZFHorzTPath_*_nosemiz subcode for the full list.
%
% Written TEST-FIRST (2026-10-30), before Epstein-Zin existed on the FHorz transition path. Until then
% slowOLG solved an Epstein-Zin model as if it were vNM, without a word, and fastOLG errored. Run against
% that code the slowOLG checks of doPart(1) fail and the first fastOLG solve stops the run.
%
% Two warnings are expected output rather than failures:
%  - "Epstein-Zin preferences do not make much sense without any shocks", from the stationary solves of
%    the no-shock models (figs 1-2 and cross test (1)).
%  - the vfoptions.WarmGlowBequestsFn / survival probability warnings, on the warm-glow legs.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): figs 1-8
% doPart(2): cross-tests
% doPart(3): Epstein-Zin parameters on the ParamPath (no figures)
% doPart(4): feasibility conventions: -Inf states and zero returns, on the path (no figures)
% doPart(5): refusals (no figures)
% doPart(6): local search, the value fn step: seed and restricted sweeps against the unrestricted EZ sweep (no figures)
% doPart(7): local search inside a general eqm transition path, against the same solve without (no figures)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzTPathEZTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzTPathEZTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzTPathEZTestsdiary.txt
fprintf('CoreFHorzTPathEZTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreFHorzTPathTests_Setup/')
addpath('./EZ_ReturnFns/')

addpath('./withEpsteinZinPreferences_subcodes/')
addpath('./withEpsteinZinPreferences_subcodes/CrossTests/')
% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorzTPath_setup

% Epstein-Zin preference parameters (the EZ vfoptions themselves are set inside each subcode), as in CoreFHorzEZTests
Params.ezgamma=3; % consumption-units risk aversion (vfoptions.EZriskaversion when EZutils=0; must be >1)
Params.ezphi=0.5; % consumption-units elasticity of intertemporal substitution (vfoptions.EZeis; cannot be 1)
Params.ezrisk=3; % utility-units additional risk aversion (vfoptions.EZriskaversion when EZutils=1)
Params.ezsigma=2; % curvature of the positiveUtils/negativeUtils utility fns (>1; the positiveUtils
% family uses a (1+x) shift so it stays strictly positive despite ezsigma>1)

% Survival-probability and warm-glow parameters, as in CoreFHorzEZTests
Params.sj=[linspace(1,0.6,N_j-1),0]; % declining survival, sj(N_j)=0: warm-glow active at every age AND the terminal-age convention is exercised
Params.wg1=2; % De Nardi warm-glow: strength of the bequest motive (theta)
Params.wg2=1; % De Nardi warm-glow: luxury-good shifter (kappa; bequests are luxury goods)
Params.wg3=Params.ezsigma; % De Nardi warm-glow: curvature, set equal to the utility curvature (utility-units cases)
Params.ezmrisk=5; % EZmortalityriskaversion (mortality risk aversion, distinct from the within-period risk aversion)

% vfoptions.exoticpreferences='EpsteinZin';
% Case 1: vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
% Case 2: vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
% Case 3: vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';

% Note: for Epstein-Zin, V_final is just the V of the final stationary solve (ValueFnIter_Case1_FHorz with
% the same vfoptions). The subcodes build it from the stationary solve, never from zeros: under Epstein-Zin
% a continuation of exactly zero is a value (a catastrophe), not 'nothing after the path'.

%% ===== doPart(1): figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): figs 1-8 =====\n')
    %% without d, without z, without e, without semiz
    figure_c=1;
    output=EZFHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, without e, without semiz
    figure_c=2;
    output=EZFHorzTPath_d_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, without e, without semiz
    figure_c=3;
    output=EZFHorzTPath_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e, without semiz
    figure_c=4;
    output=EZFHorzTPath_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, with e, without semiz
    figure_c=5;
    output=EZFHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, with e, without semiz
    figure_c=6;
    output=EZFHorzTPath_d_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, with e, without semiz
    figure_c=7;
    output=EZFHorzTPath_nod_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, with e, without semiz
    figure_c=8;
    output=EZFHorzTPath_d_z_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzTPathEZTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): figs 1-8

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    % The same model with different shock set-ups, on the moving path, plus the Epstein-Zin-specific
    % (z,e) against the joint markov (one certainty equivalent over (zprime,eprime), not two nested ones)
    output=EZFHorzTPath_CrossTests_nod_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

    output=EZFHorzTPath_CrossTests_d_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): cross-tests

%% ===== doPart(3): Epstein-Zin parameters on the ParamPath =====
if doPart(3)==1
    fprintf('\n===== doPart(3): Epstein-Zin parameters on the ParamPath =====\n')
    % Risk aversion, the EIS and the survival probabilities held on the ParamPath at new values, and a step
    % in risk aversion and in the survival probabilities part way along the path. See the header of the subcode.
    output=EZFHorzTPath_parampath(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline);
end % doPart(3): Epstein-Zin parameters on the ParamPath

%% ===== doPart(4): feasibility conventions: -Inf states and zero returns, on the path =====
if doPart(4)==1
    fprintf('\n===== doPart(4): feasibility conventions: -Inf states and zero returns, on the path =====\n')
    % The models of CoreFHorzEZTests doPart(8), on the moving path. See the header of the subcode.
    output=EZFHorzTPath_feasibility(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(4): feasibility conventions

%% ===== doPart(5): refusals =====
if doPart(5)==1
    fprintf('\n===== doPart(5): refusals =====\n')
    % What Epstein-Zin on the FHorz transition path does not do is refused by name, never silently solved
    % as something else. Checked as a set. See the header of the subcode.
    output=EZFHorzTPath_refusals(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline);
end % doPart(5): refusals

%% ===== doPart(6): local search, the value fn step =====
if doPart(6)==1
    fprintf('\n===== doPart(6): local search, the value fn step =====\n')
    % The EZ mirror of CoreFHorzTPathTests doPart(9): every layout, shock case, d/GI corner and EZ case, with
    % survival probabilities and the warm-glow. See the header of the subcode.
    output=EZFHorzTPath_localsearch(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params);
end % doPart(6): local search, the value fn step

%% ===== doPart(7): local search inside a general eqm transition path =====
if doPart(7)==1
    fprintf('\n===== doPart(7): local search inside a general eqm transition path =====\n')
    % The EZ mirror of CoreFHorzTPathTests doPart(10). See the header of the subcode.
    output=EZFHorzTPath_localsearch_shooting(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end % doPart(7): local search inside a general eqm transition path


%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreFHorzTPathEZTestsdiary.txt')

diary off

%% THINGS NOT CHECKED
% semiz and experienceasset with Epstein-Zin on a transition path (refused, doPart(5))
% Epstein-Zin with permanent types on a transition path (to come: stage 5)
% the general eqm solves with Epstein-Zin beyond doPart(7) (to come: stage 4)
