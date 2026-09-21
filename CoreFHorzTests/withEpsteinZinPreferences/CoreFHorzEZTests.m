% Implement lots of tests of the core VFI Toolkit FHorz commands, under EPSTEIN-ZIN preferences
% with/without d
% with/without z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% with/without semiz
%
% NOTE: only ONE-endogenous-state models for now. The two-endogenous-state (with2A) tests
% (mirroring QH figs 17-32) WILL BE ADDED LATER.
%
% This is the EZ mirror of CoreFHorzTests.m/CoreFHorzQHTests.m. Where QH runs each model twice
% (Naive top, Sophisticated bottom), EZ runs each model three times:
%   Case 1: consumption-units (traditional Epstein-Zin, vfoptions.EZutils=0)
%   Case 2: utility-units with positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units with negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% Each subcode then finishes with the EZ special tests (which play the role of the QH beta0=1 tests):
%   (i)   consumption-units: gamma=1/phi collapses EZ to standard vNM (Policy identical; V related
%         by V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)))
%   (ii)  utility-units: EZriskaversion=0 collapses EZ to standard vNM (positive and negative cases)
%   (iii) the same collapses under the grid interpolation layer
%   (iv)  vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
%
% The bank was written test-first, but every feature it waited on now exists: EZ runs with
% divide-and-conquer, with the grid interpolation layer (incl. DC+GI), and with semiz, and
% ValueFnFromPolicy_FHorz has an EpsteinZin branch that is exercised both with and without GI.
% Nothing here is expected to error any more.
%
% The GI layer interpolates the transformed continuation EV=E[(ezc4*V)^ezc5], so the gamma=1/phi
% collapse tests are exact under GI too.
%
% Two warnings are expected output rather than failures:
%  - "Epstein-Zin preferences do not make much sense without any shocks": Epstein-Zin without any
%    shocks is allowed (the certainty-equivalent is the identity), but every such solve says so.
%    Fires for the no-shock figures and the no-shock cross-test legs.
%  - the vfoptions.WarmGlowBequestsFn warning, on the warm-glow legs.
%
% Each subcode draws Case 1 as figure_c, Case 2 as 100+figure_c, and Case 3 as 200+figure_c.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): figs 1-8
% doPart(2): Brute-force oracle: brute-force oracle
% doPart(3): cross-tests
% doPart(4): cross-tests 2
% doPart(5): figs 9-16
% doPart(6): cross-tests 3
% doPart(7): cross-tests 4
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
if exist('../TestOutput/CoreFHorzEZTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzEZTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzEZTestsdiary.txt
fprintf('CoreFHorzEZTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreFHorzTests_Setup/')
addpath('./EZ_ReturnFns/')

addpath('./withEpsteinZinPreferences_subcodes/')
addpath('./withEpsteinZinPreferences_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

% Epstein-Zin preference parameters (the EZ vfoptions themselves are set inside each subcode)
Params.ezgamma=3; % consumption-units risk aversion (vfoptions.EZriskaversion when EZutils=0; must be >1)
Params.ezphi=0.5; % consumption-units elasticity of intertemporal substitution (vfoptions.EZeis; cannot be 1)
Params.ezrisk=3; % utility-units additional risk aversion (vfoptions.EZriskaversion when EZutils=1)
Params.ezsigma=2; % curvature of the positiveUtils/negativeUtils utility fns (>1; the positiveUtils
% family uses a (1+x) shift so it stays strictly positive despite ezsigma>1)

% Survival-probability and warm-glow parameters (special tests (vi)-(ix))
Params.sj=[linspace(1,0.6,N_j-1),0]; % declining survival, sj(N_j)=0: warm-glow active at every age AND the terminal-age convention is exercised
Params.oneminussj=1-Params.sj; % age-dependent ReturnFn parameter for the combined warm-glow reference wrappers
Params.wg1=2; % De Nardi warm-glow: strength of the bequest motive (theta)
Params.wg2=1; % De Nardi warm-glow: luxury-good shifter (kappa; bequests are luxury goods) (see LifeCycleModel12 of IntroToLifeCycleModels)
Params.wg3=Params.ezsigma; % De Nardi warm-glow: curvature, set equal to the utility curvature (utility-units cases; the cons-units warm-glow fn is in consumption units and gets its curvature from the EZ preferences)
Params.ezmrisk=5; % EZmortalityriskaversion for special test (viii) (mortality risk aversion, distinct from the within-period risk aversion)

% vfoptions.exoticpreferences='EpsteinZin';
% Case 1: vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
% Case 2: vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
% Case 3: vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';

%% The only functions worth testing are the value fn ones, as after you have Policy everything else is anyway ignoring the Epstein-Zin preferences

%% ===== doPart(1): figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): figs 1-8 =====\n')
    %% without d, without z, without e, without semiz
    % No shocks: every EZ solve prints the (expected) no-shocks warning
    figure_c=1;
    output=EZFHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d, without z, without e, without semiz
    % No shocks: every EZ solve prints the (expected) no-shocks warning
    figure_c=2;
    output=EZFHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% without d, with z, without e, without semiz
    figure_c=3;
    output=EZFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d, with z, without e, without semiz
    figure_c=4;
    output=EZFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% without d, without z, with e, without semiz
    figure_c=5;
    output=EZFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d, without z, with e, without semiz
    figure_c=6;
    output=EZFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% without d, with z, with e, without semiz
    figure_c=7;
    output=EZFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d, with z, with e, without semiz
    figure_c=8;
    output=EZFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)
end % doPart(1): figs 1-8

%% ===== doPart(2): Brute-force oracle: brute-force oracle =====
if doPart(2)==1
    fprintf('\n===== doPart(2): Brute-force oracle: brute-force oracle =====\n')
    %% Brute-force oracle: an independent implementation of the EZ recursion
    % Solves a tiny model (n_a=5, n_z=2, N_j=3, no d) twice: by explicit loops written directly
    % from the Epstein-Zin FORMULAS (gamma/phi/risk enter directly; no ezc coefficients, no
    % toolkit subcodes — the reference shares NOTHING with the toolkit implementation), and
    % through ValueFnIter_Case1_FHorz with the bank's standard EZ vfoptions. All three EZ cases,
    % plus sj+warm-glow variants that independently validate the inside-the-root terminal
    % warm-glow convention (Kraft-Munk-Weiss). This is the only test class that can catch
    % 'toolkit and reference both wrong in the same way', since every other check in this bank
    % compares two toolkit code paths against each other.
    output=EZFHorz_bruteforce_oracle(Params,DiscountFactorParamNames);
end % doPart(2): Brute-force oracle: brute-force oracle

%% ===== doPart(3): cross-tests =====
if doPart(3)==1
    fprintf('\n===== doPart(3): cross-tests =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid
    % Note: for EZ these are NOT redundant (unlike for QH): they check that an iid shock enters the
    % same single certainty-equivalent whether it lives in z or in e.

    output=EZFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=EZFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(3): cross-tests

%% ===== doPart(4): cross-tests 2 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): cross-tests 2 =====\n')
    %% CrossTests3 (EZ-specific, no QH analog): markov z plus iid e, versus the same two shocks
    % as a single joint (two-variable) markov. THE key test that the codes take one CE over the joint
    % distribution of (zprime,eprime) rather than nesting CEs (which only coincides when gamma=1/phi).

    output=EZFHorz_CrossTests3_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=EZFHorz_CrossTests3_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): cross-tests 2















%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run second half independent of first half

%% That is all the without semiz, now with semiz
% From here on, it is the eight with semiz
% From here on, use n_d_semiz and d_grid_semiz as the inputs (instead of n_d and d_grid)

% d1 is a decision variable that is not in the SemiExoStateFn

% TEST-FIRST: no EZ SemiExo solvers exist yet; ALL of figs 9-16 and the semiz cross-tests are
% expected to error at the first solve until they are implemented.

addpath('../CoreFHorzTests_Setup/')
addpath('./EZ_ReturnFns/')

addpath('./withEpsteinZinPreferences_subcodes/')
addpath('./withEpsteinZinPreferences_subcodes/CrossTests/')

addpath('./withEpsteinZinPreferences_subcodes/Semiz_subcodes/')
addpath('./EZ_ReturnFns/Semiz_ReturnFns/')
% Uses the same setup, which already had a semi-exogenous state, just that it wasn't used.
CoreFHorz_setup

% For models without d1, use:
% n_d2_semiz and d2_grid_semiz (as n_d and d_grid)
% For models with d1, use:
% n_d_semiz and d_grid_semiz (as n_d and d_grid)

% Epstein-Zin preference parameters (see notes in first half)
Params.ezgamma=3;
Params.ezphi=0.5;
Params.ezrisk=3;
Params.ezsigma=2;

% Survival-probability and warm-glow parameters (special tests (vi)-(ix))
Params.sj=[linspace(1,0.6,N_j-1),0]; % declining survival, sj(N_j)=0: warm-glow active at every age AND the terminal-age convention is exercised
Params.oneminussj=1-Params.sj; % age-dependent ReturnFn parameter for the combined warm-glow reference wrappers
Params.wg1=2; % De Nardi warm-glow: strength of the bequest motive (theta)
Params.wg2=1; % De Nardi warm-glow: luxury-good shifter (kappa; bequests are luxury goods) (see LifeCycleModel12 of IntroToLifeCycleModels)
Params.wg3=Params.ezsigma; % De Nardi warm-glow: curvature, set equal to the utility curvature (utility-units cases; the cons-units warm-glow fn is in consumption units and gets its curvature from the EZ preferences)
Params.ezmrisk=5; % EZmortalityriskaversion for special test (viii) (mortality risk aversion, distinct from the within-period risk aversion)

%% ===== doPart(5): figs 9-16 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): figs 9-16 =====\n')
    %% without d1, without z, without e, with semiz
    figure_c=9;
    output=EZFHorz_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d1, without z, without e, with semiz
    figure_c=10;
    output=EZFHorz_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% without d1, with z, without e, with semiz
    figure_c=11;
    output=EZFHorz_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d1, with z, without e, with semiz
    figure_c=12;
    output=EZFHorz_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% without d1, without z, with e, with semiz
    figure_c=13;
    output=EZFHorz_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d1, without z, with e, with semiz
    figure_c=14;
    output=EZFHorz_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% without d1, with z, with e, with semiz
    figure_c=15;
    output=EZFHorz_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

    %% with d1, with z, with e, with semiz
    figure_c=16;
    output=EZFHorz_d1_z_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
    exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
    exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)
end % doPart(5): figs 9-16

%% ===== doPart(6): cross-tests 3 =====
if doPart(6)==1
    fprintf('\n===== doPart(6): cross-tests 3 =====\n')
    %% Now some cross-tests, things like setting up a markov that is actually just an iid, make sure we get same result as just doing iid

    output=EZFHorz_CrossTests_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=EZFHorz_CrossTests_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(6): cross-tests 3

%% ===== doPart(7): cross-tests 4 =====
if doPart(7)==1
    fprintf('\n===== doPart(7): cross-tests 4 =====\n')
    %% Now some further cross-tests, using a semi-exo that is really just a markov

    output=EZFHorz_CrossTests2_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=EZFHorz_CrossTests2_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(7): cross-tests 4

%% Done! (Remember: the two-endogenous-state (with2A) EZ tests are still to be added.)
%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreFHorzEZTestsdiary.txt')

diary off
