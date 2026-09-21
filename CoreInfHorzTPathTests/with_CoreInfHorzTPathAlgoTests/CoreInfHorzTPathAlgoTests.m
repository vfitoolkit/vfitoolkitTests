% CoreInfHorzTPathAlgoTests
% Compares the transition-path algorithms against each other on four models: with and without a
% decision variable, crossed with one and two endogenous states. Each model solves the same
% null-reform general eqm path three times, from the same bumped starting guess:
%   (A) transpathoptions.GEnewprice=1, quasi-Newton on the whole price path with Broyden updates
%   (B) transpathoptions.GEnewprice=3, the shooting algorithm
%   (C) transpathoptions.GEnewprice=3 with the additionalfactor ramp
%   (F) transpathoptions.GEnewprice=2, Anderson acceleration of the same shooting map as (B)
% Plus, on the model without d only, a check that the two ways of building the FullJacobian agree.
% All of them must reach toleranceGEcondns and must agree on where the equilibrium path is, and
% Anderson at memory=0 (which is the shooting map untouched) must reproduce (B) exactly. Runtimes
% are reported throughout: the brute-force Jacobian costs (T-1)*nPrices full path solves, so (A) is
% expected to be much the slowest, while an Anderson iteration costs the same as a shooting one.

if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreInfHorzTPathAlgoTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzTPathAlgoTestsdiary.txt') % otherwise diary just appends to the previous run
end
%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): one endogenous state (figs 1-2)
% doPart(2): cross-test: stacked vs joint d_grid
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1];

diary ./TestOutput/CoreInfHorzTPathAlgoTestsdiary.txt
fprintf('CoreInfHorzTPathAlgoTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreInfHorzTPathAlgoTests_subcodes/')
addpath('../CoreInfHorzTPathTests_Setup/')
addpath('../CoreInfHorzTPath_ReturnFns/')
CoreInfHorzTPath_setup

%% ===== doPart(1): one endogenous state (figs 1-2) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): one endogenous state (figs 1-2) =====\n')
    %% One endogenous state, without d
    figure_c=1;
    output=CoreInfHorzTPathAlgo_nod_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);

    %% One endogenous state, with d
    figure_c=2;
    output=CoreInfHorzTPathAlgo_d_z_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);
end % doPart(1): one endogenous state (figs 1-2)

% %% Two endogenous states, without d
% figure_c=3;
% output=CoreInfHorzTPathAlgo_nod_z_noe_nosemiz_with2A(T,PricePath_2A,ParamPath_2A,n_d_2A,n_a_2A,n_a_2A_big,n_z_2A,d_grid_2A,a_grid_2A,a_grid_2A_big,z_grid_2A,pi_z_2A,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_2A,n_a_2A_GE,d_grid_2A,a_grid_2A_GE,figure_c);
% 
% %% Two endogenous states, with d
% figure_c=4;
% output=CoreInfHorzTPathAlgo_d_z_noe_nosemiz_with2A(T,PricePath_2A,ParamPath_2A,n_d_2A,n_a_2A,n_a_2A_big,n_z_2A,d_grid_2A,a_grid_2A,a_grid_2A_big,z_grid_2A,pi_z_2A,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_2A,n_a_2A_GE,d_grid_2A,a_grid_2A_GE,figure_c);

%% Diagnostic: is the Jacobian ill-conditioned because of the discrete d, or intrinsically?
% One Jacobian per (n_d_GE, epsprice) cell at maxiter=1, rather than four full Newton solves.
% output=CoreInfHorzTPathAlgo_jacobiandiagnostic(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c);

%% ===== doPart(2): cross-test: stacked vs joint d_grid =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-test: stacked vs joint d_grid =====\n')
    %% Cross test: a stacked-column d_grid against a joint grid with the never-chosen rows dropped
    % Its own small two-decision-variable model, so it does not depend on the four models above.
    output=CoreInfHorzTPathAlgo_crosstest_jointdgrid(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(2): cross-test: stacked vs joint d_grid

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzTPathAlgoTestsdiary.txt')

diary off
