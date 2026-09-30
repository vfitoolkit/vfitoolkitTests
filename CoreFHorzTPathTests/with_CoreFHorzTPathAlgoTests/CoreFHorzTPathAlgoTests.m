% CoreFHorzTPathAlgoTests
% Compares the FHorz transition-path algorithms against each other on two life-cycle models, with and
% without a decision variable. The FHorz counterpart of CoreInfHorzTPathAlgoTests. Each model solves the
% same null-reform general eqm path, from the same bumped starting guess, with
%   (A) transpathoptions.GEnewprice=1, quasi-Newton on the whole price path with Broyden updates, for
%       every Jacobianmethod (LudwigPath, LudwigStationary, FullJacobian, LudwigSSJ) crossed with every
%       BroydenRegularisation (minimumnorm, TikhonovRegularisation, stepcap)
%   (B) transpathoptions.GEnewprice=3, the shooting algorithm
%   (C) transpathoptions.GEnewprice=3 with the additionalfactor ramp
%   (F) transpathoptions.GEnewprice=2, Anderson acceleration of the same shooting map as (B)
% All of them must reach toleranceGEcondns and must agree on where the equilibrium path is, and Anderson
% at memory=0 (which is the shooting map untouched) must reproduce (B) exactly. Runtimes are reported.
%
% The checks particular to FHorz, or to how these algorithms were brought to it:
%   (D) the order of the howtoupdate rows does not matter
%   (H) transpathoptions.updatepert=0 and =1 build the same path
%   (I) fastOLG=0 and fastOLG=1 agree, for shooting and for the fake-news Jacobian
%   (G) the triangular FullJacobian equals the plain one (model without d only)
%   (J) the fake-news Jacobian (LudwigSSJ) against the brute-force FullJacobian (model without d only)
%   (E) quasi-Newton started halfway to the shooting solution (model with d only)
%
% The models are those of CoreFHorzTPath_nod_z_noe_nosemiz and CoreFHorzTPath_d_z_noe_nosemiz, from the
% same setup, but with a Cobb-Douglas firm supplying both r and w, and their own stationary general eqm.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): model without d (fig 1)
% doPart(2): model with d (fig 2)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathAlgoTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathAlgoTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathAlgoTestsdiary.txt
fprintf('CoreFHorzTPathAlgoTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTPathAlgoTests_subcodes/')
addpath('../CoreFHorzTPathTests_Setup/')
addpath('../CoreFHorzTPath_ReturnFns/')
% Setup so that use the same d,a,z and parameters as CoreFHorzTPathTests
CoreFHorzTPath_setup

%% ===== doPart(1): model without d (fig 1) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): model without d (fig 1) =====\n')
    figure_c=1;
    output=CoreFHorzTPathAlgo_nod_z_noe_nosemiz(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): model without d (fig 1)

%% ===== doPart(2): model with d (fig 2) =====
if doPart(2)==1
    fprintf('\n===== doPart(2): model with d (fig 2) =====\n')
    figure_c=2;
    output=CoreFHorzTPathAlgo_d_z_noe_nosemiz(T,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(2): model with d (fig 2)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTPathAlgoTestsdiary.txt')

diary off
