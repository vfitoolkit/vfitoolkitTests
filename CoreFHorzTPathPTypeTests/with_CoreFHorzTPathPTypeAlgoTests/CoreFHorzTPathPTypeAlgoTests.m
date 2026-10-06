% CoreFHorzTPathPTypeAlgoTests
% Compares the FHorz transition-path algorithms with permanent types against each other. The PType
% counterpart of CoreFHorzTPathAlgoTests. With permanent types there are two algorithms:
%   transpathoptions.GEnewprice=3, the shooting algorithm (with and without the additionalfactor ramp)
%   transpathoptions.GEnewprice=2, Anderson acceleration of the same shooting map
% (GEnewprice=1, quasi-Newton, is not implemented with permanent types.)
%
% doPart(1) and (2) are the two models of CoreFHorzTPathAlgoTests (without and with a decision
% variable), each with two permanent types that differ in kappa_j and share one Cobb-Douglas firm.
% Each solves the same null-reform general eqm path from the same bumped starting guess with shooting,
% shooting+ramp and five Anderson configurations. Every one must reach toleranceGEcondns (printed as a
% flag, 'should be one', so a run that stops short fails) and they must agree on the path. Exact checks:
%   (B) Anderson at memory=0 is the shooting algorithm
%   (C) N_i=2 identical types take exactly the iterates of the solver without permanent types, for
%       shooting and for Anderson
%   (D) the order of the howtoupdate rows does not matter
%   (I) fastOLG=0 and fastOLG=1 agree, for shooting and for Anderson
% doPart(3) is two 'islands', each ptype with its own firm, r and w, and general eqm conditions by
% ptype (transpathoptions.GEptype): shooting is separable by island, so must equal two solo solves
% exactly at a fixed number of iterations; Anderson at memory=0 must equal shooting; and both must
% converge to the solo eqm paths.
%
% The models use the setup of CoreFHorzTPathTests.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): model without d (fig 1)
% doPart(2): model with d (fig 2)
% doPart(3): GEptype islands (fig 3)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathPTypeAlgoTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathPTypeAlgoTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathPTypeAlgoTestsdiary.txt
fprintf('CoreFHorzTPathPTypeAlgoTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTPathPTypeAlgoTests_subcodes/')
addpath('../../CoreFHorzTPathTests/CoreFHorzTPathTests_Setup/')
addpath('../../CoreFHorzTPathTests/CoreFHorzTPath_ReturnFns/')
% Setup so that use the same d,a,z and parameters as CoreFHorzTPathTests
CoreFHorzTPath_setup

%% ===== doPart(1): model without d (fig 1) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): model without d (fig 1) =====\n')
    figure_c=1;
    output=CoreFHorzTPathPTypeAlgo_nod_z_noe_nosemiz(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathPTypeAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): model without d (fig 1)

%% ===== doPart(2): model with d (fig 2) =====
if doPart(2)==1
    fprintf('\n===== doPart(2): model with d (fig 2) =====\n')
    figure_c=2;
    output=CoreFHorzTPathPTypeAlgo_d_z_noe_nosemiz(T,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathPTypeAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(2): model with d (fig 2)

%% ===== doPart(3): GEptype islands (fig 3) =====
if doPart(3)==1
    fprintf('\n===== doPart(3): GEptype islands (fig 3) =====\n')
    figure_c=3;
    output=CoreFHorzTPathPTypeAlgo_GEptype_islands(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathPTypeAlgoTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): GEptype islands (fig 3)

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTPathPTypeAlgoTestsdiary.txt')

diary off
