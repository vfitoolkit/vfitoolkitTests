% PType tests for the InfHorz TPath commands (ValueFnOnTransPath_InfHorz_PType,
% AgentDistOnTransPath_InfHorz_PType, EvalFnOnTransPath_AggVars_InfHorz_PType,
% TransitionPath_InfHorz_PType). The InfHorz counterpart of CoreFHorzTPathPTypeTests.
%
% Every check has an exact oracle: the commands without permanent types, either on the same model
% (N_i=1, identical types) or one solo solve per ptype (types that differ, and types whose general eqm
% conditions are conditional on ptype, transpathoptions.GEptype, so that they do not interact).
% The return functions and the setup are those of CoreInfHorzTPathTests.
%
%   1: N_i numeric vs Names_i cell: the standalone commands, and one general eqm iteration with
%      FnsToEvaluate given per ptype (a structure keyed by the ptype names).
%   2: N_i=2 identical types vs no-PType, all 16 combinations of d (d1), z, e and semiz: the standalone commands on a
%      fixed path, and the general eqm path (shooting, maxiter=2) for N_i=1 and N_i=2. The FnsToEvaluate
%      'earnings' reads w, which is on the ParamPath, so the in-path aggregates must see the period-t
%      parameters.
%   3: two different types (sigma differs, as a structure) vs two solo solves: the standalone commands,
%      and the general eqm conditions the solver returns against those evaluated by hand on the aggregates
%      of the standalone commands (maxiter=1, at the INPUT path: when maxiter runs out the returned price
%      path has had one more update). Then a ParamPath parameter by ptype (as a structure and as an
%      N_i-by-T matrix) and a PricePath price by ptype.
%   4: Anderson acceleration (GEnewprice=2) vs no-PType, all 16 combinations, N_i=1 and N_i=2 identical.
%   5: four decoupled types (GEptype) that differ in n_a, d, z, z_grid, e and divide-and-conquer, with an
%      interest rate per ptype given as a structure with its fields out of Names_i order, plus a common
%      general eqm condition, vs four solo solves. GEcondnPath read by name.
%   6: shocks by type: one ptype whose z_grid varies over the path and one with a constant z_grid of a
%      different size, in both orders of Names_i, vs solo solves.
%   7: _tminus1 and _tplus1 (prices, a parameter, an aggregate) with initial values by ptype, GEptype,
%      vs two solo solves; and EvalFnOnTransPath_AggVars_InfHorz_PType with the same.
%   8: constraints on the prices (constrainpositive, constrainAtoB, constraintstepfraction): N_i=1 and identical types vs
%      no-PType for shooting and Anderson (Type-II and Type-I) against a binding bound; a price by ptype under a
%      binding bound (two islands) vs two solo constrained solves; refusals by name.
%
% Anderson only mixes once it has a history (iteration 1 is a plain shooting step, and by default so is
% iteration 2, anderson.warmup=2), so every Anderson check runs anderson.warmup=0 and maxiter=3, and
% prints that Anderson then differs from shooting ('should NOT be zero'). Where the general eqm condition
% is by ptype (GEptype) Anderson is not separable by ptype, so there anderson.memory=0, which is the
% shooting algorithm, is checked against shooting.
% The algorithms themselves (shooting vs Anderson, to convergence) are compared in
% with_CoreInfHorzTPathPTypeAlgoTests.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): N_i numeric vs Names_i cell (no figures)
% doPart(2): N_i=2 identical types vs no-PType, shooting (no figures)
% doPart(3): two different types vs two solo solves (no figures)
% doPart(4): Anderson vs no-PType (no figures)
% doPart(5): GE with decoupled types (GEptype) vs solo solves (no figures)
% doPart(6): shocks by type, z_grid varying over the path (no figures)
% doPart(7): GE with _tminus1/_tplus1 by ptype vs solo solves (no figures)
% doPart(8): constraints on the prices (no figures)
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreInfHorzTPathPTypeTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzTPathPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreInfHorzTPathPTypeTestsdiary.txt
fprintf('CoreInfHorzTPathPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreInfHorzTPathPTypeTests_subcodes/')
addpath('./CoreInfHorzTPathPTypeTests_Setup/')

% Reuse the return functions and the setup of CoreInfHorzTPathTests
addpath('../CoreInfHorzTPathTests/CoreInfHorzTPath_ReturnFns/')
addpath('../CoreInfHorzTPathTests/CoreInfHorzTPathTests_Setup/')

% Setup: builds on CoreInfHorzTPath_setup and adds PType-specific pieces
CoreInfHorzTPathPType_setup

combos={'nod_noz_noe','d_noz_noe','nod_z_noe','d_z_noe','nod_noz_e','d_noz_e','nod_z_e','d_z_e'};
% the semiz combos: without d1 the model has only the d2 that drives semiz (n_d2_semiz, d2_grid_semiz), with d1 both (n_d_semiz, d_grid_semiz)
combos_semiz={'nod1_noz_noe_semiz','d1_noz_noe_semiz','nod1_z_noe_semiz','d1_z_noe_semiz','nod1_noz_e_semiz','d1_noz_e_semiz','nod1_z_e_semiz','d1_z_e_semiz'};

%% ===== doPart(1): N_i numeric vs Names_i cell =====
if doPart(1)==1
    fprintf('\n===== doPart(1): N_i numeric vs Names_i cell =====\n')
    output=CoreInfHorzTPathPType_NivsNames(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline);
end % doPart(1): N_i numeric vs Names_i cell

%% ===== doPart(2): N_i=2 identical types vs no-PType, shooting =====
if doPart(2)==1
    fprintf('\n===== doPart(2): N_i=2 identical types vs no-PType, shooting =====\n')
    for cc=1:length(combos)
        output=CoreInfHorzTPathPType_NiAllIdentical_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,transpathoptionsbaseline,combos{cc});
    end
    for cc=1:length(combos_semiz)
        if strcmp(combos_semiz{cc}(1:3),'nod')
            n_d_semiz_cc=n_d2_semiz; d_grid_semiz_cc=d2_grid_semiz;
        else
            n_d_semiz_cc=n_d_semiz; d_grid_semiz_cc=d_grid_semiz;
        end
        output=CoreInfHorzTPathPType_NiAllIdentical_vsNoPType(T,PricePath,ParamPath,n_d_semiz_cc,n_a,n_z,d_grid_semiz_cc,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,transpathoptionsbaseline,combos_semiz{cc});
    end
end % doPart(2): N_i=2 identical types vs no-PType, shooting

%% ===== doPart(3): two different types vs two solo solves =====
if doPart(3)==1
    fprintf('\n===== doPart(3): two different types vs two solo solves =====\n')
    output=CoreInfHorzTPathPType_NiTwoTypes_vsTwoSolves(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(3): two different types vs two solo solves

%% ===== doPart(4): Anderson vs no-PType =====
if doPart(4)==1
    fprintf('\n===== doPart(4): Anderson vs no-PType =====\n')
    for cc=1:length(combos)
        output=CoreInfHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,combos{cc});
    end
    for cc=1:length(combos_semiz)
        if strcmp(combos_semiz{cc}(1:3),'nod')
            n_d_semiz_cc=n_d2_semiz; d_grid_semiz_cc=d2_grid_semiz;
        else
            n_d_semiz_cc=n_d_semiz; d_grid_semiz_cc=d_grid_semiz;
        end
        output=CoreInfHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz_cc,n_a,n_z,d_grid_semiz_cc,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline_semiz,simoptionsbaseline_semiz,combos_semiz{cc});
    end
end % doPart(4): Anderson vs no-PType

%% ===== doPart(5): GE with decoupled types (GEptype) vs solo solves =====
if doPart(5)==1
    fprintf('\n===== doPart(5): GE with decoupled types (GEptype) vs solo solves =====\n')
    output=CoreInfHorzTPathPType_GE_DecoupledTypes(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline);
end % doPart(5): GE with decoupled types (GEptype) vs solo solves

%% ===== doPart(6): shocks by type, z_grid varying over the path =====
if doPart(6)==1
    fprintf('\n===== doPart(6): shocks by type, z_grid varying over the path =====\n')
    output=CoreInfHorzTPathPType_ShocksByType(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,{'zvar','zfix'});
    output=CoreInfHorzTPathPType_ShocksByType(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,{'zfix','zvar'});
end % doPart(6): shocks by type, z_grid varying over the path

%% ===== doPart(7): GE with _tminus1/_tplus1 by ptype vs solo solves =====
if doPart(7)==1
    fprintf('\n===== doPart(7): GE with _tminus1/_tplus1 by ptype vs solo solves =====\n')
    output=CoreInfHorzTPathPType_GE_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(7): GE with _tminus1/_tplus1 by ptype vs solo solves

%% ===== doPart(8): constraints on the prices =====
if doPart(8)==1
    fprintf('\n===== doPart(8): constraints on the prices =====\n')
    output=CoreInfHorzTPathPType_constraints(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end % doPart(8): constraints on the prices

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzTPathPTypeTestsdiary.txt')

diary off
