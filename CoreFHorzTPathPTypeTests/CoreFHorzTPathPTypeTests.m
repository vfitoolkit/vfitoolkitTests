% PType tests for the FHorz TPath commands (ValueFnOnTransPath_Case1_FHorz_PType,
% AgentDistOnTransPath_Case1_FHorz_PType, EvalFnOnTransPath_AggVars_Case1_FHorz_PType,
% TransitionPath_Case1_FHorz_PType).
%
% Tests 1-4 use only z (no e, no semiz) and rely on existing return functions
% from CoreFHorzTPathTests/CoreFHorzTPath_ReturnFns/.
%
% The ShockTests test at the end uses all 4 (z,e) combinations across 4 PTypes
% and must be set up via Names_i + per-type structures, because the n_z,
% z_grid, pi_z and the e pieces of vfoptions/simoptions differ across types.
% (semiz is not yet part of the TPath tests, so it is not included here either.)
%
% The non-PType TPath commands themselves (fastOLG, divide-and-conquer, grid
% interpolation, lowmemory) are tested in CoreFHorzTPathTests; these tests
% focus on the PType layer.
%
% Note: per-type N_j is not allowed for TPath (unlike the core FHorz PType
% commands), so there is no PerTypeNj test here.
%
% Parts 4-9 test the general eqm transition path with permanent types (TransitionPath_Case1_FHorz_PType),
% which solves each ptype with the same single path iteration as TransitionPath_Case1_FHorz and only
% evaluates the general eqm conditions across ptypes. Each runs the shooting algorithm for maxiter=2,
% which returns the price path after one update together with the general eqm conditions on it, and
% each has an exact oracle:
%   4-7: N_i=1, and N_i=2 identical types, against the solver without permanent types, for all 16
%        combinations of d, z, e and semiz, fastOLG 1 and 0; and with two different types, the general
%        eqm conditions the solver returns against those evaluated by hand on the aggregates of the
%        standalone PType commands (so the in-path aggregates must see the period-t parameters).
%   8:   four types that differ in shocks, z_grid, n_a, d, divide-and-conquer, age weights, survival and
%        the age-1 distribution path, with general eqm conditions by ptype (GEptype), against four solo
%        solves (so each ptype must use its own shocks on the path).
%   9:   _tminus1 and _tplus1 (prices, a parameter, an aggregate) with initial values by ptype, GEptype,
%        against two solo solves.
%   10:  age weights (mewj), survival (sj) and the age-1 distribution, constant and as paths, common to all
%        ptypes and differing by ptype, against solo solves (GEptype) and their ptweight-weighted sum; the
%        standalone PType commands with an age-weights path by ptype; and a general eqm path with bequests
%        through a demographic transition that differs by ptype (shooting vs solo, Anderson vs shooting).
% The algorithms themselves (shooting vs Anderson) are compared in with_CoreFHorzTPathPTypeAlgoTests.
%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): N_i numeric vs Names_i cell
% doPart(2): N_i all types identical vs one solve
% doPart(3): figs 1-2
% doPart(4): GE vs no-PType, without semiz (fastOLG=1)
% doPart(5): GE vs no-PType, without semiz (fastOLG=0)
% doPart(6): GE vs no-PType, with semiz (fastOLG=1)
% doPart(7): GE vs no-PType, with semiz (fastOLG=0)
% doPart(8): GE with decoupled types (GEptype) vs solo solves
% doPart(9): GE with _tminus1/_tplus1 by ptype vs solo solves
% doPart(10): age weights and demographics by ptype (no figures)
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzTPathPTypeTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzTPathPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzTPathPTypeTestsdiary.txt
fprintf('CoreFHorzTPathPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))


addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzTPathPTypeTests_subcodes/')
addpath('./CoreFHorzTPathPTypeTests_subcodes/ShockTests/')
addpath('./CoreFHorzTPathPTypeTests_Setup/')
addpath('./CoreFHorzTPathPType_ReturnFns/')

% Reuse return functions and the existing setup from CoreFHorzTPathTests
addpath('../CoreFHorzTPathTests/CoreFHorzTPath_ReturnFns/')
addpath('../CoreFHorzTPathTests/CoreFHorzTPath_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzTPathTests/CoreFHorzTPathTests_Setup/')

% Setup: builds on CoreFHorzTPath_setup and adds PType-specific pieces
CoreFHorzTPathPType_setup

%% ===== doPart(1): N_i numeric vs Names_i cell =====
if doPart(1)==1
    fprintf('\n===== doPart(1): N_i numeric vs Names_i cell =====\n')
    %% 1. N_i (numeric) vs Names_i (cell) — should give identical output
    output=CoreFHorzTPathPType_NivsNames(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline);
end % doPart(1): N_i numeric vs Names_i cell

%% ===== doPart(2): N_i all types identical vs one solve =====
if doPart(2)==1
    fprintf('\n===== doPart(2): N_i all types identical vs one solve =====\n')
    %% 2. N_i with all types identical ≡ a single no-PType solve
    output=CoreFHorzTPathPType_NiAllIdentical_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,transpathoptionsbaseline);
end % doPart(2): N_i all types identical vs one solve

%% ===== doPart(3): figs 1-2 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): figs 1-2 =====\n')
    %% 3. N_i with two distinct types ≡ two no-PType solves stacked
    figure_c=1;
    output=CoreFHorzTPathPType_NiTwoTypes_vsTwoSolves(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,transpathoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathPTypeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% 4. GE transition path with maxiter=1: identical types ≡ no-PType (fastOLG=1 and 0)
    output=CoreFHorzTPathPType_GETPath_1iter(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);

    %% 5. Two identical types: single vs per-type FnsToEvaluate — should give identical output
    % (AllStats and AgeConditionalStats on the transition path are not yet
    % implemented for PType; commented-out blocks in the subcode await them.)
    output=CoreFHorzTPathPType_PerTypeFnsToEvaluate(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline);

    %% ShockTests: one PType solve with 4 different (z,e) combos
    % Requires Names_i (per-type structures) because n_z/z_grid/pi_z and the e
    % pieces of vfoptions/simoptions differ across types.
    output=CoreFHorzTPathPType_ShockTests_4types(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);

    %% 7. General eqm by ptype (GEptype): bequests left equal bequests received, by ptype
    % Two types differ by a fixed effect fe in earnings; bequests stay within
    % type, so the PType+GEptype solve must equal two independent one-type GE
    % solves. Runs the full GE machinery (stationary endpoints + shooting
    % transition); the transition is run twice, fastOLG=1 then fastOLG=0.
    % This is the slow test of the suite.
    figure_c=2;
    output=CoreFHorzTPathPType_GEptype_Bequests(T,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzTPathPTypeTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): figs 1-2

%% ===== doPart(4): GE vs no-PType, without semiz (fastOLG=1) =====
if doPart(4)==1
    fprintf('\n===== doPart(4): GE vs no-PType, without semiz (fastOLG=1) =====\n')
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_noz_noe_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_noz_noe_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_z_noe_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_z_noe_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_noz_e_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_noz_e_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_z_e_nosemiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_z_e_nosemiz',1);
end % doPart(4): GE vs no-PType, without semiz (fastOLG=1)

%% ===== doPart(5): GE vs no-PType, without semiz (fastOLG=0) =====
if doPart(5)==1
    fprintf('\n===== doPart(5): GE vs no-PType, without semiz (fastOLG=0) =====\n')
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_noz_noe_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_noz_noe_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_z_noe_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_z_noe_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_noz_e_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_noz_e_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod_z_e_nosemiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d_z_e_nosemiz',0);
end % doPart(5): GE vs no-PType, without semiz (fastOLG=0)

%% ===== doPart(6): GE vs no-PType, with semiz (fastOLG=1) =====
% Models without d1 use n_d2_semiz and d2_grid_semiz, models with d1 use n_d_semiz and d_grid_semiz
if doPart(6)==1
    fprintf('\n===== doPart(6): GE vs no-PType, with semiz (fastOLG=1) =====\n')
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_noz_noe_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_noz_noe_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_z_noe_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_z_noe_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_noz_e_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_noz_e_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_z_e_semiz',1);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_z_e_semiz',1);
end % doPart(6): GE vs no-PType, with semiz (fastOLG=1)

%% ===== doPart(7): GE vs no-PType, with semiz (fastOLG=0) =====
if doPart(7)==1
    fprintf('\n===== doPart(7): GE vs no-PType, with semiz (fastOLG=0) =====\n')
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_noz_noe_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_noz_noe_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_z_noe_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_z_noe_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_noz_e_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_noz_e_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'nod1_z_e_semiz',0);
    output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d_semiz,n_a,n_z,N_j,d_grid_semiz,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,'d1_z_e_semiz',0);
end % doPart(7): GE vs no-PType, with semiz (fastOLG=0)

%% ===== doPart(8): GE with decoupled types (GEptype) vs solo solves =====
if doPart(8)==1
    fprintf('\n===== doPart(8): GE with decoupled types (GEptype) vs solo solves =====\n')
    output=CoreFHorzTPathPType_GE_DecoupledTypes(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,1);
    output=CoreFHorzTPathPType_GE_DecoupledTypes(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,0);
end % doPart(8): GE with decoupled types (GEptype) vs solo solves

%% ===== doPart(9): GE with _tminus1/_tplus1 by ptype vs solo solves =====
if doPart(9)==1
    fprintf('\n===== doPart(9): GE with _tminus1/_tplus1 by ptype vs solo solves =====\n')
    output=CoreFHorzTPathPType_GE_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,1);
    output=CoreFHorzTPathPType_GE_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,0);
end % doPart(9): GE with _tminus1/_tplus1 by ptype vs solo solves

%% ===== doPart(10): age weights and demographics by ptype (no figures) =====
if doPart(10)==1
    fprintf('\n===== doPart(10): age weights and demographics by ptype (no figures) =====\n')
    output=CoreFHorzTPathPType_Demographics(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,1);
    output=CoreFHorzTPathPType_Demographics(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,0);
    output=CoreFHorzTPathPType_Demographics_GEbequests(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,transpathoptionsbaseline);
end % doPart(10): age weights and demographics by ptype (no figures)

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzTPathPTypeTestsdiary.txt')

diary off
