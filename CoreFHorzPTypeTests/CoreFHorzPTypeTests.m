% PType tests for the Core FHorz commands.
%
% Every check compares a PType solve against something known to be the same model: solves of
% each type without PType, or a model without PType in which z encodes the type. The doPart list
% below says what each part covers. Return functions come from CoreFHorzTests/CoreFHorz_ReturnFns/.
%
% Not covered: per-type solver options (divideandconquer, gridinterplayer for one type only),
% per-type n_a/a_grid, ptypestorecpu=1, N_i>=10 auto-naming.
%
% The ShockTests (part 7) use all 8 (z,e,semiz) combinations across 8 PTypes and must be set up
% via Names_i + per-type structures, because the n_z, z_grid, pi_z and vfoptions pieces differ
% across types.
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): N_i numeric vs Names_i cell
% doPart(2): N_i vs no-PType with identity z
% doPart(3): N_i all types identical vs one solve
% doPart(4): N_i two distinct types vs stacked solves
% doPart(5): jequaloneDist accepted three ways
% doPart(6): per-type N_j, different lifespans
% doPart(7): ShockTests, 8 shock combos
% doPart(8): AggVars, ValuesOnGrid, LifeCycleProfiles, PolicyInd2Val, SimPanel vs stacked solves
% doPart(9): SimPanel exactly, no z
% doPart(10): every grouped statistic, via z-identity
% doPart(11): a ptype of zero mass
% doPart(12): standard decision variable d
% doPart(13): z and e grids per type, struct vs trailing dim, age-dependent
% doPart(14): ExogShockFn and EiidShockFn per type
% doPart(15): DiscountFactorParamNames per type
% doPart(16): ValuesOnGrid with per-type n_z
% doPart(17): LifeCycleProfiles with conditional restrictions, agegroupings of one and of several ages
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzPTypeTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzPTypeTestsdiary.txt
fprintf('CoreFHorzPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzPTypeTests_subcodes/')
addpath('./CoreFHorzPTypeTests_subcodes/ShockTests/')
addpath('./CoreFHorzPTypeTests_Setup/')

% Reuse return functions and the existing setup from CoreFHorzTests
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzTests/CoreFHorzTests_Setup/')

%% Setup: builds on CoreFHorz_setup and adds PType-specific pieces
CoreFHorzPType_setup

%% ===== doPart(1): N_i numeric vs Names_i cell =====
if doPart(1)==1
    fprintf('\n===== doPart(1): N_i numeric vs Names_i cell =====\n')
    %% 1. N_i (numeric) vs Names_i (cell) — should give identical output
    output=CoreFHorzPType_NivsNames(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i);
end % doPart(1): N_i numeric vs Names_i cell

%% ===== doPart(2): N_i vs no-PType with identity z =====
if doPart(2)==1
    fprintf('\n===== doPart(2): N_i vs no-PType with identity z =====\n')
    %% 2. N_i vs no-PType with z having identity transitions (z encodes the type)
    output=CoreFHorzPType_NivsZidentity(n_d,n_a,N_j,d_grid,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);
end % doPart(2): N_i vs no-PType with identity z

%% ===== doPart(3): N_i all types identical vs one solve =====
if doPart(3)==1
    fprintf('\n===== doPart(3): N_i all types identical vs one solve =====\n')
    %% 3. N_i with all types identical ≡ a single no-PType solve
    output=CoreFHorzPType_NiAllIdentical_vsNoPType(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);
end % doPart(3): N_i all types identical vs one solve

%% ===== doPart(4): N_i two distinct types vs stacked solves =====
if doPart(4)==1
    fprintf('\n===== doPart(4): N_i two distinct types vs stacked solves =====\n')
    %% 4. N_i with two distinct types ≡ two no-PType solves stacked
    output=CoreFHorzPType_NiTwoTypes_vsTwoSolves(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(4): N_i two distinct types vs stacked solves

%% ===== doPart(5): jequaloneDist accepted three ways =====
if doPart(5)==1
    fprintf('\n===== doPart(5): jequaloneDist accepted three ways =====\n')
    %% 5. jequaloneDist accepted three ways: [n_a,n_z], [n_a,n_z,N_i], struct
    output=CoreFHorzPType_jequaloneDist_3ways(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i);
end % doPart(5): jequaloneDist accepted three ways

%% ===== doPart(6): per-type N_j, different lifespans =====
if doPart(6)==1
    fprintf('\n===== doPart(6): per-type N_j, different lifespans =====\n')
    %% 6. Per-type N_j (different lifespans) — PType as struct ≡ stacked solo solves
    output=CoreFHorzPType_PerTypeNj(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames);
end % doPart(6): per-type N_j, different lifespans

%% ===== doPart(7): ShockTests, 8 shock combos =====
if doPart(7)==1
    fprintf('\n===== doPart(7): ShockTests, 8 shock combos =====\n')
    %% ShockTests: one PType solve with 8 different (z,e,semiz) combos
    % Requires Names_i (per-type structures) because n_z/z_grid/pi_z and the e/semiz
    % pieces of vfoptions differ across types.
    output=CoreFHorzPType_ShockTests_8types(n_d,n_a,n_z,n_d_semiz,d_grid_semiz,n_d2_semiz,d2_grid_semiz,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(7): ShockTests, 8 shock combos

%% ===== doPart(8): AggVars, ValuesOnGrid, LifeCycleProfiles, PolicyInd2Val, SimPanel vs stacked solves =====
if doPart(8)==1
    fprintf('\n===== doPart(8): AggVars, ValuesOnGrid, LifeCycleProfiles, PolicyInd2Val, SimPanel vs stacked solves =====\n')
    %% 8. The PType commands downstream of StationaryDist, two distinct types vs two solves stacked
    output=CoreFHorzPType_DownstreamCmds(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(8): AggVars, ValuesOnGrid, LifeCycleProfiles, PolicyInd2Val, SimPanel vs stacked solves

%% ===== doPart(9): SimPanel exactly, no z =====
if doPart(9)==1
    fprintf('\n===== doPart(9): SimPanel exactly, no z =====\n')
    %% 9. With no shocks and a point-mass jequaloneDist the panel is deterministic, so it can be checked exactly
    output=CoreFHorzPType_SimPanelNoz(n_a,N_j,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,N_i);
end % doPart(9): SimPanel exactly, no z

%% ===== doPart(10): every grouped statistic, via z-identity =====
if doPart(10)==1
    fprintf('\n===== doPart(10): every grouped statistic, via z-identity =====\n')
    %% 10. Grouped AllStats/LifeCycleProfiles equal the no-PType stats when z encodes the type (as part 2);
    % also conditionalrestrictions, agegroupings, lowmemory, groupptypesforstats and whichstats
    output=CoreFHorzPType_ZidentityStats(n_a,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);
end % doPart(10): every grouped statistic, via z-identity

%% ===== doPart(11): a ptype of zero mass =====
if doPart(11)==1
    fprintf('\n===== doPart(11): a ptype of zero mass =====\n')
    %% 11. ptypeweight 0 on one type: grouped outputs must equal those of the other type alone
    output=CoreFHorzPType_ZeroWeight(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(11): a ptype of zero mass

%% ===== doPart(12): standard decision variable d =====
if doPart(12)==1
    fprintf('\n===== doPart(12): standard decision variable d =====\n')
    %% 12. Two distinct types with a standard d vs two solves stacked, every PType command
    output=CoreFHorzPType_WithD(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(12): standard decision variable d

%% ===== doPart(13): z and e grids per type, struct vs trailing dim, age-dependent =====
if doPart(13)==1
    fprintf('\n===== doPart(13): z and e grids per type, struct vs trailing dim, age-dependent =====\n')
    %% 13. z_grid/pi_z/e_grid/pi_e per type as struct and as a trailing N_i dimension, age-independent and age-dependent
    output=CoreFHorzPType_ShockGridForms(n_a,n_z,vfoptionsbaseline.n_e,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(13): z and e grids per type, struct vs trailing dim, age-dependent

%% ===== doPart(14): ExogShockFn and EiidShockFn per type =====
if doPart(14)==1
    fprintf('\n===== doPart(14): ExogShockFn and EiidShockFn per type =====\n')
    %% 14. One shock fn with per-type parameters, and a struct of shock fns, vs explicit grids
    output=CoreFHorzPType_ShockFns(n_a,n_z,vfoptionsbaseline.n_e,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(14): ExogShockFn and EiidShockFn per type

%% ===== doPart(15): DiscountFactorParamNames per type =====
if doPart(15)==1
    fprintf('\n===== doPart(15): DiscountFactorParamNames per type =====\n')
    %% 15. DiscountFactorParamNames as a struct keyed by Names_i
    output=CoreFHorzPType_DiscountFactorStruct(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(15): DiscountFactorParamNames per type

%% ===== doPart(16): ValuesOnGrid with per-type n_z =====
if doPart(16)==1
    fprintf('\n===== doPart(16): ValuesOnGrid with per-type n_z =====\n')
    %% 16. ValuesOnGrid when n_z differs by ptype (this used to error on prod() of a struct n_z)
    output=CoreFHorzPType_ValuesOnGridPerTypeShocks(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(16): ValuesOnGrid with per-type n_z

%% ===== doPart(17): LifeCycleProfiles with conditional restrictions, agegroupings of one and of several ages =====
if doPart(17)==1
    fprintf('\n===== doPart(17): LifeCycleProfiles with conditional restrictions, agegroupings of one and of several ages =====\n')
    %% 17. Restricted LifeCycleProfiles PType (per-ptype and grouped) vs solo Case1 solves, vs E[X*1_R]/E[1_R], and vs AllStats PType;
    % unequal ptweights and age weights, and a restriction whose mass differs by age and ptype
    output=CoreFHorzPType_RestrictedLifeCycle(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);
end % doPart(17): LifeCycleProfiles with conditional restrictions, agegroupings of one and of several ages

% All looks good!

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzPTypeTestsdiary.txt')

diary off
