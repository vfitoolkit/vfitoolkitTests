% PType tests for the InfHorz commands (the InfHorz counterpart of CoreFHorzPTypeTests; the parts loosely follow its parts).
%
% Every check compares a PType solve against something known to be the same model: separate single-type solves of each type, a
% model without PType where the type is a second exogenous state with the identity as its transition (z-identity), or a brute
% force from the definitions of the statistics. The doPart list below says what each part covers. Return functions and the
% shared setup (grids, parameters) come from CoreInfHorzTests.
%
% Not here, compared with CoreFHorzPTypeTests: there is no jequaloneDist, no N_j (so no per-type lifespans), no age-conditional
% command, and no iid e (ValueFnIter_InfHorz does not support e yet).
%
% No figures are drawn anywhere in this bank, so only the diary is saved.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): N_i numeric vs Names_i cell
% doPart(2): every grouped statistic of AllStats and AggVars, via z-identity
% doPart(3): N_i all types identical vs one solve
% doPart(4): two distinct types vs separate single-type solves, without and with d (V, Policy, StationaryDist, PolicyInd2Val, AggVars, AllStats)
% doPart(5): a ptype of zero mass
% doPart(6): the exogenous shock per type: z_grid/pi_z per type, n_z per type, ExogShockFn with a per-type parameter, ExogShockFn per type
% doPart(7): the discount factor per type: a per-type parameter, DiscountFactorParamNames per type
% doPart(8): no exogenous shock, FnsToEvaluate with parameters
% doPart(9): conditional restrictions in AllStats PType, every statistic per type and grouped vs brute force
% doPart(10): ValuesOnGrid PType (new command)
% doPart(11): CrossSectionCovarCorr PType (new command): per type vs single-type, grouped vs pooled, with conditional restrictions
% doPart(12): AutoCorrTransProbs PType (new command): per type vs single-type, grouped vs pooled, with timehorizons and conditional restrictions
% doPart(13): SimPanelValues PType (new command)
% doPart(14): grouped CrossSection and AutoCorr PType, via z-identity
% doPart(15): simoptions.whichcombos (and per-combination whichstats) in AllStats, AggVars, AutoCorr and CrossSection PType, with the type dimension
% Parts 10-14 test the four commands new 2026-10-07 (written test first).
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreInfHorzPTypeTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzPTypeTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreInfHorzPTypeTestsdiary.txt
fprintf('CoreInfHorzPTypeTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, CoreStructLeaves, shared by the Core banks
addpath('./CoreInfHorzPTypeTests_subcodes/')
addpath('./CoreInfHorzPTypeTests_Setup/')

% Reuse return functions and the existing setup from CoreInfHorzTests
addpath('../CoreInfHorzTests/CoreInfHorz_ReturnFns/')
addpath('../CoreInfHorzTests/CoreInfHorzTests_Setup/')

%% Setup: the same d, a, z and parameters as CoreInfHorzTests, plus the PType pieces
CoreInfHorz_setup
CoreInfHorzPType_setup

%% ===== doPart(1): N_i numeric vs Names_i cell =====
if doPart(1)==1
    fprintf('\n===== doPart(1): N_i numeric vs Names_i cell =====\n')
    %% 1. Names_i as a number or as a cell of the same names: identical V, Policy, StationaryDist, AggVars, AllStats
    output=CoreInfHorzPType_NamesForms(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(1): N_i numeric vs Names_i cell

%% ===== doPart(2): every grouped statistic of AllStats and AggVars, via z-identity =====
if doPart(2)==1
    fprintf('\n===== doPart(2): every grouped statistic of AllStats and AggVars, via z-identity =====\n')
    %% 2. The PType model is the single-type model with the type as a second exogenous state (identity transition): grouped
    % AllStats (every statistic) and AggVars equal its single-type stats; each type's AllStats equals it restricted to the type
    output=CoreInfHorzPType_Zidentity(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(2): every grouped statistic of AllStats and AggVars, via z-identity

%% ===== doPart(3): N_i all types identical vs one solve =====
if doPart(3)==1
    fprintf('\n===== doPart(3): N_i all types identical vs one solve =====\n')
    %% 3. Two identical types: each equals the single-type solve, and the grouped AllStats equals its AllStats, every leaf
    output=CoreInfHorzPType_TwoTypes(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT,1);
end % doPart(3): N_i all types identical vs one solve

%% ===== doPart(4): two distinct types vs separate single-type solves =====
if doPart(4)==1
    fprintf('\n===== doPart(4): two distinct types vs separate single-type solves =====\n')
    %% 4. Two types (a per-type wage), without and with d: per type V, Policy, StationaryDist, PolicyInd2Val, AggVars and every leaf
    % of AllStats equal the single-type solves; grouped AggVars Mean and AllStats Mean/StdDeviation/Variance equal their pooling
    output=CoreInfHorzPType_TwoTypes(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT,0);
end % doPart(4): two distinct types vs separate single-type solves

%% ===== doPart(5): a ptype of zero mass =====
if doPart(5)==1
    fprintf('\n===== doPart(5): a ptype of zero mass =====\n')
    %% 5. ptweight 0 on one type: grouped AllStats (every leaf) and AggVars equal those of the other type alone
    output=CoreInfHorzPType_ZeroWeight(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(5): a ptype of zero mass

%% ===== doPart(6): the exogenous shock per type =====
if doPart(6)==1
    fprintf('\n===== doPart(6): the exogenous shock per type =====\n')
    %% 6. z_grid/pi_z per type, n_z per type (different sizes), one ExogShockFn with a per-type parameter, ExogShockFn per type:
    % per type vs single-type solves, grouped Mean/StdDeviation vs pooled
    output=CoreInfHorzPType_ShocksPerType(n_a,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(6): the exogenous shock per type

%% ===== doPart(7): the discount factor per type =====
if doPart(7)==1
    fprintf('\n===== doPart(7): the discount factor per type =====\n')
    %% 7. A per-type discount factor parameter, and DiscountFactorParamNames per type: per type vs single-type solves
    output=CoreInfHorzPType_DiscountFactor(n_a,n_z,a_grid,z_grid,pi_z,Params,PTypeDistParamNames,simoptionsPT);
end % doPart(7): the discount factor per type

%% ===== doPart(8): no exogenous shock =====
if doPart(8)==1
    fprintf('\n===== doPart(8): no exogenous shock =====\n')
    %% 8. n_z=0, without and with d, FnsToEvaluate with a per-type and a common parameter: per type vs single-type solves
    output=CoreInfHorzPType_NoZ(n_d,n_a,d_grid,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(8): no exogenous shock

%% ===== doPart(9): conditional restrictions in AllStats PType =====
if doPart(9)==1
    fprintf('\n===== doPart(9): conditional restrictions in AllStats PType =====\n')
    %% 9. EvalFnOnAgentDist_AllStats_InfHorz_PType with simoptions.conditionalrestrictions: every statistic, per type and grouped,
    % every page vs a brute force; restricted sample masses; per type vs the single-type command; the always page (two models, a
    % per-type wage in the return function, a FnsToEvaluate and a restriction; a zero-weight type)
    output=CoreInfHorzPType_CondlRest(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames);
end % doPart(9): conditional restrictions in AllStats PType

%% ===== doPart(10): ValuesOnGrid PType =====
if doPart(10)==1
    fprintf('\n===== doPart(10): ValuesOnGrid PType =====\n')
    %% 10. EvalFnOnAgentDist_ValuesOnGrid_InfHorz_PType: per type vs single-type, with d and z, and with n_z per type
    output=CoreInfHorzPType_ValuesOnGrid(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(10): ValuesOnGrid PType

%% ===== doPart(11): CrossSectionCovarCorr PType =====
if doPart(11)==1
    fprintf('\n===== doPart(11): CrossSectionCovarCorr PType =====\n')
    %% 11. EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType: per type vs single-type (every leaf, every restriction), grouped
    % vs the pooled single-type outputs, restricted sample masses, groupptypesforstats=0; two models, one with a zero-weight type
    output=CoreInfHorzPType_CrossSection(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(11): CrossSectionCovarCorr PType

%% ===== doPart(12): AutoCorrTransProbs PType =====
if doPart(12)==1
    fprintf('\n===== doPart(12): AutoCorrTransProbs PType =====\n')
    %% 12. EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType: per type vs single-type (every leaf, every restriction, TransitionProbs
    % included), grouped (unrestricted and restricted pair outputs, horizons 1,2,3) vs the pooled single-type outputs, restricted
    % sample masses, groupptypesforstats=0; two models, one with a zero-weight type
    output=CoreInfHorzPType_AutoCorr(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(12): AutoCorrTransProbs PType

%% ===== doPart(13): SimPanelValues PType =====
if doPart(13)==1
    fprintf('\n===== doPart(13): SimPanelValues PType =====\n')
    %% 13. SimPanelValues_InfHorz_PType: sizes, no NaN, simulations per type (none of a zero-weight type), the InitialDist as one
    % distribution or per type; loosely, the panel means vs AllStats
    output=CoreInfHorzPType_SimPanel(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(13): SimPanelValues PType

%% ===== doPart(14): grouped CrossSection and AutoCorr PType, via z-identity =====
if doPart(14)==1
    fprintf('\n===== doPart(14): grouped CrossSection and AutoCorr PType, via z-identity =====\n')
    %% 14. The grouped (and per-type, and restricted) outputs of the two new PType correlation commands equal the single-type
    % commands on the (z, type) model of part 2
    output=CoreInfHorzPType_ZidentityCorr(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(14): grouped CrossSection and AutoCorr PType, via z-identity

%% ===== doPart(15): whichcombos =====
if doPart(15)==1
    fprintf('\n===== doPart(15): whichcombos =====\n')
    %% 15. simoptions.whichcombos (new 2026-10-07) in EvalFnOnAgentDist_AllStats_InfHorz_PType (with per-combination whichstats),
    % EvalFnOnAgentDist_AggVars_InfHorz_PType, EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType and
    % EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType: a selective run vs the full run leaf by leaf, per type and grouped,
    % including a type-specific selection (the trailing type dimension; a grouped slot forces every type); broadcast inputs vs their
    % expansions; invalid inputs error. The InfHorz counterpart of CoreFHorzPTypeTests part 21.
    output=CoreInfHorzPType_WhichCombos(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT);
end % doPart(15): whichcombos

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzPTypeTestsdiary.txt')

diary off
