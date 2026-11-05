% Tests of the VFI Toolkit uncertain transition path commands, FHorz (UncertainTransitionPath_FHorz and its companions;
% see VFI_Toolkit_Docs/vfitoolkitTests/UncertainTransitionPath_FHorz_proposal.md).
% The aggregate shock S: T1,...,T2 are the periods in which agents face uncertainty, in each S_{t+1} is drawn from
% pi_S(S_t,:); S_t=S_initial up to T1, frozen from T2+1 on; prices depend on the history of S, so everything lives on
% a TREE of histories (UncertainTransPathTree, shared with the InfHorz commands).
%
% This is a sub-bank of CoreFHorzTPathTests and reuses its setup and return functions. It is built in stages (F0-F4 of
% the proposal); the parts that exist so far are below.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): the tree with parameters that are vectors (age-dependent) and the per-S form ParamPath.(name).S1, .S2, ...
%            (stage F0, no solves)
% doPart(2): a single history is a deterministic transition (N_S=1, a sure switch with an age-dependent per-S
%            parameter, a 3-D pi_S, fastOLG=0 vs 1) (stage F1, as are parts 3-7)
% doPart(3): identical regimes, every history is the deterministic path: nod z (fastOLG 1 and 0), d z, nod noz,
%            nod z e, GI, DC
% doPart(4): the tree built by hand (absorbing S, one branch per switch date), cohorts aging across every branch
% doPart(5): local consistency on a full support tree with d
% doPart(6): demographics that depend on S (age weights and newborns per S): the children of a node start apart
% doPart(7): refusals, as a set, and the warning on a square (N_j=T) entry
% doPart(8): general eqm on the tree by shooting, UncertainTransitionPath_FHorz (stage F2): a sure switch and identical
%            regimes against TransitionPath_Case1_FHorz, an absorbing tree checked by hand, a restart, its refusals
% doPart(9): Anderson acceleration and price constraints on the tree (stage F3)
% doPart(10): _tplus1 and _tminus1 (in expectation over the children at a branching node), and intermediateEqns (F3)
% doPart(11): a stochastic discount factor, name_tplus1 in DiscountFactorParamNames, given prices and in general eqm (F3)
%
% Parts are independent: the setup and the addpaths sit OUTSIDE the if-blocks.
doPart=[1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzUncertainTPathTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzUncertainTPathTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzUncertainTPathTestsdiary.txt
fprintf('CoreFHorzUncertainTPathTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreFHorzTPathTests_Setup/')
addpath('../CoreFHorzTPath_ReturnFns/')
addpath('./withUncertainTransitionPath_subcodes/')

% The parent bank's setup: grids, Params, N_j, T, PricePath.r, ParamPath.w
CoreFHorzTPath_setup

% The regimes. S=1 is the parent bank's path; S=2 has a lower interest rate and wage (think of a repealed reform).
rS=[PricePath.r; PricePath.r-0.005]; % 2-by-T, a price that depends on S
wS=[ParamPath.w; 0.9*ParamPath.w];   % 2-by-T, a parameter that depends on S

%% ===== doPart(1): the tree with vector and per-S parameters =====
if doPart(1)==1
    fprintf('\n===== doPart(1): the tree with vector and per-S parameters =====\n')
    output=UncertainTPathFHorz_TreeParams(T,N_j);
end

%% ===== doPart(2): a single history is a deterministic transition =====
if doPart(2)==1
    fprintf('\n===== doPart(2): a single history is a deterministic transition =====\n')
    output=UncertainTPathFHorz_SingleHistory(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(3): identical regimes =====
if doPart(3)==1
    fprintf('\n===== doPart(3): identical regimes, every history is the deterministic path =====\n')
    output=UncertainTPathFHorz_IdenticalRegimes(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline);
end

%% ===== doPart(4): the tree built by hand =====
if doPart(4)==1
    fprintf('\n===== doPart(4): the tree built by hand =====\n')
    output=UncertainTPathFHorz_HandBuiltTree(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(5): local consistency =====
if doPart(5)==1
    fprintf('\n===== doPart(5): local consistency on a full support tree with d =====\n')
    output=UncertainTPathFHorz_LocalConsistency(T,rS,wS,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(6): demographics that depend on S =====
if doPart(6)==1
    fprintf('\n===== doPart(6): demographics that depend on S =====\n')
    output=UncertainTPathFHorz_Demographics(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(7): refusals =====
if doPart(7)==1
    fprintf('\n===== doPart(7): refusals, and the warning on a square entry =====\n')
    output=UncertainTPathFHorz_Refusals(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline);
end

%% ===== doPart(8): general eqm on the tree =====
if doPart(8)==1
    fprintf('\n===== doPart(8): general eqm on the tree, by shooting =====\n')
    output=UncertainTPathFHorz_GE(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(9): Anderson and price constraints =====
if doPart(9)==1
    fprintf('\n===== doPart(9): Anderson acceleration and price constraints on the tree =====\n')
    output=UncertainTPathFHorz_AndersonConstraints(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(10): _tplus1, _tminus1 and intermediateEqns =====
if doPart(10)==1
    fprintf('\n===== doPart(10): _tplus1, _tminus1 and intermediateEqns on the tree =====\n')
    output=UncertainTPathFHorz_tplus1IntEqns(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% ===== doPart(11): a stochastic discount factor =====
if doPart(11)==1
    fprintf('\n===== doPart(11): a stochastic discount factor on the tree =====\n')
    output=UncertainTPathFHorz_SDF(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline);
end

%% One verdict for the whole run
CoreSummary('../TestOutput/CoreFHorzUncertainTPathTestsdiary.txt')

diary off
