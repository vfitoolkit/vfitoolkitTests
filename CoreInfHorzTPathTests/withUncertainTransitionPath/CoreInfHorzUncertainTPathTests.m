% Tests of the VFI Toolkit uncertain transition path commands, InfHorz.
% The aggregate shock S: T1,...,T2 are the periods in which agents face uncertainty, in each S_{t+1} is drawn from
% pi_S(S_t,:); S_t=S_initial up to T1, frozen from T2+1 on;
% prices depend on the history of S, so everything lives on a TREE of histories (see UncertainTransPathTree).
%
% This is a sub-bank of CoreInfHorzTPathTests and reuses its setup and return functions.
%
% THE ORACLES, strongest first
%   - a tree with a single history IS a deterministic transition: ValueFnOnTransPath_InfHorz etc. on the same path
%   - identical regimes: two S with the same parameters, full support, so every history is the deterministic path.
%     Probabilities are DYADIC (0.5, 0.25) so that the expectation over identical children is exactly the child's
%     value in floating point; with 0.3/0.7 it is not, and a ULP could flip a Policy tie and fail an identity.
%   - a tree built by hand from the deterministic commands, the way the ABCP2024 replication's stopgap does it
%   - local consistency of the tree's own answer: each branching node is one step from the expectation of its children,
%     children of one node start from one distribution, every distribution has mass one

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% doPart(1): the tree itself: counts, probabilities, numbering, segments, prices on the tree (no solves)
% doPart(2): a single history is a deterministic transition (N_S=1, and a sure switch of S, with 2-D and 3-D pi_S)
% doPart(3): identical regimes, every history is the deterministic path: nod z, d z, nod noz, nod z e, GI, DC
% doPart(4): the tree built by hand (absorbing S, one branch per switch date, as in ABCP2024)
% doPart(5): local consistency on a full support tree with d
% doPart(6): refusals, each by its message
% doPart(7): UncertainTransitionPath_InfHorz, general eqm on the tree by shooting
% doPart(8): UncertainTransitionPath_InfHorz with Anderson acceleration and with price constraints
% doPart(9): _tplus1 and _tminus1 on the tree: at a branching node a condition that uses _tplus1 holds in expectation
% doPart(10): a stochastic discount factor, name_tplus1 in DiscountFactorParamNames, read at the child
%
% Parts are independent: the setup and the addpaths sit OUTSIDE the if-blocks.
doPart=[1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreInfHorzUncertainTPathTestsdiary.txt','file')
    delete('../TestOutput/CoreInfHorzUncertainTPathTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreInfHorzUncertainTPathTestsdiary.txt
fprintf('CoreInfHorzUncertainTPathTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreInfHorzTPathTests_Setup/')
addpath('../CoreInfHorzTPath_ReturnFns/')
addpath('./withUncertainTransitionPath_subcodes/')

% The parent bank's setup: grids, Params, T=100, PricePath.r, ParamPath.w, the e process
CoreInfHorzTPath_setup

% The regimes. S=1 is the parent bank's path; S=2 has a lower wage and interest rate (think of a repealed reform).
% Used by doParts 2, 4 and 5 (doPart 3 uses identical regimes).
rS=[PricePath.r; PricePath.r-0.005]; % 2-by-T, a price that depends on S
wS=[ParamPath.w; 0.9*ParamPath.w];   % 2-by-T, a parameter that depends on S

%% ===== doPart(1): the tree itself =====
if doPart(1)==1
    fprintf('\n===== doPart(1): the tree itself =====\n')
    output=UncertainTPath_TreeStructure(T);
end

%% ===== doPart(2): a single history is a deterministic transition =====
if doPart(2)==1
    fprintf('\n===== doPart(2): a single history is a deterministic transition =====\n')
    output=UncertainTPath_SingleHistory(T,rS,wS,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(3): identical regimes =====
if doPart(3)==1
    fprintf('\n===== doPart(3): identical regimes, every history is the deterministic path =====\n')
    output=UncertainTPath_IdenticalRegimes(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline);
end

%% ===== doPart(4): the tree built by hand =====
if doPart(4)==1
    fprintf('\n===== doPart(4): the tree built by hand =====\n')
    output=UncertainTPath_HandBuiltTree(T,rS,wS,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(5): local consistency =====
if doPart(5)==1
    fprintf('\n===== doPart(5): local consistency on a full support tree with d =====\n')
    output=UncertainTPath_LocalConsistency(T,rS,wS,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(6): refusals =====
if doPart(6)==1
    fprintf('\n===== doPart(6): refusals =====\n')
    output=UncertainTPath_Refusals(T,rS,wS,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(7): general eqm on the tree =====
if doPart(7)==1
    fprintf('\n===== doPart(7): general eqm on the tree, by shooting =====\n')
    % The general eqm grid of the parent bank (its n_a=101 grid makes K a step function of the prices)
    output=UncertainTPath_GE(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(8): Anderson and price constraints =====
if doPart(8)==1
    fprintf('\n===== doPart(8): Anderson acceleration and price constraints on the tree =====\n')
    output=UncertainTPath_AndersonConstraints(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(9): _tplus1 and _tminus1 =====
if doPart(9)==1
    fprintf('\n===== doPart(9): _tplus1 and _tminus1 on the tree =====\n')
    output=UncertainTPath_tplus1(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% ===== doPart(10): a stochastic discount factor =====
if doPart(10)==1
    fprintf('\n===== doPart(10): a stochastic discount factor, name_tplus1 in DiscountFactorParamNames =====\n')
    output=UncertainTPath_SDF(T,rS,wS,n_a,n_z,a_grid,z_grid,pi_z,n_a_GE,a_grid_GE,Params,DiscountFactorParamNames,transpathoptionsbaseline);
end

%% One verdict for the whole run
CoreSummary('../TestOutput/CoreInfHorzUncertainTPathTestsdiary.txt')

diary off
