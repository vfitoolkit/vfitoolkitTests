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
%
% Parts are independent: the setup and the addpaths sit OUTSIDE the if-blocks.
doPart=[1];

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

% The parent bank's setup: grids, Params, N_j, T
CoreFHorzTPath_setup

%% ===== doPart(1): the tree with vector and per-S parameters =====
if doPart(1)==1
    fprintf('\n===== doPart(1): the tree with vector and per-S parameters =====\n')
    output=UncertainTPathFHorz_TreeParams(T,N_j);
end

%% One verdict for the whole run
CoreSummary('../TestOutput/CoreFHorzUncertainTPathTestsdiary.txt')

diary off
