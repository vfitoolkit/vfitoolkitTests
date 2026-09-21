% Implement core tests of the VFI Toolkit InfHorz commands with an INHERITANCE ASSET.
%
% An inheritance asset is one whose next-period value is determined by a decision d2 together
% with the realized shock transition (z,zprime), and not by this period's asset. So the return
% function has no aprime input (Case-2 style), and vfoptions.aprimeFn takes (d2,z,zprime).
%
% WHAT THE TOOLKIT SUPPORTS, AND THEREFORE WHAT THIS BANK COVERS
% ValueFnIter_InfHorz_InheritAsset has exactly ONE live branch: d1 present, no a1, z present,
% no e. The other eight combinations are 'Have not yet implemented' errors. There is no tier
% cube to sweep, so this is a small bank whose value comes from instruments rather than from
% combinations:
%   fig 1  the supported tier, plus an independent recomputation of V from Policy
%   fig 2  cross-test: a de-risked inheritance asset must equal a plain one-asset model
%          (no figure) shape guards: the eight unsupported combinations must still error
%
% NOT tested here, because the toolkit does not have it:
%   - divide-and-conquer: not usable for InfHorz at all
%   - grid interpolation layer: ValueFnIter_InfHorz_InheritAsset has no gridinterplayer branch
%   - Howards: the Howards block in ValueFnIter_InfHorz_InheritAsset_noa1_raw is commented out,
%     so vfoptions.howards has no effect and a howards on/off check would be vacuous
%   - ValueFnFromPolicy: ValueFnFromPolicy_InfHorz has no inheritance-asset (or aprimeFn)
%     handling at all, so the bank recomputes V from Policy itself in fig 1
%   - FHorz: there is no FHorz inheritance asset
%
% WHY THIS BANK EXISTS
% The InfHorz inheritance asset had no automated coverage of any kind. It was written as a
% sibling of the experience-asset family and carries the same shapes, so it inherits their bugs
% without ever having been exercised. One such bug was found and fixed on 2026-09-10 (e0aff864:
% a2primeProbs was zeroed in place on an array built outside the iteration loop, so the answer
% depended on the starting guess); the V0-independence check in fig 1 is its regression test.


%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): supported tier: d1 noa1 z noe (fig 1)
% doPart(2): cross-test vs plain one-asset model (fig 2)
% doPart(3): shape guards
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
if exist('./TestOutput/CoreInfHorzInheritAssetTestsdiary.txt','file')
    delete('./TestOutput/CoreInfHorzInheritAssetTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreInfHorzInheritAssetTestsdiary.txt
fprintf('CoreInfHorzInheritAssetTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreInfHorzInheritAssetTests_subcodes/')
addpath('./CoreInfHorzInheritAssetTests_subcodes/CrossTests/')
addpath('./CoreInfHorzInheritAssetTests_Setup/')
addpath('./CoreInfHorzInheritAsset_ReturnFns/')
% Setup so that use the same d,a,z in all the models
CoreInfHorzInheritAsset_setup

%% ===== doPart(1): supported tier: d1 noa1 z noe (fig 1) =====
if doPart(1)==1
    fprintf('\n===== doPart(1): supported tier: d1 noa1 z noe (fig 1) =====\n')
    %% The supported tier: d1, no a1, z, no e
    figure_c=1;
    output=CoreInfHorzInheritAsset_d1_noa1_z_noe(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzInheritAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): supported tier: d1 noa1 z noe (fig 1)

%% ===== doPart(2): cross-test vs plain one-asset model (fig 2) =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-test vs plain one-asset model (fig 2) =====\n')
    %% Cross-test: de-risked inheritance asset vs plain one-asset model
    figure_c=2;
    output=CoreInfHorzInheritAsset_CrossTests_vsplain(n_d,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreInfHorzInheritAssetTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(2): cross-test vs plain one-asset model (fig 2)

%% ===== doPart(3): shape guards =====
if doPart(3)==1
    fprintf('\n===== doPart(3): shape guards =====\n')
    %% Shape guards (no figure)
    output=CoreInfHorzInheritAsset_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline);
end % doPart(3): shape guards

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreInfHorzInheritAssetTestsdiary.txt')

diary off
