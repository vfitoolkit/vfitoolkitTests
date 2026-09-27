% Implement tests of the core VFI Toolkit InfHorz commands, with ENTRY AND EXIT.
% with/without d
% with/without a semi-endogenous shock
% endogenousexit=1 (pure endogenous exit) and endogenousexit=2 (a three-way exit mixture)
% including general equilibrium with endogenous entry
%
% This is an entry-exit sub-bank of CoreInfHorzTests. Unlike the quasi-hyperbolic sub-bank it
% does NOT reuse the parent's setup or return functions: entry and exit is a firm problem
% (a is previous employment, exit means shutting down and paying firing costs), so it brings its
% own. The model shape follows Hopenhayn-Rogerson (1993).
%
% GPU ONLY. Of the eleven entry-exit raws, six are the parallel=0 and parallel=1 CPU paths and
% are deliberately out of scope. The guards cross-test probes two of them so the diary at least
% records their status. The five covered here are:
%   fig 1  EndogExit_nod_Par2_raw
%   fig 2  EndogExit_Par2_raw
%   fig 3  EndogExit_SemiEndog_nod_Par2_raw
%   fig 4  EndogExit2_nod_Par2_raw and EndogExit2_Par2_raw
%   fig 5  all of the above again, many times over, through the general eqm root-finder
%
% WHY THIS BANK EXISTS
% Entry and exit had no automated coverage of any kind. It also contains ~32 sites of a defect
% swept out of the rest of the toolkit on 2026-09-10 (issue #111): a zero weight multiplying a
% -Inf value gives NaN, which the next iteration's isnan-to-zero then launders into a
% continuation value of zero. Here the weight is the exit decision
%     ExitPolicy=((ReturnToExitMatrix-Vtemp)>0);
%     VKron=ExitPolicy.*ReturnToExitMatrix+(1-ExitPolicy).*Vtemp;
% which is a bare logical, so it is ALWAYS exactly 0 or 1 -- worse than the interpolation sites
% swept elsewhere, where extreme weights were only occasional. The model in this bank is built
% to have both a region where staying is infeasible and a region where exit is forbidden, so
% both branches of that weight meet a -Inf. See withEntryExit_Setup/CoreInfHorzEntryExit_setup.m.
%
% This bank is intended to be run BEFORE those sites are edited, so that the diary is a genuine
% pre-fix baseline. Expect the NaN census lines to be non-zero on that first run.
%
% NOT tested here:
%   - the CPU raws (parallel=0, parallel=1): out of scope, see above
%   - grid interpolation and divide-and-conquer: no entry-exit path supports either
%   - lowmemory: the dispatcher rejects it outright (probed in the guards cross-test)
%   - transition paths with entry-exit: a separate bank's business


%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): endogenousexit figs 1-4   - OFF by default, see below
% doPart(2): general equilibrium, entry and exit (fig 5)
% doPart(3): interface guards and out-of-scope probes
%
% doPart(1) IS OFF, and that is the pre-existing state of this bank rather than a new decision.
% Figs 1-4 were all green on 2026-09-16 (0 NaN everywhere, every cross-test exact), so they were
% skipped while fig 5's general equilibrium is being diagnosed. That was done here with a local
% dofigs1to4 flag, written for the same reason doPart exists ("a flag rather than commented-out
% lines, so restoring is one character and cannot be half-done"); the flag has been folded into
% doPart so the file has only one mechanism. Set doPart(1) to 1 to restore them.
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[0,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreInfHorzEntryExitTestsdiary.txt','file')
    delete('../TestOutput/CoreInfHorzEntryExitTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreInfHorzEntryExitTestsdiary.txt
fprintf('CoreInfHorzEntryExitTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

%%
addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./withEntryExit_Setup/')
addpath('./withEntryExit_ReturnFns/')
addpath('./withEntryExit_subcodes/')
addpath('./withEntryExit_subcodes/CrossTests/')
% Setup so that use the same d,a,z in all the models
CoreInfHorzEntryExit_setup

%% ===== doPart(1): endogenousexit figs 1-4 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): endogenousexit figs 1-4 =====\n')
    %% endogenousexit=1, without d
    figure_c=1;
    output=EEInfHorz_nod(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEntryExitTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% endogenousexit=1, with d
    figure_c=2;
    output=EEInfHorz_d(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEntryExitTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% endogenousexit=1, semi-endogenous shock, without d
    figure_c=3;
    output=EEInfHorz_semiendog(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEntryExitTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% endogenousexit=2, with and without d
    figure_c=4;
    output=EEInfHorz_exit2(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEntryExitTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): endogenousexit figs 1-4

%% ===== doPart(2): general equilibrium, entry and exit (fig 5) =====
if doPart(2)==1
    fprintf('\n===== doPart(2): general equilibrium, entry and exit (fig 5) =====\n')
    %% General equilibrium with endogenous entry and endogenous exit
    figure_c=5;
    output=EEInfHorz_generaleqm(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreInfHorzEntryExitTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(2): general equilibrium, entry and exit (fig 5)

%% ===== doPart(3): interface guards and out-of-scope probes =====
if doPart(3)==1
    fprintf('\n===== doPart(3): interface guards and out-of-scope probes =====\n')
    %% Interface guards and out-of-scope probes (no figure)
    output=EEInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(3): interface guards and out-of-scope probes

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreInfHorzEntryExitTestsdiary.txt')

diary off
