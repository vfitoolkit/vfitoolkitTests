% Implement tests of the core VFI Toolkit FHorz with ExpAssetsemiz commands
% experienceassetsemiz: aprime depends on (d2,a2,semiz), so semiz is always present
% (semiz drives the experience asset). semiz is the semi-exogenous state, so the
% semiz decision d3 is also always present.
% with/without d1
% with/without ordinary z
% with/without e
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate)
%
% Tier ordering follows CoreFHorzExpAssetTests and CoreFHorzExpAssetzTests:
%   noa1 (figs 1-8):    the experienceassetsemiz a2 is the only endogenous state
%   withA1 (figs 9-16): a1 (standard endogenous state) alongside the experienceassetsemiz
%   with2A1 (figs 17-24): two standard endogenous assets, a = [a1_1, a1_2 (binary), a2]

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): noa1 figs 1-8
% doPart(2): noa1 cross-tests
% doPart(3): withA1 figs 9-16
% doPart(4): withA1 cross-tests
% doPart(5): With TWO standard endogenous figs 17-24
% doPart(6): With TWO standard endogenous cross-tests
% doPart(7): with2A2 noa1 figs 25-32          (TWO experience assets, no standard asset)
% doPart(8): with2A2 WITH a1 figs 33-40       (a1 alongside the two experience assets)
% doPart(9): with2A2 cross-tests
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
%
% NOTE: fig 24 stays commented out INSIDE part 5 - it out-of-memories on this machine. doPart is for
% whole tiers; an individual call that cannot run stays commented.
%
% This run: ALL SIX PARTS, as the l_a2==1 REGRESSION GUARD for the first stage of the with2A2
% (two-experience-asset) work on this family. 19 toolkit files changed on 2026-10-14 and every one
% of them is on the ONE-experience-asset path this bank exercises:
%   both aprimeFn builders (CreateExperienceAssetsemizFnMatrix 186->348,
%     CreateaprimePolicyExperienceAssetsemiz 180->309)
%   5 aprimeFnParamNames arity sites (the (l_a2>=2) whicha term)
%   72 a2_grid -> a2_gridvals argument switches across the 8 noa1 raws
%   all 4 downstream: StationaryDist_FHorz_ExpAssetsemiz, SimPanelIndexes_FHorz_ExpAssetsemiz,
%     ValueFnFromPolicy_FHorz_ExpAssetsemiz{,_GI}
% Every change is meant to be inert at l_a2==1, and most are provably so (the arity terms add
% (l_a2>=2)=0; CreateGridvals(n,grid,1) returns the grid itself for a 1-D grid, verified; the new
% code sits inside l_a2==2 arms). The exception worth running for: StationaryDist and
% SimPanelIndexes had their l_a2==1 path RESTRUCTURED, not just wrapped - an l_a==1/2/3 index chain
% replaced by an l_a1-based corner loop - so the projection gate cannot certify those two and only
% a run can. The with2A2 tier itself does NOT exist yet (the 40 raws are still to port), so there
% is nothing new to test here; this run is purely about not having broken what already worked.
%
% Baseline to compare against, from the 2026-10-02 run of parts 5-6: 327 checks, and the ONE
% reported failure was 'Divide-and-conquer (DC2A), this should be zero: 1.000e+00' at part 5 - the
% known benign exact-tie artefact (Policy off by one index while V sits at the ULP floor), closed
% 2026-09-19. Expect that same single failure and nothing else. Parts 1-4 have not run since at
% least 2026-10-02, so their counts are the figure to record for next time.
% NOTE (2026-10-14): on this machine parts 1-4 and parts 5-6 must be run as TWO separate MATLAB
% invocations. All six in one process killed MATLAB partway through part 5 (no error message,
% after figs 17-23 and 258 of its 323 checks); parts 5-6 alone then completed with 327 checks and
% statistics identical to the 2026-10-02 run, so it is a process memory ceiling, not a defect.
% Part 5 is the heaviest tier here - two standard endogenous assets PLUS the experience asset,
% three endogenous states - and its fig 24 is already commented out for the same reason.
% So: run this as-is for parts 1-4, then set doPart=[0,0,0,0,1,1,0,0,0] and run it again.
% doPart must have ONE ENTRY PER PART - the driver reads doPart(1) to doPart(9), so a short
% vector runs the parts it covers and then dies on an index out of bounds at the first part it
% does not, AFTER all that work and BEFORE CoreSummary, leaving a diary with no summary.
%
% This run: the l_a2==1 regression guard for with2A2 STAGE 2 - the l_a2==2 interpolation arms in
% the 40 in-scope ValueFnIter raws (228 sites, plus 80 wrapped aprime-construction regions).
% Unlike stage 1, every one of those 40 files projects BYTE-IDENTICAL at l_a2==1: the port only
% wraps, it restructures nothing. So this run is cheap insurance rather than the only evidence -
% but it is worth doing separately from the with2A2 bank tier that comes next, because a failure
% now is unambiguously a stage-2 regression, whereas a failure after the tier lands could be
% either. Baseline: parts 1-4 gave 804 checks 0 failed, parts 5-6 gave 327 checks with the single
% known-benign DC2A tie, on the 2026-10-14 stage-1 run.
% Parts 1-4 came back 566 checks 0 failed on the stage-2 guard (16:45), matching stage 1 exactly;
% largest identity 1.863e-09, largest "close to zero" 2.002e-02. Now the second invocation, for
% parts 5-6: part 5 is the tier with TWO standard endogenous assets, i.e. the with-a1 raws at
% l_a1==2, which are among the 40 files stage 2 touched, so it is the part most worth confirming.
%
% This run: parts 7-9 again, the with2A2 TIER, after fixing the two failure classes the
% 2026-10-16 09:35 run turned up (570 checks, 10 failed; part 7 noa1 was already 154/0).
%   (a) 6 x the non-vacuity check in part 8, at 1.4e-03 and 3.4e-03. Not a defect: the withA1
%       tier has a small infeasible corner (a1=0 and a2_1=0), the comparison beside it still ran
%       over every finite entry and printed 0.000e+00, and a -Inf facing a finite value would
%       have shown as Inf. The check was simply written as an identity when what it must reject
%       is V being -Inf essentially EVERYWHERE (the noa1 pathology, which printed 1.000e+00).
%       Re-worded to 'close to zero' in all 12 subcodes that carry it, so closetol=0.1 separates
%       a 1e-3 corner from a 1.0 pathology with two orders of margin either side.
%   (b) 4 x cross test 7 (inert a2_1), V Inf and Policy 1.000e+02. A defect in the TEST: the
%       ordinary with2A2 return function reads earnings off a2_1 and gives a2_2 only the pension,
%       so making a2_1 the inert dimension made V depend on the inert asset and be flat in the
%       live one - the two models were not comparable. Every other family in this programme
%       carries a dedicated _2A2swap return function for exactly this test; this family was the
%       one that never got it ported. Added ReturnFn_ExpAssetsemiz_CrossTest_2A2swap_{nod1,d1}
%       (earnings off a2_2) and wired cross test 7 to it. Cross tests 6 and 8 passed clean and
%       are untouched - 6 passing while 7 failed is what localised this.
% Expect 24 checks in part 9, not 22: two of cross test 7's checks printed 'Inf', and until the
% 2026-10-16 CoreSummary patch a check whose value was Inf or NaN was SKIPPED rather than failed,
% so the part's own count was short by two and the loudest two failures were invisible.
% Parts 7-9 came back 574 checks 0 failed (10:19): 154 noa1, 394 withA1, 26 cross-tests, largest
% identity 4.657e-10 (ValueFnFromPolicy) and largest 'close to zero' 1.635e-02. Part 9 is 26 and
% not the 24 predicted, because the 'reduces to one experience asset' line carries TWO values
% (V and Policy): cross tests 6 and 7 are 5 checks per variant, cross test 8 is 3. So the old
% count of 22 was short by FOUR, not two - the two 'flat: Inf' lines went dark entirely, and in
% each 'V Inf, Policy 1.000e+02' line the V value was dropped while the Policy value still
% counted. The CoreSummary patch recovered exactly those four.
%
% This run: parts 5-6, the OUTSTANDING half of the l_a2==1 regression guard for stage 2. Parts 1-4
% were confirmed green against the current toolkit state (566 checks 0 failed, 2026-10-15 16:45;
% all 44 changed toolkit files date from 09:50 that morning and none has been touched since -
% everything changed on 2026-10-16 was test-bank or CoreSummary, not toolkit). Parts 5-6 under
% stage 2 were never run: the note above says 'now the second invocation' and stops there. They
% have to be a separate MATLAB invocation from parts 1-4 on this machine (see the 2026-10-14 NOTE).
% Baseline: 327 checks with the single known-benign DC2A policy tie, from the 2026-10-02 run.
% Part 5 is the tier with TWO standard endogenous assets, i.e. the with-a1 raws at l_a1==2, which
% are among the 40 files stage 2 touched, so it is the part most worth confirming. fig 24 stays
% commented out inside part 5 (it out-of-memories on this machine).
doPart=[0,0,0,0,1,1,0,0,0];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt','file')
    delete('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt
fprintf('CoreFHorzExpAssetsemizTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('./CoreFHorzExpAssetsemizTests_subcodes/')
addpath('./CoreFHorzExpAssetsemizTests_Setup/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/')
addpath('./CoreFHorzExpAssetsemizTests_subcodes/CrossTests/')


%% Setup so that use the same d,a,semiz,z,e in all the models that use them
CoreFHorzExpAssetsemiz_setup


%% ================= noa1 (figs 1-8): the experience asset a2 is the ONLY endogenous state =================
% TEST-FIRST: these error when run, as the ExpAssetsemiz noa1 raws do not exist in the toolkit yet
% (ValueFnIter_FHorz_ExpAssetsemiz errors on N_a1==0).
addpath('./CoreFHorzExpAssetsemizTests_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/Noa1_ReturnFns/')

%% ===== doPart(1): noa1 figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): noa1 figs 1-8 =====\n')
    %% without d1, without z, without e, noa1
    figure_c=1;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e, noa1
    figure_c=2;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e, noa1
    figure_c=3;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e, noa1
    figure_c=4;
    output=CoreFHorzExpAssetsemiz_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e, noa1
    figure_c=5;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e, noa1
    figure_c=6;
    output=CoreFHorzExpAssetsemiz_d1_z_noe_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e, noa1
    figure_c=7;
    output=CoreFHorzExpAssetsemiz_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e, noa1
    figure_c=8;
    output=CoreFHorzExpAssetsemiz_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): noa1 figs 1-8

%% ===== doPart(2): noa1 cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): noa1 cross-tests =====\n')
    %% noa1 cross-tests (TEST-FIRST: error when run, same reason as above)
    % CrossTest 4: noa1 vs withA1 model where a1 (n_a1=1) is ignored (should match bit-exact)
    output=CoreFHorzExpAssetsemiz_CrossTests4_nod1_noz_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests4_d1_noz_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,0,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests4_nod1_z_e_noa1(n_d_withoutd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests4_d1_z_e_noa1(n_d_withd1,n_a_justexpasset,n_a_justexpasset,n_z,N_j,d_grid_withd1,a_grid_justexpasset,a_grid_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): noa1 cross-tests


%% ================= withA1 (figs 9-16): a1 alongside the experienceassetsemiz =================

%% ===== doPart(3): withA1 figs 9-16 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): withA1 figs 9-16 =====\n')
    %% without d1, without z, without e
    figure_c=9;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e
    n_a_notsobig=[301,13]; % To avoid out-of-memory errors
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3; % to test Grid Interpolation (same grid, just more points)
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=10;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=11;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e(n_d_withoutd1,n_a,n_a_notsobig,0,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e
    n_a_notsobig=[201,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=12;
    output=CoreFHorzExpAssetsemiz_d1_noz_e(n_d_withd1,n_a,n_a_notsobig,0,N_j,d_grid_withd1,a_grid,a_grid_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=13;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=14;
    output=CoreFHorzExpAssetsemiz_d1_z_noe(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e
    n_a_notsobig=[301,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=15;
    output=CoreFHorzExpAssetsemiz_nod1_z_e(n_d_withoutd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e
    n_a_notsobig=[201,13];
    a1_grid_notsobig=5*linspace(0,1,n_a_notsobig(1))'.^3;
    a_grid_notsobig=[a1_grid_notsobig;a2_grid];

    figure_c=16;
    output=CoreFHorzExpAssetsemiz_d1_z_e(n_d_withd1,n_a,n_a_notsobig,n_z,N_j,d_grid_withd1,a_grid,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): withA1 figs 9-16


%% ===== doPart(4): withA1 cross-tests =====
if doPart(4)==1
    fprintf('\n===== doPart(4): withA1 cross-tests =====\n')
    %% withA1 cross-tests
    % CrossTest 1: experienceassetsemiz with a degenerate (d3-invariant) semiz vs experienceassetz (should match)
    output=CoreFHorzExpAssetsemiz_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % CrossTest 2: 'fake' experienceassetsemiz that ignores semiz vs plain experienceasset+semiz (should match)
    output=CoreFHorzExpAssetsemiz_CrossTests2_nod1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests2_d1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    % CrossTest 3: experienceassetsemiz+z vs experienceassetz with combined z=[semiz,z] (pins bothz ordering)
    output=CoreFHorzExpAssetsemiz_CrossTests3_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests3_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): withA1 cross-tests


%% ================= With TWO standard endogenous assets (figs 17-24) =================
% a = [a1_1, a1_2 (binary), a2 (experienceassetsemiz)]
% TEST-FIRST: these error when run, as the ExpAssetsemiz DC2A/GI2A/DC2A_GI2A raws do not
% exist in the toolkit yet.
addpath('./CoreFHorzExpAssetsemizTests_subcodes/With2A1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/With2A1_ReturnFns/')

% a1main is kept modest here because the binary second asset doubles the a-grid.
% n_a_2A1=[a1, a1_2, a2] and a_grid_2A1 come from CoreFHorzExpAssetsemiz_setup; a1_2 is a genuine
% multi-point second standard asset, so the subcodes take n_a/a_grid as given and build nothing.
n_a_2A1_notsobig=[151,n_a1_2,n_a_justexpasset];
a1_grid_2A1_notsobig=5*linspace(0,1,n_a_2A1_notsobig(1))'.^3;
a_grid_2A1_notsobig=[a1_grid_2A1_notsobig;a1_2_grid;a2_grid];

%% ===== doPart(5): With TWO standard endogenous figs 17-24 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): With TWO standard endogenous figs 17-24 =====\n')
    %% with2A1, without d1, without z, without e
    figure_c=17;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, without z, without e
    figure_c=18;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, without d1, without z, with e
    figure_c=19;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, without z, with e
    figure_c=20;
    output=CoreFHorzExpAssetsemiz_d1_noz_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,0,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, without d1, with z, without e
    figure_c=21;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, with z, without e
    figure_c=22;
    output=CoreFHorzExpAssetsemiz_d1_z_noe_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, without d1, with z, with e
    figure_c=23;
    output=CoreFHorzExpAssetsemiz_nod1_z_e_with2A1(n_d_withoutd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A1, with d1, with z, with e
    figure_c=24;
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
    % output=CoreFHorzExpAssetsemiz_d1_z_e_with2A1(n_d_withd1,n_a_2A1,n_a_2A1_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A1,a_grid_2A1_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    % exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    % I CANNOT RUN THIS AS IT JUST OUT-OF-MEMORY ERRORS
end % doPart(5): With TWO standard endogenous figs 17-24

%% ===== doPart(6): With TWO standard endogenous cross-tests =====
if doPart(6)==1
    fprintf('\n===== doPart(6): With TWO standard endogenous cross-tests =====\n')
    %% with2A1 cross-tests (TEST-FIRST: error when run, same reason as above)
    % CrossTest 5: a degenerate second standard asset a1_2 (single point {0}) reduces the
    % two-standard-asset model back to the with-a1 model (i.e. DC2A/GI2A vs the ordinary solvers)
    output=CoreFHorzExpAssetsemiz_CrossTests5_nod1_with2A1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests5_d1_with2A1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(6): With TWO standard endogenous cross-tests


%% ================= with2A2: TWO EXPERIENCE ASSETS (figs 25-40) =================
% vfoptions.experienceassetsemiz is the integer COUNT of a2 dimensions (2 here), not a flag.
% Grid layout a = [a1, a2_1, a2_2], and the noa1 tier is a = [a2_1, a2_2]. BOTH experience assets
% are driven by the semi-exogenous state. The grids come from the setup; aprimeFn is
% vfoptionsbaseline.aprimeFn_2A2, which takes 'whicha' AFTER semiz:
% aprimeFn(d2,a2_1,a2_2,semiz,whicha,params...). Each subcode sets experienceassetsemiz=2 itself.
addpath('./CoreFHorzExpAssetsemizTests_subcodes/With2A2_subcodes/')
addpath('./CoreFHorzExpAssetsemizTests_subcodes/With2A2_subcodes/Noa1_subcodes/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/With2A2_ReturnFns/')
addpath('./CoreFHorzExpAssetsemiz_ReturnFns/With2A2_ReturnFns/Noa1_ReturnFns/')

% The grid-interpolation comparison needs the 'notsobig' a1 grid to be strictly FINER in a1 than
% n_a_2A2 (=[51,n_a2_1,n_a2_2]), otherwise the GI-vs-dense check is vacuous. Unlike the
% experienceassetU family this one carries no u dimension, so the panel array is N_u times
% cheaper at the same grid: the heaviest variant (z+e) at 151 is ~203 MB per array, under the
% ~237 MB that ran on this machine for ExpAssetU. If a variant out-of-memories, trim ITS a1 grid
% here rather than the a2 grids -- the a2 dims are what this tier exists to exercise.
n_a_2A2_notsobig=[151,n_a2_1,n_a2_2];
a_grid_2A2_notsobig=[5*linspace(0,1,n_a_2A2_notsobig(1))'.^3;a2_1_grid;a2_2_grid];

%% ===== doPart(7): with2A2 noa1 figs 25-32 =====
if doPart(7)==1
    fprintf('\n===== doPart(7): with2A2 noa1 figs 25-32 =====\n')
    %% with2A2 noa1, nod1, noz, noe
    figure_c=25;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,0,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, d1, noz, noe
    figure_c=26;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,0,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, nod1, noz, e
    figure_c=27;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,0,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, d1, noz, e
    figure_c=28;
    output=CoreFHorzExpAssetsemiz_d1_noz_e_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,0,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, nod1, z, noe
    figure_c=29;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, d1, z, noe
    figure_c=30;
    output=CoreFHorzExpAssetsemiz_d1_z_noe_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, nod1, z, e
    figure_c=31;
    output=CoreFHorzExpAssetsemiz_nod1_z_e_noa1_with2A2(n_d_withoutd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withoutd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 noa1, d1, z, e
    figure_c=32;
    output=CoreFHorzExpAssetsemiz_d1_z_e_noa1_with2A2(n_d_withd1,n_a_2A2_justexpasset,n_a_2A2_justexpasset,n_z,N_j,d_grid_withd1,a_grid_2A2_justexpasset,a_grid_2A2_justexpasset,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(7): with2A2 noa1 figs 25-32

%% ===== doPart(8): with2A2 WITH a1 figs 33-40 =====
if doPart(8)==1
    fprintf('\n===== doPart(8): with2A2 WITH a1 figs 33-40 =====\n')
    %% with2A2 with a1, nod1, noz, noe
    figure_c=33;
    output=CoreFHorzExpAssetsemiz_nod1_noz_noe_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, d1, noz, noe
    figure_c=34;
    output=CoreFHorzExpAssetsemiz_d1_noz_noe_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,0,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, nod1, noz, e
    figure_c=35;
    output=CoreFHorzExpAssetsemiz_nod1_noz_e_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,0,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, d1, noz, e
    figure_c=36;
    output=CoreFHorzExpAssetsemiz_d1_noz_e_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,0,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, nod1, z, noe
    figure_c=37;
    output=CoreFHorzExpAssetsemiz_nod1_z_noe_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, d1, z, noe
    figure_c=38;
    output=CoreFHorzExpAssetsemiz_d1_z_noe_with2A2(n_d_withd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, nod1, z, e
    figure_c=39;
    output=CoreFHorzExpAssetsemiz_nod1_z_e_with2A2(n_d_withoutd1,n_a_2A2,n_a_2A2_notsobig,n_z,N_j,d_grid_withoutd1,a_grid_2A2,a_grid_2A2_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with2A2 with a1, d1, z, e
    % This is the heaviest of the 16: d1 and z and e on top of semiz and BOTH
    % experience assets. On the 2026-10-15 run it ran the GPU out of memory inside
    % the VFI (entireRHS_d3 in ValueFnIter_FHorz_ExpAssetsemiz_e_raw), so it gets its
    % own smaller MAIN a1 grid. Trim a1, never the a2 grids -- the a2 dims are what
    % this tier exists to exercise.
    n_a_2A2_heavy=[21,n_a2_1,n_a2_2];
    a_grid_2A2_heavy=[5*linspace(0,1,n_a_2A2_heavy(1))'.^3;a2_1_grid;a2_2_grid];
    n_a_2A2_heavy_notsobig=[61,n_a2_1,n_a2_2];
    a_grid_2A2_heavy_notsobig=[5*linspace(0,1,n_a_2A2_heavy_notsobig(1))'.^3;a2_1_grid;a2_2_grid];
    figure_c=40;
    output=CoreFHorzExpAssetsemiz_d1_z_e_with2A2(n_d_withd1,n_a_2A2_heavy,n_a_2A2_heavy_notsobig,n_z,N_j,d_grid_withd1,a_grid_2A2_heavy,a_grid_2A2_heavy_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/CoreFHorzExpAssetsemizTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

end % doPart(8): with2A2 WITH a1 figs 33-40

%% ===== doPart(9): with2A2 cross-tests =====
if doPart(9)==1
    fprintf('\n===== doPart(9): with2A2 cross-tests =====\n')
    % CrossTests 6 and 7 make each experience asset inert in turn and require EXACT reduction to
    % the one-experience-asset model. Run together they pin down the ORDERING of the two a2 dims:
    % either alone would still pass if the nested interpolation had them transposed.
    % CrossTests 8 keeps both assets alive with identical grids and laws and a symmetric return
    % function, so it is the only test that reaches the region where all four corners of the a2
    % lattice carry weight.
    % NOTE the grids: these take the ONE-experience-asset n_a/a_grid and build the SECOND
    % a2 dimension themselves (an inert one, or a symmetric twin). Handing them n_a_2A2
    % would make Model A a 3-dim-n_a model with a 2-dim n_a grid -- matching the committed
    % ExpAssetU cross-tests, which are passed n_a/n_a_big for the same reason.
    output=CoreFHorzExpAssetsemiz_CrossTests6_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests6_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests7_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests7_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests8_nod1_with2A2(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
    output=CoreFHorzExpAssetsemiz_CrossTests8_d1_with2A2(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(9): with2A2 cross-tests

%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('./TestOutput/CoreFHorzExpAssetsemizTestsdiary.txt')

diary off
