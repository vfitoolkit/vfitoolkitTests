% Tests of the VFI Toolkit discretization commands.
% Organised by the process being discretized; see description.txt
%
% Currently built: P0 to P7.

%% Which parts to run
% One entry per part, P0 first: doPart(1) is P0, doPart(2) is P1, and so on up to doPart(11) for
% P10. Set an entry to zero to skip that part. That is for building - while a later part is being
% written there is no reason to rerun the earlier ones that already pass, and P5 at M=3 alone takes
% over a minute.
%
% Skipping a part does NOT renumber the figures. Each part owns a fixed block of ten figure
% numbers - P0 gets 1 to 10, P1 gets 11 to 20, and so on - so Fig31 is P3's first figure whatever
% doPart says, and a figure saved by a previous full run is never overwritten by a different one.
% The gaps that leaves in the numbering are deliberate; no part currently uses more than six.
%
% Parts are independent: each one runs its own setup, and no part reads another part's output. So
% any subset can be run, in any combination.
%
% ALL ON. Every part is green, but only in pieces: P2/P4/P7 on one run, P1/P3/P5 on another,
% P6/P8/P9 on a third. This is the first run of the whole bank since P8 and P9 were built, and the
% first since the discretization commands were taken off the Statistics Toolbox - so it is the run
% that produces a single number for the bank rather than a set of numbers that have to be argued
% about. doPart(11) is set too; P10 was never scoped and the block below says so rather than
% silently doing nothing.
%
% WHAT TO EXPECT. P3 carries two standing reds - its worst-over-the-sweep variance and skewness
% bars, both dominated by the coarsest grid in the sweep, where the paired finest-grid checks pass
% comfortably. Everything else passed on its most recent run. Anything beyond those two is new.
%
% WHAT THIS RUN COVERS THAT NO EARLIER ONE DID. The erfc rewrite of every normcdf/normpdf/norminv
% site, the paren strip that followed it, B30 and B31 in the life-cycle gaussian-mixture commands,
% B32 and B33 in discretizeLifeCycleVAR1_Tauchen, the new e_grid option on both IIDNormal commands,
% and the new Tauchen_q=[] default on discretizeIIDNormal_Tauchen. Those have each been run, but
% never all together, and never alongside the parts that call them indirectly.
doPart=[1,1,1,1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into TestOutput as they are created)
if ~exist('./TestOutput','dir')
    mkdir('./TestOutput')
end
if exist('./TestOutput/DiscretizationMethodTestsdiary.txt','file')
    delete('./TestOutput/DiscretizationMethodTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ./TestOutput/DiscretizationMethodTestsdiary.txt

addpath('./DiscretizationMethodTests_Setup/')
addpath('./DiscretizationMethodTests_subcodes/P0_Instruments/')
addpath('./DiscretizationMethodTests_subcodes/P1_IIDNormal/')
addpath('./DiscretizationMethodTests_subcodes/P2_AR1Normal/')
addpath('./DiscretizationMethodTests_subcodes/P3_AR1GM/')
addpath('./DiscretizationMethodTests_subcodes/P4_AR1SV/')
addpath('./DiscretizationMethodTests_subcodes/P5_VAR1/')
addpath('./DiscretizationMethodTests_subcodes/P6_LCAR1/')
addpath('./DiscretizationMethodTests_subcodes/P7_LCAR1GM/')
addpath('./DiscretizationMethodTests_subcodes/P8_LCVAR1/')
addpath('./DiscretizationMethodTests_subcodes/P9_NormalOnGrid/')

output=struct();

%% ===================================================================================
%% P0: validating the instruments
%% Everything later in this bank measures discretizations by running them through
%% MarkovChainMoments(). So P0 comes first: closed-form checks that are seconds to run,
%% then the simulation-based validation that actually licenses the rest of the bank.
%% ===================================================================================
if doPart(1)==1
    DiscSetup_P0

    figure_c=1; % P0 owns figure numbers 1 to 10, and every part owns a fixed block of ten
    output.P0.MarkovChainMoments       = DiscP0_MarkovChainMoments(calibP0,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P0.MarkovChainMoments_FHorz = DiscP0_MarkovChainMoments_FHorz(calibP0,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P0.vs_panel       = DiscP0_MarkovChainMoments_vs_panel(calibP0,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P0.vs_panel_FHorz = DiscP0_MarkovChainMoments_FHorz_vs_panel(calibP0,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    DiscP0_runtimetable(output.P0);
end % doPart(1), which is P0

%% ===================================================================================
%% P1: iid normal
%% Native: discretizeIIDNormal_Tauchen, discretizeIIDNormal_TanakaToda (these own the
%% option sweep). Also run: the two discretizeIID_* originals at defaults (they are
%% currently exact copies, and are to be generalized to other distributions later), and
%% the AR(1) methods at rho=0, which discretize this same process.
%% ===================================================================================
if doPart(2)==1
    DiscSetup_IIDNormal

    figure_c=11; % P1 owns figure numbers 11 to 20

    output.P1.IIDNormal_Tauchen    = DiscP1_IIDNormal_Tauchen(calibIID,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P1.IIDNormal_TanakaToda = DiscP1_IIDNormal_TanakaToda(calibIID,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P1.IID_Tauchen          = DiscP1_IID_Tauchen(calibIID,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P1.IID_TanakaToda       = DiscP1_IID_TanakaToda(calibIID,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P1.AR1methods_rho0      = DiscP1_AR1methods_rho0(calibIID,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P1.compare              = DiscP1_compare(calibIID,znums,output.P1,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P1.crosstests           = DiscP1_crosstests(calibIID,znums);
    output.P1.downstream           = DiscP1_downstream(calibIID);
    DiscP1_runtimetable(output.P1);
end % doPart(2), which is P1

%% ===================================================================================
%% P2: AR(1) with normal innovations
%% Native: discretizeAR1_Tauchen, _Rouwenhorst, _TauchenHussey, _FarmerToda (these own
%% the option sweeps). Three calibrations: moderate persistence, high persistence, and
%% one with mew~=0, which is what pins the grid-centring convention. Plus Floden (2008)
%% Table 1 as a published external oracle, and the reductions from the VAR(1),
%% gaussian-mixture and stochastic-volatility commands to this same process.
%% ===================================================================================
if doPart(3)==1
    DiscSetup_AR1Normal

    figure_c=21; % P2 owns figure numbers 21 to 30

    output.P2.Tauchen        = DiscP2_AR1_Tauchen(calibAR1,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.Rouwenhorst    = DiscP2_AR1_Rouwenhorst(calibAR1,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.TauchenHussey  = DiscP2_AR1_TauchenHussey(calibAR1,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.FarmerToda     = DiscP2_AR1_FarmerToda(calibAR1,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.compare        = DiscP2_compare(calibAR1,znums,output.P2,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.floden         = DiscP2_floden(calibAR1,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.gridwidth        = DiscP2_gridwidth(calibAR1,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P2.crosstests     = DiscP2_crosstests(calibAR1,znums);
    output.P2.downstream     = DiscP2_downstream(calibAR1);
    DiscP2_runtimetable(output.P2,calibAR1,znums);
end % doPart(3), which is P2

%% ===================================================================================
%% P3: AR(1) with gaussian-mixture innovations
%% Native: discretizeAR1wGM_FarmerToda, and discretizeAR1wGM_Tauchen which DOES NOT
%% EXIST YET - this block is written test-first, and it will ERROR at that subcode and
%% take the rest of the run with it until the command is written. That is deliberate.
%% The gaussian AR(1) methods are also run, fitted to the mixture's first two moments,
%% as misspecified comparators: that comparison is the point of the block.
%% ===================================================================================
if doPart(4)==1
    DiscSetup_AR1GM

    figure_c=31; % P3 owns figure numbers 31 to 40

    output.P3.AR1wGM_FarmerToda = DiscP3_AR1wGM_FarmerToda(calibGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P3.AR1wGM_Tauchen    = DiscP3_AR1wGM_Tauchen(calibGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P3.gaussianmethods   = DiscP3_gaussianmethods(calibGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P3.compare           = DiscP3_compare(calibGM,znums,output.P3,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P3.gridwidth        = DiscP3_gridwidth(calibGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P3.crosstests        = DiscP3_crosstests(calibGM,znums);
    output.P3.downstream        = DiscP3_downstream(calibGM);
    DiscP3_runtimetable(output.P3,calibGM,znums);
end % doPart(4), which is P3

%% ===================================================================================
%% P4: AR(1) with log-AR(1) stochastic volatility
%% The bank's first TWO-DIMENSIONAL command: the output is a stacked grid [x_grid;z_grid]
%% with a joint transition matrix on (x,z), x varying fastest. So the sweep is nested over
%% xnum and znum independently - with equal sizes an error in the stacking is invisible.
%% Excess kurtosis is the moment stochastic volatility exists to produce, and it has an
%% exact closed form here (a convergent series), so it is the block's main accuracy target.
%% ===================================================================================
if doPart(5)==1
    DiscSetup_AR1SV

    figure_c=41; % P4 owns figure numbers 41 to 50

    output.P4.AR1wSV_FarmerToda = DiscP4_AR1wSV_FarmerToda(calibSV,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P4.AR1wSV_Tauchen    = DiscP4_AR1wSV_Tauchen(calibSV,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P4.comparator        = DiscP4_comparator(calibSV,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P4.compare           = DiscP4_compare(calibSV,znums,output.P4,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P4.gridwidth        = DiscP4_gridwidth(calibSV,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P4.crosstests        = DiscP4_crosstests(calibSV,znums);
    output.P4.downstream        = DiscP4_downstream(calibSV);
    DiscP4_runtimetable(output.P4,calibSV,znums);
end % doPart(5), which is P4

%% ===================================================================================
%% P5: VAR(1) with normally distributed innovations
%% The first block whose two commands disagree on the OUTPUT CONVENTION: VAR1_Tauchen
%% returns a stacked sum(znum)-by-1 grid indexed variable-1-fastest, VAR1_FarmerToda a
%% joint (znum^M)-by-M grid indexed last-variable-fastest. Both conventions are
%% established directly in the cross-tests rather than assumed, and the downstream solve
%% uses a calibration whose two variables have means of opposite sign so that a swap
%% cannot hide. The moment that matters is the CROSS-covariance: it is the only one a
%% pair of independent AR(1)s cannot reach, and it is why these commands exist.
%% ===================================================================================
if doPart(6)==1
    DiscSetup_VAR1

    figure_c=51; % P5 owns figure numbers 51 to 60

    output.P5.VAR1_Tauchen      = DiscP5_VAR1_Tauchen(calibVAR,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P5.VAR1_FarmerToda   = DiscP5_VAR1_FarmerToda(calibVAR,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P5.comparator        = DiscP5_comparator(calibVAR,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P5.compare           = DiscP5_compare(calibVAR,znums,output.P5,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P5.gridwidth        = DiscP5_gridwidth(calibVAR,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P5.crosstests        = DiscP5_crosstests(calibVAR,znums);
    output.P5.downstream        = DiscP5_downstream(calibVAR);
    DiscP5_runtimetable(output.P5,calibVAR,znums);
end % doPart(6), which is P5

%% ===================================================================================
%% P6: life-cycle (age-dependent) AR(1) with normal innovations
%% The first block whose every result is an AGE PROFILE, and the first to test two
%% outputs the bank has never touched: jequaloneDistz and otheroutputs.sigma_z. It is
%% also the first to check the J-1 shape of pi_z_J. Its oracle is the frozen-life-cycle
%% identity - constant parameters plus a stationary initial condition means the
%% age-dependent method must reproduce the stationary method it extends, exactly, at
%% every age - which is an oracle that is not an approximation of anything.
%% ===================================================================================
if doPart(7)==1
    DiscSetup_LCAR1

    figure_c=61; % P6 owns figure numbers 61 to 70

    output.P6.LCAR1_KFTT        = DiscP6_LCAR1_KFTT(calibLC,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P6.LCAR1_FGP         = DiscP6_LCAR1_FGP(calibLC,znums,znumsEven,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P6.LCAR1_FGPTauchen  = DiscP6_LCAR1_FGPTauchen(calibLC,znums,znumsEven,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P6.comparator        = DiscP6_comparator(calibLC,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P6.compare           = DiscP6_compare(calibLC,znums,output.P6,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P6.gridwidth        = DiscP6_gridwidth(calibLC,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P6.initialconditions = DiscP6_initialconditions(calibLC);
    output.P6.crosstests        = DiscP6_crosstests(calibLC,znums);
    output.P6.downstream        = DiscP6_downstream(calibLC);
    DiscP6_runtimetable(output.P6,calibLC,znums);
end % doPart(7), which is P6

%% ===================================================================================
%% P7: life-cycle (age-dependent) AR(1) with gaussian-mixture innovations
%% The block with an exact age profile for all FOUR moments, from the cumulant recursion
%% kn(z_j)=rho(j)^n*kn(z_{j-1})+kn(e_j). Skewness and excess kurtosis are the reason
%% gaussian mixtures are used, so they are what this block measures; mean and variance
%% are hygiene. Before P7 discretizeLifeCycleAR1wGM_KFTT appeared in exactly one subcode
%% of the whole bank, and only for the B27 regression, so its own behaviour had never
%% been measured. discretizeLifeCycleAR1wGM_Tauchen was written alongside this block.
%% The calibration is deliberately NOT mean zero, which is what lets the cross-tests see
%% the two grid-centring conventions in this family diverge - P3 could not.
%% ===================================================================================
if doPart(8)==1
    DiscSetup_LCAR1GM

    figure_c=71; % P7 owns figure numbers 71 to 80

    output.P7.LCAR1wGM_KFTT     = DiscP7_LCAR1wGM_KFTT(calibLCGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P7.LCAR1wGM_Tauchen  = DiscP7_LCAR1wGM_Tauchen(calibLCGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P7.comparator        = DiscP7_comparator(calibLCGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P7.compare           = DiscP7_compare(calibLCGM,znums,output.P7,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P7.gridwidth         = DiscP7_gridwidth(calibLCGM,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P7.zeromean          = DiscP7_zeromean(calibLCGM);
    output.P7.initialconditions = DiscP7_initialconditions(calibLCGM);
    output.P7.crosstests        = DiscP7_crosstests(calibLCGM,znums);
    output.P7.downstream        = DiscP7_downstream(calibLCGM);
    DiscP7_runtimetable(output.P7,calibLCGM,znums);
end % doPart(8), which is P7

%% ===================================================================================
%% P8: LIFE-CYCLE VAR(1), NORMAL INNOVATIONS
%% Z(j) = Mew(:,j) + Rho(:,:,j)*Z(j-1) + e(j), e(j) ~ N(0,SigmaSq(:,:,j))
%%
%% ONE COMMAND, and until this block it was the last user-facing command in
%% DiscretizationMethods with no coverage whatsoever - the string
%% discretizeLifeCycleVAR1_Tauchen did not appear anywhere in the test tree. It is also
%% the command that carried B7, B11 and B14, all three found by reading rather than by
%% running, so nothing in it had ever executed under a check.
%%
%% With no sibling method there is no accuracy race, so the weight falls on the analytic
%% age profiles from the matrix recursion and on the cross-tests, where
%% discretizeVAR1_Tauchen and discretizeLifeCycleAR1_FellaGallipoliPanTauchen act as
%% independently written implementations of special cases this command must reproduce.
%% The diagonal-calibration cross-test is the important one: P5 found the two stationary
%% VAR commands order their joint index OPPOSITELY, and a life-cycle command that picked
%% the wrong convention would still return a valid stochastic matrix with correct
%% marginals. Only reconstructing the joint from its parts catches that.
%% ===================================================================================
if doPart(9)==1
    DiscSetup_LCVAR1

    figure_c=81; % P8 owns figure numbers 81 to 90

    output.P8.LCVAR1_Tauchen = DiscP8_LCVAR1_Tauchen(calibLCV,znums,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P8.gridwidth      = DiscP8_gridwidth(calibLCV,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    DiscP8_crosstests(calibLCV,9);
    output.P8.downstream     = DiscP8_downstream(calibLCV);
end % doPart(9), which is P8

%% ===================================================================================
%% P9: A NORMAL DISTRIBUTION ON A GIVEN GRID
%%
%% NOT A PROCESS CLASS. There is no dynamics here - no transition matrix, no
%% persistence, no ages - just how to allocate probabilities to grid points you were
%% handed rather than got to choose. It sits at the end of the process blocks because
%% that is the reading order, not because it is one of them.
%%
%% WHAT MAKES IT A BLOCK. Three commands can now be asked the same question on the same
%% grid: MVNormal_ProbabilitiesOnGrid natively, and the two IIDNormal commands through an
%% e_grid option ADDED FOR THIS BLOCK. Tauchen and MVNormal integrate the density over
%% each bin, so they must agree to machine precision - that is an identity. Tanaka-Toda
%% instead solves a maximum entropy problem matching moments on the grid, so it is
%% deliberately different and should win on the moments and lose on the density. Neither
%% is right; it depends on what the grid is for.
%%
%% TWO THINGS HERE HAD NEVER RUN. MVNormal_ProbabilitiesOnGrid was called from exactly
%% two places in the whole bank before this, both at M=1, both guarding B12 - everything
%% multivariate about it was unexercised. And B13 widened its guard from l_z>=5 to l_z>5,
%% making the l_z==5 branch reachable for the first time.
%% ===================================================================================
if doPart(10)==1
    DiscSetup_NormalOnGrid

    figure_c=91; % P9 owns figure numbers 91 to 100

    output.P9.MVNormal = DiscP9_MVNormal(calibNG,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;

    output.P9.e_grid   = DiscP9_e_grid(calibNG,figure_c);
    exportgraphics(figure(figure_c),['./TestOutput/DiscretizationMethodTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
    figure_c=figure_c+1;
end % doPart(10), which is P9

%% P10 was never scoped. The proposal's block table stops at P9, and doPart has an eleventh entry
%% only so the vector does not have to change if something is ever added. It currently does nothing.
if doPart(11)==1
    fprintf('\nNote: doPart(11) refers to P10, which was never scoped, so it did nothing. \n')
end

%% One verdict for the whole run
% The bank prints its results in four different formats and there are around 2600 of them, so this
% reads the diary back and says plainly whether the run passed. See DiscSummary.m for why that is
% not something a reader should be left to do by scanning.
DiscSummary('./TestOutput/DiscretizationMethodTestsdiary.txt')

%%
diary off
