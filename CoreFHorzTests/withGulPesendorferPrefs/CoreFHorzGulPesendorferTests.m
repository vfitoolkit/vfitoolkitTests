% Implement lots of tests of the core VFI Toolkit FHorz commands, under GUL-PESENDORFER
% preferences (temptation and self-control: a temptation fn v alongside the return fn u, with
%   V_j = max_{d,a'} [ u + v + beta*E V_{j+1} ] - max_{d,a'} v
% so Policy maximizes the tempted objective and V nets off the self-control cost)
% with/without d
% with/without z and/or e (Gul-Pesendorfer is fine with no shocks at all, so all 8 combos)
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate; the z&e models also get lowmemory=2)
%
% This is the Gul-Pesendorfer mirror of CoreFHorzTests.m. The temptation fn is declared as
% vfoptions.temptationFn, with the same input signature as the ReturnFn (its parameters may
% differ; here they deliberately do: the temptation is over consumption only). Baseline
% temptation: v = lambdaGP*u_c(c) + shiftGP with lambdaGP=0.1, shiftGP=0 (a la Krusell,
% Kuruscu & Smith); shiftGP exists purely so the cross tests can exercise the
% constant-cancellation and shift-invariance properties.
%
% STATE (2026-09-07, see GulPesendorfer_FHorz_proposal.md in the toolkit repo): the whole
% bank is GPU-green. Figs 1-24 + the four earlier cross-test files under toolkit e9d6d8aa /
% bank 3681a51 (746 zero-checks), and the semiz-with2A section (figs 25-32 + its cross
% tests; 266 checks) under the GP semiz-2A wave (24 raws: SemiExo_DC2A/GI2A/DC2A_GI2A,
% plus 2A branches in the GP SemiExo dispatchers and FromPolicy support). Memory note
% throughout: GP holds the return matrix and its temptation twin simultaneously (~2x
% core), hence the reduced big-grids on the largest cases.

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): figs 1-8
% doPart(2): cross-tests
% doPart(3): figs 9-16
% doPart(4): cross-tests 2
% doPart(5): figs 17-24
% doPart(6): cross-tests 3
% doPart(7): figs 25-32
% doPart(8): cross-tests 4
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1,1,1,1,1,1,1];

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt
fprintf('CoreFHorzGulPesendorferTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreFHorzTests_Setup/')
addpath('../CoreFHorz_ReturnFns/')

addpath('./withGulPesendorferPrefs_subcodes/')
addpath('./withGulPesendorferPrefs_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

%% Gul-Pesendorfer: the temptation parameters
% The temptation fns live in withGulPesendorferPrefs_subcodes/ (GPTemptationFn_*), one per
% case signature; each subcode declares its own vfoptions.temptationFn (mirroring how the
% ReturnFn is declared per-subcode)
Params.lambdaGP=0.1; % temptation strength: v = lambdaGP*u_c(c) + shiftGP
Params.shiftGP=0; % constant shift of the temptation utility; should never change anything

% vfoptions.exoticpreferences='GulPesendorfer' is set inside each subcode

%% The only functions worth testing are the value fn ones, as after you have Policy everything else is anyway ignoring the temptation

%% Gul-Pesendorfer without declaring vfoptions.temptationFn is deliberately an error
vfoptionstemp.exoticpreferences='GulPesendorfer';
ReturnFn_none=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
try
    [Vtemp,Policytemp]=ValueFnIter_Case1_FHorz(0,n_a,0,N_j,[],a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionstemp);
    fprintf('GulPesendorfer without temptationFn: FAIL, this should error and did not \n')
catch ME
    if contains(ME.message,'temptationFn')
        fprintf('GulPesendorfer without temptationFn errors as intended :) \n')
    else
        fprintf('GulPesendorfer without temptationFn: FAIL, errored but with the wrong message: %s \n',ME.message)
    end
end
clear vfoptionstemp

%% ===== doPart(1): figs 1-8 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): figs 1-8 =====\n')
    %% without d, with z, without e
    figure_c=1;
    output=GPFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e
    figure_c=2;
    output=GPFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, with e
    figure_c=3;
    output=GPFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, with e
    figure_c=4;
    output=GPFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, with e
    figure_c=5;
    output=GPFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, with e
    figure_c=6;
    output=GPFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, without e (Gul-Pesendorfer works with no shocks at all)
    figure_c=7;
    output=GPFHorz_nod_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, without e
    figure_c=8;
    output=GPFHorz_d_noz_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(1): figs 1-8

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    %% The cross tests (see the comments at the top of the cross-test subcodes for what they cover)
    output=GPFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=GPFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(2): cross-tests


%% That is all the without semiz, now with semiz
% d1 is a decision variable that is not in the SemiExoStateFn; d2 (binary search effort) drives
% the semiz (employed/not-employed) transitions. Uses the same setup (which already carried the
% semi-exogenous state, it just wasn't used by the nosemiz cases).
addpath('../CoreFHorz_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(3): figs 9-16 =====
if doPart(3)==1
    fprintf('\n===== doPart(3): figs 9-16 =====\n')
    %% without d1, without z, without e (semiz)
    figure_c=9;
    output=GPFHorz_nod1_noz_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e (semiz)
    figure_c=10;
    output=GPFHorz_d1_noz_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e (semiz)
    figure_c=11;
    output=GPFHorz_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e (semiz)
    figure_c=12;
    output=GPFHorz_d1_z_noe_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e (semiz)
    figure_c=13;
    output=GPFHorz_nod1_noz_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e (semiz)
    figure_c=14;
    output=GPFHorz_d1_noz_e_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e (semiz)
    figure_c=15;
    output=GPFHorz_nod1_z_e_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e (semiz)
    % GP holds the return matrix AND its temptation twin simultaneously (~2x the core solver's
    % memory), so the d1+z+e semiz case gets a not-so-big grid for its big-grid moments block
    % (the 1001-point grid OOMs the GPU here; same remedy as the with2A d-variants below)
    n_a_semiznotsobig=501;
    a_grid_semiznotsobig=5*linspace(0,1,n_a_semiznotsobig)'.^3;
    figure_c=16;
    output=GPFHorz_d1_z_e_semiz(n_d_semiz,n_a,n_a_semiznotsobig,n_z,N_j,d_grid_semiz,a_grid,a_grid_semiznotsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(3): figs 9-16

%% ===== doPart(4): cross-tests 2 =====
if doPart(4)==1
    fprintf('\n===== doPart(4): cross-tests 2 =====\n')
    %% The semiz cross tests (see the comments at the top of the cross-test subcodes for what they cover)
    output=GPFHorz_CrossTests_nod1_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=GPFHorz_CrossTests_d1_semiz(n_d_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(4): cross-tests 2

%% with2A: TWO standard endogenous states (triggers the DC2A/GI2A/DC2A_GI2A code paths), under Gul-Pesendorfer
% Mirror of the QH/Ambiguity banks' with2A sections, reusing the same With2A ReturnFns as the
% exponential suite (all 8 shock combos, since Gul-Pesendorfer works with no shocks).

addpath('./withGulPesendorferPrefs_subcodes/With2A_subcodes/')
addpath('./withGulPesendorferPrefs_subcodes/With2A_subcodes/CrossTests/')
addpath('../CoreFHorz_ReturnFns/With2A_ReturnFns/')

% Redefine the asset grid to two endogenous states for this section
n_a_2A=[n_a,4];
n_a_2A_big=[n_a_big,4];
a2_grid_2A=[0;1;2;3];
a_grid_2A=[a_grid; a2_grid_2A];
a_grid_2A_big=[a_grid_big; a2_grid_2A];
n_a_notsobig=[501,4]; % to test Grid Interpolation without OOM in the d variants
a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
Params.phi1=3; % second endo-state preference params
Params.phi2=0.1;

%% ===== doPart(5): figs 17-24 =====
if doPart(5)==1
    fprintf('\n===== doPart(5): figs 17-24 =====\n')
    %% without d, with z, without e (with2A)
    figure_c=17;
    output=GPFHorz_nod_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, without e (with2A)
    figure_c=18;
    output=GPFHorz_d_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, with e (with2A)
    figure_c=19;
    output=GPFHorz_nod_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, with e (with2A)
    figure_c=20;
    output=GPFHorz_d_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, with z, with e (with2A)
    % notsobig here too (GP's temptation twin doubles matrix memory; the nod z+e 2A big grid
    % is the one nod case in the same OOM range as the d-variants)
    figure_c=21;
    output=GPFHorz_nod_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, with z, with e (with2A)
    % even notsobig OOMs here under GP (d x z x e x 2A is the biggest case, and the temptation
    % twin makes it three ~4GB matrices at once), so this one case gets a smaller grid again
    n_a_smaller2A=[301,4];
    a_grid_smaller2A=[5*linspace(0,1,n_a_smaller2A(1))'.^3; a2_grid_2A];
    figure_c=22;
    output=GPFHorz_d_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_smaller2A,n_z,N_j,d_grid,a_grid_2A,a_grid_smaller2A,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d, without z, without e (with2A)
    figure_c=23;
    output=GPFHorz_nod_noz_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d, without z, without e (with2A)
    figure_c=24;
    output=GPFHorz_d_noz_noe_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(5): figs 17-24

%% ===== doPart(6): cross-tests 3 =====
if doPart(6)==1
    fprintf('\n===== doPart(6): cross-tests 3 =====\n')
    %% The with2A cross tests (see the comments at the top of the cross-test subcodes for what they cover)
    output=GPFHorz_CrossTests_nod_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=GPFHorz_CrossTests_d_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(6): cross-tests 3

%% That is all the without semiz (with2A), now with semiz (with2A)
% The cross of the two tiers: two standard endogenous states AND a semi-exogenous state.
% Memory note: several cases get reduced big-grids (tighter than the QH bank's, since GP
% holds the return matrix and its temptation twin simultaneously, ~2x the core memory).
addpath('./withGulPesendorferPrefs_subcodes/With2A_subcodes/Semiz_subcodes/')
addpath('./withGulPesendorferPrefs_subcodes/With2A_subcodes/Semiz_subcodes/CrossTests/')
addpath('../CoreFHorz_ReturnFns/With2A_ReturnFns/Semiz_ReturnFns/')

%% ===== doPart(7): figs 25-32 =====
if doPart(7)==1
    fprintf('\n===== doPart(7): figs 25-32 =====\n')
    %% without d1, without z, without e (semiz, with2A)
    figure_c=25;
    output=GPFHorz_nod1_noz_noe_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, without e (semiz, with2A)
    figure_c=26;
    n_a_notsobig=[501,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_d1_noz_noe_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, without e (semiz, with2A)
    figure_c=27;
    n_a_notsobig=[501,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_nod1_z_noe_semiz_with2A(n_d2_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, without e (semiz, with2A)
    figure_c=28;
    n_a_notsobig=[301,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_d1_z_noe_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, without z, with e (semiz, with2A)
    figure_c=29;
    n_a_notsobig=[501,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_nod1_noz_e_semiz_with2A(n_d2_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, without z, with e (semiz, with2A)
    figure_c=30;
    n_a_notsobig=[401,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_d1_noz_e_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% without d1, with z, with e (semiz, with2A)
    figure_c=31;
    n_a_notsobig=[401,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_nod1_z_e_semiz_with2A(n_d2_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

    %% with d1, with z, with e (semiz, with2A)
    figure_c=32;
    n_a_notsobig=[201,4];
    a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
    output=GPFHorz_d1_z_e_semiz_with2A(n_d_semiz,n_a_2A,n_a_notsobig,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzGulPesendorferTests_Fig',num2str(figure_c),'.png'],'Resolution',150)
end % doPart(7): figs 25-32

%% ===== doPart(8): cross-tests 4 =====
if doPart(8)==1
    fprintf('\n===== doPart(8): cross-tests 4 =====\n')
    %% The semiz-with2A cross tests (see the comments at the top of the cross-test subcodes for what they cover)
    output=GPFHorz_CrossTests_nod1_semiz_with2A(n_d2_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d2_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

    output=GPFHorz_CrossTests_d1_semiz_with2A(n_d_semiz,n_a_2A,n_a_2A_big,n_z,N_j,d_grid_semiz,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
end % doPart(8): cross-tests 4

%% Done
%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreFHorzGulPesendorferTestsdiary.txt')

diary off
