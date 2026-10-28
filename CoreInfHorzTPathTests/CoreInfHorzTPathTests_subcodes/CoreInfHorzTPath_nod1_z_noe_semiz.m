function output=CoreInfHorzTPath_nod1_z_noe_semiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c)
% without d1, with z, without e, with a semi-exogenous state semiz (fig 11).
% WRITTEN 2026-10-27, TEST-FIRST, against the InfHorz transition path semiz support that is being built
% in the toolkit right now. Built from the live CoreInfHorzTPath_nod_z_noe_nosemiz (fig 3), the subcode with the same
% z/e structure, and adapted to semiz the way the stationary CoreInfHorzTests semiz subcodes adapt
% their nosemiz siblings.
%
% THE MODEL. As in CoreInfHorzTests doPart(9)-(11): semiz is employment (n_semiz=2, semiz_grid=[0;1], 1
% is employed), its transitions pi_semiz(semiz,semiz',d2) depend on a binary search effort d2, the LAST
% decision variable (CoreInfHorzTPathSetup_SemiExoStateFn: d2=1 makes both states ten times more
% persistent, and costs searcheffortcost in utility). The unemployed get uempbenefit.
% n_d is n_d2_semiz and d_grid is d2_grid_semiz: the only decision variable is d2.
%
% USER-FACING INTERFACE (the shared conventions of the InfHorz TPath semiz build): semiz goes in through
% vfoptions/simoptions (n_semiz, semiz_grid, SemiExoStateFn, l_dsemiz), exactly as for the stationary
% commands, and simoptions also carries d_grid. Exogenous states are ordered (semiz,z,e), semiz fastest;
% ReturnFn and FnsToEvaluate take (d,aprime,a,semiz,z,e,...). VPath is [n_a,n_semiz,(n_z),(n_e),T], PolicyPath
% [rows,n_a,n_semiz,(n_z),(n_e),T] with rows d1 (if any), then d2, then aprime (then L2 and L2flag with GI).
%
% WHAT IS DIFFERENT FROM THE NOSEMIZ SIBLING, AND WHY:
%  - The stationary InfHorz semiz solver does NOT support divide-and-conquer, so V_final, Policy_final and
%    the initial agent dist come from a non-DC stationary solve (as they do in the sibling anyway), and the
%    stationary general eqm solves use GI without DC; the PATH supports DC, so every with/without-DC path
%    comparison of the sibling is kept.
%  - lowmemory=2 is refused with semiz by the single-step raws (lowmemory>=2 errors, with or without e),
%    so only lowmemory=1 is checked, and the lowmemory=2 lines of the nosemiz e siblings are dropped.
%    lowmemory=1 loops over (semiz,z) on the plain tier, but on the divide-and-conquer tier without e it
%    is accepted and solved as lowmemory=0, so that one pair is identical by construction (kept, cheap).
%  - The plain and DC1 raws maximise jointly over (d1,aprime) within each d2 (d1 is not refined out).
%  - A path-varying semiz block, against FHorz: probfindjob on the ParamPath, so pi_semiz moves along the
%    path. Plus the same path given as an explicit 4-D pi_semiz [n_semiz,n_semiz,n_d2,T], which must
%    give the identical answer.
%  - The AutoCorrTransProbs dense-kernel oracle carries the pi_semiz(semiz,semiz',d2) factor, with d2 read
%    off the policy at each state.
%  - General eqm: the closed-economy firm of the z subcodes (both r and w general eqm prices), in every
%    semiz subcode including the noz noe ones. Unlike the deterministic nosemiz noz noe model, the semiz
%    model has idiosyncratic income risk (losing the job), so households hold precautionary savings and
%    K(r) rises continuously from near zero toward infinity as r approaches 1/beta-1: the stationary
%    general eqm is well posed. The top of a_grid_GE must not bind, and that is checked.
%  - The no-GI stationary general eqm is started from the GI answer, as in the with-d nosemiz subcodes,
%    in every semiz subcode: d2 is a discrete choice even without d1, so N moves in steps with the prices.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% The agent distribution is compared against itself after T periods of a constant path, with a
% 1e-6 bar. StationaryDist_InfHorz defaults simoptions.tolerance to 1e-6 on the L-Infinity norm of
% the distribution, so at the default the reference object carries exactly as much error as the bar
% allows and the check can only fail: under a fixed policy the forward iteration moves TOWARD the
% true fixed point, so what it measures is the starting distribution's own error. Converging the
% stationary dist three orders tighter than the bar makes the check test the path iteration
% instead. maxit is 1e6 and warns if missed, so this cannot fail silently.
simoptions.tolerance=1e-9;
% semi-exogenous state, passed via vfoptions/simoptions (simoptions also needs d_grid)
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptions.l_dsemiz=vfoptionsbaseline.l_dsemiz;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.l_dsemiz=simoptionsbaseline.l_dsemiz;
simoptions.d_grid=d_grid;
n_semiz=vfoptions.n_semiz;
semiz_grid=vfoptions.semiz_grid;
n_d2=n_d(end); % d2 is the last decision variable
d2_grid=d_grid(end-n_d2+1:end);

ReturnFn=@(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w) w*z*semiz;
FnsToEvaluate.nextassets=@(d2,aprime,a,semiz,z) aprime; % aprime-dependent, so it tests the policy decode directly
FnsToEvaluate.searcheffort=@(d2,aprime,a,semiz,z) d2; % the decision that drives semiz, so it tests the d2 decode directly
FnsToEvaluate.semizval=@(d2,aprime,a,semiz,z) semiz; % the semi-exogenous state itself (the employment rate), moves in the path-varying semiz block
FnsToEvaluate.zval=@(d2,aprime,a,semiz,z) z; % the markov shock itself, so the path of its mean/std deviation tracks the z process


%% Period-0 VFI: gives the final-step V/Policy (used as both V_final for TPath and the steady state to compare against)
% No divide-and-conquer in any stationary solve: the stationary InfHorz semiz solver does not support it.
vfoptions1=vfoptions;
simoptions1=simoptions;
[V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% GI variant of Policy_final shape (extra L2 channels)
vfoptions1_GI=vfoptions;
vfoptions1_GI.gridinterplayer=1;
vfoptions1_GI.ngridinterp=5;
simoptions1_GI=simoptions;
simoptions1_GI.gridinterplayer=vfoptions1_GI.gridinterplayer;
simoptions1_GI.ngridinterp=vfoptions1_GI.ngridinterp;
[V_final_GI,Policy_final_GI]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_GI);

% Big-grid versions
[V_final_big,Policy_final_big]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
[V_final_big_GI,Policy_final_big_GI]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_GI);

% Stationary dist under period-0 prices (used as initial dist for TPath)
AgentDist_initial=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z,pi_z,simoptions1,Params,[]);
AgentDist_initial_big=StationaryDist_InfHorz(Policy_final_big,n_d,n_a_big,n_z,pi_z,simoptions1,Params,[]);

%% With and without divide-and-conquer
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);

vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[VPath2,PolicyPath2]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

% lowmemory should give same answer (lowmemory=2 is refused with semiz, so only lowmemory=1; see the header)
% The vfoptions2 (divide-and-conquer, no GI) pair below is IDENTICAL BY CONSTRUCTION: without an e, the
% DC1 semiz raws accept lowmemory=1 but solve it as lowmemory=0. Kept because it is cheap and checks the
% option is accepted; it is NOT evidence that a lowmemory loop works (the vfoptions1 pair is).
vfoptions1.lowmemory=1;
[VPath1B,PolicyPath1B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[VPath2B,PolicyPath2B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
vfoptions2.lowmemory=0;

fprintf('low memory, this should be zero: %.3e \n',max(abs(VPath1(:)-VPath1B(:))))
fprintf('low memory, this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath1B(:))))
fprintf('low memory, this should be zero: %.3e \n',max(abs(VPath2(:)-VPath2B(:))))
fprintf('low memory, this should be zero: %.3e \n',max(abs(PolicyPath2(:)-PolicyPath2B(:))))

%% ValueFnFromPolicy along the path: recomputing V from the policy must reproduce VPath
% ValueFnFromPolicyOnTransPath_InfHorz evaluates the Bellman at the given PolicyPath instead of
% maximising, so it must return exactly the VPath that produced that policy. (Divide-and-conquer
% is irrelevant here -- there is no maximisation -- so this covers the no-GI tier as a whole.)
% With semiz the policy's d2 row also picks the pi_semiz(.,.,d2) slice of the expectation, so this
% checks that the d2 index is read back the way the maximisation wrote it.
VfromPolicyPath=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPath2,V_final,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions2);
fprintf('ValueFnFromPolicy along TPath (no GI), this should be zero: %.3e \n',max(abs(VfromPolicyPath(:)-VPath2(:))))
clear VfromPolicyPath

%% Cross-check against an FHorz solve (with DC, no GI)
% The TPath value fn iteration is a finite-horizon problem in disguise: encode the price/param
% path as age-dependent parameters, hand the terminal V in via vfoptions.V_Jplus1, and the two
% must agree. Timing: ValueFnOnTransPath_InfHorz sets VPath(:,:,T)=V_final and then, for
% t=T-1,...,1, computes V_t from params_t with continuation V_{t+1}; row T of the paths is never
% read. ValueFnIter_Case1_FHorz with V_Jplus1 computes V_{N_j} from params_{N_j} with
% continuation V_Jplus1. So the FHorz model has N_j=T-1 periods, its age-dependent parameters
% are rows 1,...,T-1 of the paths, and V_FHorz(:,:,j) must equal VPath(:,:,j).
% This is the only check in this file that can catch an off-by-one in the path indexing: every
% other comparison here either runs a constant path (where t and t+-1 are indistinguishable) or
% moves both sides of the comparison together.
% Done at the divide-and-conquer tier (cheaper); VPath2 equals VPath1 exactly, so this covers
% the no-DC tier too.
%
% FHorz with semiz and V_Jplus1 (SemiExogShockSetup_FHorz): pi_semiz_J is built from the SemiExoStateFn
% at the age-jj parameters with N_jFH slices (not N_jFH-1), the last being the transition from period
% N_jFH into the V_Jplus1 period; V_Jplus1 is reshaped to [N_a,N_semiz*N_z(,N_e)], so V_final goes in
% as it comes out of the stationary solve. Here probfindjob is a scalar, so every slice is the pi_semiz
% the path uses. The FHorz semiz solvers support divide-and-conquer, so vfoptions2 is used as is.
N_jFH=T-1;
ParamsFH=Params;
ParamsFH.r=PricePath.r(1:N_jFH);
ParamsFH.w=ParamPath.w(1:N_jFH);

vfoptionsFH=vfoptions2;
vfoptionsFH.V_Jplus1=V_final;
[V_FH,Policy_FH]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFH,DiscountFactorParamNames,[],vfoptionsFH);

% Reshape to (everything else)-by-periods, so these lines are agnostic to the state dimensions
temp1=reshape(VPath2,[],T); temp2=reshape(V_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, V: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyPath2,[],T); temp2=reshape(Policy_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, Policy: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath2,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions2),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, PolicyVals: %.3e \n',max(abs(dev(:))))

% --- carry the same construction through the agent dist and the downstream statistics ---
% The FHorz stationary dist bakes in the age weights (StationaryDist_FHorz_Iteration_raw ends
% with .*AgeWeights), whereas the TPath dist has mass one in every period. So use uniform age
% weights and divide each age slice back out by its own mass. LifeCycleProfiles renormalises
% within each age group itself, so the weights cancel there and need no undoing.
% whichstats must be pinned on both sides: EvalFnOnTransPath_AllStats_InfHorz defaults to
% ones(7,1) but LifeCycleProfiles_FHorz_Case1 defaults to [1,1,1,2,1,2,1], and entries 4 and 6
% take a different path through StatsFromWeightedGrid.
% The FnsToEvaluate use w, which lives on the ParamPath, so these lines additionally check that
% the evaluation code indexes the path the same way FHorz indexes age-dependent parameters --
% a different code path from the value fn iteration checked above.
ParamsFH.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFH.mewjFH(end)=1-sum(ParamsFH.mewjFH(1:end-1)); % sum must be one to within 1e-15
simoptionsFH=simoptions2;
simoptionsFH.npoints=100;
simoptionsFH.nquantiles=20;
simoptionsFH.whichstats=ones(7,1);

AgentDistPathFH=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath2, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptionsFH);
StatDist_FH=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FH,n_d,n_a,n_z,N_jFH,pi_z,ParamsFH,simoptionsFH);
temp1=reshape(AgentDistPathFH,[],T); temp2=reshape(StatDist_FH,[],N_jFH);
temp2=temp2./sum(temp2,1); % undo the age weighting, so each age carries mass one like each TPath period
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AgentDist: %.3e \n',max(abs(dev(:))))

AggVarsPathFH=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsFH);
AllStatsPathFH=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsFH);
AgeCondStats_FH=LifeCycleProfiles_FHorz_Case1(StatDist_FH,Policy_FH,FnsToEvaluate,ParamsFH,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFH);

for fnname=fieldnames(FnsToEvaluate)'
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AggVars %s Mean: %.3e \n',fnname{1},max(abs(dev(:))))
end
% AllStats: every field the two commands have in common (periods/ages are the trailing dim in both)
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AllStats assets %s: %.3e \n',statname{1},max(abs(dev(:))))
end
% Kept for the non-vacuity check of the path-varying semiz block below: the employment rate along this
% (fixed semiz process) path.
semizMeanFixed=gather(AggVarsPathFH.semizval.Mean);

clear V_FH Policy_FH ParamsFH vfoptionsFH temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

%% Cross-check against an FHorz solve, with a shock process that varies along the path
% Same construction as the block above, but now the exogenous process itself moves: z_grid and
% pi_z are period dependent rather than fixed. This is the bank's only coverage of
% transpathoptions.zpathtrivial=0 alongside a semiz, and it is checked the same way, against an FHorz
% solve of the identical problem, which is what makes it an exact comparison rather than a plausibility
% check. The semiz process is held fixed here (it moves in its own block below).
%
% Timing. discretizeLifeCycleAR1_KFTT returns pi_z_J with slice j the transition from period j to
% j+1 (it shifts to that convention at the end of that file), which is also what
% ValueFnOnTransPath_InfHorz uses (V_t is computed from pi_z_T(:,:,t), see the backward loop in
% TransitionPath_InfHorz_substeps_Step1_ValueFnIter) and what ValueFnIter_Case1_FHorz documents.
% So the two line up slice for slice, exactly as the price and param paths do.
%
% Why Jz=T+1: KFTT returns J-1 transition matrices for J ages, since there is no age J+1 to
% transition to. The TPath wants T of them (slice T is never read, as VPath(:,:,T)=V_final is
% assigned rather than computed, but the setup still checks the size), and FHorz with V_Jplus1 wants
% N_jFH=T-1. So ask for T+1 ages, take all T transition slices for the TPath and the first T-1 for
% FHorz. The grids come back with T+1 columns, one more than the TPath needs, so they are trimmed.
%
% Note V_final is the stationary V under the ORIGINAL z_grid, whereas the KFTT grid at t=T holds
% different values. Economically odd, but irrelevant to the comparison: both sides are handed the
% identical array and index it by z-index, so the check stays exact.
Jz=T+1;
kfttmew=zeros(1,Jz);
kfttrho=linspace(0.95,0.75,Jz);  % persistence falls along the path
kfttsigma=0.03*linspace(1,2,Jz); % innovation std deviation doubles along the path
kfttoptions.nMoments=2; % the default of 4 cannot be hit with this few grid points, and the
                        % unmet-moments warning is just noise here: any valid markov chain that
                        % actually varies by period serves this test
[z_grid_KFTT,pi_z_KFTT]=discretizeLifeCycleAR1_KFTT(kfttmew,kfttrho,kfttsigma,n_z,Jz,kfttoptions);
z_grid_KFTT=exp(z_grid_KFTT); % same exp() convention as the baseline z_grid
z_grid_KFTT=z_grid_KFTT(:,1:T); % KFTT returns Jz=T+1 columns; the TPath wants T

[VPathZ,PolicyPathZ]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid_KFTT, pi_z_KFTT, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
% ValueFnFromPolicy with the z process varying along the path. Until 2026-10-28 the z-only branch of
% ValueFnFromPolicyOnTransPath_InfHorz read a fixed z_grid and pi_z, so here it would have evaluated every period at
% the period-1 process; the fixed-process check above could never see that.
VfromPolicyPathZ=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPathZ,V_final,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid_KFTT, pi_z_KFTT, ReturnFn, Params, DiscountFactorParamNames, vfoptions2);
fprintf('ValueFnFromPolicy along TPath with z varying along the path, this should be zero: %.3e \n',max(abs(VfromPolicyPathZ(:)-VPathZ(:))))
clear VfromPolicyPathZ

ParamsFHZ=Params;
ParamsFHZ.r=PricePath.r(1:N_jFH);
ParamsFHZ.w=ParamPath.w(1:N_jFH);
vfoptionsFHZ=vfoptions2;
vfoptionsFHZ.V_Jplus1=V_final;
[V_FHZ,Policy_FHZ]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid_KFTT(:,1:N_jFH),pi_z_KFTT(:,:,1:N_jFH),ReturnFn,ParamsFHZ,DiscountFactorParamNames,[],vfoptionsFHZ);

temp1=reshape(VPathZ,[],T); temp2=reshape(V_FHZ,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying z), this should be zero, V: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyPathZ,[],T); temp2=reshape(Policy_FHZ,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying z), this should be zero, Policy: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPathZ,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FHZ,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions2),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying z), this should be zero, PolicyVals: %.3e \n',max(abs(dev(:))))

% Same age-weight handling as the fixed-process block above
ParamsFHZ.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFHZ.mewjFH(end)=1-sum(ParamsFHZ.mewjFH(1:end-1));
simoptionsFHZ=simoptions2;
simoptionsFHZ.npoints=100;
simoptionsFHZ.nquantiles=20;
simoptionsFHZ.whichstats=ones(7,1);

AgentDistPathZ=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPathZ, n_d, n_a, n_z, pi_z_KFTT, T, Params, transpathoptionsbaseline, simoptionsFHZ);
StatDist_FHZ=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FHZ,n_d,n_a,n_z,N_jFH,pi_z_KFTT(:,:,1:N_jFH),ParamsFHZ,simoptionsFHZ);
temp1=reshape(AgentDistPathZ,[],T); temp2=reshape(StatDist_FHZ,[],N_jFH);
temp2=temp2./sum(temp2,1); % undo the age weighting, so each age carries mass one like each TPath period
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying z), this should be zero, AgentDist: %.3e \n',max(abs(dev(:))))

AggVarsPathZ=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathZ, PolicyPathZ, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid_KFTT, transpathoptionsbaseline, simoptionsFHZ);
AllStatsPathZ=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathZ, PolicyPathZ, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid_KFTT, transpathoptionsbaseline, simoptionsFHZ);
AgeCondStats_FHZ=LifeCycleProfiles_FHorz_Case1(StatDist_FHZ,Policy_FHZ,FnsToEvaluate,ParamsFHZ,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid_KFTT(:,1:N_jFH),simoptionsFHZ);

for fnname=fieldnames(FnsToEvaluate)'
    dev=AggVarsPathZ.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FHZ.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (path-varying z), this should be zero, AggVars %s Mean: %.3e \n',fnname{1},max(abs(dev(:))))
end
% The shock FnsToEvaluate, which is what this block is really about: with the process moving, the
% mean and the std deviation of the shock both move, and the two solvers must agree on both.
for fnname={'zval','semizval'}
    for statname={'Mean','StdDeviation'}
        temp1=AllStatsPathZ.(fnname{1}).(statname{1}); temp2=AgeCondStats_FHZ.(fnname{1}).(statname{1});
        dev=temp1(:,1:N_jFH)-temp2;
        fprintf('TPath vs FHorz (path-varying z), this should be zero, %s %s: %.3e \n',fnname{1},statname{1},max(abs(dev(:))))
    end
end
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathZ.assets.(statname{1}); temp2=AgeCondStats_FHZ.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (path-varying z), this should be zero, AllStats assets %s: %.3e \n',statname{1},max(abs(dev(:))))
end

% Non-vacuity. If the time-varying z_grid/pi_z were silently ignored, everything above would still
% agree -- both sides would just be solving the fixed-process model -- and the block would prove
% nothing. Here the KFTT innovation std deviation doubles along the path, so the shock's own std
% deviation must visibly move. Printed in the inverted form, which CoreSummary counts and fails if
% the range is not well above zero.
temp1=AllStatsPathZ.zval.StdDeviation;
fprintf('Path-varying z bit, range across periods of zval StdDeviation, this should NOT be zero: %.3e \n',max(temp1(:))-min(temp1(:)))

clear Jz kfttoptions kfttmew kfttrho kfttsigma z_grid_KFTT pi_z_KFTT VPathZ PolicyPathZ ParamsFHZ vfoptionsFHZ V_FHZ Policy_FHZ simoptionsFHZ AgentDistPathZ StatDist_FHZ AggVarsPathZ AllStatsPathZ AgeCondStats_FHZ temp1 temp2 dev

%% Cross-check against an FHorz solve, with a semiz transition that varies along the path
% The semiz counterpart of the path-varying z block, and the bank's only coverage of
% transpathoptions.semizpathtrivial=0. The SemiExoStateFn takes probfindjob as a parameter, so putting
% probfindjob on the ParamPath makes pi_semiz period dependent: SemiExogShockSetup_InfHorz_TPath builds
% pi_semiz_T [n_semiz,n_semiz,n_d2,T], slice t from probfindjob(t). Checked against an FHorz solve of the
% identical problem, where probfindjob is an age-dependent parameter.
%
% TIMING. pi_semiz_T(:,:,:,t) is the transition from t to t+1: V_t takes its expectation over semiz' with
% slice t, and the agent dist from t to t+1 uses slice t. FHorz: pi_semiz_J(:,:,:,j) is the transition
% from age j to j+1, built from the age-j parameters (and with V_Jplus1 there are N_jFH slices, the last
% the transition into the V_Jplus1 period). So TPath period t is FHorz age j=t, exactly as for the price
% and param paths, and probfindjob(1:N_jFH) is what FHorz gets.
%
% NOT SYMMETRIC IN TIME, on purpose. probfindjob rises monotonically, 0.5 to 0.9, over the whole path, so
% an implementation that reads slice t+1 (or t-1) where slice t is meant moves every job-finding
% probability by 0.4/(T-1) per period -- 0.04 of it under d2=1 -- and the V, Policy and AgentDist
% comparisons catch it. A path that only moved near the end could let an off-by-one through.
%
% V_final and AgentDist_initial are the stationary ones under the scalar Params.probfindjob; as in the z
% block, both sides are handed the identical arrays, so the check stays exact.
probfindjob_T=linspace(0.5,0.9,T);
ParamPathS=ParamPath;
ParamPathS.probfindjob=probfindjob_T;
[VPathS,PolicyPathS]=ValueFnOnTransPath_InfHorz(PricePath, ParamPathS, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);

ParamsFHS=Params;
ParamsFHS.r=PricePath.r(1:N_jFH);
ParamsFHS.w=ParamPath.w(1:N_jFH);
ParamsFHS.probfindjob=probfindjob_T(1:N_jFH); % age-dependent: FHorz builds slice j from probfindjob(j)
vfoptionsFHS=vfoptions2;
vfoptionsFHS.V_Jplus1=V_final;
[V_FHS,Policy_FHS]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFHS,DiscountFactorParamNames,[],vfoptionsFHS);

temp1=reshape(VPathS,[],T); temp2=reshape(V_FHS,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying semiz), this should be zero, V: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyPathS,[],T); temp2=reshape(Policy_FHS,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying semiz), this should be zero, Policy: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPathS,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FHS,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions2),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying semiz), this should be zero, PolicyVals: %.3e \n',max(abs(dev(:))))

% Same age-weight handling as the fixed-process block above
ParamsFHS.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFHS.mewjFH(end)=1-sum(ParamsFHS.mewjFH(1:end-1));
simoptionsFHS=simoptions2;
simoptionsFHS.npoints=100;
simoptionsFHS.nquantiles=20;
simoptionsFHS.whichstats=ones(7,1);

AgentDistPathS=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPathS, PolicyPathS, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptionsFHS);
StatDist_FHS=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FHS,n_d,n_a,n_z,N_jFH,pi_z,ParamsFHS,simoptionsFHS);
temp1=reshape(AgentDistPathS,[],T); temp2=reshape(StatDist_FHS,[],N_jFH);
temp2=temp2./sum(temp2,1); % undo the age weighting, so each age carries mass one like each TPath period
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying semiz), this should be zero, AgentDist: %.3e \n',max(abs(dev(:))))

AggVarsPathS=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathS, PolicyPathS, PricePath, ParamPathS, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsFHS);
AllStatsPathS=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathS, PolicyPathS, PricePath, ParamPathS, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsFHS);
AgeCondStats_FHS=LifeCycleProfiles_FHorz_Case1(StatDist_FHS,Policy_FHS,FnsToEvaluate,ParamsFHS,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFHS);

for fnname=fieldnames(FnsToEvaluate)'
    dev=AggVarsPathS.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FHS.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (path-varying semiz), this should be zero, AggVars %s Mean: %.3e \n',fnname{1},max(abs(dev(:))))
end
for fnname={'semizval','searcheffort'}
    for statname={'Mean','StdDeviation'}
        temp1=AllStatsPathS.(fnname{1}).(statname{1}); temp2=AgeCondStats_FHS.(fnname{1}).(statname{1});
        dev=temp1(:,1:N_jFH)-temp2;
        fprintf('TPath vs FHorz (path-varying semiz), this should be zero, %s %s: %.3e \n',fnname{1},statname{1},max(abs(dev(:))))
    end
end
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathS.assets.(statname{1}); temp2=AgeCondStats_FHS.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (path-varying semiz), this should be zero, AllStats assets %s: %.3e \n',statname{1},max(abs(dev(:))))
end

% Non-vacuity. If the probfindjob on the ParamPath were silently ignored (pi_semiz built once, from the
% scalar Params.probfindjob), everything above would still agree -- both sides would just be solving the
% fixed-process model -- and the block would prove nothing. A range-across-periods check is NOT enough
% here: the employment rate moves along the baseline path anyway, because r and w move and so does the
% choice of d2. So compare against the same path with the fixed process (VPath2, and the employment rate
% kept from the fixed-process FHorz block above): they must differ. Inverted form, counted.
fprintf('Path-varying semiz bit, VPath vs the fixed-semiz-process VPath, max abs difference, this should NOT be zero: %.3e \n',max(abs(VPathS(:)-VPath2(:))))
fprintf('Path-varying semiz bit, employment rate (semizval Mean) vs the fixed-semiz-process path, max abs difference, this should NOT be zero: %.3e \n',max(abs(gather(AggVarsPathS.semizval.Mean(:))-semizMeanFixed(:))))

% The same path again, with pi_semiz given explicitly as [n_semiz,n_semiz,n_d2,T] instead of through
% the SemiExoStateFn and a ParamPath parameter (the other route to semizpathtrivial=0). Built here by
% loops from the setup's own SemiExoStateFn, so it is the identical process, and the ParamPath is the
% ORIGINAL one (no probfindjob): VPath, PolicyPath and AgentDistPath must be identical to the above.
pi_semiz_T=zeros(n_semiz,n_semiz,n_d2,T);
for tt=1:T
    for s1=1:n_semiz
        for s2=1:n_semiz
            for d2c=1:n_d2
                pi_semiz_T(s1,s2,d2c,tt)=CoreInfHorzTPathSetup_SemiExoStateFn(semiz_grid(s1),semiz_grid(s2),d2_grid(d2c),probfindjob_T(tt),Params.problosejob);
            end
        end
    end
end
vfoptionsS4=rmfield(vfoptions2,'SemiExoStateFn');
vfoptionsS4.pi_semiz=pi_semiz_T;
simoptionsS4=rmfield(simoptions2,'SemiExoStateFn');
simoptionsS4.pi_semiz=pi_semiz_T;
[VPathS4,PolicyPathS4]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsS4);
AgentDistPathS4=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPathS4, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptionsS4);
fprintf('Path-varying semiz, explicit 4-D pi_semiz vs SemiExoStateFn with probfindjob on the ParamPath, this should be zero, V: %.3e \n',max(abs(VPathS4(:)-VPathS(:))))
fprintf('Path-varying semiz, explicit 4-D pi_semiz vs SemiExoStateFn with probfindjob on the ParamPath, this should be zero, Policy: %.3e \n',max(abs(PolicyPathS4(:)-PolicyPathS(:))))
fprintf('Path-varying semiz, explicit 4-D pi_semiz vs SemiExoStateFn with probfindjob on the ParamPath, this should be zero, AgentDist: %.3e \n',max(abs(AgentDistPathS4(:)-AgentDistPathS(:))))

clear probfindjob_T ParamPathS VPathS PolicyPathS ParamsFHS vfoptionsFHS V_FHS Policy_FHS simoptionsFHS AgentDistPathS StatDist_FHS AggVarsPathS AllStatsPathS AgeCondStats_FHS semizMeanFixed pi_semiz_T vfoptionsS4 simoptionsS4 VPathS4 PolicyPathS4 AgentDistPathS4 temp1 temp2 dev

%%
clear VPath1 VPath2 PolicyPath1 PolicyPath2 VPath1B VPath2B PolicyPath1B PolicyPath2B

%% Solve with grid-interpolation. With and without divide-and-conquer
vfoptions3=vfoptions1_GI;
simoptions3=simoptions1_GI;
[VPath3,PolicyPath3]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);

vfoptions4=vfoptions3;
vfoptions4.divideandconquer=1;
simoptions4=simoptions3;
[VPath4,PolicyPath4]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);

fprintf('Divide-and-conquer (with GI), this should be zero: %.3e \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI), this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

% lowmemory should give same answer
vfoptions3.lowmemory=1;
[VPath3B,PolicyPath3B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[VPath4B,PolicyPath4B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
vfoptions4.lowmemory=0;

fprintf('low memory, this should be zero: %.3e \n',max(abs(VPath3(:)-VPath3B(:))))
fprintf('low memory, this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath3B(:))))
fprintf('low memory, this should be zero: %.3e \n',max(abs(VPath4(:)-VPath4B(:))))
fprintf('low memory, this should be zero: %.3e \n',max(abs(PolicyPath4(:)-PolicyPath4B(:))))

%% ValueFnFromPolicy along the path, grid-interpolation tier
VfromPolicyPath=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPath4,V_final_GI,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions4);
fprintf('ValueFnFromPolicy along TPath (with GI), this should be zero: %.3e \n',max(abs(VfromPolicyPath(:)-VPath4(:))))
clear VfromPolicyPath

%% Cross-check against an FHorz solve (with DC and GI)
% Same construction as the no-GI version above, against the grid-interpolation tier. Not
% redundant with it: GI genuinely changes the solution, so it needs its own comparison.
% THE RAW POLICY INDEXES ARE NOT COMPARED at this tier, only PolicyVals: under GI the indexes can
% legitimately differ, because (L1,L2) and (L1+1,L2-(ngridinterp+1)) encode the same aprime when the
% optimum sits on a coarse grid point, and the two implementations need not break that tie the same
% way. PolicyVals is immune to it and is the one that has to be zero.
N_jFH=T-1;
ParamsFH=Params;
ParamsFH.r=PricePath.r(1:N_jFH);
ParamsFH.w=ParamPath.w(1:N_jFH);

vfoptionsFH_GI=vfoptions4;
vfoptionsFH_GI.V_Jplus1=V_final_GI;
[V_FH,Policy_FH]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFH,DiscountFactorParamNames,[],vfoptionsFH_GI);

temp1=reshape(VPath4,[],T); temp2=reshape(V_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, V: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath4,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions4),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, PolicyVals: %.3e \n',max(abs(dev(:))))

% --- carry the same construction through the agent dist and the downstream statistics ---
% (same age-weight and whichstats handling as the no-GI block above)
ParamsFH.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFH.mewjFH(end)=1-sum(ParamsFH.mewjFH(1:end-1)); % sum must be one to within 1e-15
simoptionsFH=simoptions4;
simoptionsFH.npoints=100;
simoptionsFH.nquantiles=20;
simoptionsFH.whichstats=ones(7,1);

AgentDistPathFH=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath4, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptionsFH);
StatDist_FH=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FH,n_d,n_a,n_z,N_jFH,pi_z,ParamsFH,simoptionsFH);
temp1=reshape(AgentDistPathFH,[],T); temp2=reshape(StatDist_FH,[],N_jFH);
temp2=temp2./sum(temp2,1); % undo the age weighting, so each age carries mass one like each TPath period
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, AgentDist: %.3e \n',max(abs(dev(:))))

AggVarsPathFH=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath4, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsFH);
AllStatsPathFH=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath4, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsFH);
AgeCondStats_FH=LifeCycleProfiles_FHorz_Case1(StatDist_FH,Policy_FH,FnsToEvaluate,ParamsFH,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFH);

for fnname=fieldnames(FnsToEvaluate)'
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (with DC and GI), this should be zero, AggVars %s Mean: %.3e \n',fnname{1},max(abs(dev(:))))
end
% AllStats: every field the two commands have in common (periods/ages are the trailing dim in both)
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (with DC and GI), this should be zero, AllStats assets %s: %.3e \n',statname{1},max(abs(dev(:))))
end

clear V_FH Policy_FH ParamsFH vfoptionsFH_GI temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

%%
clear VPath3 VPath4 PolicyPath3 PolicyPath4 VPath3B VPath4B PolicyPath3B PolicyPath4B

%% Big a_grid: moments along path should be close with/without grid interp
[VPath2b,PolicyPath2b]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big, Params, n_d, n_a_big, n_z, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPath2=AgentDistOnTransPath_InfHorz(AgentDist_initial_big, PricePath, ParamPath, PolicyPath2b, n_d, n_a_big, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, d_grid, a_grid_big, z_grid, transpathoptionsbaseline, simoptions2);

[VPath4b,PolicyPath4b]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_big_GI, Policy_final_big_GI, Params, n_d, n_a_big, n_z, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_InfHorz(AgentDist_initial_big, PricePath, ParamPath, PolicyPath4b, n_d, n_a_big, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, d_grid, a_grid_big, z_grid, transpathoptionsbaseline, simoptions4);

%% SimPanel along the path: per-period panel mean should reproduce the AggVars along the path
% (Monte Carlo simulation, so this is a roughly-equal check, not machine precision)
% Done twice: off the no-GI policy path, and off the GI policy path. The GI one matters
% because with gridinterplayer=1 the Policy carries the extra L2 and L2flag channels that
% SimPanelValues_TransPath_InfHorz has to strip. With semiz the panel also has to draw semiz' from
% pi_semiz(semiz,.,d2) with d2 read off the policy, which is what the employment-rate line checks.
% Each panel is compared against the AggVars built from the SAME policy path.
simoptionsSP=simoptions2;
simoptionsSP.numbersims=10^4;
simoptionsSP.simperiods=T;
SimPanelTPath=SimPanelValues_TransPath_InfHorz(PolicyPath2b, PricePath, ParamPath, T, AgentDist_initial_big, n_d, n_a_big, n_z, pi_z, d_grid, a_grid_big, z_grid, FnsToEvaluate, Params, transpathoptionsbaseline, simoptionsSP);

simoptionsSP_GI=simoptions4;
simoptionsSP_GI.numbersims=10^4;
simoptionsSP_GI.simperiods=T;
SimPanelTPath_GI=SimPanelValues_TransPath_InfHorz(PolicyPath4b, PricePath, ParamPath, T, AgentDist_initial_big, n_d, n_a_big, n_z, pi_z, d_grid, a_grid_big, z_grid, FnsToEvaluate, Params, transpathoptionsbaseline, simoptionsSP_GI);

fprintf('SimPanel along TPath: per-period panel mean should roughly match AggVarsPath (Monte Carlo) \n')
dev=abs(mean(SimPanelTPath.earnings,2)'-AggVarsPath2.earnings.Mean);
fprintf('SimPanel (no GI),   earnings:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.earnings.Mean)))
dev=abs(mean(SimPanelTPath.assets,2)'-AggVarsPath2.assets.Mean);
fprintf('SimPanel (no GI),   assets:     max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.assets.Mean)))
dev=abs(mean(SimPanelTPath.nextassets,2)'-AggVarsPath2.nextassets.Mean);
fprintf('SimPanel (no GI),   nextassets: max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.nextassets.Mean)))
dev=abs(mean(SimPanelTPath.semizval,2)'-AggVarsPath2.semizval.Mean);
fprintf('SimPanel (no GI),   semizval:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.semizval.Mean)))
dev=abs(mean(SimPanelTPath_GI.earnings,2)'-AggVarsPath4.earnings.Mean);
fprintf('SimPanel (with GI), earnings:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.earnings.Mean)))
dev=abs(mean(SimPanelTPath_GI.assets,2)'-AggVarsPath4.assets.Mean);
fprintf('SimPanel (with GI), assets:     max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.assets.Mean)))
dev=abs(mean(SimPanelTPath_GI.nextassets,2)'-AggVarsPath4.nextassets.Mean);
fprintf('SimPanel (with GI), nextassets: max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.nextassets.Mean)))
dev=abs(mean(SimPanelTPath_GI.semizval,2)'-AggVarsPath4.semizval.Mean);
fprintf('SimPanel (with GI), semizval:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.semizval.Mean)))
[AggVarsPath2.earnings.Mean; mean(SimPanelTPath.earnings,2)']
[AggVarsPath2.assets.Mean;   mean(SimPanelTPath.assets,2)']

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
% The DISTRIBUTION comparison is reported as a DIAGNOSTIC rather than a check: grid interpolation
% changes which of the lowest nodes a pile-up sits on, and a_grid_big is cubed, so the spacing down
% there is ~1e-5 (see the nosemiz subcodes). The moments are the claim, and are checked.
% [] absorbs the exogenous states, whatever they are.
AD2=reshape(gather(AgentDistPath2),prod(n_a_big),[],T);
AD4=reshape(gather(AgentDistPath4),prod(n_a_big),[],T);
ADcdf2=cumsum(sum(AD2,2),1); ADcdf4=cumsum(sum(AD4,2),1);
[ADdevmax,ADdevwhere]=max(abs(AgentDistPath2(:)-AgentDistPath4(:)));
AggVarNamesGI=fieldnames(AggVarsPath2);
devmoment=0;
for ff=1:length(AggVarNamesGI)
    m2=AggVarsPath2.(AggVarNamesGI{ff}).Mean; m4=AggVarsPath4.(AggVarNamesGI{ff}).Mean;
    devmoment=max(devmoment,max(abs(m2-m4)./max(abs(m4),1e-12)));
end
fprintf('AgentDist along TPath with/without grid interp, max relative difference in the AggVars means, this should be close to zero: %.3e \n',gather(devmoment))
fprintf('   [diagnostics, not checks: pointwise density difference %.3e at flat index %i of %i; asset \n',gather(ADdevmax),gather(ADdevwhere),numel(AgentDistPath2))
fprintf('    CDF difference %.3e; total mass difference per period %.3e] \n',max(abs(ADcdf2(:)-ADcdf4(:))),max(abs(squeeze(sum(sum(AD2,1),2))-squeeze(sum(sum(AD4,1),2)))))
clear AD2 AD4 ADcdf2 ADcdf4 ADdevmax ADdevwhere AggVarNamesGI devmoment m2 m4
[AggVarsPath2.earnings.Mean; AggVarsPath4.earnings.Mean]
[AggVarsPath2.assets.Mean; AggVarsPath4.assets.Mean]

%% AutoCorrTransProbs along the path
% EvalFnOnTransPath_AutoCorrTransProbs_InfHorz returns Mean and StdDeviation (exact references in
% AggVars and AllStats), AutoCovariance, AutoCorrelation, and optionally a transition probability
% matrix between the values of a FnToEvaluate. AutoCovariance and TransitionProbs are checked against a
% dense kernel built here (the command pushes the measures and never forms the transition matrix, so
% this is its only oracle). Run on the small grid: the dense kernel is (N_a*N_semiz*N_z*N_e)^2 per period.
%
% NOT GUARDED. The semiz support on this command is being built test-first: an error here is the toolkit
% not yet covering this model, not a test defect.
simoptionsAC=simoptions2;
simoptionsAC.whichstats=ones(7,1); % pin, so AllStats and this command are asked for the same things
% assets only: nextassets is the chosen aprime set, whose number of distinct values changes
% along the path, and EvalFnOnTransPath_AutoCorrTransProbs_InfHorz preallocates
% TransitionProbs at a fixed size from the tt=2 slice, so it errors when that count moves.
simoptionsAC.transprobs={'assets'};
% PolicyPath2 was cleared further up, so re-solve the plain path here, on the small grid.
[~,PolicyPathAC]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPathAC=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPathAC, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptionsAC);
AggVarsPathAC=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathAC, PolicyPathAC, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsAC);
AllStatsPathAC=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathAC, PolicyPathAC, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsAC);
CorrTransProbsPathAC=EvalFnOnTransPath_AutoCorrTransProbs_InfHorz(FnsToEvaluate, AgentDistPathAC, PolicyPathAC, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, transpathoptionsbaseline, simoptionsAC);

for fnname=fieldnames(FnsToEvaluate)'
    dev=CorrTransProbsPathAC.(fnname{1}).Mean-AggVarsPathAC.(fnname{1}).Mean;
    fprintf('AutoCorrTransProbs vs AggVars along TPath, this should be zero, Mean %s: %.3e \n',fnname{1},max(abs(dev)))
    dev=CorrTransProbsPathAC.(fnname{1}).StdDeviation-AllStatsPathAC.(fnname{1}).StdDeviation;
    fprintf('AutoCorrTransProbs vs AllStats along TPath, this should be zero, StdDeviation %s: %.3e \n',fnname{1},max(abs(dev)))
end

% searcheffort (d2) can be the same for everyone in some period, and then its std deviation is exactly
% zero and AutoCorrelation is 0/0, which is NaN, and CoreSummary rightly fails a printed NaN. So the
% identity is only checked where it is defined, and otherwise a non-check line says why it was not.
for fnname=fieldnames(FnsToEvaluate)'
    sdthis=CorrTransProbsPathAC.(fnname{1}).StdDeviation(2:T);
    sdlag=CorrTransProbsPathAC.(fnname{1}).StdDeviation(1:T-1);
    if all(sdlag.*sdthis>0)
        dev=CorrTransProbsPathAC.(fnname{1}).AutoCorrelation(1:T-1)-CorrTransProbsPathAC.(fnname{1}).AutoCovariance(1:T-1)./(sdlag.*sdthis);
        fprintf('AutoCorrTransProbs internal, corr=cov/(sd*sd), this should be zero, %s: %.3e \n',fnname{1},max(abs(dev)))
        fprintf('AutoCorrTransProbs internal, max abs AutoCorrelation (should be at most one), %s: %2.8f \n',fnname{1},max(abs(CorrTransProbsPathAC.(fnname{1}).AutoCorrelation(1:T-1))))
    else
        fprintf('AutoCorrTransProbs internal, corr=cov/(sd*sd), %s: not checked, the std deviation is zero in %i of the %i period pairs, where AutoCorrelation is 0/0 \n',fnname{1},sum(~(sdlag.*sdthis>0)),T-1)
    end
end

% The command only fills TransitionProbs in when the set of distinct values is the same in
% consecutive periods, so check the field is there before using it.
for fnname={'assets'}
    if isfield(CorrTransProbsPathAC.(fnname{1}),'TransitionProbs')
        TPac=CorrTransProbsPathAC.(fnname{1}).TransitionProbs;
        fprintf('AutoCorrTransProbs TransitionProbs %s: size %s, smallest entry (should not be negative): %2.10f \n',fnname{1},mat2str(size(TPac)),min(TPac(:)))
        fprintf('AutoCorrTransProbs TransitionProbs %s, row sums should be one, max abs deviation: %2.10f \n',fnname{1},max(abs(sum(TPac,2)-1),[],'all'))
        fprintf('AutoCorrTransProbs TransitionProbs %s, column sums, max abs deviation from one: %2.10f \n',fnname{1},max(abs(sum(TPac,1)-1),[],'all'))
    else
        fprintf('AutoCorrTransProbs TransitionProbs %s: not returned, the distinct values are not the same in consecutive periods \n',fnname{1})
    end
end

% Dense-kernel check of the push:
% T_lag(s,s')=1{a'=aprime_{t-1}(s)} * pi_semiz(semiz,semiz',d2_{t-1}(s)) * pi_z(z,z') * pi_e(e'), built by loops
% from PolicyPathAC(:,...,tt-1), with d2_{t-1}(s) the d2 index in the policy at state s (row length(n_d),
% the last decision row; aprime is row length(n_d)+1). J=dist_{t-1}.*T_lag is the joint distribution of
% (state at t-1, state at t). AutoCovariance(tt-1) is the sum of J.*((x-mu_{t-1})(x-mu_t)'), and
% TransitionProbs(:,:,tt-1) is J summed over the value bins divided by the origin bin's mass (NaN where
% the origin bin has no mass, as the command returns). State s=(a,semiz,z,e) in column-major order, a
% fastest, then semiz, z, e (missing ones have size one here). pi_semiz is built here from the setup's
% SemiExoStateFn, independently of the toolkit, and is constant on this path (probfindjob is a scalar);
% on a path-varying pi_semiz it would be slice tt-1, the transition from tt-1 to tt.
% The semizval AutoCovariance is the line that checks the pi_semiz(.,.,d2) factor of the push, i.e. that
% the agent dist moves semiz with the d2 the agent chose at t-1 and not some other d2.
N_aAC=prod(n_a);
N_semizAC=prod(n_semiz);
N_zAC=prod(n_z);
pi_zAC=gather(pi_z);
N_eAC=1; % no e
pi_eAC=1;
pi_semizAC=zeros(N_semizAC,N_semizAC,n_d2);
for s1=1:N_semizAC
    for s2=1:N_semizAC
        for d2c=1:n_d2
            pi_semizAC(s1,s2,d2c)=CoreInfHorzTPathSetup_SemiExoStateFn(semiz_grid(s1),semiz_grid(s2),d2_grid(d2c),Params.probfindjob,Params.problosejob);
        end
    end
end
N_sAC=N_aAC*N_semizAC*N_zAC*N_eAC;
PolicyPathAC_cpu=reshape(gather(PolicyPathAC),[size(PolicyPathAC,1),N_aAC,N_semizAC,N_zAC,N_eAC,T]);
AgentDistPathAC_cpu=reshape(gather(AgentDistPathAC),[N_sAC,T]);
xa=gpuArray(repmat(gather(a_grid(:)),N_semizAC*N_zAC*N_eAC,1)); % assets
xs=gpuArray(repmat(repelem(gather(semiz_grid(:)),N_aAC,1),N_zAC*N_eAC,1)); % semizval
xz=gpuArray(repmat(repelem(gather(z_grid(:)),N_aAC*N_semizAC,1),N_eAC,1)); % zval
[~,~,idxa]=unique(gather(xa));
Sa=gpuArray(full(sparse(1:N_sAC,idxa,1,N_sAC,max(idxa)))); % bin indicator for assets
dev_cov_a=0; dev_cov_s=0; dev_cov_z=0; dev_cov_e=0; dev_tp=0; nanmismatch=0;
for tt=2:T
    Tlag=zeros(N_sAC,N_sAC);
    for ee=1:N_eAC
        for zz=1:N_zAC
            for ss=1:N_semizAC
                rows=(1:N_aAC)'+N_aAC*(ss-1)+N_aAC*N_semizAC*(zz-1)+N_aAC*N_semizAC*N_zAC*(ee-1);
                ap=reshape(PolicyPathAC_cpu(length(n_d)+1,:,ss,zz,ee,tt-1),[N_aAC,1]); % the aprime index
                d2ind=reshape(PolicyPathAC_cpu(length(n_d),:,ss,zz,ee,tt-1),[N_aAC,1]); % the d2 index
                for ee2=1:N_eAC
                    for zz2=1:N_zAC
                        for ss2=1:N_semizAC
                            cols=ap+N_aAC*(ss2-1)+N_aAC*N_semizAC*(zz2-1)+N_aAC*N_semizAC*N_zAC*(ee2-1);
                            probs=reshape(pi_semizAC(ss,ss2,d2ind),[N_aAC,1])*pi_zAC(zz,zz2)*pi_eAC(ee2);
                            Tlag(rows+(cols-1)*N_sAC)=Tlag(rows+(cols-1)*N_sAC)+probs;
                        end
                    end
                end
            end
        end
    end
    Tlag=gpuArray(Tlag);
    dlag=gpuArray(AgentDistPathAC_cpu(:,tt-1)); dnow=gpuArray(AgentDistPathAC_cpu(:,tt));
    J=dlag.*Tlag;
    mua_lag=sum(dlag.*xa); mua=sum(dnow.*xa);
    mus_lag=sum(dlag.*xs); mus=sum(dnow.*xs);
    dev_cov_a=max(dev_cov_a,abs(gather(CorrTransProbsPathAC.assets.AutoCovariance(tt-1))-gather(sum(sum(J.*((xa-mua_lag)*(xa-mua)'))))));
    dev_cov_s=max(dev_cov_s,abs(gather(CorrTransProbsPathAC.semizval.AutoCovariance(tt-1))-gather(sum(sum(J.*((xs-mus_lag)*(xs-mus)'))))));
    muz_lag=sum(dlag.*xz); muz=sum(dnow.*xz);
    dev_cov_z=max(dev_cov_z,abs(gather(CorrTransProbsPathAC.zval.AutoCovariance(tt-1))-gather(sum(sum(J.*((xz-muz_lag)*(xz-muz)'))))));
    if isfield(CorrTransProbsPathAC.assets,'TransitionProbs')
        P_bf=gather((Sa'*J*Sa)./(Sa'*dlag));
        P=gather(CorrTransProbsPathAC.assets.TransitionProbs(:,:,tt-1));
        dd=abs(P-P_bf); dd=dd(~isnan(dd));
        dev_tp=max(dev_tp,max([dd(:);0])); nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(P_bf(:)));
    end
end
fprintf('AutoCorrTransProbs along TPath, assets AutoCovariance vs dense kernel (max over periods), this should be zero: %.3e \n',dev_cov_a)
fprintf('AutoCorrTransProbs along TPath, semizval AutoCovariance vs dense kernel (max over periods), this should be zero: %.3e \n',dev_cov_s)
fprintf('AutoCorrTransProbs along TPath, zval AutoCovariance vs dense kernel (max over periods), this should be zero: %.3e \n',dev_cov_z)
if isfield(CorrTransProbsPathAC.assets,'TransitionProbs')
    fprintf('AutoCorrTransProbs along TPath, assets TransitionProbs vs binned dense kernel (max over periods), this should be zero: %.3e, NaN-pattern mismatches %i \n',dev_tp,nanmismatch)
end
clear simoptionsAC PolicyPathAC AgentDistPathAC AggVarsPathAC AllStatsPathAC CorrTransProbsPathAC TPac sdthis sdlag dev
clear N_aAC N_semizAC N_zAC N_eAC pi_zAC pi_eAC pi_semizAC N_sAC PolicyPathAC_cpu AgentDistPathAC_cpu xa xs xz xe idxa Sa Tlag dlag dnow J P_bf P dd rows cols ap d2ind probs dev_cov_a dev_cov_s dev_cov_z dev_cov_e dev_tp nanmismatch

%% Plots
fig=figure(figure_c);
subplot(3,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
title('Earnings Mean'); legend('1','2')
subplot(3,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean'); legend('1','2')
subplot(3,1,3); plot(1:1:T,AggVarsPath2.semizval.Mean, 1:1:T,AggVarsPath4.semizval.Mean)
title('Employment rate (semiz Mean)'); legend('1','2')

clear VPath2b VPath4b AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% Constant path -- VPath should equal repmat of stationary V, similarly Policy and AgentDist
% Done three ways: (1) no DC, no GI; (2) with DC, no GI; (3) with DC and GI.
% Both the raw Policy (indexes) and the PolicyVals (the actual d/aprime values) are compared on the
% two no-GI tiers. ON THE GI TIER ONLY PolicyVals IS COMPARED: the Policy INDEXES can legitimately
% differ by ngridinterp+1, because (L1,L2) and (L1+1,L2-(ngridinterp+1)) encode the same aprime when
% the optimum sits exactly on a coarse grid point. PolicyVals is immune to this.
% The stationary references are solved WITHOUT divide-and-conquer (not supported for the stationary semiz
% solve), so tiers (2) and (3) also check the path's DC against the stationary solver's full search.
% T is the last dimension of every path object, so V_final(:) repeated T times is the constant path.
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.sigma=Params.sigma*ones(1,T);

% (1) without divide-and-conquer, without grid interpolation
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions1);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions1);

Vfin_rep=repmat(V_final(:),1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, V: %.3e \n',max(abs(VPath1(:)-Vfin_rep(:))))
Polfin_rep=repmat(Policy_final(:),1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, Policy: %.3e \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
PolicyValsfin=PolicyInd2Val_InfHorz(Policy_final,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyValsfin_rep=repmat(PolicyValsfin(:),1,T); % T is the last dimension of PolicyValsPath
fprintf('Constant TPath (no DC, no GI), this should be zero, PolicyVals: %.3e \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
AD_rep=repmat(AgentDist_initial(:),1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, AgentDist (note: tolerance=1e-6): %.3e \n',max(abs(AgentDistPath1(:)-AD_rep(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1

% (2) with divide-and-conquer, without grid interpolation (same targets as (1))
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions2);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2);

fprintf('Constant TPath (with DC, no GI), this should be zero, V: %.3e \n',max(abs(VPath1(:)-Vfin_rep(:))))
fprintf('Constant TPath (with DC, no GI), this should be zero, Policy: %.3e \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
fprintf('Constant TPath (with DC, no GI), this should be zero, PolicyVals: %.3e \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
fprintf('Constant TPath (with DC, no GI), this should be zero, AgentDist (note: tolerance=1e-6): %.3e \n',max(abs(AgentDistPath1(:)-AD_rep(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1

% (3) with divide-and-conquer, with grid interpolation
% (uses the GI terminal V/Policy, and the stationary dist computed with GI)
AgentDist_initial_GI=StationaryDist_InfHorz(Policy_final_GI,n_d,n_a,n_z,pi_z,simoptions4,Params,[]);
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial_GI, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions4);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4);

Vfin_rep_GI=repmat(V_final_GI(:),1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, V: %.3e \n',max(abs(VPath1(:)-Vfin_rep_GI(:))))
% No raw Policy index comparison on the GI tier (two legal spellings of a point on a coarse node; see above).
PolicyValsfin_GI=PolicyInd2Val_InfHorz(Policy_final_GI,n_d,n_a,n_z,d_grid,a_grid,vfoptions4);
PolicyValsfin_GI_rep=repmat(PolicyValsfin_GI(:),1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, PolicyVals: %.3e \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_GI_rep(:))))
AD_rep_GI=repmat(AgentDist_initial_GI(:),1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, AgentDist (note: tolerance=1e-6): %.3e \n',max(abs(AgentDistPath1(:)-AD_rep_GI(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1
%% General equilibrium: solve the stationary GE first, then run a null-reform transition path
% The only part of this file that exercises the GE machinery itself (the Anderson algorithm, the
% price update rule, the convergence test). Everything else either runs a constant path through
% the non-GE commands, or calls TransitionPath_InfHorz with a dummy GE eqn and maxiter=1, so
% nothing ever updates.
%
% A Cobb-Douglas firm supplies BOTH prices from its first order conditions:
%     r = firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha)) - firmdelta   (MPK less depreciation)
%     w = (1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha))                (MPL)
% so w is a general eqm price here, not a ParamPath entry as it is elsewhere in this file.
% K is aggregate assets and N aggregate effective labour, N=E[z*semiz]: only the employed work.
%
% WELL POSED, AND SO THE CLOSED ECONOMY OF THE z SUBCODES, IN EVERY SEMIZ SUBCODE (including the noz
% noe ones, where the nosemiz sibling has to use a small open economy instead). The nosemiz noz noe model
% has no risk at all, so with beta*(1+r)<1 everyone runs assets down to a=0, K=0, and the closed-economy
% MPK is infinite. Here semiz IS a shock: an employed household loses the job (and drops to uempbenefit)
% with positive probability whatever d2 it chooses, so there is idiosyncratic income risk, a precautionary
% motive, and K(r)>0 rising continuously toward infinity as r approaches 1/beta-1 -- the standard
% Aiyagari setting, with a unique crossing of the firm's capital demand. The one way it could still fail
% here is the top of a_grid_GE binding (5 is the grid maximum), which would cap K artificially; that is
% checked directly below, as the mass on the top grid point.
%
% Both solvers are given the SAME residual-form GeneralEqmEqns, and the path uses GEnewprice=2 (Anderson),
% as in the nosemiz subcodes. firm* names are used because the 2A entrepreneur ReturnFn already uses alpha
% and delta with its own (Kitao) meanings.
%
% The stationary InfHorz semiz solver has no divide-and-conquer, so every STATIONARY solve below uses the
% plain or the GI tier, while the transition paths on the DC+GI tier do use divide-and-conquer.
% simoptions.d_grid has to be the GE block's own d grid.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

simoptions1GE=simoptions1;
simoptions1GE.d_grid=d_grid_GE;

FnsToEvaluateGE.K=@(d2,aprime,a,semiz,z) a;
FnsToEvaluateGE.N=@(d2,aprime,a,semiz,z) z*semiz;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

% THE NO-GI STATIONARY SOLVE IS STARTED FROM THE GI ANSWER, as in the with-d nosemiz subcodes: on coarse
% grids K=E[a] and N are step functions of the prices (here d2 is a discrete choice whether or not there
% is a d1, so N=E[z*semiz] moves in steps as households switch search effort), the general eqm objective is
% piecewise-constant, and fminsearch can stall where it starts. HeteroAgentStationaryEqm_InfHorz takes
% its initial guess from Parameters.r and Parameters.w. The GI tier is GI WITHOUT divide-and-conquer.
vfoptionsGEstart=vfoptions1;
vfoptionsGEstart.gridinterplayer=1;
vfoptionsGEstart.ngridinterp=5;
simoptionsGEstart=simoptions1GE;
simoptionsGEstart.gridinterplayer=1;
simoptionsGEstart.ngridinterp=5;
[p_eqmstart,~]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptionsGEstart, vfoptionsGEstart);
ParamsGEstart=Params;
ParamsGEstart.r=p_eqmstart.r;
ParamsGEstart.w=p_eqmstart.w;
fprintf('Stationary GE: the GI tier is solved first to give the no-GI solve a starting guess, r=%2.8f, w=%2.8f \n',p_eqmstart.r,p_eqmstart.w)
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGEstart, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions1GE, vfoptions1);
fprintf('Stationary GE: r=%2.8f, w=%2.8f \n',p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)
fprintf('Stationary GE, r below 1/beta-1 (precautionary savings, the well-posed side), this should be one: %i \n',p_eqm.r<1/Params.beta-1)

% The path needs the EQUILIBRIUM V_final and initial dist, else the null reform is not exact
ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions1);
AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d_GE,n_a_GE,n_z,pi_z,simoptions1GE,ParamsGE,[]);
% Well-posedness on this grid: the top asset grid point must hold (almost) no mass, else K is being capped
% by the grid rather than chosen. The dist is [n_a_GE,(exogenous states)], so the top point is row n_a_GE.
temp1=reshape(gather(AgentDist_initialGE),prod(n_a_GE),[]);
fprintf('Stationary GE, mass on the top point of a_grid_GE (the grid must not bind), this should be close to zero: %.3e \n',sum(temp1(end,:)))
clear temp1

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=transpathoptionsbaseline;
transpathoptionsGE.maxiter=25; % the no-change solve starts at the answer, so it needs few iterations
transpathoptionsGE.verbose=0;
% GEnewprice=2: Anderson acceleration (see the nosemiz subcodes). howtoupdate columns are {GEcondn name,
% price name, add, factor}; both residuals are (price - its firm FOC), so a POSITIVE residual means the
% price is too high: add=0 for both. Rows in the order of the PricePath fields (r then w).
transpathoptionsGE.GEnewprice=2;
transpathoptionsGE.GEnewprice2.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

% (i) Start from the equilibrium path. The update formulae return the same prices, so the solver
% must leave the path alone.
PricePathGEOut=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE, simoptions1GE, vfoptions1, []);
% The CapitalMarket residual is in units of r and is how far r has to move to clear the market, so the
% ratio below should be about one (a diagnostic, not a counted check; see the nosemiz subcodes).
fprintf('Null-reform GE from the equilibrium path, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathGEOut.r-p_eqm.r))/abs(GEcondns.CapitalMarket))
fprintf('Null-reform GE from the equilibrium path, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathGEOut.r-p_eqm.r)),max(abs(PricePathGEOut.w-p_eqm.w)))

% (ii) Start from a bumped guess: both prices bumped by between 1%% and 5%%, by a different amount
% in every period, over the first two-thirds of the path only (r rising 1->5%%, w falling 5->1%%, so
% the two have different profiles). The solver must pull both back to the equilibrium path.
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+bumpr);
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+bumpw);
transpathoptionsGE.maxiter=100; % the bumped solve has to actually travel, so give it more iterations
PricePathBumpedOut=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE, simoptions1GE, vfoptions1, []);
fprintf('Null-reform GE bumped guess, max initial deviation: r %2.8f, w %2.8f \n',max(abs(PricePathBumped.r-p_eqm.r)),max(abs(PricePathBumped.w-p_eqm.w)))
fprintf('Null-reform GE from a bumped guess, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathBumpedOut.r-p_eqm.r))/abs(GEcondns.CapitalMarket))
fprintf('Null-reform GE from a bumped guess, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathBumpedOut.r-p_eqm.r)),max(abs(PricePathBumpedOut.w-p_eqm.w)))

%% The same general equilibrium tests again, on the divide-and-conquer + grid-interpolation tier
% Everything above runs on the plain tier. This repeat changes ONLY the solution method -- same
% FnsToEvaluateGE, same GeneralEqmEqnsGE, same grids, same howtoupdate, same bump. The stationary solves
% are GI WITHOUT divide-and-conquer (vfoptionsGEstat_GI; the stationary semiz solver has no DC), the
% transition paths are DC+GI (vfoptionsGE_GI).
vfoptionsGEstat_GI=vfoptions1;
vfoptionsGEstat_GI.gridinterplayer=1;
vfoptionsGEstat_GI.ngridinterp=5;
vfoptionsGE_GI=vfoptionsGEstat_GI;
vfoptionsGE_GI.divideandconquer=1;
simoptionsGE_GI=simoptions1GE;
simoptionsGE_GI.gridinterplayer=1;
simoptionsGE_GI.ngridinterp=5;

[p_eqm_GI,GEcondns_GI]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptionsGE_GI, vfoptionsGEstat_GI);
fprintf('Stationary GE (with GI): r=%2.8f, w=%2.8f \n',p_eqm_GI.r,p_eqm_GI.w)
fprintf('Stationary GE conditions (with GI), these should be close to zero: CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns_GI.CapitalMarket,GEcondns_GI.LabourMarket)
% THE TWO TIERS MUST AGREE ON THE EQUILIBRIUM, as a counted check. Not an identity: the two tiers
% discretise differently, so a correct pair still differs by the discretisation gap. The close-to-zero
% bar of 0.1 passes that comfortably and still catches a stall.
fprintf('Stationary GE, the two tiers must agree on r, this should be close to zero: %.3e \n',abs(p_eqm_GI.r-p_eqm.r))
fprintf('Stationary GE, the two tiers must agree on w, this should be close to zero: %.3e \n',abs(p_eqm_GI.w-p_eqm.w))

ParamsGE_GI=Params;
ParamsGE_GI.r=p_eqm_GI.r;
ParamsGE_GI.w=p_eqm_GI.w;
[V_finalGE_GI,Policy_finalGE_GI]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE_GI,DiscountFactorParamNames,[],vfoptionsGEstat_GI);
AgentDist_initialGE_GI=StationaryDist_InfHorz(Policy_finalGE_GI,n_d_GE,n_a_GE,n_z,pi_z,simoptionsGE_GI,ParamsGE_GI,[]);

PricePathGE_GI.r=p_eqm_GI.r*ones(1,T);
PricePathGE_GI.w=p_eqm_GI.w*ones(1,T);

transpathoptionsGE.maxiter=25; % as above: the no-change solve starts at the answer
PricePathGEOut_GI=TransitionPath_InfHorz(PricePathGE_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('Null-reform GE (with DC and GI) from the equilibrium path, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathGEOut_GI.r-p_eqm_GI.r))/abs(GEcondns_GI.CapitalMarket))
fprintf('Null-reform GE (with DC and GI) from the equilibrium path, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathGEOut_GI.r-p_eqm_GI.r)),max(abs(PricePathGEOut_GI.w-p_eqm_GI.w)))

% Same bump profile as the plain tier (bumpr/bumpw/Tbump are reused)
PricePathBumped_GI=PricePathGE_GI;
PricePathBumped_GI.r(1:Tbump)=PricePathGE_GI.r(1:Tbump).*(1+bumpr);
PricePathBumped_GI.w(1:Tbump)=PricePathGE_GI.w(1:Tbump).*(1+bumpw);
transpathoptionsGE.maxiter=100; % as above: the bumped solve has to actually travel
PricePathBumpedOut_GI=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('Null-reform GE (with DC and GI) bumped guess, max initial deviation: r %2.8f, w %2.8f \n',max(abs(PricePathBumped_GI.r-p_eqm_GI.r)),max(abs(PricePathBumped_GI.w-p_eqm_GI.w)))
fprintf('Null-reform GE (with DC and GI) from a bumped guess, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathBumpedOut_GI.r-p_eqm_GI.r))/abs(GEcondns_GI.CapitalMarket))
fprintf('Null-reform GE (with DC and GI) from a bumped guess, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathBumpedOut_GI.r-p_eqm_GI.r)),max(abs(PricePathBumpedOut_GI.w-p_eqm_GI.w)))

%% Does the convergence criterion actually decide, or does maxiter always cut the loop off first?
% Every other GE solve in this file runs to a maxiter it reaches, so none of them ever exercise
% transpathoptions.toleranceGEcondns. Here the cap is generous and Anderson is what gets there quickly.
% (In the nod nosemiz subcodes too; the with-d ones leave it out, and so does the d1 semiz one.)
transpathoptionsGE_conv=transpathoptionsGE;
transpathoptionsGE_conv.maxiter=250;
[PricePathConv,GEcondnPathConv]=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE_conv, simoptionsGE_GI, vfoptionsGE_GI);
% GEcondnPathConv holds the general eqm eqns along the returned path, one row per eqn. This rebuilds
% the exact number the solver's loop tests: scalarize within each period by multiGEcriterion=1 (the
% root-sum-of-squares, with the default weights of one), then take the L-Infinity norm over time.
GEcondnDist=max(sqrt(GEcondnPathConv.CapitalMarket.^2+GEcondnPathConv.LabourMarket.^2));
fprintf('Convergence criterion, GE condn distance of the returned path: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',GEcondnDist,1e-4)

% The same solve with one more iteration allowed. If the loop stopped because the general eqm
% conditions got small then the extra iteration is never taken and the two paths are identical.
transpathoptionsGE_conv.maxiter=251;
PricePathConv_plus1=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE_conv, simoptionsGE_GI, vfoptionsGE_GI);
fprintf('Convergence criterion decided rather than maxiter, one more iteration changes nothing, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathConv_plus1.r-PricePathConv.r)),max(abs(PricePathConv_plus1.w-PricePathConv.w)))
clear transpathoptionsGE_conv PricePathConv PricePathConv_plus1 GEcondnPathConv GEcondnDist

clear FnsToEvaluateGE vfoptionsGEstat_GI vfoptionsGE_GI simoptionsGE_GI simoptions1GE vfoptionsGEstart simoptionsGEstart p_eqmstart ParamsGEstart p_eqm_GI GEcondns_GI ParamsGE_GI V_finalGE_GI Policy_finalGE_GI AgentDist_initialGE_GI PricePathGE_GI PricePathGEOut_GI PricePathBumped_GI PricePathBumpedOut_GI GeneralEqmEqnsGE heteroagentoptionsGE p_eqm GEcondns ParamsGE V_finalGE Policy_finalGE AgentDist_initialGE PricePathGE ParamPathGE transpathoptionsGE PricePathGEOut PricePathBumped PricePathBumpedOut bumpr bumpw Tbump

%% Run the GE transition path, but with transpathoptions.maxiter=1, so it ends after one iteration -- just a shape-check (the core is tested elsewhere)
transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(earnings) 0;
% factor=0 means the dummy eqn never moves the price: this only tests that the returned objects have the right shapes.
transpathoptions.GEnewprice=2;
transpathoptions.GEnewprice2.howtoupdate={'dummy','r',0,0};

PricePath2=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);

% Big grid + DC
PricePath3A=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, n_d, n_a_big, n_z, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions2, vfoptions2, []);

% Big grid + DC + GI
PricePath3B=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_big_GI, AgentDist_initial_big, n_d, n_a_big, n_z, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions4, vfoptions4, []);

fprintf('One iter of TPath, big-grid with/without GI, this should be close to zero: %.3e \n',max(abs(PricePath3A.r-PricePath3B.r)))

%%
output=struct();

end
