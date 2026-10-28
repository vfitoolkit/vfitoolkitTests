function output=CoreInfHorzTPath_nod_noz_e_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c)
% without d, without z, with e (fig 5).
% REGENERATED 2026-10-27 from the live CoreInfHorzTPath_nod_z_noe_nosemiz (fig 3), adapted to
% this model. The six dark doPart(1) subcodes had fallen behind the two live z ones (tolerances,
% the deleted spelling-unsafe GI Policy-index comparisons, the FHorz cross-checks including the
% path-varying shock block, the stationary-GE two-tier invariant wired as a real check, the
% convergence-criterion block),
% so each was rebuilt from its live sibling rather than patched. Support for no z with an iid e on
% the InfHorz transition path commands is being built in the toolkit test-first, so this subcode
% is written against the intended interface: "no z" is n_z=0, z_grid=[], pi_z=[], passed
% positionally as in the stationary CoreInfHorzTests noz subcodes, and the iid e is
% vfoptions/simoptions .n_e, .e_grid, .pi_e, as in the stationary e subcodes.

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
n_d=0; d_grid=[];
n_d_GE=0; d_grid_GE=[]; % nod: no decision variable in the GE block either
n_z=0; z_grid=[]; pi_z=[]; % noz: no markov z (the inputs stay positional, as empty)
% iid shocks
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(aprime,a,e) a;
FnsToEvaluate.earnings=@(aprime,a,e,w) w*e;
FnsToEvaluate.nextassets=@(aprime,a,e) aprime; % aprime-dependent, so it tests the policy decode directly
FnsToEvaluate.eval=@(aprime,a,e) e; % the iid shock itself (fixed except in the path-varying e block, where its mean/std deviation must move)


%% Period-0 VFI: gives the final-step V/Policy (used as both V_final for TPath and the steady state to compare against)
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

% With an iid e and no z: lowmemory=1 loops over e, and lowmemory=2 (the FHorz convention: loop over z
% and e) has only e to loop over, so it should take the same route as lowmemory=1. The lowmemory=2
% lines below are kept to check that the no-z case does not break lowmemory=2, not as separate evidence.
% lowmemory should give same answer
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

% lowmemory 2 should give same answer
vfoptions1.lowmemory=2;
[VPath1C,PolicyPath1C]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=2;
[VPath2C,PolicyPath2C]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
vfoptions2.lowmemory=0;

fprintf('low memory 2, this should be zero: %.3e \n',max(abs(VPath1(:)-VPath1C(:))))
fprintf('low memory 2, this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath1C(:))))
fprintf('low memory 2, this should be zero: %.3e \n',max(abs(VPath2(:)-VPath2C(:))))
fprintf('low memory 2, this should be zero: %.3e \n',max(abs(PolicyPath2(:)-PolicyPath2C(:))))

%% ValueFnFromPolicy along the path: recomputing V from the policy must reproduce VPath
% ValueFnFromPolicyOnTransPath_InfHorz evaluates the Bellman at the given PolicyPath instead of
% maximising, so it must return exactly the VPath that produced that policy. (Divide-and-conquer
% is irrelevant here -- there is no maximisation -- so this covers the no-GI tier as a whole.)
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
N_jFH=T-1;
ParamsFH=Params;
ParamsFH.r=PricePath.r(1:N_jFH);
ParamsFH.w=ParamPath.w(1:N_jFH);

% vfoptions2 carries n_e, e_grid and pi_e, so the FHorz side solves with the same iid e. pi_e is
% [n_e,1], which ExogShockSetup_FHorz broadcasts to the N_jFH+1 columns that V_Jplus1 needs (the last
% column is the e of the V_Jplus1 period). The timing is set out at the path-varying e block below.
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

for fnname={'assets','earnings','nextassets','eval'}
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AggVars %s Mean: %.3e \n',fnname{1},max(abs(dev(:))))
end
% AllStats: every field the two commands have in common (periods/ages are the trailing dim in both)
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AllStats assets %s: %.3e \n',statname{1},max(abs(dev(:))))
end

clear V_FH Policy_FH ParamsFH vfoptionsFH temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

%% Cross-check against an FHorz solve, with a shock process that varies along the path
% The path-varying z part is skipped for this subcode: there is no markov z here, so
% transpathoptions.zpathtrivial has nothing to vary. The iid e does vary, in the next block.
fprintf('TPath vs FHorz (path-varying z): skipped, this model has no markov z \n')

%% Cross-check against an FHorz solve, with an iid e process that varies along the path
% The iid-e counterpart of the path-varying z block: e_grid and pi_e are period dependent, given to the
% transition path commands as [n_e,T] in vfoptions/simoptions. This is the bank's only coverage of
% transpathoptions.epathtrivial=0 on these commands against an independent solver (doPart(8) varies
% only e_grid, and checks an EiidShockFn against the same explicit grid, never pi_e and never against
% FHorz), and it is checked the same way as the z block, against an FHorz solve of the identical problem.
%
% TIMING ("Reading A"): pi_e(:,t) is the distribution of the e realized IN period t. On the path, V_t
% takes its expectation over e_{t+1} with pi_e_T(:,t+1), and the agent dist from t to t+1 draws e_{t+1}
% from pi_e_T(:,t+1). FHorz uses the same convention: pi_e_J(:,j) is the distribution of e at age j,
% and the age-j value fn uses pi_e_J(:,j+1). So TPath period t is FHorz age j=t, exactly as for the
% price and param paths. With V_Jplus1 the FHorz side needs N_jFH+1=T columns of pi_e, i.e.
% pi_e_T(:,1:T), the last being the distribution of e in the V_Jplus1 period; e_grid has one column
% per age, e_grid_T(:,1:N_jFH). Column 1 of pi_e is read by neither side: the period-1 e is whatever
% AgentDist_initial says.
%
% NOT SYMMETRIC IN TIME, on purpose. The weight on the highest e point rises monotonically along the
% whole path, and the e values themselves drift up, so an implementation that reads column t where t+1
% is meant (or t+1 where t is meant) solves a visibly different problem and the comparison catches it.
% A pi_e path symmetric about T/2, or one that only moved near the end, could let an off-by-one through.
%
% V_final is the stationary V under the ORIGINAL e process; as in the z block, both sides are handed the
% identical array and index it by e-index, so the check stays exact.
lambdae=0.5*linspace(0,1,T); % weight moved onto the highest e point, rising monotonically from 0 to 1/2
pi_e_T=(1-lambdae).*vfoptionsbaseline.pi_e+lambdae.*[zeros(vfoptionsbaseline.n_e-1,1);1]; % [n_e,T], each column sums to one
e_grid_T=vfoptionsbaseline.e_grid.*linspace(0.9,1.1,T); % [n_e,T], the e values drift up along the path

vfoptionsE=vfoptions2;
vfoptionsE.e_grid=e_grid_T;
vfoptionsE.pi_e=pi_e_T;
[VPathE,PolicyPathE]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);

ParamsFHE=Params;
ParamsFHE.r=PricePath.r(1:N_jFH);
ParamsFHE.w=ParamPath.w(1:N_jFH);
vfoptionsFHE=vfoptions2;
vfoptionsFHE.e_grid=e_grid_T(:,1:N_jFH); % one column per age
vfoptionsFHE.pi_e=pi_e_T(:,1:N_jFH+1); % N_jFH+1 columns: the last is the distribution of e in the V_Jplus1 period
vfoptionsFHE.V_Jplus1=V_final;
[V_FHE,Policy_FHE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFHE,DiscountFactorParamNames,[],vfoptionsFHE);

temp1=reshape(VPathE,[],T); temp2=reshape(V_FHE,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying e), this should be zero, V: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyPathE,[],T); temp2=reshape(Policy_FHE,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying e), this should be zero, Policy: %.3e \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPathE,n_d,n_a,n_z,T,d_grid,a_grid,vfoptionsE),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FHE,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptionsFHE),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying e), this should be zero, PolicyVals: %.3e \n',max(abs(dev(:))))

% Same age-weight handling as the fixed-process block above
ParamsFHE.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFHE.mewjFH(end)=1-sum(ParamsFHE.mewjFH(1:end-1));
simoptionsE=simoptions2;
simoptionsE.e_grid=e_grid_T;
simoptionsE.pi_e=pi_e_T;
simoptionsE.npoints=100;
simoptionsE.nquantiles=20;
simoptionsE.whichstats=ones(7,1);
simoptionsFHE=simoptions2;
simoptionsFHE.e_grid=e_grid_T(:,1:N_jFH);
simoptionsFHE.pi_e=pi_e_T(:,1:N_jFH); % the agent dist has no V_Jplus1 period, so one column per age
simoptionsFHE.npoints=100;
simoptionsFHE.nquantiles=20;
simoptionsFHE.whichstats=ones(7,1);

AgentDistPathE=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPathE, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptionsE);
StatDist_FHE=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FHE,n_d,n_a,n_z,N_jFH,pi_z,ParamsFHE,simoptionsFHE);
temp1=reshape(AgentDistPathE,[],T); temp2=reshape(StatDist_FHE,[],N_jFH);
temp2=temp2./sum(temp2,1); % undo the age weighting, so each age carries mass one like each TPath period
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (path-varying e), this should be zero, AgentDist: %.3e \n',max(abs(dev(:))))

AggVarsPathE=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathE, PolicyPathE, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsE);
AllStatsPathE=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathE, PolicyPathE, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptionsE);
AgeCondStats_FHE=LifeCycleProfiles_FHorz_Case1(StatDist_FHE,Policy_FHE,FnsToEvaluate,ParamsFHE,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFHE);

for fnname=fieldnames(FnsToEvaluate)'
    dev=AggVarsPathE.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FHE.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (path-varying e), this should be zero, AggVars %s Mean: %.3e \n',fnname{1},max(abs(dev(:))))
end
% The e shock FnsToEvaluate, which is what this block is really about: with the e process moving, the
% mean and the std deviation of e both move, and the two solvers must agree on both.
for fnname={'eval'}
    for statname={'Mean','StdDeviation'}
        temp1=AllStatsPathE.(fnname{1}).(statname{1}); temp2=AgeCondStats_FHE.(fnname{1}).(statname{1});
        dev=temp1(:,1:N_jFH)-temp2;
        fprintf('TPath vs FHorz (path-varying e), this should be zero, %s %s: %.3e \n',fnname{1},statname{1},max(abs(dev(:))))
    end
end
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathE.assets.(statname{1}); temp2=AgeCondStats_FHE.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (path-varying e), this should be zero, AllStats assets %s: %.3e \n',statname{1},max(abs(dev(:))))
end

% Non-vacuity. If the [n_e,T] e_grid/pi_e were silently ignored (or read as the period-1 column
% throughout), everything above would still agree -- both sides would just be solving the fixed-process
% model -- and the block would prove nothing. Here the weight on the highest e rises along the path and
% the e values drift up, so the e shock's own mean and std deviation must visibly move. Printed in the
% inverted form, which CoreSummary counts and fails if the range is not well above zero.
for fnname={'eval'}
    for statname={'Mean','StdDeviation'}
        temp1=AllStatsPathE.(fnname{1}).(statname{1});
        fprintf('Path-varying e bit, range across periods of %s %s, this should NOT be zero: %.3e \n',fnname{1},statname{1},max(temp1(:))-min(temp1(:)))
    end
end

clear lambdae pi_e_T e_grid_T vfoptionsE VPathE PolicyPathE ParamsFHE vfoptionsFHE V_FHE Policy_FHE simoptionsE simoptionsFHE AgentDistPathE StatDist_FHE AggVarsPathE AllStatsPathE AgeCondStats_FHE temp1 temp2 dev

%%
clear VPath1 VPath2 PolicyPath1 PolicyPath2 VPath1B VPath2B PolicyPath1B PolicyPath2B VPath1C VPath2C PolicyPath1C PolicyPath2C

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

% With an iid e and no z: lowmemory=1 loops over e, and lowmemory=2 (the FHorz convention: loop over z
% and e) has only e to loop over, so it should take the same route as lowmemory=1. The lowmemory=2
% lines below are kept to check that the no-z case does not break lowmemory=2, not as separate evidence.
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

% lowmemory 2 should give same answer
vfoptions3.lowmemory=2;
[VPath3C,PolicyPath3C]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=2;
[VPath4C,PolicyPath4C]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
vfoptions4.lowmemory=0;

fprintf('low memory 2, this should be zero: %.3e \n',max(abs(VPath3(:)-VPath3C(:))))
fprintf('low memory 2, this should be zero: %.3e \n',max(abs(PolicyPath3(:)-PolicyPath3C(:))))
fprintf('low memory 2, this should be zero: %.3e \n',max(abs(VPath4(:)-VPath4C(:))))
fprintf('low memory 2, this should be zero: %.3e \n',max(abs(PolicyPath4(:)-PolicyPath4C(:))))

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
% way. PolicyVals is immune to it and is the one that has to be zero. See the note at the check.
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
temp1=reshape(PolicyPath4,[],T); temp2=reshape(Policy_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
% DELETED, same reason as the other raw Policy index comparisons under the interpolation layer, and
% this one PASSED -- the two spellings coincided here. A wrong check that passes is worse than a
% wrong check that fails, because nothing ever prompts you to look at it. The PolicyVals comparison
% immediately below is the real check and is spelling-proof by construction.
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath4,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions4),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, PolicyVals: %.3e \n',max(abs(dev(:))))

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

for fnname={'assets','earnings','nextassets','eval'}
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
clear VPath3 VPath4 PolicyPath3 PolicyPath4 VPath3B VPath4B PolicyPath3B PolicyPath4B VPath3C VPath4C PolicyPath3C PolicyPath4C

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
% SimPanelValues_TransPath_InfHorz has to strip (l_daprime=size(PolicyPath,1)-2*gridinterplayer).
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
dev=abs(mean(SimPanelTPath_GI.earnings,2)'-AggVarsPath4.earnings.Mean);
fprintf('SimPanel (with GI), earnings:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.earnings.Mean)))
dev=abs(mean(SimPanelTPath_GI.assets,2)'-AggVarsPath4.assets.Mean);
fprintf('SimPanel (with GI), assets:     max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.assets.Mean)))
dev=abs(mean(SimPanelTPath_GI.nextassets,2)'-AggVarsPath4.nextassets.Mean);
fprintf('SimPanel (with GI), nextassets: max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.nextassets.Mean)))
[AggVarsPath2.earnings.Mean; mean(SimPanelTPath.earnings,2)']
[AggVarsPath2.assets.Mean;   mean(SimPanelTPath.assets,2)']

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
% The heading above promises the same MOMENTS, and AggVarsPath2/AggVarsPath4 are already computed
% further up -- until now compared only against their own SimPanels, never against each other. So
% compare them, which is the claim actually being made.
%
% The DISTRIBUTION comparison is reported below as a DIAGNOSTIC rather than a check. Grid
% interpolation changes which of the lowest nodes a pile-up sits on, and a_grid_big is cubed, so the
% spacing down there is ~1e-5: a large difference in mass-by-node is a negligible difference in
% assets. Measured: the max sits at a-index 29 of 1001 (a=1.1e-4), and the cumulative difference is
% LARGER than the pointwise one, which rules out the adjacent-node-splitting explanation -- the mass
% genuinely sits on different bottom nodes. No fix removes that; it is what GI does.
% [] absorbs the exogenous states, whatever they are: the z subcodes write prod(n_z) there, which is 0
% with no z, and with an e the dist has a further dimension, so that form would error in this subcode.
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
% EvalFnOnTransPath_AutoCorrTransProbs_InfHorz is the one TPath evaluation command the rest of this
% bank never touches. It returns Mean and StdDeviation (which it has to compute anyway as
% intermediates to the correlation), AutoCovariance, AutoCorrelation, and optionally a transition
% probability matrix between the values of a FnToEvaluate.
%
% Mean and StdDeviation have exact references: the same two objects come out of
% EvalFnOnTransPath_AggVars_InfHorz and EvalFnOnTransPath_AllStats_InfHorz, so three different
% commands must agree to machine precision given the same agent distribution and policy path.
%
% AutoCovariance and AutoCorrelation have no independent reference here, so they are only checked
% for internal consistency: the command sets Corr=Covar/(stddev_lag*stddev), and a correlation
% cannot exceed one in absolute value. Both are written into entries 1:T-1 (entry tt-1 pairs
% period tt with period tt-1) with the last entry left as zero, which is why the checks stop at T-1.
%
% TransitionProbs is asked for on assets and nextassets. It is meant to be a transition matrix
% between consecutive periods' values, so each row should sum to one. Run on the small grid: the
% command builds the full (N_a*N_z)-by-(N_a*N_z) transition matrix, and TransitionProbs is
% (number of distinct values)^2-by-(T-1), both of which grow fast in n_a.
%
% NOT GUARDED HERE. The z subcodes skip this block for the noz and e models, because the command
% had no branch for N_z=0 and no e support. That support is being built test-first, so in this
% subcode the block runs: an error here is the toolkit not yet covering this model, not a test defect.
simoptionsAC=simoptions2;
simoptionsAC.whichstats=ones(7,1); % pin, so AllStats and this command are asked for the same things
% assets only: nextassets is the chosen aprime set, whose number of distinct values changes
% along the path, and EvalFnOnTransPath_AutoCorrTransProbs_InfHorz preallocates
% TransitionProbs at a fixed size from the tt=2 slice, so it errors when that count moves.
simoptionsAC.transprobs={'assets'};
% PolicyPath2 was cleared further up, so re-solve the plain path here. Deliberately the small
% grid, not n_a_big: the command builds the full (N_a*N_z)-by-(N_a*N_z) transition matrix, and
% TransitionProbs is (number of distinct values)^2-by-(T-1), both of which grow fast in n_a.
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

for fnname=fieldnames(FnsToEvaluate)'
    sdthis=CorrTransProbsPathAC.(fnname{1}).StdDeviation(2:T);
    sdlag=CorrTransProbsPathAC.(fnname{1}).StdDeviation(1:T-1);
    dev=CorrTransProbsPathAC.(fnname{1}).AutoCorrelation(1:T-1)-CorrTransProbsPathAC.(fnname{1}).AutoCovariance(1:T-1)./(sdlag.*sdthis);
    fprintf('AutoCorrTransProbs internal, corr=cov/(sd*sd), this should be zero, %s: %.3e \n',fnname{1},max(abs(dev)))
    fprintf('AutoCorrTransProbs internal, max abs AutoCorrelation (should be at most one), %s: %2.8f \n',fnname{1},max(abs(CorrTransProbsPathAC.(fnname{1}).AutoCorrelation(1:T-1))))
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
clear simoptionsAC PolicyPathAC AgentDistPathAC AggVarsPathAC AllStatsPathAC CorrTransProbsPathAC TPac sdthis sdlag dev

%% Plots
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
title('Earnings Mean'); legend('1','2')
subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean'); legend('1','2')

clear VPath2b VPath4b AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% Constant path -- VPath should equal repmat of stationary V, similarly Policy and AgentDist
% Done three ways: (1) no DC, no GI; (2) with DC, no GI; (3) with DC and GI.
% Both the raw Policy (indexes) and the PolicyVals (the actual d/aprime values) are compared on the
% two no-GI tiers. ON THE GI TIER ONLY PolicyVals IS COMPARED: the Policy INDEXES can legitimately
% differ by ngridinterp+1, because (L1,L2) and (L1+1,L2-(ngridinterp+1)) encode the same aprime when
% the optimum sits exactly on a coarse grid point. PolicyVals is immune to this.
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.sigma=Params.sigma*ones(1,T);

% (1) without divide-and-conquer, without grid interpolation
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions1);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions1);

Vfin_rep=repmat(V_final,1,1,1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, V: %.3e \n',max(abs(VPath1(:)-Vfin_rep(:))))
Polfin_rep=repmat(Policy_final,1,1,1,1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, Policy: %.3e \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
PolicyValsfin=PolicyInd2Val_InfHorz(Policy_final,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyValsfin_rep=repmat(PolicyValsfin(:),1,T); % T is the last dimension of PolicyValsPath
fprintf('Constant TPath (no DC, no GI), this should be zero, PolicyVals: %.3e \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
AD_rep=repmat(AgentDist_initial,1,1,1,T);
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

Vfin_rep_GI=repmat(V_final_GI,1,1,1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, V: %.3e \n',max(abs(VPath1(:)-Vfin_rep_GI(:))))
% Polfin_rep_GI is gone with the raw Policy index comparison it was the only input to.
% THE RAW POLICY INDEX COMPARISON IS DELETED, not fixed: under the interpolation layer a point on a
% coarse node has TWO legal spellings, (P1,L2=1) and (P1-1,L2=1+ngridinterp+1), so a raw channel
% difference of exactly 1+ngridinterp is the same aprime named two ways. That is what this reported --
% 6.000e+00 at ngridinterp=5 -- while the PolicyVals comparison beside it read 4e-16 or exactly zero.
% The values agreeing to machine precision IS the check; the index comparison was testing the
% encoding, and testing it wrongly. See the header of the LS1_GI1 raws for why neither spelling is
% canonical.
PolicyValsfin_GI=PolicyInd2Val_InfHorz(Policy_final_GI,n_d,n_a,n_z,d_grid,a_grid,vfoptions4);
PolicyValsfin_GI_rep=repmat(PolicyValsfin_GI(:),1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, PolicyVals: %.3e \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_GI_rep(:))))
AD_rep_GI=repmat(AgentDist_initial_GI,1,1,1,T);
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
% K is aggregate assets and N aggregate labour supply.
%
% Both solvers are given the SAME residual-form GeneralEqmEqns. That works because the path uses
% transpathoptions.GEnewprice=2 rather than the default 1: GEnewprice=1 treats the eqns as price
% UPDATING FORMULAE (updatePricePathNew_TPath_tt sets PricePathNew_tt=p_i directly), whereas
% GEnewprice=2 treats them as residuals and updates via howtoupdate, which is the same convention
% HeteroAgentStationaryEqm_InfHorz uses. One definition of the firm, not two.
%
% firm* names are used because the 2A entrepreneur ReturnFn already uses alpha and delta with its
% own (Kitao) meanings.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,e) a;
FnsToEvaluateGE.N=@(aprime,a,e) e;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions1, vfoptions1);
fprintf('Stationary GE: r=%2.8f, w=%2.8f \n',p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)

% The path needs the EQUILIBRIUM V_final and initial dist, else the null reform is not exact
ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions1);
AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d_GE,n_a_GE,n_z,pi_z,simoptions1,ParamsGE,[]);

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=transpathoptionsbaseline;
transpathoptionsGE.maxiter=25; % the no-change solve starts at the answer, so it needs few iterations
transpathoptionsGE.verbose=0;
% GEnewprice=2: Anderson acceleration. The base fixed-point map is the shooting one -- treat the
% GeneralEqmEqns as residuals and update each price by a damped fraction of its own residual -- and
% Anderson combines the last few iterates of that map to take a quasi-Newton-like step without any
% derivatives. It is the best of the three algorithms on these models, which is why this bank uses
% it throughout; the ALGORITHMS are compared in the with_CoreInfHorzTPathAlgoTests sub-bank, not
% here. howtoupdate columns are {GEcondn name, price name, add, factor}, and
% updatePricePathNew_TPath_tt does new = old + add*factor*residual - (1-add)*factor*residual.
% Both residuals are of the form (price - its firm FOC), so a POSITIVE residual means the price is
% too high and must come down: hence add=0 (subtract) for both. The rows have to be in the same
% order as the PricePath fields (r then w), because add and factor are applied to the prices
% positionally by row. setupGEnewprice3_shooting has a reorder meant to lift that requirement,
% but it matches on the wrong column of howtoupdate and so never fires.
transpathoptionsGE.GEnewprice=2;
transpathoptionsGE.GEnewprice2.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

% (i) Start from the equilibrium path. The update formulae return the same prices, so the solver
% must leave the path alone.
PricePathGEOut=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE, simoptions1, vfoptions1, []);
% The CapitalMarket eqn is r-(MPK-delta), so its residual is already in units of r, and is exactly how
% far r has to move to clear the market. The path moves r there, so the deviation from p_eqm.r equals
% that residual and the ratio below should be about one. A ratio near one says the path is right and
% p_eqm is off by the stationary solver's own residual; a ratio far from one says the path is wrong.
% w has no such identity (it moves only because K moves, with a model-dependent elasticity), so w is
% reported as a diagnostic rather than checked.
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
PricePathBumpedOut=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE, simoptions1, vfoptions1, []);
fprintf('Null-reform GE bumped guess, max initial deviation: r %2.8f, w %2.8f \n',max(abs(PricePathBumped.r-p_eqm.r)),max(abs(PricePathBumped.w-p_eqm.w)))
% Same ratio, from a different starting point: two starts landing on the same r is what says the path
% has converged rather than been cut off by maxiter.
fprintf('Null-reform GE from a bumped guess, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathBumpedOut.r-p_eqm.r))/abs(GEcondns.CapitalMarket))
fprintf('Null-reform GE from a bumped guess, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathBumpedOut.r-p_eqm.r)),max(abs(PricePathBumpedOut.w-p_eqm.w)))

%% The same general equilibrium tests again, on the divide-and-conquer + grid-interpolation tier
% Everything above runs on the plain tier (vfoptions1 is a bare struct). This repeat changes ONLY
% the solution method -- same FnsToEvaluateGE, same GeneralEqmEqnsGE, same grids, same
% howtoupdate, same bump -- so any difference is attributable to DC+GI alone. Two things come out
% of it: whether GI smooths the GE objective enough to tighten the residuals (interpolating the
% asset policy off the coarse grid does far more for that than adding grid points), and a
% cross-check that the two tiers agree on r_eqm and w_eqm.
vfoptionsGE_GI=vfoptions1;
vfoptionsGE_GI.divideandconquer=1;
vfoptionsGE_GI.gridinterplayer=1;
vfoptionsGE_GI.ngridinterp=5;
simoptionsGE_GI=simoptions1;
simoptionsGE_GI.gridinterplayer=1;
simoptionsGE_GI.ngridinterp=5;

[p_eqm_GI,GEcondns_GI]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptionsGE_GI, vfoptionsGE_GI);
fprintf('Stationary GE (with DC and GI): r=%2.8f, w=%2.8f \n',p_eqm_GI.r,p_eqm_GI.w)
fprintf('Stationary GE conditions (with DC and GI), these should be close to zero: CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns_GI.CapitalMarket,GEcondns_GI.LabourMarket)
% THE TWO TIERS MUST AGREE ON THE EQUILIBRIUM, and this is now a CHECK and not a print. It used to
% read 'the two tiers should agree', which is not one of the four phrases CoreSummary counts, so the
% number was computed and printed on every run and never read by anything. That is a worse failure
% than a missing test: the invariant is stated in the code and measured every run, and only the
% enforcement was absent, so a grep for the counted phrases finds nothing wrong.
% FIVE OF THIS BANK'S TEN SUBCODES HAVE A d, and in those the no-GI stationary solve stalls rather
% than converging: K=E[a] and N=E[d*z] are step functions of the prices on coarse grids, so the
% general eqm objective is piecewise-constant and fminsearch stops where it starts. On the d_z model
% that stall is at w=0.68066 against the GI tier's w=0.43767, so this line has been printing a
% difference of about 0.24 beneath a sentence asserting it should be zero.
% The fix, measured in CoreInfHorzTPath_localsearch_shooting: solve WITH GI first, then re-solve the
% no-GI tier STARTING FROM that answer. HeteroAgentStationaryEqm_InfHorz takes its initial guess from
% Parameters.r and Parameters.w. That moved the d_z model to w=0.44304 at a residual 90x smaller.
% NOT AN IDENTITY. The two tiers discretise differently, so a correct pair still differs by the
% discretisation gap -- about 5e-03 in w on the with-d model. The close-to-zero bar of 0.1 passes
% that comfortably and still catches a stall an order of magnitude above it.
fprintf('Stationary GE, the two tiers must agree on r, this should be close to zero: %.3e \n',abs(p_eqm_GI.r-p_eqm.r))
fprintf('Stationary GE, the two tiers must agree on w, this should be close to zero: %.3e \n',abs(p_eqm_GI.w-p_eqm.w))

ParamsGE_GI=Params;
ParamsGE_GI.r=p_eqm_GI.r;
ParamsGE_GI.w=p_eqm_GI.w;
[V_finalGE_GI,Policy_finalGE_GI]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE_GI,DiscountFactorParamNames,[],vfoptionsGE_GI);
AgentDist_initialGE_GI=StationaryDist_InfHorz(Policy_finalGE_GI,n_d_GE,n_a_GE,n_z,pi_z,simoptionsGE_GI,ParamsGE_GI,[]);

PricePathGE_GI.r=p_eqm_GI.r*ones(1,T);
PricePathGE_GI.w=p_eqm_GI.w*ones(1,T);

transpathoptionsGE.maxiter=25; % as above: the no-change solve starts at the answer
PricePathGEOut_GI=TransitionPath_InfHorz(PricePathGE_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('Null-reform GE (with DC and GI) from the equilibrium path, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathGEOut_GI.r-p_eqm_GI.r))/abs(GEcondns_GI.CapitalMarket))
fprintf('Null-reform GE (with DC and GI) from the equilibrium path, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathGEOut_GI.r-p_eqm_GI.r)),max(abs(PricePathGEOut_GI.w-p_eqm_GI.w)))

% Same bump profile as the plain tier (bumpr/bumpw/Tbump are reused, so the two tiers face an
% identical perturbation and the results are directly comparable)
PricePathBumped_GI=PricePathGE_GI;
PricePathBumped_GI.r(1:Tbump)=PricePathGE_GI.r(1:Tbump).*(1+bumpr);
PricePathBumped_GI.w(1:Tbump)=PricePathGE_GI.w(1:Tbump).*(1+bumpw);
transpathoptionsGE.maxiter=100; % as above: the bumped solve has to actually travel
PricePathBumpedOut_GI=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('Null-reform GE (with DC and GI) bumped guess, max initial deviation: r %2.8f, w %2.8f \n',max(abs(PricePathBumped_GI.r-p_eqm_GI.r)),max(abs(PricePathBumped_GI.w-p_eqm_GI.w)))
fprintf('Null-reform GE (with DC and GI) from a bumped guess, deviation in r over the stationary CapitalMarket residual, this should be about one: %2.4f \n',max(abs(PricePathBumpedOut_GI.r-p_eqm_GI.r))/abs(GEcondns_GI.CapitalMarket))
fprintf('Null-reform GE (with DC and GI) from a bumped guess, deviation from p_eqm, r: %2.10f, w: %2.10f \n',max(abs(PricePathBumpedOut_GI.r-p_eqm_GI.r)),max(abs(PricePathBumpedOut_GI.w-p_eqm_GI.w)))

%% Does the convergence criterion actually decide, or does maxiter always cut the loop off first?
% Every other GE solve in this bank runs to a maxiter it reaches, so none of them ever exercise
% transpathoptions.toleranceGEcondns. Here the cap is generous and Anderson is what gets there
% quickly. This used to lean on GEnewprice=3's additionalfactor ramp to converge inside the cap;
% that option is shooting-only and now lives in the with_CoreInfHorzTPathAlgoTests sub-bank.
transpathoptionsGE_conv=transpathoptionsGE;
transpathoptionsGE_conv.maxiter=250;
[PricePathConv,GEcondnPathConv]=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, transpathoptionsGE_conv, simoptionsGE_GI, vfoptionsGE_GI);
% GEcondnPathConv holds the general eqm eqns along the returned path, one row per eqn. This rebuilds
% the exact number the solver's loop tests: scalarize within each period by multiGEcriterion=1 (the
% root-sum-of-squares, with the default weights of one), then take the L-Infinity norm over time.
GEcondnDist=max(sqrt(GEcondnPathConv.CapitalMarket.^2+GEcondnPathConv.LabourMarket.^2));
fprintf('Convergence criterion, GE condn distance of the returned path: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',GEcondnDist,1e-4)

% The same solve with one more iteration allowed. If the loop stopped because the general eqm
% conditions got small then the extra iteration is never taken and the two paths are identical; if it
% stopped because it ran out of iterations then that extra one moves the prices. Far cheaper than
% raising the cap by a large factor, and just as decisive.
transpathoptionsGE_conv.maxiter=251;
PricePathConv_plus1=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, transpathoptionsGE_conv, simoptionsGE_GI, vfoptionsGE_GI);
fprintf('Convergence criterion decided rather than maxiter, one more iteration changes nothing, this should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathConv_plus1.r-PricePathConv.r)),max(abs(PricePathConv_plus1.w-PricePathConv.w)))
clear transpathoptionsGE_conv PricePathConv PricePathConv_plus1 GEcondnPathConv GEcondnDist

clear FnsToEvaluateGE vfoptionsGE_GI simoptionsGE_GI p_eqm_GI GEcondns_GI ParamsGE_GI V_finalGE_GI Policy_finalGE_GI AgentDist_initialGE_GI PricePathGE_GI PricePathGEOut_GI PricePathBumped_GI PricePathBumpedOut_GI GeneralEqmEqnsGE heteroagentoptionsGE p_eqm GEcondns ParamsGE V_finalGE Policy_finalGE AgentDist_initialGE PricePathGE ParamPathGE transpathoptionsGE PricePathGEOut PricePathBumped PricePathBumpedOut bumpr bumpw Tbump

%% Run the GE transition path, but with transpathoptions.maxiter=1, so it ends after one iteration -- just a shape-check (the core is tested elsewhere)
transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(earnings) 0;
% transpathoptions.GEnewprice has no default while the Newton options are being built, so it has to be
% set. factor=0 means the dummy eqn never moves the price, which is what this shape-check wants: it is
% only testing that the returned objects have the right shapes.
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
