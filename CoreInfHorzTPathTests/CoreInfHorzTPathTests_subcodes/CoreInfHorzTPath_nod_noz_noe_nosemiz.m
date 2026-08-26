function output=CoreInfHorzTPath_nod_noz_noe_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
n_d_GE=0; d_grid_GE=[]; % nod: no decision variable in the GE block either
n_z=0; z_grid=[]; pi_z=[];

ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(aprime,a) a;
FnsToEvaluate.earnings=@(aprime,a,w) w;
FnsToEvaluate.nextassets=@(aprime,a) aprime; % aprime-dependent, so it tests the policy decode directly


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

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

%% ValueFnFromPolicy along the path: recomputing V from the policy must reproduce VPath
% ValueFnFromPolicyOnTransPath_InfHorz evaluates the Bellman at the given PolicyPath instead of
% maximising, so it must return exactly the VPath that produced that policy. (Divide-and-conquer
% is irrelevant here -- there is no maximisation -- so this covers the no-GI tier as a whole.)
VfromPolicyPath=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPath2,V_final,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions2);
fprintf('ValueFnFromPolicy along TPath (no GI), this should be zero: %2.10f \n',max(abs(VfromPolicyPath(:)-VPath2(:))))
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

vfoptionsFH=vfoptions2;
vfoptionsFH.V_Jplus1=V_final;
[V_FH,Policy_FH]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFH,DiscountFactorParamNames,[],vfoptionsFH);

% Reshape to (everything else)-by-periods, so these lines are agnostic to the state dimensions
temp1=reshape(VPath2,[],T); temp2=reshape(V_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, V: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyPath2,[],T); temp2=reshape(Policy_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, Policy: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath2,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions2),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, PolicyVals: %2.10f \n',max(abs(dev(:))))

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
fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AgentDist: %2.10f \n',max(abs(dev(:))))

AggVarsPathFH=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AllStatsPathFH=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AgeCondStats_FH=LifeCycleProfiles_FHorz_Case1(StatDist_FH,Policy_FH,FnsToEvaluate,ParamsFH,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFH);

for fnname={'assets','earnings','nextassets'}
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AggVars %s Mean: %2.10f \n',fnname{1},max(abs(dev(:))))
end
% AllStats: every field the two commands have in common (periods/ages are the trailing dim in both)
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (with DC, no GI), this should be zero, AllStats assets %s: %2.10f \n',statname{1},max(abs(dev(:))))
end

clear V_FH Policy_FH ParamsFH vfoptionsFH temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

%% Cross-check against an FHorz solve, with a shock process that varies along the path
% Skipped for this subcode: there is no markov z here, so transpathoptions.zpathtrivial has
% nothing to vary. (epathtrivial, the iid-e analogue, is a separate flag and is not covered.)
fprintf('TPath vs FHorz (path-varying z): skipped, this model has no markov z \n')

%%
clear VPath1 VPath2 PolicyPath1 PolicyPath2

%% Solve with grid-interpolation. With and without divide-and-conquer
vfoptions3=vfoptions1_GI;
simoptions3=simoptions1_GI;
[VPath3,PolicyPath3]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);

vfoptions4=vfoptions3;
vfoptions4.divideandconquer=1;
simoptions4=simoptions3;
[VPath4,PolicyPath4]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);

fprintf('Divide-and-conquer (with GI), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('Divide-and-conquer (with GI), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

%% ValueFnFromPolicy along the path, grid-interpolation tier
VfromPolicyPath=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPath4,V_final_GI,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions4);
fprintf('ValueFnFromPolicy along TPath (with GI), this should be zero: %2.10f \n',max(abs(VfromPolicyPath(:)-VPath4(:))))
clear VfromPolicyPath

%% Cross-check against an FHorz solve (with DC and GI)
% Same construction as the no-GI version above, against the grid-interpolation tier. Not
% redundant with it: GI genuinely changes the solution, so it needs its own comparison.
% Note on the raw Policy line: under GI the indexes can legitimately differ, because (L1,L2) and
% (L1+1,L2-(ngridinterp+1)) encode the same aprime when the optimum sits on a coarse grid point,
% and the two implementations need not break that tie the same way. The PolicyVals line is
% immune to it, and is the one that has to be zero.
N_jFH=T-1;
ParamsFH=Params;
ParamsFH.r=PricePath.r(1:N_jFH);
ParamsFH.w=ParamPath.w(1:N_jFH);

vfoptionsFH_GI=vfoptions4;
vfoptionsFH_GI.V_Jplus1=V_final_GI;
[V_FH,Policy_FH]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFH,DiscountFactorParamNames,[],vfoptionsFH_GI);

temp1=reshape(VPath4,[],T); temp2=reshape(V_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, V: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyPath4,[],T); temp2=reshape(Policy_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, Policy: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath4,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions4),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('TPath vs FHorz (with DC and GI), this should be zero, PolicyVals: %2.10f \n',max(abs(dev(:))))

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
fprintf('TPath vs FHorz (with DC and GI), this should be zero, AgentDist: %2.10f \n',max(abs(dev(:))))

AggVarsPathFH=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath4, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AllStatsPathFH=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath4, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AgeCondStats_FH=LifeCycleProfiles_FHorz_Case1(StatDist_FH,Policy_FH,FnsToEvaluate,ParamsFH,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFH);

for fnname={'assets','earnings','nextassets'}
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('TPath vs FHorz (with DC and GI), this should be zero, AggVars %s Mean: %2.10f \n',fnname{1},max(abs(dev(:))))
end
% AllStats: every field the two commands have in common (periods/ages are the trailing dim in both)
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('TPath vs FHorz (with DC and GI), this should be zero, AllStats assets %s: %2.10f \n',statname{1},max(abs(dev(:))))
end

clear V_FH Policy_FH ParamsFH vfoptionsFH_GI temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

%%
clear VPath3 VPath4 PolicyPath3 PolicyPath4

%% Big a_grid: moments along path should be close with/without grid interp
[VPath2b,PolicyPath2b]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big, Params, n_d, n_a_big, n_z, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPath2=AgentDistOnTransPath_InfHorz(AgentDist_initial_big, PricePath, ParamPath, PolicyPath2b, n_d, n_a_big, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, d_grid, a_grid_big, z_grid, simoptions2);

[VPath4b,PolicyPath4b]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_big_GI, Policy_final_big_GI, Params, n_d, n_a_big, n_z, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_InfHorz(AgentDist_initial_big, PricePath, ParamPath, PolicyPath4b, n_d, n_a_big, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, d_grid, a_grid_big, z_grid, simoptions4);

%% SimPanel along the path: per-period panel mean should reproduce the AggVars along the path
% (Monte Carlo simulation, so this is a roughly-equal check, not machine precision)
% Done twice: off the no-GI policy path, and off the GI policy path. The GI one matters
% because with gridinterplayer=1 the Policy carries the extra L2 and L2flag channels that
% SimPanelValues_TransPath_InfHorz has to strip (l_daprime=size(PolicyPath,1)-2*gridinterplayer).
% Each panel is compared against the AggVars built from the SAME policy path.
simoptionsSP=simoptions2;
simoptionsSP.numbersims=10^4;
simoptionsSP.simperiods=T;
SimPanelTPath=SimPanelValues_TransPath_InfHorz(PolicyPath2b, PricePath, ParamPath, T, AgentDist_initial_big, n_d, n_a_big, n_z, pi_z, d_grid, a_grid_big, z_grid, FnsToEvaluate, Params, simoptionsSP);

simoptionsSP_GI=simoptions4;
simoptionsSP_GI.numbersims=10^4;
simoptionsSP_GI.simperiods=T;
SimPanelTPath_GI=SimPanelValues_TransPath_InfHorz(PolicyPath4b, PricePath, ParamPath, T, AgentDist_initial_big, n_d, n_a_big, n_z, pi_z, d_grid, a_grid_big, z_grid, FnsToEvaluate, Params, simoptionsSP_GI);

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
fprintf('AgentDist along TPath with/without grid interp, this should be close to zero: %2.8f \n',max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
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
% Guarded: the command has no branch for N_z=0, and no e support at all (it defaults simoptions.n_e
% but never reads it), so the noz and e subcodes would either error or silently return the wrong
% thing. Remove the guard once the toolkit covers those cases.
if prod(n_z)>0 && ~isfield(simoptions,'n_e')
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
    AggVarsPathAC=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathAC, PolicyPathAC, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsAC);
    AllStatsPathAC=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathAC, PolicyPathAC, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsAC);
    CorrTransProbsPathAC=EvalFnOnTransPath_AutoCorrTransProbs_InfHorz(FnsToEvaluate, AgentDistPathAC, PolicyPathAC, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, pi_z, simoptionsAC);

    for fnname=fieldnames(FnsToEvaluate)'
        dev=CorrTransProbsPathAC.(fnname{1}).Mean-AggVarsPathAC.(fnname{1}).Mean;
        fprintf('AutoCorrTransProbs vs AggVars along TPath, this should be zero, Mean %s: %2.10f \n',fnname{1},max(abs(dev)))
        dev=CorrTransProbsPathAC.(fnname{1}).StdDeviation-AllStatsPathAC.(fnname{1}).StdDeviation;
        fprintf('AutoCorrTransProbs vs AllStats along TPath, this should be zero, StdDeviation %s: %2.10f \n',fnname{1},max(abs(dev)))
    end

    for fnname=fieldnames(FnsToEvaluate)'
        sdthis=CorrTransProbsPathAC.(fnname{1}).StdDeviation(2:T);
        sdlag=CorrTransProbsPathAC.(fnname{1}).StdDeviation(1:T-1);
        dev=CorrTransProbsPathAC.(fnname{1}).AutoCorrelation(1:T-1)-CorrTransProbsPathAC.(fnname{1}).AutoCovariance(1:T-1)./(sdlag.*sdthis);
        fprintf('AutoCorrTransProbs internal, corr=cov/(sd*sd), this should be zero, %s: %2.10f \n',fnname{1},max(abs(dev)))
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
else
    fprintf('AutoCorrTransProbs along TPath: skipped, EvalFnOnTransPath_AutoCorrTransProbs_InfHorz has no N_z=0 branch and no e support \n')
end

%% Plots
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath2.earnings.Mean, 1:1:T,AggVarsPath4.earnings.Mean)
title('Earnings Mean'); legend('1','2')
subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean'); legend('1','2')

clear VPath2b VPath4b AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4

%% Constant path -- VPath should equal repmat of stationary V, similarly Policy and AgentDist
% Done three ways: (1) no DC, no GI; (2) with DC, no GI; (3) with DC and GI.
% Both the raw Policy (indexes) and the PolicyVals (the actual d/aprime values) are compared.
% Note: with the grid interpolation layer the Policy INDEXES can legitimately differ by
% ngridinterp+1, because (L1,L2) and (L1+1,L2-(ngridinterp+1)) encode the same aprime when the
% optimum sits exactly on a coarse grid point. The PolicyVals comparison is immune to this.
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.sigma=Params.sigma*ones(1,T);

% (1) without divide-and-conquer, without grid interpolation
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions1);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions1);

Vfin_rep=repmat(V_final,1,1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, V: %2.8f \n',max(abs(VPath1(:)-Vfin_rep(:))))
Polfin_rep=repmat(Policy_final,1,1,1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
PolicyValsfin=PolicyInd2Val_InfHorz(Policy_final,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyValsfin_rep=repmat(PolicyValsfin(:),1,T); % T is the last dimension of PolicyValsPath
fprintf('Constant TPath (no DC, no GI), this should be zero, PolicyVals: %2.8f \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
AD_rep=repmat(AgentDist_initial,1,1,T);
fprintf('Constant TPath (no DC, no GI), this should be zero, AgentDist (note: tolerance=1e-6): %2.8f \n',max(abs(AgentDistPath1(:)-AD_rep(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1

% (2) with divide-and-conquer, without grid interpolation (same targets as (1))
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions2);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2);

fprintf('Constant TPath (with DC, no GI), this should be zero, V: %2.8f \n',max(abs(VPath1(:)-Vfin_rep(:))))
fprintf('Constant TPath (with DC, no GI), this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
fprintf('Constant TPath (with DC, no GI), this should be zero, PolicyVals: %2.8f \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
fprintf('Constant TPath (with DC, no GI), this should be zero, AgentDist (note: tolerance=1e-6): %2.8f \n',max(abs(AgentDistPath1(:)-AD_rep(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1

% (3) with divide-and-conquer, with grid interpolation
% (uses the GI terminal V/Policy, and the stationary dist computed with GI)
AgentDist_initial_GI=StationaryDist_InfHorz(Policy_final_GI,n_d,n_a,n_z,pi_z,simoptions4,Params,[]);
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial_GI, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, transpathoptionsbaseline, simoptions4);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4);

Vfin_rep_GI=repmat(V_final_GI,1,1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, V: %2.8f \n',max(abs(VPath1(:)-Vfin_rep_GI(:))))
Polfin_rep_GI=repmat(Policy_final_GI,1,1,1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Polfin_rep_GI(:))))
PolicyValsfin_GI=PolicyInd2Val_InfHorz(Policy_final_GI,n_d,n_a,n_z,d_grid,a_grid,vfoptions4);
PolicyValsfin_GI_rep=repmat(PolicyValsfin_GI(:),1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, PolicyVals: %2.8f \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_GI_rep(:))))
AD_rep_GI=repmat(AgentDist_initial_GI,1,1,T);
fprintf('Constant TPath (with DC and GI), this should be zero, AgentDist (note: tolerance=1e-6): %2.8f \n',max(abs(AgentDistPath1(:)-AD_rep_GI(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1
%% General equilibrium: solve the stationary GE first, then run a null-reform transition path
% The only part of this file that exercises the GE machinery itself (the shooting algorithm, the
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
% transpathoptions.GEnewprice=3 rather than the default 1: GEnewprice=1 treats the eqns as price
% UPDATING FORMULAE (updatePricePathNew_TPath_tt sets PricePathNew_tt=p_i directly), whereas
% GEnewprice=3 treats them as residuals and updates via howtoupdate, which is the same convention
% HeteroAgentStationaryEqm_InfHorz uses. One definition of the firm, not two.
%
% firm* names are used because the 2A entrepreneur ReturnFn already uses alpha and delta with its
% own (Kitao) meanings.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a) a;
FnsToEvaluateGE.N=@(aprime,a) 1;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions1, vfoptions1);
fprintf('Stationary GE: r=%2.8f, w=%2.8f \n',p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%2.10f, LabourMarket=%2.10f \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)

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
% GEnewprice=3: treat the GeneralEqmEqns as residuals and update each price by a damped fraction
% of its own residual. howtoupdate columns are {GEcondn name, price name, add, factor}, and
% updatePricePathNew_TPath_tt does new = old + add*factor*residual - (1-add)*factor*residual.
% Both residuals are of the form (price - its firm FOC), so a POSITIVE residual means the price is
% too high and must come down: hence add=0 (subtract) for both. Rows are given in the same order
% as the PricePath fields (r then w), which is the order setupGEnewprice3_shooting expects.
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

% (i) Start from the equilibrium path. The update formulae return the same prices, so the solver
% must leave the path alone.
PricePathGEOut=TransitionPath_InfHorz(PricePathGE, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE, simoptions1, vfoptions1, []);
fprintf('Null-reform GE from the equilibrium path, this should be zero, r: %2.10f \n',max(abs(PricePathGEOut.r-p_eqm.r)))
fprintf('Null-reform GE from the equilibrium path, this should be zero, w: %2.10f \n',max(abs(PricePathGEOut.w-p_eqm.w)))

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
fprintf('Null-reform GE from a bumped guess, should converge back toward zero, r: %2.10f \n',max(abs(PricePathBumpedOut.r-p_eqm.r)))
fprintf('Null-reform GE from a bumped guess, should converge back toward zero, w: %2.10f \n',max(abs(PricePathBumpedOut.w-p_eqm.w)))

%% GEnewprice3 additional factors: ramp the update factors over the shooting iterations
% howtoupdate can be given with 6 columns instead of 4, the extra two being f_add and t_add: the
% factor in column 4 is used as it is on iteration 1, f_add times it from iteration t_add on, and
% linearly in between. The point is to start cautiously and speed up once the path has settled.
%
% (a) is the check that matters: no additionalfactor, and additionalfactor=[1,2], must be bit
% identical, since f_add=1 makes the ramp the identity whatever t_add is. Run from the bumped guess so the prices
% actually move, and for only a few iterations, because a difference in the update rule shows up on
% the first one and there is no reason to pay for a hundred.
transpathoptionsGE_AF=transpathoptionsGE;
transpathoptionsGE_AF.maxiter=5;
PricePathAF_noAF=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE_AF, simoptions1, vfoptions1, []);

transpathoptionsGE_AF.GEnewprice3.additionalfactor=[1,1,2];
PricePathAF_withAF=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsGE_AF, simoptions1, vfoptions1, []);
fprintf('additionalfactor no-op (f_add=1), unset vs [1,1,2], this should be zero, r: %2.10f \n',max(abs(PricePathAF_noAF.r-PricePathAF_withAF.r)))
fprintf('additionalfactor no-op (f_add=1), unset vs [1,1,2], this should be zero, w: %2.10f \n',max(abs(PricePathAF_noAF.w-PricePathAF_withAF.w)))

% The ramp itself is measured in the DC+GI pass below, not here: this plain-tier p_eqm is not a
% reliable equilibrium (see the GEcondns printed above), so 'distance from p_eqm' would not mean
% what it says. The no-op check above is unaffected, being a bit-equality test.
clear transpathoptionsGE_AF PricePathAF_noAF PricePathAF_withAF

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
fprintf('Stationary GE conditions (with DC and GI), these should be close to zero: CapitalMarket=%2.10f, LabourMarket=%2.10f \n',GEcondns_GI.CapitalMarket,GEcondns_GI.LabourMarket)
fprintf('Stationary GE, the two tiers should agree, difference in r: %2.8f, in w: %2.8f \n',abs(p_eqm_GI.r-p_eqm.r),abs(p_eqm_GI.w-p_eqm.w))

ParamsGE_GI=Params;
ParamsGE_GI.r=p_eqm_GI.r;
ParamsGE_GI.w=p_eqm_GI.w;
[V_finalGE_GI,Policy_finalGE_GI]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE_GI,DiscountFactorParamNames,[],vfoptionsGE_GI);
AgentDist_initialGE_GI=StationaryDist_InfHorz(Policy_finalGE_GI,n_d_GE,n_a_GE,n_z,pi_z,simoptionsGE_GI,ParamsGE_GI,[]);

PricePathGE_GI.r=p_eqm_GI.r*ones(1,T);
PricePathGE_GI.w=p_eqm_GI.w*ones(1,T);

transpathoptionsGE.maxiter=25; % as above: the no-change solve starts at the answer
PricePathGEOut_GI=TransitionPath_InfHorz(PricePathGE_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('Null-reform GE (with DC and GI) from the equilibrium path, this should be zero, r: %2.10f \n',max(abs(PricePathGEOut_GI.r-p_eqm_GI.r)))
fprintf('Null-reform GE (with DC and GI) from the equilibrium path, this should be zero, w: %2.10f \n',max(abs(PricePathGEOut_GI.w-p_eqm_GI.w)))

% Same bump profile as the plain tier (bumpr/bumpw/Tbump are reused, so the two tiers face an
% identical perturbation and the results are directly comparable)
PricePathBumped_GI=PricePathGE_GI;
PricePathBumped_GI.r(1:Tbump)=PricePathGE_GI.r(1:Tbump).*(1+bumpr);
PricePathBumped_GI.w(1:Tbump)=PricePathGE_GI.w(1:Tbump).*(1+bumpw);
transpathoptionsGE.maxiter=100; % as above: the bumped solve has to actually travel
PricePathBumpedOut_GI=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('Null-reform GE (with DC and GI) bumped guess, max initial deviation: r %2.8f, w %2.8f \n',max(abs(PricePathBumped_GI.r-p_eqm_GI.r)),max(abs(PricePathBumped_GI.w-p_eqm_GI.w)))
fprintf('Null-reform GE (with DC and GI) from a bumped guess, should converge back toward zero, r: %2.10f \n',max(abs(PricePathBumpedOut_GI.r-p_eqm_GI.r)))
fprintf('Null-reform GE (with DC and GI) from a bumped guess, should converge back toward zero, w: %2.10f \n',max(abs(PricePathBumpedOut_GI.w-p_eqm_GI.w)))

% A genuine additionalfactor ramp, measured on the DC+GI tier because that is where the stationary
% solve is sound (the plain tier's GE conditions are far enough from zero that distance from its
% p_eqm is not a meaningful measure of convergence). There is no exact reference for what a ramp
% should give, so this is reported rather than asserted: it says whether ramping the factors up got
% closer to the equilibrium in the same number of iterations as the un-ramped solve just above
% (PricePathBumpedOut_GI: same starting guess, same maxiter=100). Convergence speed is a property of
% the model, not a correctness condition.
transpathoptionsGE_AF=transpathoptionsGE;
transpathoptionsGE_AF.maxiter=100;
transpathoptionsGE_AF.GEnewprice3.additionalfactor=[2,10,50]; % hold at factor for 10 iterations, then ramp to 2x by iteration 50
PricePathAF_ramp=TransitionPath_InfHorz(PricePathBumped_GI, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE_GI, DiscountFactorParamNames, transpathoptionsGE_AF, simoptionsGE_GI, vfoptionsGE_GI, []);
fprintf('additionalfactor ramp (f_add=2 over iterations 10 to 50) vs no ramp, deviation from eqm, r: %2.10f vs %2.10f \n',max(abs(PricePathAF_ramp.r-p_eqm_GI.r)),max(abs(PricePathBumpedOut_GI.r-p_eqm_GI.r)))
fprintf('additionalfactor ramp (f_add=2 over iterations 10 to 50) vs no ramp, deviation from eqm, w: %2.10f vs %2.10f \n',max(abs(PricePathAF_ramp.w-p_eqm_GI.w)),max(abs(PricePathBumpedOut_GI.w-p_eqm_GI.w)))
clear transpathoptionsGE_AF PricePathAF_ramp

clear FnsToEvaluateGE vfoptionsGE_GI simoptionsGE_GI p_eqm_GI GEcondns_GI ParamsGE_GI V_finalGE_GI Policy_finalGE_GI AgentDist_initialGE_GI PricePathGE_GI PricePathGEOut_GI PricePathBumped_GI PricePathBumpedOut_GI GeneralEqmEqnsGE heteroagentoptionsGE p_eqm GEcondns ParamsGE V_finalGE Policy_finalGE AgentDist_initialGE PricePathGE ParamPathGE transpathoptionsGE PricePathGEOut PricePathBumped PricePathBumpedOut bumpr bumpw Tbump

%% Run the GE transition path with transpathoptions.maxiter=1 -- shape check only
transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(earnings) 0;

PricePath2=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);

% Big grid + DC
PricePath3A=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, n_d, n_a_big, n_z, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions2, vfoptions2, []);

% Big grid + DC + GI
PricePath3B=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_big_GI, AgentDist_initial_big, n_d, n_a_big, n_z, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions4, vfoptions4, []);

fprintf('One iter of TPath, big-grid with/without GI, this should be close to zero: %2.8f \n',max(abs(PricePath3A.r-PricePath3B.r)))

%%
output=struct();

end
