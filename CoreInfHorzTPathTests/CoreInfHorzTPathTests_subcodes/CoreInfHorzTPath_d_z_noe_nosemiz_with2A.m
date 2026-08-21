function output=CoreInfHorzTPath_d_z_noe_nosemiz_with2A(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c)
% Two-endogenous-state InfHorz TPath test, mirroring the fig-4 (d, z) battery.
% Model: simplified Kitao (2008) with Bruggemann (2021) endogenous labor supply, no taxes --
% the same model CoreInfHorzTests uses for its d 2A subcode.
%   Endogenous states: asset, occupation (0=worker, 1=entrepreneur)
%   Decision variable:  l (labor supply; entrepreneurs must supply exactly lbar)
%   Exogenous markov z: eta (labor productivity), theta (entrepreneurial ability)
% n_a(2)=2 is what triggers the DC2A/GI2A code paths; the decision variable additionally
% routes through the Refine variants of them.
%
% Two things differ from the 1A subcodes and are worth knowing before reading the checks:
%
% (1) V is -Inf at infeasible states (an entrepreneur with near-zero assets has no feasible
%     positive consumption). -Inf minus -Inf is NaN, and max() omits NaN, so the comparisons
%     below silently ignore states where both sides agree on -Inf -- which is what we want. A
%     state where only ONE side is -Inf gives +-Inf, so max() reports Inf and the check still
%     fails loudly. So no special handling is needed, but these lines are not doing quite what
%     they appear to at a glance.
%
% (2) The plain tier (no DC, no GI) works for 2A as it stands: pure discretization is generic
%     in N_a=prod(n_a), so it does not care how many endogenous states there are. The
%     GI-without-DC tier (vfoptions3) does NOT. ValueFnIter_InfHorz_TPath_SingleStep_GI1_nod_raw
%     assumes N_a==length(a_grid), which is false for two endogenous states: here N_a=202 while
%     the stacked a_grid has 103 entries, so its interp1 on line 17 errors. The dispatcher sends
%     gridinterplayer==1 with divideandconquer==0 to that 1A-only raw with no length(n_a) check.
%     A GI2A-without-DC raw is the missing piece -- the stationary solver already has one
%     (ValueFnIter_InfHorz_preGI2A_* / _postGI2A_*), so the gap is in the TPath tier only.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

vfoptions.verbose_advice=0; % 2-endo GI would otherwise sound the postGI advice on every solve

% Terminal parameters are the end of the path (same convention as the 1A subcodes)
Params.r=PricePath.r(T);
Params.w=ParamPath.w(T);

ReturnFn=@(l,aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi,xi,sigma2,lbar) ...
    ReturnFn_d_z_noe_nosemiz_with2A(l,aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi,xi,sigma2,lbar);

% Setup some FnsToEvaluate (functions of (l,aprime,eprime,a,e,eta,theta) -- decision variable first)
FnsToEvaluate.assets=@(l,aprime,eprime,a,e,eta,theta) a;
FnsToEvaluate.entrepreneur=@(l,aprime,eprime,a,e,eta,theta) e; % fraction who are entrepreneurs
FnsToEvaluate.nextassets=@(l,aprime,eprime,a,e,eta,theta) aprime; % aprime-dependent, so it tests the policy decode directly
FnsToEvaluate.laborsupply=@(l,aprime,eprime,a,e,eta,theta) l; % d-dependent, so it tests the decision-variable decode


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

fprintf('with2A: Divide-and-conquer, this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('with2A: Divide-and-conquer, this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))

% lowmemory should give same answer
vfoptions1.lowmemory=1;
[VPath1B,PolicyPath1B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[VPath2B,PolicyPath2B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
vfoptions2.lowmemory=0;

fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath1B(:))))
fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath1B(:))))
fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(VPath2(:)-VPath2B(:))))
fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(PolicyPath2(:)-PolicyPath2B(:))))

%% ValueFnFromPolicy along the path: recomputing V from the policy must reproduce VPath
% ValueFnFromPolicyOnTransPath_InfHorz evaluates the Bellman at the given PolicyPath instead of
% maximising, so it must return exactly the VPath that produced that policy. (Divide-and-conquer
% is irrelevant here -- there is no maximisation -- so this covers the no-GI tier as a whole.)
VfromPolicyPath=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPath2,V_final,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions2);
fprintf('ValueFnFromPolicy along TPath (no GI), this should be zero: %2.10f \n',max(abs(VfromPolicyPath(:)-VPath2(:))))
clear VfromPolicyPath

%% Cross-check against an FHorz solve (with DC, no GI)
% See the 1A subcodes for the full explanation of the mapping. In short: the TPath value fn
% iteration is a finite-horizon problem in disguise, so encode the price/param path as
% age-dependent parameters, pass the terminal V via vfoptions.V_Jplus1, and the FHorz model has
% N_j=T-1 periods with V_FHorz(:,:,:,j) equal to VPath(:,:,:,j).
N_jFH=T-1;
ParamsFH=Params;
ParamsFH.r=PricePath.r(1:N_jFH);
ParamsFH.w=ParamPath.w(1:N_jFH);

vfoptionsFH=vfoptions2;
vfoptionsFH.V_Jplus1=V_final;
[V_FH,Policy_FH]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFH,DiscountFactorParamNames,[],vfoptionsFH);

temp1=reshape(VPath2,[],T); temp2=reshape(V_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC, no GI), this should be zero, V: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyPath2,[],T); temp2=reshape(Policy_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC, no GI), this should be zero, Policy: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath2,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions2),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC, no GI), this should be zero, PolicyVals: %2.10f \n',max(abs(dev(:))))

% --- carry the same construction through the agent dist and the downstream statistics ---
ParamsFH.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFH.mewjFH(end)=1-sum(ParamsFH.mewjFH(1:end-1)); % sum must be one to within 1e-15
simoptionsFH=simoptions2;
simoptionsFH.npoints=100;
simoptionsFH.nquantiles=20;
simoptionsFH.whichstats=ones(7,1); % the two commands ship different defaults; pin both

AgentDistPathFH=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath2, n_d, n_a, n_z, pi_z, T, Params, simoptionsFH);
StatDist_FH=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FH,n_d,n_a,n_z,N_jFH,pi_z,ParamsFH,simoptionsFH);
temp1=reshape(AgentDistPathFH,[],T); temp2=reshape(StatDist_FH,[],N_jFH);
temp2=temp2./sum(temp2,1); % undo the age weighting, so each age carries mass one like each TPath period
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC, no GI), this should be zero, AgentDist: %2.10f \n',max(abs(dev(:))))

AggVarsPathFH=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AllStatsPathFH=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AgeCondStats_FH=LifeCycleProfiles_FHorz_Case1(StatDist_FH,Policy_FH,FnsToEvaluate,ParamsFH,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFH);

for fnname={'assets','entrepreneur','nextassets','laborsupply'}
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('with2A: TPath vs FHorz (with DC, no GI), this should be zero, AggVars %s Mean: %2.10f \n',fnname{1},max(abs(dev(:))))
end
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('with2A: TPath vs FHorz (with DC, no GI), this should be zero, AllStats assets %s: %2.10f \n',statname{1},max(abs(dev(:))))
end

clear V_FH Policy_FH ParamsFH vfoptionsFH temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

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

fprintf('with2A: Divide-and-conquer (with GI), this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath4(:))))
fprintf('with2A: Divide-and-conquer (with GI), this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath4(:))))

% lowmemory should give same answer
vfoptions3.lowmemory=1;
[VPath3B,PolicyPath3B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions3);
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[VPath4B,PolicyPath4B]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
vfoptions4.lowmemory=0;

fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(VPath3(:)-VPath3B(:))))
fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(PolicyPath3(:)-PolicyPath3B(:))))
fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(VPath4(:)-VPath4B(:))))
fprintf('with2A: low memory, this should be zero: %2.8f \n',max(abs(PolicyPath4(:)-PolicyPath4B(:))))

%% ValueFnFromPolicy along the path, grid-interpolation tier
VfromPolicyPath=ValueFnFromPolicyOnTransPath_InfHorz(PolicyPath4,V_final_GI,ParamPath,PricePath,T,n_d,n_a,n_z,d_grid,a_grid,z_grid, pi_z, ReturnFn, Params, DiscountFactorParamNames, vfoptions4);
fprintf('ValueFnFromPolicy along TPath (with GI), this should be zero: %2.10f \n',max(abs(VfromPolicyPath(:)-VPath4(:))))
clear VfromPolicyPath

%% Cross-check against an FHorz solve (with DC and GI)
N_jFH=T-1;
ParamsFH=Params;
ParamsFH.r=PricePath.r(1:N_jFH);
ParamsFH.w=ParamPath.w(1:N_jFH);

vfoptionsFH_GI=vfoptions4;
vfoptionsFH_GI.V_Jplus1=V_final_GI;
[V_FH,Policy_FH]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsFH,DiscountFactorParamNames,[],vfoptionsFH_GI);

temp1=reshape(VPath4,[],T); temp2=reshape(V_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC and GI), this should be zero, V: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyPath4,[],T); temp2=reshape(Policy_FH,[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC and GI), this should be zero, Policy: %2.10f \n',max(abs(dev(:))))
temp1=reshape(PolicyInd2Val_InfHorz_TPath(PolicyPath4,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4),[],T);
temp2=reshape(PolicyInd2Val_FHorz(Policy_FH,n_d,n_a,n_z,N_jFH,d_grid,a_grid,vfoptions4),[],N_jFH);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC and GI), this should be zero, PolicyVals: %2.10f \n',max(abs(dev(:))))

% --- carry the same construction through the agent dist and the downstream statistics ---
ParamsFH.mewjFH=ones(1,N_jFH)/N_jFH;
ParamsFH.mewjFH(end)=1-sum(ParamsFH.mewjFH(1:end-1));
simoptionsFH=simoptions4;
simoptionsFH.npoints=100;
simoptionsFH.nquantiles=20;
simoptionsFH.whichstats=ones(7,1);

AgentDistPathFH=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath4, n_d, n_a, n_z, pi_z, T, Params, simoptionsFH);
StatDist_FH=StationaryDist_FHorz_Case1(AgentDist_initial,{'mewjFH'},Policy_FH,n_d,n_a,n_z,N_jFH,pi_z,ParamsFH,simoptionsFH);
temp1=reshape(AgentDistPathFH,[],T); temp2=reshape(StatDist_FH,[],N_jFH);
temp2=temp2./sum(temp2,1);
dev=temp1(:,1:N_jFH)-temp2;
fprintf('with2A: TPath vs FHorz (with DC and GI), this should be zero, AgentDist: %2.10f \n',max(abs(dev(:))))

AggVarsPathFH=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath4, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AllStatsPathFH=EvalFnOnTransPath_AllStats_InfHorz(FnsToEvaluate, AgentDistPathFH, PolicyPath4, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, simoptionsFH);
AgeCondStats_FH=LifeCycleProfiles_FHorz_Case1(StatDist_FH,Policy_FH,FnsToEvaluate,ParamsFH,[],n_d,n_a,n_z,N_jFH,d_grid,a_grid,z_grid,simoptionsFH);

for fnname={'assets','entrepreneur','nextassets','laborsupply'}
    dev=AggVarsPathFH.(fnname{1}).Mean(1:N_jFH)-AgeCondStats_FH.(fnname{1}).Mean;
    fprintf('with2A: TPath vs FHorz (with DC and GI), this should be zero, AggVars %s Mean: %2.10f \n',fnname{1},max(abs(dev(:))))
end
for statname={'Mean','Median','StdDeviation','Variance','Gini','LorenzCurve','QuantileCutoffs','QuantileMeans'}
    temp1=AllStatsPathFH.assets.(statname{1}); temp2=AgeCondStats_FH.assets.(statname{1});
    dev=temp1(:,1:N_jFH)-temp2;
    fprintf('with2A: TPath vs FHorz (with DC and GI), this should be zero, AllStats assets %s: %2.10f \n',statname{1},max(abs(dev(:))))
end

clear V_FH Policy_FH ParamsFH vfoptionsFH_GI temp1 temp2 dev AgentDistPathFH StatDist_FH AggVarsPathFH AllStatsPathFH AgeCondStats_FH simoptionsFH

%%
clear VPath3 VPath4 PolicyPath3 PolicyPath4 VPath3B VPath4B PolicyPath3B PolicyPath4B

%% Big asset grid: moments along path should be close with/without grid interp
[VPath2b,PolicyPath2b]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_big, Policy_final_big, Params, n_d, n_a_big, n_z, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPath2=AgentDistOnTransPath_InfHorz(AgentDist_initial_big, PricePath, ParamPath, PolicyPath2b, n_d, n_a_big, n_z, pi_z, T, Params, simoptions2);
AggVarsPath2=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath2, PolicyPath2b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, d_grid, a_grid_big, z_grid, simoptions2);

[VPath4b,PolicyPath4b]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_big_GI, Policy_final_big_GI, Params, n_d, n_a_big, n_z, d_grid, a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath4=AgentDistOnTransPath_InfHorz(AgentDist_initial_big, PricePath, ParamPath, PolicyPath4b, n_d, n_a_big, n_z, pi_z, T, Params, simoptions4);
AggVarsPath4=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath4, PolicyPath4b, PricePath, ParamPath, Params, T, n_d, n_a_big, n_z, d_grid, a_grid_big, z_grid, simoptions4);

%% SimPanel along the path: per-period panel mean should reproduce the AggVars along the path
% (Monte Carlo simulation, so this is a roughly-equal check, not machine precision)
% Done twice: off the no-GI policy path, and off the GI policy path.
simoptionsSP=simoptions2;
simoptionsSP.numbersims=10^4;
simoptionsSP.simperiods=T;
SimPanelTPath=SimPanelValues_TransPath_InfHorz(PolicyPath2b, PricePath, ParamPath, T, AgentDist_initial_big, n_d, n_a_big, n_z, pi_z, d_grid, a_grid_big, z_grid, FnsToEvaluate, Params, simoptionsSP);

simoptionsSP_GI=simoptions4;
simoptionsSP_GI.numbersims=10^4;
simoptionsSP_GI.simperiods=T;
SimPanelTPath_GI=SimPanelValues_TransPath_InfHorz(PolicyPath4b, PricePath, ParamPath, T, AgentDist_initial_big, n_d, n_a_big, n_z, pi_z, d_grid, a_grid_big, z_grid, FnsToEvaluate, Params, simoptionsSP_GI);

fprintf('with2A: SimPanel along TPath: per-period panel mean should roughly match AggVarsPath (Monte Carlo) \n')
dev=abs(mean(SimPanelTPath.entrepreneur,2)'-AggVarsPath2.entrepreneur.Mean);
fprintf('with2A: SimPanel (no GI),   entrepreneur: max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.entrepreneur.Mean)))
dev=abs(mean(SimPanelTPath.assets,2)'-AggVarsPath2.assets.Mean);
fprintf('with2A: SimPanel (no GI),   assets:       max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.assets.Mean)))
dev=abs(mean(SimPanelTPath.nextassets,2)'-AggVarsPath2.nextassets.Mean);
fprintf('with2A: SimPanel (no GI),   nextassets:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.nextassets.Mean)))
dev=abs(mean(SimPanelTPath.laborsupply,2)'-AggVarsPath2.laborsupply.Mean);
fprintf('with2A: SimPanel (no GI),   laborsupply:  max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath2.laborsupply.Mean)))
dev=abs(mean(SimPanelTPath_GI.entrepreneur,2)'-AggVarsPath4.entrepreneur.Mean);
fprintf('with2A: SimPanel (with GI), entrepreneur: max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.entrepreneur.Mean)))
dev=abs(mean(SimPanelTPath_GI.assets,2)'-AggVarsPath4.assets.Mean);
fprintf('with2A: SimPanel (with GI), assets:       max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.assets.Mean)))
dev=abs(mean(SimPanelTPath_GI.nextassets,2)'-AggVarsPath4.nextassets.Mean);
fprintf('with2A: SimPanel (with GI), nextassets:   max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.nextassets.Mean)))
dev=abs(mean(SimPanelTPath_GI.laborsupply,2)'-AggVarsPath4.laborsupply.Mean);
fprintf('with2A: SimPanel (with GI), laborsupply:  max abs deviation: %2.8f, max abs percentage deviation: %2.4f%% \n',max(dev),100*max(dev./abs(AggVarsPath4.laborsupply.Mean)))

fprintf('with2A: With/without grid interp, should get much the same moments (for big asset grid) \n')
fprintf('with2A: AgentDist along TPath with/without grid interp, this should be close to zero: %2.8f \n',max(abs(AgentDistPath2(:)-AgentDistPath4(:))))
[AggVarsPath2.entrepreneur.Mean; AggVarsPath4.entrepreneur.Mean]
[AggVarsPath2.assets.Mean; AggVarsPath4.assets.Mean]

%% Plots
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,AggVarsPath2.entrepreneur.Mean, 1:1:T,AggVarsPath4.entrepreneur.Mean)
title('Entrepreneur Share'); legend('1','2')
subplot(2,1,2); plot(1:1:T,AggVarsPath2.assets.Mean, 1:1:T,AggVarsPath4.assets.Mean)
title('Assets Mean'); legend('1','2')

clear VPath2b VPath4b AggVarsPath4 PolicyPath4b AgentDistPath2 AgentDistPath4 SimPanelTPath SimPanelTPath_GI

%% Constant path -- VPath should equal repmat of stationary V, similarly Policy and AgentDist
% Done three ways: (1) no DC, no GI; (2) with DC, no GI; (3) with DC and GI.
% As in the 1A subcodes, with the grid interpolation layer the Policy INDEXES can legitimately
% differ by ngridinterp+1, because (L1,L2) and (L1+1,L2-(ngridinterp+1)) encode the same aprime
% when the optimum sits exactly on a coarse grid point. The PolicyVals comparison is immune.
PricePathConstant.r=Params.r*ones(1,T);
ParamPathConstant.w=Params.w*ones(1,T);

% (1) without divide-and-conquer, without grid interpolation
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions1);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, simoptions1);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions1);

Vfin_rep=repmat(V_final,1,1,1,1,T);
fprintf('with2A: Constant TPath (no DC, no GI), this should be zero, V: %2.8f \n',max(abs(VPath1(:)-Vfin_rep(:))))
Polfin_rep=repmat(Policy_final,1,1,1,1,1,T);
fprintf('with2A: Constant TPath (no DC, no GI), this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
PolicyValsfin=PolicyInd2Val_InfHorz(Policy_final,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyValsfin_rep=repmat(PolicyValsfin(:),1,T); % T is the last dimension of PolicyValsPath
fprintf('with2A: Constant TPath (no DC, no GI), this should be zero, PolicyVals: %2.8f \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
AD_rep=repmat(AgentDist_initial,1,1,1,1,T);
fprintf('with2A: Constant TPath (no DC, no GI), this should be zero, AgentDist (note: tolerance=1e-6): %2.8f \n',max(abs(AgentDistPath1(:)-AD_rep(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1

% (2) with divide-and-conquer, without grid interpolation (same targets as (1))
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final, Policy_final, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions2);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, simoptions2);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions2);

fprintf('with2A: Constant TPath (with DC, no GI), this should be zero, V: %2.8f \n',max(abs(VPath1(:)-Vfin_rep(:))))
fprintf('with2A: Constant TPath (with DC, no GI), this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Polfin_rep(:))))
fprintf('with2A: Constant TPath (with DC, no GI), this should be zero, PolicyVals: %2.8f \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_rep(:))))
fprintf('with2A: Constant TPath (with DC, no GI), this should be zero, AgentDist (note: tolerance=1e-6): %2.8f \n',max(abs(AgentDistPath1(:)-AD_rep(:))))

clear VPath1 PolicyPath1 AgentDistPath1 PolicyValsPath1

% (3) with divide-and-conquer, with grid interpolation
AgentDist_initial_GI=StationaryDist_InfHorz(Policy_final_GI,n_d,n_a,n_z,pi_z,simoptions4,Params,[]);
[VPath1,PolicyPath1]=ValueFnOnTransPath_InfHorz(PricePathConstant, ParamPathConstant, T, V_final_GI, Policy_final_GI, Params, n_d, n_a, n_z, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions4);
AgentDistPath1=AgentDistOnTransPath_InfHorz(AgentDist_initial_GI, PricePathConstant, ParamPathConstant, PolicyPath1, n_d, n_a, n_z, pi_z, T, Params, simoptions4);
PolicyValsPath1=PolicyInd2Val_InfHorz_TPath(PolicyPath1,n_d,n_a,n_z,T,d_grid,a_grid,vfoptions4);

Vfin_rep_GI=repmat(V_final_GI,1,1,1,1,T);
fprintf('with2A: Constant TPath (with DC and GI), this should be zero, V: %2.8f \n',max(abs(VPath1(:)-Vfin_rep_GI(:))))
Polfin_rep_GI=repmat(Policy_final_GI,1,1,1,1,1,T);
fprintf('with2A: Constant TPath (with DC and GI), this should be zero, Policy: %2.8f \n',max(abs(PolicyPath1(:)-Polfin_rep_GI(:))))
PolicyValsfin_GI=PolicyInd2Val_InfHorz(Policy_final_GI,n_d,n_a,n_z,d_grid,a_grid,vfoptions4);
PolicyValsfin_GI_rep=repmat(PolicyValsfin_GI(:),1,T);
fprintf('with2A: Constant TPath (with DC and GI), this should be zero, PolicyVals: %2.8f \n',max(abs(PolicyValsPath1(:)-PolicyValsfin_GI_rep(:))))
AD_rep_GI=repmat(AgentDist_initial_GI,1,1,1,1,T);
fprintf('with2A: Constant TPath (with DC and GI), this should be zero, AgentDist (note: tolerance=1e-6): %2.8f \n',max(abs(AgentDistPath1(:)-AD_rep_GI(:))))

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
% Corporate sector, alongside the entrepreneurs. Kitao/Bruggemann structure: households' total
% assets are split between capital tied up in entrepreneurs' own businesses and capital rented to
% the corporate sector, and the labour entrepreneurs hire comes out of the same pool the corporate
% sector draws on. So the corporate factor inputs are RESIDUALS, and it is those -- not total
% assets and total labour -- that enter the firm's first order conditions.
% firmdelta is set equal to the entrepreneurs' depreciation rate so there is one depreciation
% rate in the model.
Params.firmalpha=0.36;
Params.firmdelta=Params.delta;
Params.firmA=1;

FnsToEvaluateGE.A=@(l,aprime,eprime,a,e,eta,theta) a; % total household assets
FnsToEvaluateGE.K_noncorp=@(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi) EntreCapital_d_with2A(l,aprime,eprime,a,e,theta,r,w,delta,upsilon1,upsilon2,leverage,phi); % capital inside entrepreneurs' own businesses
FnsToEvaluateGE.L_workers=@(l,aprime,eprime,a,e,eta,theta) (e==0)*l*eta; % labour supplied to the market (entrepreneurs' own eta goes to their own firm)
FnsToEvaluateGE.N_hired=@(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi) EntreHiredLabour_d_with2A(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi); % labour entrepreneurs hire in

% Intermediate eqns are evaluated after the AggVars are written into Parameters and before the
% GeneralEqmEqns, with each result written back as a named parameter (see the useintermediateEqns
% block in TransitionPath_InfHorz_shooting / HeteroAgentStationaryEqm_InfHorz). That is exactly
% what is needed to form the corporate residuals once, and use them in both GE conditions.
intermediateEqnsGE.CorporateCapital=@(A,K_noncorp) A-K_noncorp;
intermediateEqnsGE.CorporateLabour=@(L_workers,N_hired) L_workers-N_hired;

GeneralEqmEqnsGE.CapitalMarket=@(r,CorporateCapital,CorporateLabour,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(CorporateCapital^(firmalpha-1))*(CorporateLabour^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,CorporateCapital,CorporateLabour,firmalpha,firmA) w-((1-firmalpha)*firmA*(CorporateCapital^firmalpha)*(CorporateLabour^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w
heteroagentoptionsGE.intermediateEqns=intermediateEqnsGE;

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
transpathoptionsGE.intermediateEqns=intermediateEqnsGE; % same corporate residuals as the stationary solve

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

clear FnsToEvaluateGE vfoptionsGE_GI simoptionsGE_GI p_eqm_GI GEcondns_GI ParamsGE_GI V_finalGE_GI Policy_finalGE_GI AgentDist_initialGE_GI PricePathGE_GI PricePathGEOut_GI PricePathBumped_GI PricePathBumpedOut_GI GeneralEqmEqnsGE intermediateEqnsGE heteroagentoptionsGE p_eqm GEcondns ParamsGE V_finalGE Policy_finalGE AgentDist_initialGE PricePathGE ParamPathGE transpathoptionsGE PricePathGEOut PricePathBumped PricePathBumpedOut bumpr bumpw Tbump

%% Run the GE transition path with transpathoptions.maxiter=1 -- shape check only
transpathoptions.maxiter=1;
GeneralEqmEqns.dummy=@(entrepreneur) 0;

PricePath2=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);

% Big grid + DC
PricePath3A=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_big, AgentDist_initial_big, n_d, n_a_big, n_z, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions2, vfoptions2, []);

% Big grid + DC + GI
PricePath3B=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_big_GI, AgentDist_initial_big, n_d, n_a_big, n_z, d_grid,a_grid_big,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions4, vfoptions4, []);

fprintf('with2A: One iter of TPath, big-grid with/without GI, this should be close to zero: %2.8f \n',max(abs(PricePath3A.r-PricePath3B.r)))

%%
output=struct();

end
