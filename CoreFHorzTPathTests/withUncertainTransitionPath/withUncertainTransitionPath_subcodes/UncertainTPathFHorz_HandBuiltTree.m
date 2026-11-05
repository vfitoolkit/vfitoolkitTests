function output=UncertainTPathFHorz_HandBuiltTree(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptions)
% The tree built BY HAND from the deterministic FHorz commands, against the uncertain FHorz commands. Absorbing S:
% pi_S=[P,1-P;0,1], S=1 until it switches to 2 for good, which can land in any period D1,...,D2 (agents face the
% uncertainty in T1=D1-1,...,T2=D2-1). So there is path 0 (never switches) and one branch per switch date s. The FHorz
% counterpart of UncertainTPath_HandBuiltTree (InfHorz). Written 2026-11-05, before the first run.
%   Branch s: periods s,...,T, the S=2 prices, a deterministic transition to the S=2 final eqm, starting from path 0's
%     agent dist in period s (the demographics do not depend on S, so that is the child's start).
%   Path 0: from D2 on a deterministic transition to the S=1 final eqm; in each period t=D2-1,...,D1-1, one period with
%     next period value fn P*V0_{t+1}+(1-P)*Vbranch(t+1)_{t+1} for every age; before that, deterministic.
% Every cohort alive at a switch date ages across it, so this checks the expectation at every age, the start of every
% branch, and the numbering, independently of the tree's code. P=0.75 so that P and 1-P are exact; the expectation by
% hand is written in the same order as the command's (S=1 child first), so V is compared as an identity.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);
NV=N_a*N_z*N_j;
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
AggNames=fieldnames(FnsToEvaluate);
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(1,ceil(n_z/2))=1;

P=0.75;
pi_S=[P,1-P;0,1];
D1=3; D2=6; % the periods the switch can land in
T1=D1-1; T2=D2-1;

V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end
Params1=Params; Params1.r=rS(1,1); Params1.w=wS(1,1);
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
Policyshape=size(Policy_final.S1); % [1,n_a,n_z,N_j]
NP=prod(Policyshape);

%% By hand: backward, the branches
VBr=cell(D2,1); PolicyBr=cell(D2,1);
for s=D1:D2
    clear PricePathB ParamPathB
    PricePathB.r=rS(2,s:T);
    ParamPathB.w=wS(2,s:T);
    [VBr{s},PolicyBr{s}]=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathB, T-s+1, V_final.S2, Policy_final.S2, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    VBr{s}=gather(reshape(VBr{s},[NV,T-s+1]));
    PolicyBr{s}=gather(reshape(PolicyBr{s},[NP,T-s+1]));
end
%% By hand: backward, path 0
V0=zeros(NV,T); Policy0=zeros(NP,T);
clear PricePathB ParamPathB
PricePathB.r=rS(1,D2:T);
ParamPathB.w=wS(1,D2:T);
[Vtemp,Policytemp]=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathB, T-D2+1, V_final.S1, Policy_final.S1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
V0(:,D2:T)=gather(reshape(Vtemp,[NV,T-D2+1]));
Policy0(:,D2:T)=gather(reshape(Policytemp,[NP,T-D2+1]));
for tt=D2-1:-1:1
    if tt+1>=D1
        Vnext=sum([V0(:,tt+1),VBr{tt+1}(:,1)].*[P,1-P],2); % the command's order and form: the S=1 child, then the S=2 child
    else
        Vnext=V0(:,tt+1);
    end
    clear PricePathB ParamPathB
    PricePathB.r=rS(1,[tt,tt]);
    ParamPathB.w=wS(1,[tt,tt]);
    [Vtemp,Policytemp]=ValueFnOnTransPath_Case1_FHorz(PricePathB, ParamPathB, 2, gpuArray(reshape(Vnext,[n_a,n_z,N_j])), reshape(Policy0(:,tt+1),Policyshape), Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    Vtemp=gather(reshape(Vtemp,[NV,2]));
    Policytemp=gather(reshape(Policytemp,[NP,2]));
    V0(:,tt)=Vtemp(:,1);
    Policy0(:,tt)=Policytemp(:,1);
end
%% By hand: forward, path 0, then each branch from path 0's dist in its switch period
clear PricePath0 ParamPath0
PricePath0.r=rS(1,:);
ParamPath0.w=wS(1,:);
AgentDist0=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath0, ParamPath0, gpuArray(reshape(Policy0,[Policyshape,T])), AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions, simoptions);
AggVars0=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, AgentDist0, gpuArray(reshape(Policy0,[Policyshape,T])), PricePath0, ParamPath0, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions, simoptions);
AgentDist0=gather(reshape(AgentDist0,[NV,T]));
AgentDistBr=cell(D2,1); AggVarsBr=cell(D2,1);
for s=D1:D2
    clear PricePathB ParamPathB
    PricePathB.r=rS(2,s:T);
    ParamPathB.w=wS(2,s:T);
    temp=AgentDistOnTransPath_Case1_FHorz(reshape(AgentDist0(:,s),[n_a,n_z,N_j]), jequaloneDist, PricePathB, ParamPathB, gpuArray(reshape(PolicyBr{s},[Policyshape,T-s+1])), AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T-s+1,Params, transpathoptions, simoptions);
    AggVarsBr{s}=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, temp, gpuArray(reshape(PolicyBr{s},[Policyshape,T-s+1])), PricePathB, ParamPathB, Params, T-s+1, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions, simoptions);
    AgentDistBr{s}=gather(reshape(temp,[NV,T-s+1]));
end

%% The uncertain commands
PricePathU.r=rS;
ParamPathU.w=wS;
[VPathU,PolicyPathU]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPathU=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PolicyPathU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
AggVarsTreeU=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluate, AgentDistPathU, PolicyPathU, PricePathU, ParamPathU, Params, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptions, simoptions);
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
nNodes=Tree.nNodes;
VU=gather(reshape(VPathU,[NV,nNodes]));
PU=gather(reshape(PolicyPathU,[NP,nNodes]));
DU=gather(reshape(AgentDistPathU,[NV,nNodes]));

%% Compare, history by history
fprintf('Uncertain FHorz hand-built tree: the tree has D2-D1+2 histories, should be one: %i \n',Tree.nHist==D2-D1+2)
hh0=find(all(Tree.historyS==1,2));
nodes0=Tree.histories(hh0,:);
dd=abs(VU(:,nodes0)-V0); dd(VU(:,nodes0)==V0)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz hand-built tree: path 0, V (every age), this should be zero: %.3e \n',max(dd,[],'all'))
fprintf('Uncertain FHorz hand-built tree: path 0, Policy, this should be zero: %.3e \n',max(abs(PU(:,nodes0)-Policy0),[],'all'))
fprintf('Uncertain FHorz hand-built tree: path 0, AgentDist, this should be zero: %.3e \n',max(abs(DU(:,nodes0)-AgentDist0),[],'all'))
aggerr=0;
for ff=1:length(AggNames)
    aggerr=max(aggerr,max(abs(AggVarsTreeU.(AggNames{ff}).Mean(nodes0)-gather(AggVars0.(AggNames{ff}).Mean(:))')));
end
fprintf('Uncertain FHorz hand-built tree: path 0, AggVars, this should be zero: %.3e \n',aggerr)
Verr=0; Perr=0; Derr=0; aggerr=0; probok=1;
for s=D1:D2
    hh=find(Tree.historyS(:,s)==2 & Tree.historyS(:,s-1)==1);
    probok=probok*(abs(Tree.historyprob(hh)-P^(s-D1)*(1-P))<10^(-15));
    nodesb=Tree.histories(hh,s:T);
    dd=abs(VU(:,nodesb)-VBr{s}); dd(VU(:,nodesb)==VBr{s})=0; dd(isnan(dd))=Inf;
    Verr=max(Verr,max(dd,[],'all'));
    Perr=max(Perr,max(abs(PU(:,nodesb)-PolicyBr{s}),[],'all'));
    Derr=max(Derr,max(abs(DU(:,nodesb)-AgentDistBr{s}),[],'all'));
    for ff=1:length(AggNames)
        aggerr=max(aggerr,max(abs(AggVarsTreeU.(AggNames{ff}).Mean(nodesb)-gather(AggVarsBr{s}.(AggNames{ff}).Mean(:))')));
    end
    probok=probok*isequal(Tree.histories(hh,1:s-1),nodes0(1:s-1));
end
fprintf('Uncertain FHorz hand-built tree: every branch, V, this should be zero: %.3e \n',Verr)
fprintf('Uncertain FHorz hand-built tree: every branch, Policy, this should be zero: %.3e \n',Perr)
fprintf('Uncertain FHorz hand-built tree: every branch, AgentDist, this should be zero: %.3e \n',Derr)
fprintf('Uncertain FHorz hand-built tree: every branch, AggVars, this should be zero: %.3e \n',aggerr)
fprintf('Uncertain FHorz hand-built tree: each branch has probability P^(s-D1)*(1-P) and shares path 0''s nodes before s, should be one: %i \n',probok)
% For scale: the possibility of the switch moves path 0 away from the deterministic S=1 path before D2
[VD,~]=ValueFnOnTransPath_Case1_FHorz(PricePath0, ParamPath0, T, V_final.S1, Policy_final.S1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VD=gather(reshape(VD,[NV,T]));
dd=abs(V0(:,1:D2-1)-VD(:,1:D2-1)); dd(V0(:,1:D2-1)==VD(:,1:D2-1))=0; dd(isnan(dd))=0;
fprintf('Uncertain FHorz hand-built tree: for scale, path 0 vs the deterministic S=1 path before D2, V, should NOT be zero: %.3e \n',max(dd,[],'all'))
dd=abs(V0(:,D2:T)-VD(:,D2:T)); dd(V0(:,D2:T)==VD(:,D2:T))=0;
fprintf('Uncertain FHorz hand-built tree: and from D2 on they are the same deterministic path, this should be zero: %.3e \n',max(dd,[],'all'))

output=struct();

end
