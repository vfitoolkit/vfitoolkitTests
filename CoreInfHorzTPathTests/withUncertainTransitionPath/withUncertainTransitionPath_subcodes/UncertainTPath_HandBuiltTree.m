function output=UncertainTPath_HandBuiltTree(T,rS,wS,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptions)
% The tree built BY HAND from the deterministic commands, the way the ABCP2024 replication's stopgap (ABCP2024_Tree.m)
% does it, against the uncertain commands. Absorbing S: pi_S=[P,1-P;0,1], S=1 until it switches to 2 for good, which
% can happen in any period T1,...,T2. So there is path 0 (never switches) and one branch per switch date s.
%   Branch s: periods s,...,T, the S=2 prices and parameters, a deterministic transition to the S=2 final eqm, starting
%     from path 0's agent dist in period s (which is fixed by path 0's policy in s-1, before S_s is known).
%   Path 0: from T2 on a deterministic transition to the S=1 final eqm; in each period t=T2-1,...,T1-1, one period
%     with next period value fn P*V0_{t+1}+(1-P)*Vbranch(t+1)_{t+1}; before that, deterministic.
% Nothing here uses the tree's numbering except to look up the histories at the end, so this checks the segments,
% the expectation at branching nodes, the starting distribution of each branch, and the numbering, independently.
% P=0.75 so that P and 1-P are exact; the hand expectation is written in the same order as the command's
% (S=1 child first), so V is compared as an identity.

vfoptions=struct();
simoptions=struct();
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.nextassets=@(aprime,a,z) aprime;
FnsToEvaluate.income=@(aprime,a,z,r,w) r*a+w*z;
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);

P=0.75;
pi_S=[P,1-P;0,1];
T1=3; T2=6;

V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end
ParamsI=Params; ParamsI.r=rS(1,1); ParamsI.w=wS(1,1);
[~,Policy_initial]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsI,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_initial,n_d,n_a,n_z,pi_z,simoptions,ParamsI,[]);

%% By hand: backward, the branches
VBr=cell(T2,1); PolicyBr=cell(T2,1);
for s=T1:T2
    clear PricePathB ParamPathB
    PricePathB.r=rS(2,s:T);
    ParamPathB.w=wS(2,s:T);
    [VBr{s},PolicyBr{s}]=ValueFnOnTransPath_InfHorz(PricePathB, ParamPathB, T-s+1, V_final.S2, Policy_final.S2, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    VBr{s}=gather(reshape(VBr{s},[N_a*N_z,T-s+1]));
    PolicyBr{s}=gather(reshape(PolicyBr{s},[1,N_a*N_z,T-s+1]));
end
%% By hand: backward, path 0
V0=zeros(N_a*N_z,T); Policy0=zeros(1,N_a*N_z,T);
clear PricePathB ParamPathB
PricePathB.r=rS(1,T2:T);
ParamPathB.w=wS(1,T2:T);
[Vtemp,Policytemp]=ValueFnOnTransPath_InfHorz(PricePathB, ParamPathB, T-T2+1, V_final.S1, Policy_final.S1, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
V0(:,T2:T)=gather(reshape(Vtemp,[N_a*N_z,T-T2+1]));
Policy0(:,:,T2:T)=gather(reshape(Policytemp,[1,N_a*N_z,T-T2+1]));
for tt=T2-1:-1:1
    if tt+1>=T1
        Vnext=P*V0(:,tt+1)+(1-P)*VBr{tt+1}(:,1); % same order as the command: the S=1 child, then the S=2 child
    else
        Vnext=V0(:,tt+1);
    end
    clear PricePathB ParamPathB
    PricePathB.r=rS(1,[tt,tt]);
    ParamPathB.w=wS(1,[tt,tt]);
    [Vtemp,Policytemp]=ValueFnOnTransPath_InfHorz(PricePathB, ParamPathB, 2, gpuArray(reshape(Vnext,[N_a,N_z])), reshape(Policy0(:,:,tt+1),[1,N_a,N_z]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    Vtemp=gather(reshape(Vtemp,[N_a*N_z,2]));
    Policytemp=gather(reshape(Policytemp,[1,N_a*N_z,2]));
    V0(:,tt)=Vtemp(:,1);
    Policy0(:,:,tt)=Policytemp(:,:,1);
end
%% By hand: forward, path 0, then each branch from path 0's dist in its switch period
clear PricePath0 ParamPath0
PricePath0.r=rS(1,:);
ParamPath0.w=wS(1,:);
AgentDist0=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath0, ParamPath0, reshape(Policy0,[1,N_a,N_z,T]), n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptions);
AgentDist0=gather(reshape(AgentDist0,[N_a*N_z,T]));
AggVars0=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, gpuArray(reshape(AgentDist0,[N_a,N_z,T])), gpuArray(reshape(Policy0,[1,N_a,N_z,T])), PricePath0, ParamPath0, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptions, simoptions);
AgentDistBr=cell(T2,1); AggVarsBr=cell(T2,1);
for s=T1:T2
    clear PricePathB ParamPathB
    PricePathB.r=rS(2,s:T);
    ParamPathB.w=wS(2,s:T);
    temp=AgentDistOnTransPath_InfHorz(reshape(AgentDist0(:,s),[N_a,N_z]), PricePathB, ParamPathB, reshape(PolicyBr{s},[1,N_a,N_z,T-s+1]), n_d,n_a,n_z,pi_z,T-s+1,Params,transpathoptions,simoptions);
    AgentDistBr{s}=gather(reshape(temp,[N_a*N_z,T-s+1]));
    AggVarsBr{s}=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, temp, gpuArray(reshape(PolicyBr{s},[1,N_a,N_z,T-s+1])), PricePathB, ParamPathB, Params, T-s+1, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptions, simoptions);
end

%% The uncertain commands
PricePathU.r=rS;
ParamPathU.w=wS;
[VPathU,PolicyPathU]=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPathU=AgentDistOnUncertainTransPath_InfHorz(AgentDist_initial, PricePathU, ParamPathU, PolicyPathU, T, T1, T2, pi_S, 1, n_d,n_a,n_z,pi_z, Params, transpathoptions, simoptions);
AggVarsTreeU=EvalFnOnUncertainTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathU, PolicyPathU, PricePathU, ParamPathU, Params, T, T1, T2, pi_S, 1, n_d,n_a,n_z, d_grid,a_grid,z_grid, transpathoptions, simoptions);
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
nNodes=Tree.nNodes;
VU=gather(reshape(VPathU,[N_a*N_z,nNodes]));
PU=gather(reshape(PolicyPathU,[1,N_a*N_z,nNodes]));
DU=gather(reshape(AgentDistPathU,[N_a*N_z,nNodes]));
AggNames=fieldnames(FnsToEvaluate);

%% Compare, history by history
fprintf('Hand-built tree: the tree has T2-T1+2 histories, should be one: %i \n',Tree.nHist==T2-T1+2)
hh0=find(all(Tree.historyS==1,2));
nodes0=Tree.histories(hh0,:);
dd=abs(VU(:,nodes0)-V0); dd(VU(:,nodes0)==V0)=0; dd(isnan(dd))=Inf;
fprintf('Hand-built tree: path 0, V, this should be zero: %.3e \n',max(dd,[],'all'))
fprintf('Hand-built tree: path 0, Policy, this should be zero: %.3e \n',max(abs(PU(:,:,nodes0)-Policy0),[],'all'))
fprintf('Hand-built tree: path 0, AgentDist, this should be zero: %.3e \n',max(abs(DU(:,nodes0)-AgentDist0),[],'all'))
aggerr=0;
for ff=1:length(AggNames)
    aggerr=max(aggerr,max(abs(AggVarsTreeU.(AggNames{ff}).Mean(nodes0)-gather(AggVars0.(AggNames{ff}).Mean(:))')));
end
fprintf('Hand-built tree: path 0, AggVars, this should be zero: %.3e \n',aggerr)
Verr=0; Perr=0; Derr=0; aggerr=0; probok=1;
for s=T1:T2
    hh=find(Tree.historyS(:,s)==2 & Tree.historyS(:,s-1)==1);
    probok=probok*(abs(Tree.historyprob(hh)-P^(s-T1)*(1-P))<10^(-15));
    nodesb=Tree.histories(hh,s:T);
    dd=abs(VU(:,nodesb)-VBr{s}); dd(VU(:,nodesb)==VBr{s})=0; dd(isnan(dd))=Inf;
    Verr=max(Verr,max(dd,[],'all'));
    Perr=max(Perr,max(abs(PU(:,:,nodesb)-PolicyBr{s}),[],'all'));
    Derr=max(Derr,max(abs(DU(:,nodesb)-AgentDistBr{s}),[],'all'));
    for ff=1:length(AggNames)
        aggerr=max(aggerr,max(abs(AggVarsTreeU.(AggNames{ff}).Mean(nodesb)-gather(AggVarsBr{s}.(AggNames{ff}).Mean(:))')));
    end
    % and before the switch, the branch IS path 0
    probok=probok*isequal(Tree.histories(hh,1:s-1),nodes0(1:s-1));
end
fprintf('Hand-built tree: every branch, V, this should be zero: %.3e \n',Verr)
fprintf('Hand-built tree: every branch, Policy, this should be zero: %.3e \n',Perr)
fprintf('Hand-built tree: every branch, AgentDist, this should be zero: %.3e \n',Derr)
fprintf('Hand-built tree: every branch, AggVars, this should be zero: %.3e \n',aggerr)
fprintf('Hand-built tree: each branch has probability P^(s-T1)*(1-P) and shares path 0''s nodes before s, should be one: %i \n',probok)
% For scale: the possibility of the switch moves path 0 away from the deterministic S=1 path before T2
[VD,~]=ValueFnOnTransPath_InfHorz(PricePath0, ParamPath0, T, V_final.S1, Policy_final.S1, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VD=gather(reshape(VD,[N_a*N_z,T]));
dd=abs(V0(:,1:T2-1)-VD(:,1:T2-1)); dd(V0(:,1:T2-1)==VD(:,1:T2-1))=0;
fprintf('Hand-built tree: for scale, path 0 vs the deterministic S=1 path before T2, V, should NOT be zero: %.3e \n',max(dd,[],'all'))
dd=abs(V0(:,T2:T)-VD(:,T2:T)); dd(V0(:,T2:T)==VD(:,T2:T))=0;
fprintf('Hand-built tree: and from T2 on they are the same deterministic path, this should be zero: %.3e \n',max(dd,[],'all'))

output=struct();

end
