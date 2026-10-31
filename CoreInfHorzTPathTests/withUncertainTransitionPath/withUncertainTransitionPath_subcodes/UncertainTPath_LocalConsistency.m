function output=UncertainTPath_LocalConsistency(T,rS,wS,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptions)
% Local consistency of the tree's own answer, on a full support tree (every S can follow every S) with d and regimes
% that differ, where nothing can be compared to a single deterministic path:
%   (1) every branching node is one period of the deterministic step from the expectation of its children's V
%   (2) in period T every node is at the final eqm of its S
%   (3) every segment that runs to T, re-solved on its own from the final eqm of its S, is the tree on that segment
%   (4) the children of a node all start from the same distribution, its distribution moved on one period by its policy
%   (5) every distribution on the tree has mass one

vfoptions=struct();
simoptions=struct();
ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
FnsToEvaluate.assets=@(d,aprime,a,z) a;
N_a=prod(n_a); N_z=prod(n_z);

pi_S=[0.5,0.5;0.5,0.5];
T1=3; T2=5;

V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end
ParamsI=Params; ParamsI.r=rS(1,1); ParamsI.w=wS(1,1);
[~,Policy_initial]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsI,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_initial,n_d,n_a,n_z,pi_z,simoptions,ParamsI,[]);

PricePathU.r=rS;
ParamPathU.w=wS;
[VPathU,PolicyPathU]=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPathU=AgentDistOnUncertainTransPath_InfHorz(AgentDist_initial, PricePathU, ParamPathU, PolicyPathU, T, T1, T2, pi_S, 1, n_d,n_a,n_z,pi_z, Params, transpathoptions, simoptions);
[Tree,PricePathNodes,ParamPathNodes]=UncertainTransPathTree(T,T1,T2,pi_S,1,PricePathU,ParamPathU);
nNodes=Tree.nNodes;
nP=size(PolicyPathU,1);
VU=gather(reshape(VPathU,[N_a*N_z,nNodes]));
PU=gather(reshape(PolicyPathU,[nP,N_a*N_z,nNodes]));
DU=gather(reshape(AgentDistPathU,[N_a*N_z,nNodes]));

%% (1) branching nodes
branching=find(Tree.nodenchild>1)';
Verr=0; Perr=0;
for kk=branching
    children=Tree.nodefirstchild(kk):(Tree.nodefirstchild(kk)+Tree.nodenchild(kk)-1);
    Vnext=zeros(N_a*N_z,1);
    for cc=children
        Vnext=Vnext+Tree.nodecondprob(cc)*VU(:,cc);
    end
    clear PricePathB ParamPathB
    PricePathB.r=PricePathNodes.r([kk,kk]);
    ParamPathB.w=ParamPathNodes.w([kk,kk]);
    [Vtemp,Policytemp]=ValueFnOnTransPath_InfHorz(PricePathB, ParamPathB, 2, gpuArray(reshape(Vnext,[N_a,N_z])), reshape(PU(:,:,children(1)),[nP,N_a,N_z]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    Vtemp=gather(reshape(Vtemp,[N_a*N_z,2]));
    Policytemp=gather(reshape(Policytemp,[nP,N_a*N_z,2]));
    dd=abs(Vtemp(:,1)-VU(:,kk)); dd(Vtemp(:,1)==VU(:,kk))=0; dd(isnan(dd))=Inf;
    Verr=max(Verr,max(dd));
    Perr=max(Perr,max(abs(Policytemp(:,:,1)-PU(:,:,kk)),[],'all'));
end
% The expectation here is summed child by child, the command sums the children in one call: the same two terms in
% the same order, but this is an identity only up to the ULP of V, so it is checked to the bank's restol
fprintf('Local consistency: %i branching nodes, each is one step from the expectation of its children, V, this should be zero: %.3e \n',length(branching),Verr)
fprintf('Local consistency: same, Policy, this should be zero: %.3e \n',Perr)

%% (2) period T is the final eqm of each S
lastnodes=find(Tree.nodedate==T)';
Verr=0;
for kk=lastnodes
    Vf=gather(reshape(V_final.(['S',num2str(Tree.nodeS(kk))]),[N_a*N_z,1]));
    dd=abs(Vf-VU(:,kk)); dd(Vf==VU(:,kk))=0; dd(isnan(dd))=Inf;
    Verr=max(Verr,max(dd));
end
fprintf('Local consistency: every node in period T is the final eqm V of its S, this should be zero: %.3e \n',Verr)

%% (3) segments that run to T, re-solved on their own
Verr=0; Perr=0; nseg3=0;
for ss=1:Tree.nSeg
    sn=Tree.segnodes{ss};
    if Tree.nodedate(sn(end))==T
        nseg3=nseg3+1;
        clear PricePathB ParamPathB
        PricePathB.r=PricePathNodes.r(sn);
        ParamPathB.w=ParamPathNodes.w(sn);
        SS=Tree.nodeS(sn(end));
        [Vtemp,Policytemp]=ValueFnOnTransPath_InfHorz(PricePathB, ParamPathB, length(sn), V_final.(['S',num2str(SS)]), Policy_final.(['S',num2str(SS)]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
        Vtemp=gather(reshape(Vtemp,[N_a*N_z,length(sn)]));
        Policytemp=gather(reshape(Policytemp,[nP,N_a*N_z,length(sn)]));
        dd=abs(Vtemp-VU(:,sn)); dd(Vtemp==VU(:,sn))=0; dd(isnan(dd))=Inf;
        Verr=max(Verr,max(dd,[],'all'));
        Perr=max(Perr,max(abs(Policytemp-PU(:,:,sn)),[],'all'));
    end
end
fprintf('Local consistency: the %i segments that run to T, re-solved on their own, V, this should be zero: %.3e \n',nseg3,Verr)
fprintf('Local consistency: same, Policy, this should be zero: %.3e \n',Perr)

%% (4) the children of a node start from its distribution moved on one period
Derr=0; Dsib=0;
for kk=find(Tree.nodenchild>0)'
    children=Tree.nodefirstchild(kk):(Tree.nodefirstchild(kk)+Tree.nodenchild(kk)-1);
    clear PricePathB ParamPathB
    PricePathB.r=PricePathNodes.r([kk,kk]);
    ParamPathB.w=ParamPathNodes.w([kk,kk]);
    temp=AgentDistOnTransPath_InfHorz(reshape(DU(:,kk),[N_a,N_z]), PricePathB, ParamPathB, reshape(PU(:,:,[kk,kk]),[nP,N_a,N_z,2]), n_d,n_a,n_z,pi_z,2,Params,transpathoptions,simoptions);
    temp=gather(reshape(temp,[N_a*N_z,2]));
    Derr=max(Derr,max(abs(DU(:,children)-temp(:,2)),[],'all'));
    Dsib=max(Dsib,max(abs(DU(:,children)-DU(:,children(1))),[],'all'));
end
fprintf('Local consistency: the children of a node all start from the same distribution, this should be zero: %.3e \n',Dsib)
fprintf('Local consistency: that is the node''s distribution moved on one period by its policy, this should be zero: %.3e \n',Derr)

%% (5) mass one
fprintf('Local consistency: every distribution on the tree has mass one, this should be zero: %.3e \n',max(abs(sum(DU,1)-1)))
% For scale: the children of a branching node have different V (the regimes differ)
kk=branching(end);
children=Tree.nodefirstchild(kk):(Tree.nodefirstchild(kk)+Tree.nodenchild(kk)-1);
dd=abs(VU(:,children(1))-VU(:,children(2))); dd(VU(:,children(1))==VU(:,children(2)))=0;
fprintf('Local consistency: for scale, the two children of a branching node have different V, should NOT be zero: %.3e \n',max(dd))

output=struct();

end
