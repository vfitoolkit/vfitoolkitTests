function output=UncertainTPathFHorz_LocalConsistency(T,rS,wS,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptions)
% Local consistency of the tree's own answer, on a full-support tree with a decision variable (labour), for the FHorz
% uncertain commands. Written 2026-11-05, before the first run.
%   (a) every branching node is one step from the expectation of its children, for every age: a 2-period
%       ValueFnOnTransPath_Case1_FHorz from sum_c p_c V_c, as an identity
%   (b) the children of one node start from the same distribution (the demographics do not depend on S), exactly
%   (c) every node's distribution has the age weights as the mass of each age
%   (d) for scale, the two children of a branching node have different V (the regimes differ)

vfoptions=struct();
simoptions=struct();
N_a=prod(n_a); N_z=prod(n_z);
NV=N_a*N_z*N_j;
ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(1,ceil(n_z/2))=1;

pi_S=[0.75,0.25;0.25,0.75]; T1=3; T2=5;
V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end
Params1=Params; Params1.r=rS(1,1); Params1.w=wS(1,1);
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
Policyshape=size(Policy_final.S1);
NP=prod(Policyshape);

PricePathU.r=rS;
ParamPathU.w=wS;
[VPathU,PolicyPathU]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDistPathU=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PolicyPathU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
[Tree,PriceNodes,ParamNodes]=UncertainTransPathTree(T,T1,T2,pi_S,1,PricePathU,ParamPathU);
nNodes=Tree.nNodes;
VU=gather(reshape(VPathU,[NV,nNodes]));
PU=gather(reshape(PolicyPathU,[NP,nNodes]));
DU=gather(reshape(AgentDistPathU,[N_a*N_z,N_j,nNodes]));
fprintf('Uncertain FHorz local consistency: %i histories, %i nodes \n',Tree.nHist,nNodes)

%% (a) and (b) at every branching node
errV=0; errstart=0; scaleV=0;
for kk=find(Tree.nodenchild>1)'
    children=Tree.nodefirstchild(kk)+(0:Tree.nodenchild(kk)-1);
    Vnext=sum(VU(:,children).*reshape(Tree.nodecondprob(children),1,[]),2); % the command's form
    clear PricePathH ParamPathH
    PricePathH.r=PriceNodes.r([kk,kk]);
    ParamPathH.w=ParamNodes.w([kk,kk]);
    Vh=ValueFnOnTransPath_Case1_FHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vnext,[n_a,n_z,N_j])), reshape(PU(:,children(1)),Policyshape), Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    Vh=gather(reshape(Vh,[NV,2]));
    dd=abs(Vh(:,1)-VU(:,kk)); dd(Vh(:,1)==VU(:,kk))=0; dd(isnan(dd))=Inf;
    errV=max(errV,max(dd));
    errstart=max(errstart,max(abs(DU(:,:,children(1))-DU(:,:,children(2))),[],'all'));
    dd=abs(VU(:,children(1))-VU(:,children(2))); dd(isnan(dd))=0;
    scaleV=max(scaleV,max(dd));
end
fprintf('Uncertain FHorz local consistency (a): every branching node is one step from sum_c p_c V_c (every age), this should be zero: %.3e \n',errV)
fprintf('Uncertain FHorz local consistency (b): the children of a node start from the same distribution, this should be zero: %.3e \n',errstart)
fprintf('Uncertain FHorz local consistency (d): for scale, the children of a branching node have V differing by, should NOT be zero: %.3e \n',scaleV)
%% (c) the mass of each age at every node is its age weight
massU=reshape(sum(DU,1),[N_j,nNodes]);
fprintf('Uncertain FHorz local consistency (c): at every node the mass of each age is its age weight, this should be zero: %.3e \n',max(abs(massU-Params.mewj(:)),[],'all'))

output=struct();

end
