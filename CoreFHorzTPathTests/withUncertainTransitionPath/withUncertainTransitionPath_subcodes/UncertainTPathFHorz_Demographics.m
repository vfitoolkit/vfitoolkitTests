function output=UncertainTPathFHorz_Demographics(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptions)
% Demographics that depend on S: age weights (mewj, given per S on the ParamPath, N_j-by-T fields) and newborns
% (jequaloneDist given per S). The policy chosen at a node is chosen before next period's S is known, so every age j>=2
% of the children of one node is the node's distribution moved on by its policy; the newborns and the age weights are
% the child's, so the children start from DIFFERENT distributions. Written 2026-11-05, before the first run.
% The prices do not depend on S here, so every node in a period has the same policy and only the demographics differ.
%   (a) a sure switch (one history) against AgentDistOnTransPath_Case1_FHorz with the mewj and jequaloneDist paths of
%       that history built by hand, as an identity
%   (b) on the absorbing tree, the mass of each age at every node is its S's age weight in its period
%   (c) the newborns (age 1) of every node, as a distribution, are its S's jequaloneDist
%   (d) the children of a branching node have the same distribution, conditional on age, at every age j>=2
%   (e) the DISCRIMINATOR: those children differ in their newborns and in their age weights
%   (f) with demographics that do not depend on S, the children of a branching node start identical, exactly

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jvec=(1:1:N_j)';

% Newborns per S: in S=2 they arrive with some assets (grid point 5)
jeq1=zeros([n_a,n_z]); jeq1(1,ceil(n_z/2))=1;
jeq2=zeros([n_a,n_z]); jeq2(5,ceil(n_z/2))=1;
jeqS.S1=jeq1; jeqS.S2=jeq2;
% Age weights per S: in S=2 the population is older, from the period the switch lands on
mewj1=Params.mewj(:);
mewj2=mewj1.*(1+0.05*jvec); mewj2=mewj2/sum(mewj2);
ParamPathU.w=[wS(1,:); wS(1,:)];
ParamPathU.mewj.S1=repmat(mewj1,1,T);
ParamPathU.mewj.S2=repmat(mewj2,1,T);
PricePathU.r=[rS(1,:); rS(1,:)];

V_final=struct(); Policy_final=struct();
ParamsT=Params; ParamsT.r=rS(1,T); ParamsT.w=wS(1,T);
[V_fin,Policy_fin]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsT,DiscountFactorParamNames,[],vfoptions);
V_final.S1=V_fin; V_final.S2=V_fin; Policy_final.S1=Policy_fin; Policy_final.S2=Policy_fin;
Params1=Params; Params1.r=rS(1,1); Params1.w=wS(1,1);
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jeq1,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);

%% (a) A sure switch to S=2 in period 2, against the deterministic command with the history's demographics by hand
pi_SA=[0,1;0,1];
[~,PA]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, pi_SA, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DA=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jeqS, PricePathU, ParamPathU, PA, AgeWeightParamNames, T, 1, 1, pi_SA, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
clear PricePathD ParamPathD
PricePathD.r=rS(1,:);
ParamPathD.w=wS(1,:);
ParamPathD.mewj=[mewj1, repmat(mewj2,1,T-1)];
jeqD=zeros([n_a,n_z,T]);
jeqD(:,:,1)=jeq1;
jeqD(:,:,2:T)=repmat(jeq2,1,1,T-1);
DD=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jeqD, PricePathD, ParamPathD, gpuArray(PA), AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions, simoptions);
fprintf('Uncertain FHorz demographics (a): a sure switch with per-S age weights and newborns is the deterministic path with them by hand, AgentDist, this should be zero: %.3e \n',max(abs(gather(DA(:))-gather(DD(:)))))

%% The absorbing tree
pi_S=[0.75,0.25;0,1]; T1=2; T2=4;
[~,PU]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DU=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jeqS, PricePathU, ParamPathU, PU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
nNodes=Tree.nNodes;
DU=gather(reshape(DU,[N_a*N_z,N_j,nNodes]));
massU=reshape(sum(DU,1),[N_j,nNodes]);
% (b) mass of each age
massHand=repmat(mewj1,1,nNodes);
massHand(:,Tree.nodeS==2)=repmat(mewj2,1,sum(Tree.nodeS==2));
fprintf('Uncertain FHorz demographics (b): at every node the mass of each age is its S''s age weight, this should be zero: %.3e \n',max(abs(massU-massHand),[],'all'))
% (c) the newborns
errnew=0;
for kk=2:nNodes % period 1 is the initial dist, which is the stationary one
    newborns=DU(:,1,kk)/sum(DU(:,1,kk));
    if Tree.nodeS(kk)==1
        errnew=max(errnew,max(abs(newborns-jeq1(:))));
    else
        errnew=max(errnew,max(abs(newborns-jeq2(:))));
    end
end
fprintf('Uncertain FHorz demographics (c): the newborns of every node are its S''s jequaloneDist, this should be zero: %.3e \n',errnew)
% (d) and (e) the two children of every branching node
errcond=0; scalenew=0; scalemass=0;
for kk=find(Tree.nodenchild>1)'
    children=Tree.nodefirstchild(kk)+(0:Tree.nodenchild(kk)-1);
    cond1=DU(:,2:N_j,children(1))./sum(DU(:,2:N_j,children(1)),1);
    cond2=DU(:,2:N_j,children(2))./sum(DU(:,2:N_j,children(2)),1);
    errcond=max(errcond,max(abs(cond1-cond2),[],'all'));
    scalenew=max(scalenew,max(abs(DU(:,1,children(1))/sum(DU(:,1,children(1)))-DU(:,1,children(2))/sum(DU(:,1,children(2))))));
    scalemass=max(scalemass,max(abs(massU(:,children(1))-massU(:,children(2)))));
end
fprintf('Uncertain FHorz demographics (d): the children of a branching node, conditional on age, are the same at every age j>=2, this should be zero: %.3e \n',errcond)
fprintf('Uncertain FHorz demographics (e): for scale, their newborns differ by, should NOT be zero: %.3e \n',scalenew)
fprintf('Uncertain FHorz demographics (e): for scale, their age weights differ by, should NOT be zero: %.3e \n',scalemass)

%% (f) Demographics that do not depend on S: the children start identical
ParamPathF.w=ParamPathU.w;
DF=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jeq1, PricePathU, ParamPathF, PU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
DF=gather(reshape(DF,[N_a*N_z*N_j,nNodes]));
errF=0;
for kk=find(Tree.nodenchild>1)'
    children=Tree.nodefirstchild(kk)+(0:Tree.nodenchild(kk)-1);
    errF=max(errF,max(abs(DF(:,children(1))-DF(:,children(2)))));
end
fprintf('Uncertain FHorz demographics (f): with demographics that do not depend on S, the children of a branching node start identical, this should be zero: %.3e \n',errF)

output=struct();

end
