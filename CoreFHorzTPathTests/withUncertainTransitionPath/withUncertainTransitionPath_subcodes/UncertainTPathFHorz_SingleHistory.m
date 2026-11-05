function output=UncertainTPathFHorz_SingleHistory(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptions)
% A tree with a single history IS a deterministic transition: ValueFnOnUncertainTransPath_FHorz,
% AgentDistOnUncertainTransPath_FHorz and EvalFnOnUncertainTransPath_AggVars_FHorz against ValueFnOnTransPath_Case1_FHorz,
% AgentDistOnTransPath_Case1_FHorz and EvalFnOnTransPath_AggVars_Case1_FHorz on the same path, as identities.
% Written 2026-11-05, before the first run.
%   (a) N_S=1
%   (b) a sure switch to S=2 in period 2, with an age-dependent parameter (kappa_j, N_j-by-T) given per S
%   (c) the same with a 3-D pi_S (all slices equal), against (b)
%   (d) fastOLG=0 against fastOLG=1 on the tree of (b)
% The model is the parent bank's nod z: households with an age profile of earnings kappa_j, retirement and a pension.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
jequaloneDist=zeros([n_a,n_z]);
jequaloneDist(1,ceil(n_z/2))=1;
AggNames=fieldnames(FnsToEvaluate);
tvec=1:1:T;

% kappa_j per S: in S=2 every working age earns 10 percent less, from period 2 on
kappa1=repmat(Params.kappa_j(:),1,T);
kappa2=0.9*kappa1;

% The initial agent dist: the stationary dist at the period 1 (S=1) prices
Params1=Params; Params1.r=rS(1,1); Params1.w=wS(1,1);
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
% The final eqm value fn of each S (at its period T prices; the S=2 one with the lower kappa_j)
V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    if SS==2
        ParamsS.kappa_j=0.9*Params.kappa_j;
    end
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end

%% (a) N_S=1
clear PricePathU ParamPathU PricePathD ParamPathD
PricePathU.r=rS(1,:);
ParamPathU.w=wS(1,:);
[VU,PU]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, 1, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DU=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PU, AgeWeightParamNames, T, 1, 1, 1, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
AU=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluate, DU, PU, PricePathU, ParamPathU, Params, T, 1, 1, 1, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptions, simoptions);
[VD,PD]=ValueFnOnTransPath_Case1_FHorz(PricePathU, ParamPathU, T, V_final.S1, Policy_final.S1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DD=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PD, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions, simoptions);
AD=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, DD, PD, PricePathU, ParamPathU, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions, simoptions);
VU=gather(VU(:)); VD=gather(VD(:));
dd=abs(VU-VD); dd(VU==VD)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz single history (a): N_S=1, V is ValueFnOnTransPath_Case1_FHorz, this should be zero: %.3e \n',max(dd))
fprintf('Uncertain FHorz single history (a): N_S=1, Policy, this should be zero: %.3e \n',max(abs(gather(PU(:))-gather(PD(:)))))
fprintf('Uncertain FHorz single history (a): N_S=1, AgentDist, this should be zero: %.3e \n',max(abs(gather(DU(:))-gather(DD(:)))))
aggerr=0;
for ff=1:length(AggNames)
    aggerr=max(aggerr,max(abs(AU.(AggNames{ff}).Mean(:)-gather(AD.(AggNames{ff}).Mean(:)))));
end
fprintf('Uncertain FHorz single history (a): N_S=1, AggVars, this should be zero: %.3e \n',aggerr)

%% (b) A sure switch to S=2 in period 2, with kappa_j per S
pi_SB=[0,1;0,1];
clear PricePathU ParamPathU PricePathD ParamPathD
PricePathU.r=rS;
ParamPathU.w=wS;
ParamPathU.kappa_j.S1=kappa1;
ParamPathU.kappa_j.S2=kappa2;
Spath=[1,2*ones(1,T-1)];
idx=Spath+2*(tvec-1);
PricePathD.r=rS(idx);
ParamPathD.w=wS(idx);
ParamPathD.kappa_j=kappa1;
ParamPathD.kappa_j(:,Spath==2)=kappa2(:,Spath==2);
[VUb,PUb]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, pi_SB, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DUb=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PUb, AgeWeightParamNames, T, 1, 1, pi_SB, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
AUb=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluate, DUb, PUb, PricePathU, ParamPathU, Params, T, 1, 1, pi_SB, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptions, simoptions);
[VDb,PDb]=ValueFnOnTransPath_Case1_FHorz(PricePathD, ParamPathD, T, V_final.S2, Policy_final.S2, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DDb=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePathD, ParamPathD, PDb, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions, simoptions);
ADb=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, DDb, PDb, PricePathD, ParamPathD, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions, simoptions);
VUbv=gather(VUb(:)); VDbv=gather(VDb(:));
dd=abs(VUbv-VDbv); dd(VUbv==VDbv)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz single history (b): a sure switch with kappa_j per S (N_j-by-T fields), V, this should be zero: %.3e \n',max(dd))
fprintf('Uncertain FHorz single history (b): Policy, this should be zero: %.3e \n',max(abs(gather(PUb(:))-gather(PDb(:)))))
fprintf('Uncertain FHorz single history (b): AgentDist, this should be zero: %.3e \n',max(abs(gather(DUb(:))-gather(DDb(:)))))
aggerr=0;
for ff=1:length(AggNames)
    aggerr=max(aggerr,max(abs(AUb.(AggNames{ff}).Mean(:)-gather(ADb.(AggNames{ff}).Mean(:)))));
end
fprintf('Uncertain FHorz single history (b): AggVars, this should be zero: %.3e \n',aggerr)
fprintf('Uncertain FHorz single history (b): for scale, earnings move with the switch (kappa_j is cut) by, should NOT be zero: %.3e \n',abs(AUb.earnings.Mean(2)-AUb.earnings.Mean(1)))

%% (c) A 3-D pi_S, all slices equal, gives the tree of (b) (T2=2: two periods of draws, both sure to be S=2)
[VUc,PUc]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 2, cat(3,pi_SB,pi_SB), 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VUcv=gather(VUc(:));
dd=abs(VUcv-VUbv); dd(VUcv==VUbv)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz single history (c): a 3-D pi_S (two equal slices) gives the single history of (b), V, this should be zero: %.3e \n',max(dd))
fprintf('Uncertain FHorz single history (c): and Policy, this should be zero: %.3e \n',max(abs(gather(PUc(:))-gather(PUb(:)))))

%% (d) fastOLG=0 against fastOLG=1 on the tree of (b)
transpathoptions0=transpathoptions; transpathoptions0.fastOLG=0;
transpathoptions1=transpathoptions; transpathoptions1.fastOLG=1;
[V0,P0]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, pi_SB, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions0, vfoptions);
[V1,P1]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, pi_SB, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions1, vfoptions);
D0=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, P0, AgeWeightParamNames, T, 1, 1, pi_SB, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions0, simoptions);
D1=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, P1, AgeWeightParamNames, T, 1, 1, pi_SB, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions1, simoptions);
V0=gather(V0(:)); V1=gather(V1(:));
dd=abs(V0-V1); dd(V0==V1)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz single history (d): fastOLG=0 vs 1 on the tree, V, this should be zero: %.3e \n',max(dd))
fprintf('Uncertain FHorz single history (d): fastOLG=0 vs 1 on the tree, Policy, this should be zero: %.3e \n',max(abs(gather(P0(:))-gather(P1(:)))))
fprintf('Uncertain FHorz single history (d): fastOLG=0 vs 1 on the tree, AgentDist, this should be zero: %.3e \n',max(abs(gather(D0(:))-gather(D1(:)))))

output=struct();

end
