function output=CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest: RiskyAsset vs degenerate ExperienceAssetu (nod1, +semiz, WITH a1).
% a=[a1,a2]: a1=standard safe asset (a1prime chosen directly), a2=risky/exp asset (set by aprimeFn).
% n_a=[n_a1,n_a2]; a_grid=[a1_grid;a2_grid]. BOTH sides carry a genuine employment semiz.
%   R: riskyasset, refine_d=[0,1,1,1].   E: experienceassetu, l_dexperienceassetu=2.
%   E aprimeFn ignores its a2 input; E ReturnFn ignores riskyshare (aprimeFn-only on side R).
% Only 2 shock blocks: (noz,noe,+semiz) and (noz,e,+semiz).

%% Side R (riskyasset)
vfoptionsR=struct();
vfoptionsR.riskyasset=1;
vfoptionsR.refine_d=[0,1,1,1];
vfoptionsR.aprimeFn=@(riskyshare,savings,u,r) aprimeFn_CoreTestRiskyAsset(riskyshare,savings,u,r);
vfoptionsR.n_u=vfoptionsbaseline.n_u; vfoptionsR.u_grid=vfoptionsbaseline.u_grid; vfoptionsR.pi_u=vfoptionsbaseline.pi_u;
vfoptionsR.n_semiz=vfoptionsbaseline.n_semiz; vfoptionsR.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptionsR.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsR=struct();
simoptionsR.riskyasset=1; simoptionsR.refine_d=vfoptionsR.refine_d; simoptionsR.aprimeFn=vfoptionsR.aprimeFn;
simoptionsR.n_u=vfoptionsR.n_u; simoptionsR.u_grid=vfoptionsR.u_grid; simoptionsR.pi_u=vfoptionsR.pi_u;
simoptionsR.n_semiz=simoptionsbaseline.n_semiz; simoptionsR.semiz_grid=simoptionsbaseline.semiz_grid; simoptionsR.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsR.d_grid=d_grid; simoptionsR.a_grid=a_grid;

%% Side E (degenerate expassetu)
vfoptionsE=struct();
vfoptionsE.experienceassetu=1;
vfoptionsE.l_dexperienceassetu=2;
vfoptionsE.aprimeFn=@(riskyshare,savings,a2,u,r) aprimeFn_CoreTestRiskyAsset(riskyshare,savings,u,r); % ignores a2
vfoptionsE.n_u=vfoptionsbaseline.n_u; vfoptionsE.u_grid=vfoptionsbaseline.u_grid; vfoptionsE.pi_u=vfoptionsbaseline.pi_u;
vfoptionsE.n_semiz=vfoptionsbaseline.n_semiz; vfoptionsE.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptionsE.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsE=struct();
simoptionsE.experienceassetu=1; simoptionsE.aprimeFn=vfoptionsE.aprimeFn;
simoptionsE.n_u=vfoptionsE.n_u; simoptionsE.u_grid=vfoptionsE.u_grid; simoptionsE.pi_u=vfoptionsE.pi_u;
simoptionsE.n_semiz=simoptionsbaseline.n_semiz; simoptionsE.semiz_grid=simoptionsbaseline.semiz_grid; simoptionsE.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsE.d_grid=d_grid; simoptionsE.a_grid=a_grid;

%% (1) noz, noe, +semiz
jequaloneDist=zeros([n_a,vfoptionsR.n_semiz],'gpuArray'); jequaloneDist(1,1,ceil(vfoptionsR.n_semiz/2))=1;
ReturnFn_R=@(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_E=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
FnsR.a2=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz) a2; FnsE.a2=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz) a2;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,0,N_j,[],Params,simoptionsR);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsR);
[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,0,N_j,[],Params,simoptionsE);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsE);
fprintf('CrossTest ExpAssetu (nod1 noz noe, +semiz withA1): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.a2.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.a2.Mean-AllStatsE.a2.Mean))

%% (2) noz, e, +semiz
vfoptionsR_e=vfoptionsR; vfoptionsR_e.n_e=vfoptionsbaseline.n_e; vfoptionsR_e.e_grid=vfoptionsbaseline.e_grid; vfoptionsR_e.pi_e=vfoptionsbaseline.pi_e;
simoptionsR_e=simoptionsR; simoptionsR_e.n_e=simoptionsbaseline.n_e; simoptionsR_e.e_grid=simoptionsbaseline.e_grid; simoptionsR_e.pi_e=simoptionsbaseline.pi_e;
vfoptionsE_e=vfoptionsE; vfoptionsE_e.n_e=vfoptionsbaseline.n_e; vfoptionsE_e.e_grid=vfoptionsbaseline.e_grid; vfoptionsE_e.pi_e=vfoptionsbaseline.pi_e;
simoptionsE_e=simoptionsE; simoptionsE_e.n_e=simoptionsbaseline.n_e; simoptionsE_e.e_grid=simoptionsbaseline.e_grid; simoptionsE_e.pi_e=simoptionsbaseline.pi_e;
jequaloneDist=zeros([n_a,vfoptionsR.n_semiz,vfoptionsbaseline.n_e],'gpuArray'); jequaloneDist(1,1,ceil(vfoptionsR.n_semiz/2),ceil(vfoptionsbaseline.n_e/2))=1;
ReturnFn_R=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_E=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
FnsR.a2=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,e) a2; FnsE.a2=@(riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,e) a2;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR_e);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,0,N_j,[],Params,simoptionsR_e);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsR_e);
[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE_e);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,0,N_j,[],Params,simoptionsE_e);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsE_e);
fprintf('CrossTest ExpAssetu (nod1 noz e, +semiz withA1): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.a2.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.a2.Mean-AllStatsE.a2.Mean))

%%
output=struct();

end
