function output=CoreFHorzRiskyAsset_CrossTests_ExpAssetu_d1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest: RiskyAsset vs a degenerate ExperienceAssetu (d1, +semiz, single asset).
% As CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1_semiz but with a d1 (labour supply h).
%   n_d=[n_d1,n_d2,n_d3,n_d4]=[h,riskyshare,savings,dsemiz].
%   R: refine_d=[1,1,1,1]        (h ReturnFn-only, riskyshare aprimeFn-only, savings both, dsemiz the semiz decision)
%   E: l_dexperienceassetu=2     (riskyshare+savings feed aprimeFn; h is d1; trailing dsemiz is the semiz decision)
%   E aprimeFn ignores a2; E ReturnFn ignores riskyshare (aprimeFn-only on side R).
% Only 2 shock blocks: (noz,noe,+semiz) and (noz,e,+semiz).

%% Side R (riskyasset)
vfoptionsR=struct();
vfoptionsR.riskyasset=1;
vfoptionsR.refine_d=[1,1,1,1];
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
jequaloneDist=zeros([n_a,vfoptionsR.n_semiz],'gpuArray'); jequaloneDist(1,ceil(vfoptionsR.n_semiz/2))=1;
ReturnFn_R=@(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_E=@(h,riskyshare,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
FnsR.assets=@(h,riskyshare,savings,dsemiz,a,semiz) a; FnsE.assets=@(h,riskyshare,savings,dsemiz,a,semiz) a;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,0,N_j,[],Params,simoptionsR);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsR);
[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,0,N_j,[],Params,simoptionsE);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsE);
fprintf('CrossTest ExpAssetu (d1 noz noe, +semiz): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.assets.Mean-AllStatsE.assets.Mean))

%% (2) noz, e, +semiz
vfoptionsR_e=vfoptionsR; vfoptionsR_e.n_e=vfoptionsbaseline.n_e; vfoptionsR_e.e_grid=vfoptionsbaseline.e_grid; vfoptionsR_e.pi_e=vfoptionsbaseline.pi_e;
simoptionsR_e=simoptionsR; simoptionsR_e.n_e=simoptionsbaseline.n_e; simoptionsR_e.e_grid=simoptionsbaseline.e_grid; simoptionsR_e.pi_e=simoptionsbaseline.pi_e;
vfoptionsE_e=vfoptionsE; vfoptionsE_e.n_e=vfoptionsbaseline.n_e; vfoptionsE_e.e_grid=vfoptionsbaseline.e_grid; vfoptionsE_e.pi_e=vfoptionsbaseline.pi_e;
simoptionsE_e=simoptionsE; simoptionsE_e.n_e=simoptionsbaseline.n_e; simoptionsE_e.e_grid=simoptionsbaseline.e_grid; simoptionsE_e.pi_e=simoptionsbaseline.pi_e;
jequaloneDist=zeros([n_a,vfoptionsR.n_semiz,vfoptionsbaseline.n_e],'gpuArray'); jequaloneDist(1,ceil(vfoptionsR.n_semiz/2),ceil(vfoptionsbaseline.n_e/2))=1;
ReturnFn_R=@(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_E=@(h,riskyshare,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
FnsR.assets=@(h,riskyshare,savings,dsemiz,a,semiz,e) a; FnsE.assets=@(h,riskyshare,savings,dsemiz,a,semiz,e) a;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR_e);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,0,N_j,[],Params,simoptionsR_e);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsR_e);
[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE_e);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,0,N_j,[],Params,simoptionsE_e);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsE_e);
fprintf('CrossTest ExpAssetu (d1 noz e, +semiz): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.assets.Mean-AllStatsE.assets.Mean))

%%
output=struct();

end
