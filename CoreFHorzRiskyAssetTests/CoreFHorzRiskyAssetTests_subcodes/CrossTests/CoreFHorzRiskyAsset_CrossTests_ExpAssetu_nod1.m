function output=CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest: RiskyAsset vs a degenerate ExperienceAssetu (nod1, nosemiz, single asset).
%
% Idea: riskyasset's aprimeFn is aprimeFn(d2,d3,u)   -- next asset from decisions and u, NOT from current a.
%       expassetu's aprimeFn is aprimeFn(d,a2,u)     -- next asset from decisions, current a2, and u.
% If we take an expassetu whose aprimeFn simply IGNORES its a2 input, it reproduces the riskyasset
% law of motion exactly. Both toolkits interpolate the continuous aprime onto a2_grid with the
% identical lower-point-probability scheme (CreateRiskyAssetFnMatrix vs CreateExperienceAssetuFnMatrix),
% so V, Policy and StationaryDist should match to machine precision.
%
% Mapping (Side R = riskyasset, Side E = degenerate expassetu):
%   decisions : n_d=[n_d2,n_d3]=[riskyshare,savings] on both sides, same d_grid, same order
%   R: refine_d=[0,1,1]           (d2=riskyshare in aprimeFn only, d3=savings in both)
%   E: l_dexperienceassetu=2      (BOTH decisions feed the aprimeFn; nod1)
%   R aprimeFn: @(riskyshare,savings,u,r)
%   E aprimeFn: @(riskyshare,savings,a2,u,r)   -- ignores a2
%   R ReturnFn: @(savings,a,...)
%   E ReturnFn: @(riskyshare,savings,a,...)     -- ignores riskyshare
%
% n_d=n_d_withoutd1; d_grid=d_grid_withoutd1. Single asset (the riskyasset), so noa1.

%% Common riskyasset (Side R) options
vfoptionsR=struct();
vfoptionsR.riskyasset=1;
vfoptionsR.refine_d=[0,1,1];
vfoptionsR.aprimeFn=@(riskyshare,savings,u,r) aprimeFn_CoreTestRiskyAsset(riskyshare,savings,u,r);
vfoptionsR.n_u=vfoptionsbaseline.n_u;
vfoptionsR.u_grid=vfoptionsbaseline.u_grid;
vfoptionsR.pi_u=vfoptionsbaseline.pi_u;
simoptionsR=struct();
simoptionsR.riskyasset=1;
simoptionsR.refine_d=vfoptionsR.refine_d;
simoptionsR.aprimeFn=vfoptionsR.aprimeFn;
simoptionsR.n_u=vfoptionsR.n_u;
simoptionsR.u_grid=vfoptionsR.u_grid;
simoptionsR.pi_u=vfoptionsR.pi_u;
simoptionsR.d_grid=d_grid;
simoptionsR.a_grid=a_grid;

%% Common degenerate expassetu (Side E) options
vfoptionsE=struct();
vfoptionsE.experienceassetu=1;
vfoptionsE.l_dexperienceassetu=2; % BOTH decisions (riskyshare,savings) feed the aprimeFn
vfoptionsE.aprimeFn=@(riskyshare,savings,a2,u,r) aprimeFn_CoreTestRiskyAsset(riskyshare,savings,u,r); % ignores a2
vfoptionsE.n_u=vfoptionsbaseline.n_u;
vfoptionsE.u_grid=vfoptionsbaseline.u_grid;
vfoptionsE.pi_u=vfoptionsbaseline.pi_u;
simoptionsE=struct();
simoptionsE.experienceassetu=1;
simoptionsE.aprimeFn=vfoptionsE.aprimeFn;
simoptionsE.n_u=vfoptionsE.n_u;
simoptionsE.u_grid=vfoptionsE.u_grid;
simoptionsE.pi_u=vfoptionsE.pi_u;
simoptionsE.d_grid=d_grid;
simoptionsE.a_grid=a_grid;

%% (1) noz, noe
jequaloneDist=zeros(n_a,1,'gpuArray'); jequaloneDist(1)=1;

ReturnFn_R=@(savings,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_E=@(riskyshare,savings,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension);
FnsR.assets=@(riskyshare,savings,a) a;
FnsE.assets=@(riskyshare,savings,a) a;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,0,N_j,[],Params,simoptionsR);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsR);

[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,0,N_j,[],Params,simoptionsE);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsE);

fprintf('CrossTest ExpAssetu (nod1 noz noe): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.assets.Mean-AllStatsE.assets.Mean))

%% (2) z, noe
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn_R=@(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_E=@(riskyshare,savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsR.assets=@(riskyshare,savings,a,z) a;
FnsE.assets=@(riskyshare,savings,a,z) a;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsR);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsR);

[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsE);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsE);

fprintf('CrossTest ExpAssetu (nod1 z noe): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.assets.Mean-AllStatsE.assets.Mean))

%% (3) noz, e
vfoptionsR_e=vfoptionsR; vfoptionsR_e.n_e=vfoptionsbaseline.n_e; vfoptionsR_e.e_grid=vfoptionsbaseline.e_grid; vfoptionsR_e.pi_e=vfoptionsbaseline.pi_e;
simoptionsR_e=simoptionsR; simoptionsR_e.n_e=simoptionsbaseline.n_e; simoptionsR_e.e_grid=simoptionsbaseline.e_grid; simoptionsR_e.pi_e=simoptionsbaseline.pi_e;
vfoptionsE_e=vfoptionsE; vfoptionsE_e.n_e=vfoptionsbaseline.n_e; vfoptionsE_e.e_grid=vfoptionsbaseline.e_grid; vfoptionsE_e.pi_e=vfoptionsbaseline.pi_e;
simoptionsE_e=simoptionsE; simoptionsE_e.n_e=simoptionsbaseline.n_e; simoptionsE_e.e_grid=simoptionsbaseline.e_grid; simoptionsE_e.pi_e=simoptionsbaseline.pi_e;

jequaloneDist=zeros([n_a,vfoptionsbaseline.n_e],'gpuArray'); jequaloneDist(1,ceil(vfoptionsbaseline.n_e/2))=1;

ReturnFn_R=@(savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_E=@(riskyshare,savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsR.assets=@(riskyshare,savings,a,e) a;
FnsE.assets=@(riskyshare,savings,a,e) a;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR_e);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,0,N_j,[],Params,simoptionsR_e);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsR_e);

[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE_e);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,0,N_j,[],Params,simoptionsE_e);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsE_e);

fprintf('CrossTest ExpAssetu (nod1 noz e): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.assets.Mean-AllStatsE.assets.Mean))

%% (4) z, e
jequaloneDist=zeros([n_a,n_z,vfoptionsbaseline.n_e],'gpuArray'); jequaloneDist(1,ceil(n_z/2),ceil(vfoptionsbaseline.n_e/2))=1;

ReturnFn_R=@(savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_E=@(riskyshare,savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsR.assets=@(riskyshare,savings,a,z,e) a;
FnsE.assets=@(riskyshare,savings,a,z,e) a;

[VR,PolicyR]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_R,Params,DiscountFactorParamNames,[],vfoptionsR_e);
StationaryDistR=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyR,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsR_e);
AllStatsR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistR,PolicyR,FnsR,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsR_e);

[VE,PolicyE]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_E,Params,DiscountFactorParamNames,[],vfoptionsE_e);
StationaryDistE=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyE,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsE_e);
AllStatsE=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDistE,PolicyE,FnsE,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsE_e);

fprintf('CrossTest ExpAssetu (nod1 z e): should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(VR(:)-VE(:))),max(abs(PolicyR(:)-PolicyE(:))),max(abs(StationaryDistR(:)-StationaryDistE(:))),abs(AllStatsR.assets.Mean-AllStatsE.assets.Mean))

%%
output=struct();

end
