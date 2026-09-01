function output=CoreFHorzRiskyAsset_CrossTests_zase_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (z-as-e) WITH a semiz background: RiskyAsset, nod1, single asset.
% Same three z-as-e equivalences as zase_nod1, but every model also carries the (employment)
% semi-exogenous state. State-dimension ordering is (a, semiz, z, e).
% Driver passes the SEMIZ n_d=[riskyshare,savings,dsemiz]; refine_d=[0,1,1,1].

% z as a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

ReturnFn_none=@(savings,dsemiz,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z=@(savings,dsemiz,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(savings,dsemiz,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e=@(savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz(savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze=@(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Riskyasset + semiz base options
vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1,1];
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions=struct();
simoptions.riskyasset=1; simoptions.refine_d=vfoptions.refine_d; simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.n_u=vfoptions.n_u; simoptions.u_grid=vfoptions.u_grid; simoptions.pi_u=vfoptions.pi_u;
simoptions.n_semiz=vfoptions.n_semiz; simoptions.semiz_grid=vfoptions.semiz_grid; simoptions.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions.d_grid=d_grid; simoptions.a_grid=a_grid;
n_semiz=vfoptions.n_semiz;

% Options carrying an e shock
vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e; vfoptions_withe.e_grid=vfoptionsbaseline.e_grid; vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptions_withe=simoptions;
simoptions_withe.n_e=simoptionsbaseline.n_e; simoptions_withe.e_grid=simoptionsbaseline.e_grid; simoptions_withe.pi_e=simoptionsbaseline.pi_e;

%% (A) A single point z (value 1, prob 1) is the same as no z
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray'); jequaloneDist_none(1,ceil(n_semiz/2))=1;

[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions);

jequaloneDist_none_z=zeros([n_a,n_semiz,1],'gpuArray'); jequaloneDist_none_z(1,ceil(n_semiz/2),1)=1;
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none_z,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions);

fprintf('CrossTest zase+semiz (A) single-point z == no z, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V0(:)-V0z(:))),max(abs(Policy0(:)-Policy0z(:))),max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% (B) A markov-z that is really an iid == the same shock done as e
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray'); jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

jequaloneDist_e=zeros([n_a,n_semiz,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_e(1,ceil(n_semiz/2),ceil(vfoptions_withe.n_e/2))=1;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_e,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_withe);

fprintf('CrossTest zase+semiz (B) iid-markov-z == e, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))),max(abs(StationaryDist1(:)-StationaryDist2(:))))

%% (C) z&e code with the 'other' a single point reproduces the z-only model
% z just a single point
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_withe);
jequaloneDist3=zeros([n_a,n_semiz,1,vfoptions_withe.n_e],'gpuArray'); jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptions_withe.n_e/2))=1;
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptions_withe);
V3=squeeze(V3); Policy3=squeeze(Policy3); StationaryDist3=squeeze(StationaryDist3);

fprintf('CrossTest zase+semiz (C) z&e with z=1, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V1(:)-V3(:))),max(abs(Policy1(:)-Policy3(:))),max(abs(StationaryDist1(:)-StationaryDist3(:))))

% e just a single point
vfoptions_ze2=vfoptions; vfoptions_ze2.n_e=1; vfoptions_ze2.e_grid=1; vfoptions_ze2.pi_e=1;
simoptions_ze2=simoptions; simoptions_ze2.n_e=1; simoptions_ze2.e_grid=1; simoptions_ze2.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray'); jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1;
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze2);
V4=squeeze(V4); Policy4=squeeze(Policy4); StationaryDist4=squeeze(StationaryDist4);

fprintf('CrossTest zase+semiz (C) z&e with e=1, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V1(:)-V4(:))),max(abs(Policy1(:)-Policy4(:))),max(abs(StationaryDist1(:)-StationaryDist4(:))))

%%
output=struct();

end
