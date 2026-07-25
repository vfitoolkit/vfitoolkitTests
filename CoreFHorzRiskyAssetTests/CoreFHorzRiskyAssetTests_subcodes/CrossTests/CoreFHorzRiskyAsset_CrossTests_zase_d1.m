function output=CoreFHorzRiskyAsset_CrossTests_zase_d1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (z-as-e): RiskyAsset, d1, nosemiz, single asset.
% As CoreFHorzRiskyAsset_CrossTests_zase_nod1 but with a d1 (labour supply h). n_d=[h,riskyshare,savings].
% Three internal equivalences (all riskyasset, all should be machine-precision zero):
%   (A) a single-point z (value 1, prob 1) gives the same as no z at all
%   (B) a markov-z that is really an iid gives the same as the same shock done as e (iid)
%   (C) solving the z&e code with the 'other' shock a single point reproduces the z-only model

% For crosstests, set z to be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in the same place in the ReturnFn

ReturnFn_none=@(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_z=@(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_e=@(h,savings,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_noz_e_nosemiz(h,savings,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_ze=@(h,savings,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

% Riskyasset base options
vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.refine_d=[1,1,1];
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions=struct();
simoptions.riskyasset=1; simoptions.refine_d=vfoptions.refine_d; simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.n_u=vfoptions.n_u; simoptions.u_grid=vfoptions.u_grid; simoptions.pi_u=vfoptions.pi_u;
simoptions.d_grid=d_grid; simoptions.a_grid=a_grid;

% Options carrying an e shock
vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e; vfoptions_withe.e_grid=vfoptionsbaseline.e_grid; vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptions_withe=simoptions;
simoptions_withe.n_e=simoptionsbaseline.n_e; simoptions_withe.e_grid=simoptionsbaseline.e_grid; simoptions_withe.pi_e=simoptionsbaseline.pi_e;

%% (A) A single point z (value 1, prob 1) is the same as no shocks
jequaloneDist_none=zeros([n_a,1],'gpuArray'); jequaloneDist_none(1,1,1)=1;

[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions);

[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions);

fprintf('CrossTest zase (A) single-point z == no z, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V0(:)-V0z(:))),max(abs(Policy0(:)-Policy0z(:))),max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% (B) A markov-z that is really an iid == the same shock done as e
jequaloneDist_z=zeros([n_a,n_z],'gpuArray'); jequaloneDist_z(1,ceil(n_z/2))=1;

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_withe);

fprintf('CrossTest zase (B) iid-markov-z == e, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))),max(abs(StationaryDist1(:)-StationaryDist2(:))))

%% (C) z&e code with the 'other' a single point reproduces the z-only model
% First, make z just a single point (so only e is active, but e is a copy of z)
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_withe);
jequaloneDist3=zeros([n_a,1,vfoptions_withe.n_e],'gpuArray'); jequaloneDist3(1,1,ceil(vfoptions_withe.n_e/2))=1;
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptions_withe);
V3=squeeze(V3); Policy3=squeeze(Policy3); StationaryDist3=squeeze(StationaryDist3);

fprintf('CrossTest zase (C) z&e with z=1, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V1(:)-V3(:))),max(abs(Policy1(:)-Policy3(:))),max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just a single point (so only z is active)
vfoptions_ze2=vfoptions; vfoptions_ze2.n_e=1; vfoptions_ze2.e_grid=1; vfoptions_ze2.pi_e=1;
simoptions_ze2=simoptions; simoptions_ze2.n_e=1; simoptions_ze2.e_grid=1; simoptions_ze2.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
jequaloneDist4=zeros([n_a,n_z,1],'gpuArray'); jequaloneDist4(1,ceil(n_z/2),1)=1;
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze2);
V4=squeeze(V4); Policy4=squeeze(Policy4); StationaryDist4=squeeze(StationaryDist4);

fprintf('CrossTest zase (C) z&e with e=1, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V1(:)-V4(:))),max(abs(Policy1(:)-Policy4(:))),max(abs(StationaryDist1(:)-StationaryDist4(:))))

%%
output=struct();

end
