function output=CoreFHorzExpAsset_CrossTests_d1_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Markov-as-iid equivalence cross-test, noa1 nosemiz d1 version.

% Set up z to copy e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

ReturnFn_none=@(d1,d2,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_noa1_nosemiz(d1,d2,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_z   =@(d1,d2,a,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_noe_noa1_nosemiz(d1,d2,a,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_e   =@(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_e_noa1_nosemiz(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_ze  =@(d1,d2,a,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_e_noa1_nosemiz(d1,d2,a,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

vfoptions=struct();
simoptions=struct();
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e;
vfoptions_withe.e_grid=vfoptionsbaseline.e_grid;
vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptions_withe=simoptions;
simoptions_withe.n_e=simoptionsbaseline.n_e;
simoptions_withe.e_grid=simoptionsbaseline.e_grid;
simoptions_withe.pi_e=simoptionsbaseline.pi_e;

%% (1) z degenerate == no shocks
jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1;

[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions);

[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions);

fprintf('Cross test (noa1+d1): z degenerate vs no shocks, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V0(:)-V0z(:))),max(abs(Policy0(:)-Policy0z(:))),max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% (2) Markov-as-iid == iid e
jequaloneDist_z=zeros([n_a,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_z/2))=1;

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_withe);

fprintf('Cross test (noa1+d1): iid-markov-z == iid-e, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))),max(abs(StationaryDist1(:)-StationaryDist2(:))))

%% (3) z+e with z degenerate
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_withe);
jequaloneDist3=zeros([n_a,1,vfoptions_withe.n_e],'gpuArray');
jequaloneDist3(1,1,ceil(vfoptions_withe.n_e/2))=1;
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptions_withe);
V3=squeeze(V3); Policy3=squeeze(Policy3); StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test (noa1+d1): z+e with z degenerate, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V1(:)-V3(:))),max(abs(Policy1(:)-Policy3(:))),max(abs(StationaryDist1(:)-StationaryDist3(:))))

%% (4) z+e with e degenerate
vfoptions_ze=vfoptions; vfoptions_ze.n_e=1; vfoptions_ze.e_grid=1; vfoptions_ze.pi_e=1;
simoptions_ze=simoptions; simoptions_ze.n_e=1; simoptions_ze.e_grid=1; simoptions_ze.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze);
jequaloneDist4=zeros([n_a,n_z,1],'gpuArray');
jequaloneDist4(1,ceil(n_z/2),1)=1;
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze);
V4=squeeze(V4); Policy4=squeeze(Policy4); StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test (noa1+d1): z+e with e degenerate, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V1(:)-V4(:))),max(abs(Policy1(:)-Policy4(:))),max(abs(StationaryDist1(:)-StationaryDist4(:))))

output=struct();

end
