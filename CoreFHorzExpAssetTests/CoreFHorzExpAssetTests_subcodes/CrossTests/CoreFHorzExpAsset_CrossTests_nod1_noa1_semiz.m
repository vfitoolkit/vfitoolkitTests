function output=CoreFHorzExpAsset_CrossTests_nod1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Markov-as-iid cross-test, noa1 + semiz, nod1. n_d=[n_d2, n_d3]; d_grid=[d2_grid; d3_grid].
% PENDING TOOLKIT SUPPORT: errors at first VFI call because ExpAsset+SemiExo+noa1 is not yet implemented.

% Set up z to copy e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

ReturnFn_none=@(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_noa1_semiz(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z   =@(d2,d3,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_noa1_semiz(d2,d3,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e   =@(d2,d3,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_noa1_semiz(d2,d3,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze  =@(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_noa1_semiz(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Common semiz + experienceasset setup
vfoptions=struct();
simoptions=struct();
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
n_semiz=vfoptionsbaseline.n_semiz;

vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e;
vfoptions_withe.e_grid=vfoptionsbaseline.e_grid;
vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptions_withe=simoptions;
simoptions_withe.n_e=simoptionsbaseline.n_e;
simoptions_withe.e_grid=simoptionsbaseline.e_grid;
simoptions_withe.pi_e=simoptionsbaseline.pi_e;

%% (1) z degenerate == no shocks (both with semiz)
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(n_semiz/2))=1;

[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions);

[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions);

fprintf('Cross test (noa1+semiz): z degenerate vs no shocks, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V0(:)-V0z(:))),max(abs(Policy0(:)-Policy0z(:))),max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% (2) Markov-as-iid == iid e (both with semiz)
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1;

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_withe);

fprintf('Cross test (noa1+semiz): iid-markov-z == iid-e, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))),max(abs(StationaryDist1(:)-StationaryDist2(:))))

output=struct();

end
