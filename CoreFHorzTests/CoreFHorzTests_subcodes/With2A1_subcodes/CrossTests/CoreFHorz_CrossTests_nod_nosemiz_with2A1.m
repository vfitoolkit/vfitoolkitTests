function output=CoreFHorzTwoEndo_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

n_d=0;
d_grid=[];

% For crosstests, set up z to just be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in same place in earnings

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
ReturnFn_none=@(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_noz_noe_nosemiz(a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_z=@(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_noe_nosemiz(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_e=@(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_noz_e_nosemiz(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_ze=@(a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_e_nosemiz(a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);

% Setup some FnsToEvaluate
% FnsToEvaluate_z.assets=@(aprime,a,z) a;
% FnsToEvaluate_z.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
% FnsToEvaluate_e.assets=@(aprime,a,e) a;
% FnsToEvaluate_e.earnings=@(aprime,a,e,w,kappa_j) w*kappa_j*e;


%% Solving with just a single points for z with value 1 and prob 1 gives us same as no shocks
jequaloneDist_none=zeros([n_a,1],'gpuArray');
jequaloneDist_none(1)=1; % no assets

vfoptions_z=struct();
simoptions_z=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_z);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptions_z);

vfoptions_z=struct();
simoptions_z=struct();
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptions_z);

fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid
% zeros assets, mid points for any shocks
jequaloneDist_z=zeros([n_a,n_z],'gpuArray');
jequaloneDist_z(1,1,ceil(n_z/2))=1; % no assets, midpoint shock

vfoptions_z=struct();
simoptions_z=struct();
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_z);

vfoptions_e=vfoptions;
simoptions_e=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_e);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptions_e);

fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist2(:))))

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
vfoptions_ze1=vfoptions;
simoptions_ze1=simoptions;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze1);
jequaloneDist3=zeros([n_a,1,vfoptions_ze1.n_e],'gpuArray');
jequaloneDist3(1,1,1,ceil(vfoptions_ze1.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptions_ze1);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1
vfoptions_ze2.n_e=1;
vfoptions_ze2.e_grid=1;
vfoptions_ze2.pi_e=1;
simoptions_ze2=vfoptions_ze2;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
jequaloneDist4=zeros([n_a,n_z,1],'gpuArray');
jequaloneDist4(1,1,ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze2);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end