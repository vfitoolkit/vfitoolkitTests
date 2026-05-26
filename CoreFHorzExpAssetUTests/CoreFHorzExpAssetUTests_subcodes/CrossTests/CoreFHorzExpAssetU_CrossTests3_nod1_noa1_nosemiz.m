function output=CoreFHorzExpAssetU_CrossTests3_nod1_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross test 3 noa1 (nod1, nosemiz): ExpAssetU+noa1 with degenerate u (n_u=1, prob 1) should equal
% ExpAsset+noa1. Both sides use the same noa1 ReturnFn, same d/a grids, and the
% a2prime formula u*(...) reduces to (...) when u==1.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.

ReturnFn=@(d2,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_noa1_nosemiz(d2,a,r,w,kappa_j,sigma,agej,Jr,pension);

aprimeFn_expasset=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

%% Initial dist (no shocks)
jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1;

%% Side A: experienceasset (noa1)
vfoptionsA=struct(); simoptionsA=struct();
vfoptionsA.experienceasset=1;
simoptionsA.experienceasset=1;
vfoptionsA.aprimeFn=aprimeFn_expasset;
simoptionsA.aprimeFn=aprimeFn_expasset;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

%% Side B: experienceassetu (noa1) with degenerate u
vfoptionsB=struct(); simoptionsB=struct();
vfoptionsB.experienceassetu=1;
simoptionsB.experienceassetu=1;
vfoptionsB.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptionsB.aprimeFn=vfoptionsB.aprimeFn;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
vfoptionsB.n_u=1;
vfoptionsB.u_grid=1;
vfoptionsB.pi_u=1;
simoptionsB.n_u=1;
simoptionsB.u_grid=1;
simoptionsB.pi_u=1;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy1,n_d,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3 (noa1+nod1): expassetu with u==1 reduces to expasset, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V0(:)-V1(:))),max(abs(Policy0(:)-Policy1(:))),max(abs(StationaryDist0(:)-StationaryDist1(:))))

output=struct();

end
