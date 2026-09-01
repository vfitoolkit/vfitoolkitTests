function output=CoreFHorzExpAsset_CrossTests3_d1_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 3 noa1 (d1 version): ExpAsset noa1+d1 with degenerate aprimeFn(d2,a2)=d2 should match a
% standard 1-endo Case1 with d1 (aprime plays the role of d2).
% n_d input = [n_d1, n_d2] (the original); n_a is scalar. d_grid = [d1_grid; d2_grid].

% Side A: standard 1-endo Case1 with d=d1, aprime as the (additional) choice.
ReturnFn_oneendo=@(d1,aprime,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_OneEndo_d1_noz_noe_noa1_nosemiz(d1,aprime,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% Side B: ExpAsset noa1 d1. d2 (= last d) plays the role of aprime.
ReturnFn_expasset=@(d1,d2,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_noa1_nosemiz(d1,d2,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

aprimeFn=@(d2,a2) d2;
n_d1_only=n_d(1);
d1_grid_only=d_grid(1:n_d1_only);

% For ExpAsset side: keep d1, replace d2 with a_grid
n_d_expasset=[n_d1_only, n_a];
d_grid_expasset=[d1_grid_only; a_grid];

jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1;

% Side A: standard Case1 with d=d1
vfoptionsA=struct();
simoptionsA=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d1_only,n_a,0,N_j,d1_grid_only,a_grid,[],[],ReturnFn_oneendo,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d1_only,n_a,0,N_j,[],Params,simoptionsA);

% Side B: ExpAsset noa1 with d1
vfoptionsB=struct();
simoptionsB=struct();
vfoptionsB.experienceasset=1;
simoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn;
simoptionsB.aprimeFn=aprimeFn;
simoptionsB.d_grid=d_grid_expasset;
simoptionsB.a_grid=a_grid;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d_expasset,n_a,0,N_j,d_grid_expasset,a_grid,[],[],ReturnFn_expasset,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy1,n_d_expasset,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3 (noa1+d1): expasset noa1 is just a standard 1-endo state, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V0(:)-V1(:))),max(abs(Policy0(:)-Policy1(:))),max(abs(StationaryDist0(:)-StationaryDist1(:))))

output=struct();

end
