function output=CoreFHorzExpAsset_CrossTests3_nod1_noa1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 3 noa1: ExpAsset noa1 with degenerate aprimeFn(d2,a2)=d2 should match a standard 1-endo Case1
% problem with the same ReturnFn (aprime plays the role of d2).
% Both have no z, no e, no semiz.

% Side A: standard 1-endo Case1 (no d). aprime is the choice (over a_grid).
ReturnFn_oneendo=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_OneEndo_nod1_noz_noe_noa1_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);

% Side B: ExpAsset noa1. d2 plays the role of aprime (via aprimeFn=@(d2,a2) d2).
ReturnFn_expasset=@(d2,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_noa1_nosemiz(d2,a,r,w,kappa_j,sigma,agej,Jr,pension);

aprimeFn=@(d2,a2) d2; % d2 IS a2prime
d_grid_expasset=a_grid; % set d2 grid = a_grid so d2 takes the same values as next-period a
n_d_expasset=n_a;

jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1;

% Side A: standard Case1 with no d, aprime as choice
vfoptionsA=struct();
simoptionsA=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(0,n_a,0,N_j,[],a_grid,[],[],ReturnFn_oneendo,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,0,n_a,0,N_j,[],Params,simoptionsA);

% Side B: ExpAsset noa1
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

fprintf('Cross test 3 (noa1): expasset noa1 is just a standard 1-endo state, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V0(:)-V1(:))),max(abs(Policy0(:)-Policy1(:))),max(abs(StationaryDist0(:)-StationaryDist1(:))))

output=struct();

end
