function output=CoreFHorzExpAsset_CrossTests3_nod1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 3 noa1+semiz nod1: ExpAsset noa1+semiz with degenerate aprimeFn(d2,a2)=d2 should match
% standard 1-endo Case1 with semiz (where aprime plays the role of d2).
% n_d=[n_d2, n_d3]; d_grid=[d2_grid; d3_grid].
% PENDING TOOLKIT SUPPORT: side B errors at VFI call (ExpAsset+SemiExo+noa1 not implemented).

n_d3_only=n_d(2); % d3 = dsemiz
d3_grid_only=d_grid(n_d(1)+1:end);

% Common semiz setup
semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

ReturnFn_oneendo=@(d3,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_OneEndo_nod1_noz_noe_noa1_semiz(d3,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_expasset=@(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_noa1_semiz(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn=@(d2,a2) d2;

n_d_expasset=[n_a, n_d3_only]; % d2 = aprime over a_grid, then d3 = dsemiz
d_grid_expasset=[a_grid; d3_grid_only];

jequaloneDist=zeros([n_a,semizopts.n_semiz],'gpuArray');
jequaloneDist(1,ceil(semizopts.n_semiz/2))=1;

% Side A: standard Case1 with semiz, n_d = n_d3 (just dsemiz; aprime is the other choice)
vfoptionsA=semizopts;
simoptionsA=semizopts;
simoptionsA.d_grid=d3_grid_only;
simoptionsA.a_grid=a_grid;
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d3_only,n_a,0,N_j,d3_grid_only,a_grid,[],[],ReturnFn_oneendo,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy0,n_d3_only,n_a,0,N_j,[],Params,simoptionsA);

% Side B: ExpAsset noa1 + semiz
vfoptionsB=semizopts;
vfoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn;
simoptionsB=semizopts;
simoptionsB.experienceasset=1;
simoptionsB.aprimeFn=aprimeFn;
simoptionsB.d_grid=d_grid_expasset;
simoptionsB.a_grid=a_grid;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d_expasset,n_a,0,N_j,d_grid_expasset,a_grid,[],[],ReturnFn_expasset,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1,n_d_expasset,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3 (noa1+semiz): expasset noa1+semiz is just a standard 1-endo+semiz state, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V0(:)-V1(:))),max(abs(Policy0(:)-Policy1(:))),max(abs(StationaryDist0(:)-StationaryDist1(:))))

output=struct();

end
