function output=CoreFHorzExpAssetU_CrossTests3_nod1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross test 3 noa1+semiz (nod1): ExpAssetU+noa1+semiz with degenerate u (n_u=1, prob 1) should equal
% ExpAsset+noa1+semiz. Both sides use the same noa1+semiz ReturnFn, same d/a/semiz grids, and the
% a2prime formula u*(...) reduces to (...) when u==1.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.
% n_d=[n_d2; n_d3]; d_grid=[d2_grid; d3_grid].

ReturnFn=@(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_noa1_semiz(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn_expasset=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

%% Initial dist (semiz only)
jequaloneDist_none=zeros([n_a,simoptionsbaseline.n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(simoptionsbaseline.n_semiz/2))=1;

%% Side A: experienceasset + semiz (noa1)
vfoptionsA=struct(); simoptionsA=struct();
vfoptionsA.experienceasset=1;
simoptionsA.experienceasset=1;
vfoptionsA.aprimeFn=aprimeFn_expasset;
simoptionsA.aprimeFn=aprimeFn_expasset;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
vfoptionsA.n_semiz=vfoptionsbaseline.n_semiz;
vfoptionsA.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsA.n_semiz=simoptionsbaseline.n_semiz;
simoptionsA.semiz_grid=simoptionsbaseline.semiz_grid;
simoptionsA.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

%% Side B: experienceassetu + semiz (noa1), with degenerate u
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
vfoptionsB.n_semiz=vfoptionsbaseline.n_semiz;
vfoptionsB.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptionsB.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsB.n_semiz=simoptionsbaseline.n_semiz;
simoptionsB.semiz_grid=simoptionsbaseline.semiz_grid;
simoptionsB.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy1,n_d,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3 (noa1+nod1+semiz): expassetu+semiz with u==1 reduces to expasset+semiz, this should be zero: V %.3e, Policy %.3e, Dist %.3e \n',max(abs(V0(:)-V1(:))),max(abs(Policy0(:)-Policy1(:))),max(abs(StationaryDist0(:)-StationaryDist1(:))))

output=struct();

end
