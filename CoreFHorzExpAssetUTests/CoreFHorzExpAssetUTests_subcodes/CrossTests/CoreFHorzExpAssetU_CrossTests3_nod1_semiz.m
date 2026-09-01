function output=CoreFHorzExpAssetU_CrossTests3_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Check that solving with experienceassetu and a degenerate u shock (n_u=1, prob 1)
% gives the same answer as solving with experienceasset (with semiz on both sides).
% Both sides use the same ReturnFn, same d/a grids, same semiz setup, and an
% a2prime formula that agrees pointwise; side A uses the canonical
% experienceasset interface, side B uses the experienceassetu interface with
% u==1 so that u*(...) equals (...).

ReturnFn=@(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% aprimeFn for experienceasset side: no u in signature.
% Matches the experienceassetu baseline u*(phi1*(1-d2)+(1-phi2)*a2) when u==1.
aprimeFn_expasset=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

%% Initial dist (semiz, no z/e)
jequaloneDist_none=zeros([n_a,simoptionsbaseline.n_semiz],'gpuArray');
jequaloneDist_none(1,1,ceil(simoptionsbaseline.n_semiz/2))=1; % no assets

%% Side A: experienceasset + semiz
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

%% Side B: experienceassetu + semiz, with degenerate u (n_u=1, prob 1)
vfoptionsB.experienceassetu=1;
simoptionsB.experienceassetu=1;
vfoptionsB.aprimeFn=vfoptionsbaseline.aprimeFn; % u*(phi1*(1-d2)+(1-phi2)*a2)
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

fprintf('Cross test 3: expassetu with u==1 reduces to expasset, this should be zero: %.3e \n',max(abs(V0(:)-V1(:))))
fprintf('Cross test 3: expassetu with u==1 reduces to expasset, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy1(:))))
fprintf('Cross test 3: expassetu with u==1 reduces to expasset, this should be zero: %.3e \n',max(abs(StationaryDist0(:)-StationaryDist1(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
