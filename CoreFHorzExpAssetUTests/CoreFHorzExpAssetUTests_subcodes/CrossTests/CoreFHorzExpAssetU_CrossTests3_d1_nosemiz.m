function output=CoreFHorzExpAssetU_CrossTests3_d1_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Check that solving with experienceassetu and a degenerate u shock (n_u=1, prob 1)
% gives the same answer as solving with experienceasset.
% Both sides use the same ReturnFn, same d/a grids, and an a2prime formula
% that agrees pointwise; side A uses the canonical experienceasset interface,
% side B uses the experienceassetu interface with u==1 so that u*(...) equals (...).

ReturnFn=@(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

% aprimeFn for experienceasset side: no u in signature.
% Matches the experienceassetu baseline u*(phi1*(1-d2)+(1-phi2)*a2) when u==1.
aprimeFn_expasset=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

%% Initial dist (no shocks)
jequaloneDist_none=zeros([n_a],'gpuArray');
jequaloneDist_none(1,1)=1; % no assets

%% Side A: experienceasset
vfoptionsA.experienceasset=1;
simoptionsA.experienceasset=1;
vfoptionsA.aprimeFn=aprimeFn_expasset;
simoptionsA.aprimeFn=aprimeFn_expasset;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

%% Side B: experienceassetu with degenerate u (n_u=1, prob 1)
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
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy1,n_d,n_a,0,N_j,1,Params,simoptionsB);

fprintf('Cross test 3: expassetu with u==1 reduces to expasset, this should be zero: %.3e \n',max(abs(V0(:)-V1(:))))
fprintf('Cross test 3: expassetu with u==1 reduces to expasset, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy1(:))))
fprintf('Cross test 3: expassetu with u==1 reduces to expasset, this should be zero: %.3e \n',max(abs(StationaryDist0(:)-StationaryDist1(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
