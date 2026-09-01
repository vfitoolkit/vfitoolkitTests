function output=CoreFHorzExpAssetsemiz_CrossTests2_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 2 for experienceassetsemiz: a 'fake' experienceassetsemiz whose aprimeFn
% ignores semiz must match plain experienceasset (with the same, genuinely
% d3-dependent, semiz present alongside). Only the asset-type machinery differs.
% Should give the same V, Policy and StationaryDist.

% n_d passed as [n_d2,n_d3] (both semiz decision d3 and expasset decision d2 present).
n_semiz=vfoptionsbaseline.n_semiz;
semiz_grid=vfoptionsbaseline.semiz_grid;

jequaloneDist=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist(1,1,ceil(n_semiz/2))=1;

ReturnFn=@(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_nod1_noz_noe(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% 'fake' experienceassetsemiz aprimeFn: takes semiz but ignores it
aprimeFn_fakesemiz=@(d2,a2,semiz,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;
% Plain experienceasset aprimeFn
aprimeFn_plain=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

%% Model A: 'fake' experienceassetsemiz that ignores semiz
vfoptionsA=struct();
vfoptionsA.n_semiz=n_semiz;
vfoptionsA.semiz_grid=semiz_grid;
vfoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptionsA.experienceassetsemiz=1;
vfoptionsA.aprimeFn=aprimeFn_fakesemiz;
simoptionsA=struct();
simoptionsA.n_semiz=n_semiz;
simoptionsA.semiz_grid=semiz_grid;
simoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsA.experienceassetsemiz=1;
simoptionsA.aprimeFn=aprimeFn_fakesemiz;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

%% Model B: plain experienceasset (with the same semiz present)
vfoptionsB=struct();
vfoptionsB.n_semiz=n_semiz;
vfoptionsB.semiz_grid=semiz_grid;
vfoptionsB.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB=struct();
simoptionsB.n_semiz=n_semiz;
simoptionsB.semiz_grid=semiz_grid;
simoptionsB.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsB.experienceasset=1;
simoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('CrossTest2 (fake-semiz-ignored experienceassetsemiz vs plain experienceasset), this should be zero: V %.3e, Policy %.3e, Dist %.3e \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
