function output=CoreFHorzExpAssete_CrossTests2_d1_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 2 for experienceassete (d1 version): 'fake' experienceassete whose aprimeFn ignores e,
% noa1 version: the experience asset a2 is the only endogenous state.
% vs plain experienceasset. Both have e present in the model.

ReturnFn=@(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_noz_e_noa1(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

aprimeFn_fakee=@(d2,a2,e,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;
aprimeFn_plain=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

% Model A: 'fake' experienceassete that ignores e
vfoptionsA=struct();
vfoptionsA.experienceassete=1;
vfoptionsA.aprimeFn=aprimeFn_fakee;
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct();
simoptionsA.experienceassete=1;
simoptionsA.aprimeFn=aprimeFn_fakee;
simoptionsA.n_e=vfoptionsA.n_e;
simoptionsA.e_grid=vfoptionsA.e_grid;
simoptionsA.pi_e=vfoptionsA.pi_e;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;

jequaloneDist=zeros([n_a,vfoptionsA.n_e],'gpuArray');
jequaloneDist(1,ceil(vfoptionsA.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

% Model B: plain experienceasset
vfoptionsB=struct();
vfoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn_plain;
vfoptionsB.n_e=vfoptionsbaseline.n_e;
vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
vfoptionsB.pi_e=vfoptionsbaseline.pi_e;
simoptionsB=struct();
simoptionsB.experienceasset=1;
simoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB.n_e=vfoptionsB.n_e;
simoptionsB.e_grid=vfoptionsB.e_grid;
simoptionsB.pi_e=vfoptionsB.pi_e;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('CrossTest2 (noa1: fake-e-ignored experienceassete vs plain experienceasset; d1), this should be zero: V %.3e, Policy %.3e, Dist %.3e \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
