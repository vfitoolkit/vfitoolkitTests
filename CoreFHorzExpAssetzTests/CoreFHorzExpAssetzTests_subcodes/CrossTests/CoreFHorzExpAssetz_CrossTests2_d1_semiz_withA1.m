function output=CoreFHorzExpAssetz_CrossTests2_d1_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 2 for experienceassetz+semiz (d1 version): 'fake' experienceassetz whose aprimeFn
% ignores z, vs plain experienceasset (both with semiz, both with z present).
% Should give same V, Policy, StationaryDist.

ReturnFn=@(d1,d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetz_d1_z_noe_semiz(d1,d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn_fakez=@(d2,a2,z,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;
aprimeFn_plain=@(d2,a2,phi1,phi2) phi1*(1-d2)+(1-phi2)*a2;

% Common semiz setup
semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

% Model A: 'fake' experienceassetz that ignores z
vfoptionsA=semizopts;
vfoptionsA.experienceassetz=1;
vfoptionsA.aprimeFn=aprimeFn_fakez;
simoptionsA=semizopts;
simoptionsA.experienceassetz=1;
simoptionsA.aprimeFn=aprimeFn_fakez;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

jequaloneDist=zeros([n_a,vfoptionsA.n_semiz,n_z],'gpuArray');
jequaloneDist(1,1,ceil(vfoptionsA.n_semiz/2),ceil(n_z/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: plain experienceasset (with semiz)
vfoptionsB=semizopts;
vfoptionsB.experienceasset=1;
vfoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB=semizopts;
simoptionsB.experienceasset=1;
simoptionsB.aprimeFn=aprimeFn_plain;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('CrossTest2+semiz (fake-z-ignored experienceassetz vs plain experienceasset; d1), this should be zero: V %.3e, Policy %.3e, Dist %.3e \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
