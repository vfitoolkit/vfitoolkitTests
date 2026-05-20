function output=CoreFHorzExpAssetze_CrossTests_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassetze: 'fake' experienceassetze whose aprimeFn ignores e (becomes effectively experienceassetz),
% vs actual experienceassetz. Both have z+e shocks. Should give same V, Policy, StationaryDist.

ReturnFn=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetze_nod1_z_e(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);

% 'fake' experienceassetze aprimeFn: ignores e (uses z only)
aprimeFn_fakee=@(d2,a2,z,e,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;
% experienceassetz aprimeFn
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;

% Model A: 'fake' experienceassetze that ignores e
vfoptionsA=struct();
vfoptionsA.experienceassetze=1;
vfoptionsA.aprimeFn=aprimeFn_fakee;
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct();
simoptionsA.experienceassetze=1;
simoptionsA.aprimeFn=aprimeFn_fakee;
simoptionsA.n_e=vfoptionsA.n_e;
simoptionsA.e_grid=vfoptionsA.e_grid;
simoptionsA.pi_e=vfoptionsA.pi_e;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;

jequaloneDist=zeros([n_a,n_z,vfoptionsA.n_e],'gpuArray');
jequaloneDist(1,1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: experienceassetz (with same e shock present, but aprimeFn doesn't depend on e)
vfoptionsB=struct();
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_z;
vfoptionsB.n_e=vfoptionsbaseline.n_e;
vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
vfoptionsB.pi_e=vfoptionsbaseline.pi_e;
simoptionsB=struct();
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_z;
simoptionsB.n_e=vfoptionsB.n_e;
simoptionsB.e_grid=vfoptionsB.e_grid;
simoptionsB.pi_e=vfoptionsB.pi_e;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('CrossTest1 (fake-e-ignored experienceassetze vs experienceassetz), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
