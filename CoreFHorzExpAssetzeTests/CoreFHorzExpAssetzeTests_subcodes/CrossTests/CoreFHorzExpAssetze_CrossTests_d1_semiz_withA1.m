function output=CoreFHorzExpAssetze_CrossTests_d1_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassetze+semiz (d1 version): 'fake' experienceassetze+semiz whose
% aprimeFn ignores e (becomes effectively experienceassetz+semiz), vs actual experienceassetz+semiz.
% Both have z+e shocks and semiz. Should give same V, Policy, StationaryDist.
%
% Pending toolkit support: requires BOTH experienceassetze+SemiExo AND experienceassetz+SemiExo
% dispatchers in ValueFnIter_Case1_FHorz (currently neither exists).

ReturnFn=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetze_d1_z_e_semiz(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn_fakee=@(d2,a2,z,e,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;

semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

vfoptionsA=semizopts;
vfoptionsA.experienceassetze=1;
vfoptionsA.aprimeFn=aprimeFn_fakee;
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=semizopts;
simoptionsA.experienceassetze=1;
simoptionsA.aprimeFn=aprimeFn_fakee;
simoptionsA.n_e=vfoptionsA.n_e;
simoptionsA.e_grid=vfoptionsA.e_grid;
simoptionsA.pi_e=vfoptionsA.pi_e;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

jequaloneDist=zeros([n_a,vfoptionsA.n_semiz,n_z,vfoptionsA.n_e],'gpuArray');
jequaloneDist(1,1,ceil(vfoptionsA.n_semiz/2),ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsB=semizopts;
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_z;
vfoptionsB.n_e=vfoptionsbaseline.n_e;
vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
vfoptionsB.pi_e=vfoptionsbaseline.pi_e;
simoptionsB=semizopts;
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_z;
simoptionsB.n_e=vfoptionsB.n_e;
simoptionsB.e_grid=vfoptionsB.e_grid;
simoptionsB.pi_e=vfoptionsB.pi_e;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('CrossTest1+semiz (fake-e-ignored experienceassetze vs experienceassetz; d1), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
