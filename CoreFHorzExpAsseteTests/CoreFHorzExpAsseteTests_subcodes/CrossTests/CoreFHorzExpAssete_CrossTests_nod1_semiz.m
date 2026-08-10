function output=CoreFHorzExpAssete_CrossTests_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassete+semiz: compare experienceassete+semiz (with normal iid e) to
% experienceassetz+semiz (with z that is iid-markov; transition matrix has all rows identical).
% Should give same V, Policy, StationaryDist.

% Override z to be iid-markov matching e setup
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

% ReturnFns: same formula on both sides; the trailing shock slot is e on side A, z on side B
ReturnFn_eside=@(d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_nod1_noz_e_semiz(d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_zside=@(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_nod1_noz_e_semiz(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn_e=@(d2,a2,e,phi1,phi2) phi1*(1-d2)*e+(1-phi2)*a2;
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;

% Common semiz setup
semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

% Model A: experienceassete+semiz with iid e
vfoptionsA=semizopts;
vfoptionsA.experienceassete=1;
vfoptionsA.aprimeFn=aprimeFn_e;
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=semizopts;
simoptionsA.experienceassete=1;
simoptionsA.aprimeFn=aprimeFn_e;
simoptionsA.n_e=vfoptionsA.n_e;
simoptionsA.e_grid=vfoptionsA.e_grid;
simoptionsA.pi_e=vfoptionsA.pi_e;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;

% Side A has no z -- jequaloneDist layout is (n_a, n_semiz, n_e)
jequaloneDist_A=zeros([n_a,vfoptionsA.n_semiz,vfoptionsA.n_e],'gpuArray');
jequaloneDist_A(1,1,ceil(vfoptionsA.n_semiz/2),ceil(vfoptionsA.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_eside,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

% Model B: experienceassetz+semiz with iid-markov z
vfoptionsB=semizopts;
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_z;
simoptionsB=semizopts;
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_z;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid;

% Side B has z (not e) -- jequaloneDist layout is (n_a, n_semiz, n_z)
jequaloneDist_B=zeros([n_a,vfoptionsB.n_semiz,n_z],'gpuArray');
jequaloneDist_B(1,1,ceil(vfoptionsB.n_semiz/2),ceil(n_z/2))=1;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_zside,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('CrossTest1+semiz (experienceassete iid-e vs experienceassetz iid-markov-z; nod1), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
