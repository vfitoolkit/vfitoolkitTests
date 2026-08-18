function output=CoreFHorzExpAssete_CrossTests_d1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassete+semiz (d1 version): compare experienceassete+semiz (iid e) to
% noa1 version: the experience asset a2 is the only endogenous state.
% experienceassetz+semiz (iid-markov z). Should give same V, Policy, StationaryDist.

% Override z to be iid-markov matching e setup
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

ReturnFn_eside=@(d1,d2,d3,a,semiz,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_d1_noz_e_noa1_semiz(d1,d2,d3,a,semiz,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_zside=@(d1,d2,d3,a,semiz,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_d1_noz_e_noa1_semiz(d1,d2,d3,a,semiz,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn_e=@(d2,a2,e,phi1,phi2) phi1*(1-d2)*e+(1-phi2)*a2;
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;

semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

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

jequaloneDist_A=zeros([n_a,vfoptionsA.n_semiz,vfoptionsA.n_e],'gpuArray');
jequaloneDist_A(1,ceil(vfoptionsA.n_semiz/2),ceil(vfoptionsA.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_eside,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

vfoptionsB=semizopts;
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_z;
simoptionsB=semizopts;
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_z;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid;

jequaloneDist_B=zeros([n_a,vfoptionsB.n_semiz,n_z],'gpuArray');
jequaloneDist_B(1,ceil(vfoptionsB.n_semiz/2),ceil(n_z/2))=1;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_zside,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('CrossTest1+semiz (noa1: experienceassete iid-e vs experienceassetz iid-markov-z; d1), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
