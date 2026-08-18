function output=CoreFHorzExpAssete_CrossTests_nod1_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassete: compare experienceassete (with normal iid e) to
% noa1 version: the experience asset a2 is the only endogenous state.
% experienceassetz (with z that is iid-markov; transition matrix has all rows identical).
% Should give same V, Policy, StationaryDist.

% Override z to be iid-markov matching e setup
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

% ReturnFns (same formula, different shock variable name)
ReturnFn_eside=@(d2,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssete_nod1_noz_e_noa1(d2,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_zside=@(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetz_nod1_z_noe_noa1(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

aprimeFn_e=@(d2,a2,e,phi1,phi2) phi1*(1-d2)*e+(1-phi2)*a2;
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;

% Model A: experienceassete with iid e (the natural setup)
vfoptionsA=struct();
vfoptionsA.experienceassete=1;
vfoptionsA.aprimeFn=aprimeFn_e;
vfoptionsA.n_e=vfoptionsbaseline.n_e;
vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct();
simoptionsA.experienceassete=1;
simoptionsA.aprimeFn=aprimeFn_e;
simoptionsA.n_e=vfoptionsA.n_e;
simoptionsA.e_grid=vfoptionsA.e_grid;
simoptionsA.pi_e=vfoptionsA.pi_e;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;

jequaloneDist=zeros([n_a,vfoptionsA.n_e],'gpuArray');
jequaloneDist(1,ceil(vfoptionsA.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_eside,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

% Model B: experienceassetz with iid-markov z (same shock values, but as Markov with identical rows)
vfoptionsB=struct();
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_z;
simoptionsB=struct();
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_z;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_zside,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('CrossTest1 (noa1: experienceassete iid-e vs experienceassetz iid-markov-z), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
