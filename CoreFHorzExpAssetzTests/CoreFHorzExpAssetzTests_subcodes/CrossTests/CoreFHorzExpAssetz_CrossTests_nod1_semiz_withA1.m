function output=CoreFHorzExpAssetz_CrossTests_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassetz+semiz: compare experienceassetz+semiz (with z that is iid-markov;
% transition matrix has all rows identical) to experienceassete+semiz (with normal iid e).
% Should give same V, Policy, StationaryDist.
%
% Pending toolkit support: requires BOTH experienceassetz+SemiExo AND experienceassete+SemiExo
% dispatchers in ValueFnIter_Case1_FHorz (currently neither exists).

% Override z to be iid-markov matching e setup
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

% ReturnFns: same formula on both sides; the trailing shock slot is z on side A, e on side B
ReturnFn_zside=@(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetz_nod1_z_noe_semiz(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_eside=@(d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetz_nod1_z_noe_semiz(d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% aprimeFns: same formula, different shock variable name
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;
aprimeFn_e=@(d2,a2,e,phi1,phi2) phi1*(1-d2)*e+(1-phi2)*a2;

% Common semiz setup
semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

% Model A: experienceassetz+semiz with iid-markov z
vfoptionsA=semizopts;
vfoptionsA.experienceassetz=1;
vfoptionsA.aprimeFn=aprimeFn_z;
simoptionsA=semizopts;
simoptionsA.experienceassetz=1;
simoptionsA.aprimeFn=aprimeFn_z;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

jequaloneDist=zeros([n_a,vfoptionsA.n_semiz,n_z],'gpuArray');
jequaloneDist(1,1,ceil(vfoptionsA.n_semiz/2),ceil(n_z/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_zside,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: experienceassete+semiz with iid e (matching the z setup above)
vfoptionsB=semizopts;
vfoptionsB.experienceassete=1;
vfoptionsB.aprimeFn=aprimeFn_e;
vfoptionsB.n_e=n_z;
vfoptionsB.e_grid=z_grid;
vfoptionsB.pi_e=pi_z(1,:)';
simoptionsB=semizopts;
simoptionsB.experienceassete=1;
simoptionsB.aprimeFn=aprimeFn_e;
simoptionsB.n_e=vfoptionsB.n_e;
simoptionsB.e_grid=vfoptionsB.e_grid;
simoptionsB.pi_e=vfoptionsB.pi_e;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

% Side B has no z — jequaloneDist needs the (n_a, n_semiz, n_e) layout
jequaloneDist_B=zeros([n_a,vfoptionsB.n_semiz,vfoptionsB.n_e],'gpuArray');
jequaloneDist_B(1,1,ceil(vfoptionsB.n_semiz/2),ceil(vfoptionsB.n_e/2))=1;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_eside,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('CrossTest1+semiz (experienceassetz iid-markov-z vs experienceassete iid-e; nod1), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
