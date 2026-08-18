function output=CoreFHorzExpAssetz_CrossTests_d1_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassetz (d1, noa1 version): compare experienceassetz (with z that is
% iid-markov) to experienceassete (with normal iid e). The experience asset a2 is the only
% endogenous state on both sides. Should give same V, Policy, StationaryDist.
%
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.
% n_d input = [n_d1, n_d2]. d_grid = [d1_grid; d2_grid].

% Override z to be iid-markov matching e setup
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

% ReturnFns
ReturnFn_zside=@(d1,d2,a,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssetz_d1_z_noe_noa1(d1,d2,a,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_eside=@(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_noz_e_noa1(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;
aprimeFn_e=@(d2,a2,e,phi1,phi2) phi1*(1-d2)*e+(1-phi2)*a2;

% Model A: experienceassetz with iid-markov z
vfoptionsA=struct();
vfoptionsA.experienceassetz=1;
vfoptionsA.aprimeFn=aprimeFn_z;
simoptionsA=struct();
simoptionsA.experienceassetz=1;
simoptionsA.aprimeFn=aprimeFn_z;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_zside,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Model B: experienceassete with iid e
vfoptionsB=struct();
vfoptionsB.experienceassete=1;
vfoptionsB.aprimeFn=aprimeFn_e;
vfoptionsB.n_e=n_z;
vfoptionsB.e_grid=z_grid;
vfoptionsB.pi_e=pi_z(1,:)';
simoptionsB=struct();
simoptionsB.experienceassete=1;
simoptionsB.aprimeFn=aprimeFn_e;
simoptionsB.n_e=vfoptionsB.n_e;
simoptionsB.e_grid=vfoptionsB.e_grid;
simoptionsB.pi_e=vfoptionsB.pi_e;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_eside,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('CrossTest1 (noa1: experienceassetz iid-markov-z vs experienceassete iid-e; d1), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
