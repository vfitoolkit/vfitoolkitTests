function output=CoreFHorzExpAssetze_CrossTests4_d1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 4 for experienceassetze+semiz (d1 version): experienceassetze+semiz where z is iid-markov
% compared to experienceassete+semiz with a 2-dim e (n_e=[n_z, n_e_orig]) where e1 takes the z values
% and e2 takes the original e values. The two models are mathematically the same problem.
%
% Pending toolkit support: requires BOTH experienceassetze+SemiExo AND experienceassete+SemiExo
% dispatchers in ValueFnIter_Case1_FHorz (currently neither exists).

pi_z_iid_row=rand(1,n_z); % random transition probabilites
pi_z_iid_row=pi_z_iid_row./sum(pi_z_iid_row); % normalize to sum to one
pi_z=repmat(pi_z_iid_row,n_z,1);

n_e_orig=vfoptionsbaseline.n_e;
e_grid_orig=vfoptionsbaseline.e_grid;
pi_e_orig=vfoptionsbaseline.pi_e;

ReturnFn_A=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetze_d1_z_e_semiz(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B=@(d1,d2,d3,a1prime,a1,a2,semiz,e1,e2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetze_d1_z_e_semiz(d1,d2,d3,a1prime,a1,a2,semiz,e1,e2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

aprimeFn_A=@(d2,a2,z,e,phi1,phi2) phi1*(1-d2)*z*e+(1-phi2)*a2;
aprimeFn_B=@(d2,a2,e1,e2,phi1,phi2) phi1*(1-d2)*e1*e2+(1-phi2)*a2;

semizopts=struct();
semizopts.n_semiz=vfoptionsbaseline.n_semiz;
semizopts.semiz_grid=vfoptionsbaseline.semiz_grid;
semizopts.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

vfoptionsA=semizopts;
vfoptionsA.experienceassetze=1;
vfoptionsA.aprimeFn=aprimeFn_A;
vfoptionsA.n_e=n_e_orig;
vfoptionsA.e_grid=e_grid_orig;
vfoptionsA.pi_e=pi_e_orig;
simoptionsA=semizopts;
simoptionsA.experienceassetze=1;
simoptionsA.aprimeFn=aprimeFn_A;
simoptionsA.n_e=vfoptionsA.n_e;
simoptionsA.e_grid=vfoptionsA.e_grid;
simoptionsA.pi_e=vfoptionsA.pi_e;
simoptionsA.d_grid=d_grid;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

jequaloneDist_A=zeros([n_a,vfoptionsA.n_semiz,n_z,n_e_orig],'gpuArray');
jequaloneDist_A(1,1,ceil(vfoptionsA.n_semiz/2),ceil(n_z/2),ceil(n_e_orig/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

n_e_B=[n_z, n_e_orig];
e_grid_B=[z_grid; e_grid_orig];
pi_e_B=kron(pi_e_orig, pi_z_iid_row');

vfoptionsB=semizopts;
vfoptionsB.experienceassete=1;
vfoptionsB.aprimeFn=aprimeFn_B;
vfoptionsB.n_e=n_e_B;
vfoptionsB.e_grid=e_grid_B;
vfoptionsB.pi_e=pi_e_B;
simoptionsB=semizopts;
simoptionsB.experienceassete=1;
simoptionsB.aprimeFn=aprimeFn_B;
simoptionsB.n_e=vfoptionsB.n_e;
simoptionsB.e_grid=vfoptionsB.e_grid;
simoptionsB.pi_e=vfoptionsB.pi_e;
simoptionsB.d_grid=d_grid;
simoptionsB.a_grid=a_grid;

jequaloneDist_B=zeros([n_a,vfoptionsB.n_semiz,n_z,n_e_orig],'gpuArray');
jequaloneDist_B(1,1,ceil(vfoptionsB.n_semiz/2),ceil(n_z/2),ceil(n_e_orig/2))=1;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('CrossTest4+semiz (experienceassetze iid-markov-z + e vs experienceassete with 2-dim e; d1), this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(Policy_A(:)-Policy_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
