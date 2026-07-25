function output=CoreFHorzExpAssetsemiz_CrossTests3_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 3 for experienceassetsemiz: the definitive test of the bothz=[semiz,z]
% index ordering (i.e. the repmat-vs-repelem expansion of the a2primeIndex).
% With BOTH semiz and an ordinary z present, and a degenerate semiz decision (n_d3=1),
% experienceassetsemiz must match experienceassetz with a COMBINED z=[semiz,z] where
% the aprimeFn uses only the first (semiz) component.
% semiz is the fast index of bothz, so the combined pi is kron(pi_z, pi_semiz) and
% the combined grid is [semiz_grid; z_grid]. Should give the same V and StationaryDist.
% A repmat<->repelem swap in the toolkit changes which semiz value pairs with which
% a2prime, so it would break V here (and only here: CrossTests 1 and 2 pass either way).

% n_d passed as [n_d2,n_d3]. Use a degenerate d3 (single value 0).
n_d2=n_d(1);
d2_grid=d_grid(1:n_d2);
n_semiz=vfoptionsbaseline.n_semiz;
semiz_grid=vfoptionsbaseline.semiz_grid;

% semiz transition (dsemiz=0 branch of the SemiExoStateFn), as a plain Markov pi
pi_semiz=[1-Params.probfindjob, Params.probfindjob; Params.problosejob, 1-Params.problosejob];

% Initial dist over (a, semiz, z) -- same layout for both models (semiz fast)
jequaloneDist=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist(1,1,ceil(n_semiz/2),ceil(n_z/2))=1;

% aprimeFns
aprimeFn_semiz=@(d2,a2,semiz,phi1,phi2) phi1*(1-d2)*semiz+(1-phi2)*a2;
aprimeFn_combined=@(d2,a2,z1,z2,phi1,phi2) phi1*(1-d2)*z1+(1-phi2)*a2; % z1 plays semiz, z2 ignored by asset

%% Model A: experienceassetsemiz with semiz (degenerate decision) AND ordinary z
n_d_A=[n_d2,1];
d_grid_A=[d2_grid; 0];
ReturnFn_A=@(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_nod1_z_noe(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

vfoptionsA=struct();
vfoptionsA.n_semiz=n_semiz;
vfoptionsA.semiz_grid=semiz_grid;
vfoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptionsA.experienceassetsemiz=1;
vfoptionsA.aprimeFn=aprimeFn_semiz;
simoptionsA=struct();
simoptionsA.n_semiz=n_semiz;
simoptionsA.semiz_grid=semiz_grid;
simoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptionsA.experienceassetsemiz=1;
simoptionsA.aprimeFn=aprimeFn_semiz;
simoptionsA.d_grid=d_grid_A;
simoptionsA.a_grid=a_grid;
simoptionsA.z_grid=z_grid;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d_A,n_a,n_z,N_j,d_grid_A,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d_A,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model B: experienceassetz with combined z=[semiz,z] (semiz is the fast index)
n_d_B=n_d2;
d_grid_B=d2_grid;
n_z_B=[n_semiz,n_z];
z_grid_B=[semiz_grid; z_grid];
pi_z_B=kron(pi_z,pi_semiz); % semiz fast -> pi_semiz is the second kron argument
ReturnFn_B=@(d2,a1prime,a1,a2,z1,z2,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit) ReturnFn_ExpAssetsemiz_CrossTest3_2z_nod1(d2,a1prime,a1,a2,z1,z2,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit);

vfoptionsB=struct();
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_combined;
simoptionsB=struct();
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_combined;
simoptionsB.d_grid=d_grid_B;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid_B;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d_B,n_a,n_z_B,N_j,d_grid_B,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d_B,n_a,n_z_B,N_j,pi_z_B,Params,simoptionsB);

fprintf('CrossTest3 (experienceassetsemiz+z vs experienceassetz combined-z; pins bothz ordering), this should be zero: V %2.8f, Dist %2.8f \n', max(abs(V_A(:)-V_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
