function output=CoreFHorzExpAssetsemiz_CrossTests_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% Cross-test 1 for experienceassetsemiz: with a degenerate semiz decision (n_d3=1,
% so the semi-exo transition does not depend on d3), semiz is just an ordinary
% Markov state. So experienceassetsemiz (aprime depends on semiz) must match
% experienceassetz (aprime depends on z), with z playing the role of semiz.
% Should give the same V and StationaryDist.
% (Policy is not compared directly: the experienceassetsemiz side carries an extra,
%  payoff-irrelevant d3 index that the experienceassetz side does not have.)

% n_d passed as [n_d2,n_d3]. Use a degenerate d3 (single value 0).
n_d2=n_d(1);
d2_grid=d_grid(1:n_d2);
n_semiz=vfoptionsbaseline.n_semiz;
semiz_grid=vfoptionsbaseline.semiz_grid;

% semiz transition (dsemiz=0 branch of the SemiExoStateFn), as a plain Markov pi
pi_semiz=[1-Params.probfindjob, Params.probfindjob; Params.problosejob, 1-Params.problosejob];

jequaloneDist=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist(1,1,ceil(n_semiz/2))=1;

% aprimeFns (same formula, different shock variable name)
aprimeFn_semiz=@(d2,a2,semiz,phi1,phi2) phi1*(1-d2)*semiz+(1-phi2)*a2;
aprimeFn_z=@(d2,a2,z,phi1,phi2) phi1*(1-d2)*z+(1-phi2)*a2;

%% Model A: experienceassetsemiz with degenerate semiz decision (no ordinary z)
n_d_A=[n_d2,1];
d_grid_A=[d2_grid; 0];
ReturnFn_A=@(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_nod1_noz_noe(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

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

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d_A,n_a,0,N_j,d_grid_A,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_A,n_d_A,n_a,0,N_j,[],Params,simoptionsA);

%% Model B: experienceassetz with z playing the role of semiz
n_d_B=n_d2;
d_grid_B=d2_grid;
n_z_B=n_semiz;
z_grid_B=semiz_grid;
pi_z_B=pi_semiz;
ReturnFn_B=@(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit) ReturnFn_ExpAssetsemiz_CrossTest_zside_nod1(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit);

vfoptionsB=struct();
vfoptionsB.experienceassetz=1;
vfoptionsB.aprimeFn=aprimeFn_z;
simoptionsB=struct();
simoptionsB.experienceassetz=1;
simoptionsB.aprimeFn=aprimeFn_z;
simoptionsB.d_grid=d_grid_B;
simoptionsB.a_grid=a_grid;
simoptionsB.z_grid=z_grid_B;

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d_B,n_a,n_z_B,N_j,d_grid_B,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_B,n_d_B,n_a,n_z_B,N_j,pi_z_B,Params,simoptionsB);

fprintf('CrossTest1 (experienceassetsemiz degenerate-semiz vs experienceassetz), this should be zero: V %.3e, Dist %.3e \n', max(abs(V_A(:)-V_B(:))), max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();

end
