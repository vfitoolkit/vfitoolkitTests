function output=CoreFHorzExpAssetU_CrossTests2_nod1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 2 noa1+semiz nod1: semiz state as a plain Markov should give same answer as ExpAsset noa1 with z.
% n_d=[n_d2, n_d3]; d_grid=[d2_grid; d3_grid].

% Pieces used by both sides
n_d_withoutsemiz=n_d(1); % drop d3
d_grid_withoutsemiz=d_grid(1:n_d(1));

Params.uempbenefit=0;
Params.searcheffortcost=0;

% Hardcoded 2-state z that matches the semiz (employed/unemployed) shock
n_z=2;
z_grid=[0.6;1.4];
pi_z=[1-Params.probfindjob, Params.probfindjob;...
      Params.problosejob, 1-Params.problosejob];
Params.z1=z_grid(1);
Params.z2=z_grid(2);

% Use the JustAMarkov semiz state fn so the semiz transitions match pi_z above
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzExpAssetSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);

vfoptions.experienceassetu=1;
simoptions.experienceassetu=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions.n_u=vfoptions.n_u;
simoptions.u_grid=vfoptions.u_grid;
simoptions.pi_u=vfoptions.pi_u;
vfoptions.n_semiz=n_z;
vfoptions.semiz_grid=z_grid;
vfoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.n_semiz=n_z;
simoptions.semiz_grid=z_grid;
simoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
n_semiz=vfoptions.n_semiz;

ReturnFn_semiz=@(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_noa1_semiz(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z=@(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_noa1_nosemiz(d2,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

%% Side A: plain markov z (no semiz)
jequaloneDist=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist(1,ceil(n_semiz/2))=1;

vfoptionsA=struct();
simoptionsA=struct();
vfoptionsA.experienceassetu=vfoptions.experienceassetu;
simoptionsA.experienceassetu=simoptions.experienceassetu;
vfoptionsA.aprimeFn=vfoptions.aprimeFn;
simoptionsA.aprimeFn=simoptions.aprimeFn;
simoptionsA.d_grid=d_grid_withoutsemiz;
simoptionsA.a_grid=simoptions.a_grid;
vfoptionsA.n_u=vfoptions.n_u;
vfoptionsA.u_grid=vfoptions.u_grid;
vfoptionsA.pi_u=vfoptions.pi_u;
simoptionsA.n_u=simoptions.n_u;
simoptionsA.u_grid=simoptions.u_grid;
simoptionsA.pi_u=simoptions.pi_u;
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(n_d_withoutsemiz,n_a,n_z,N_j,d_grid_withoutsemiz,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1A,n_d_withoutsemiz,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Side B: semiz that mimics the markov (no z)
vfoptionsB=vfoptions;
simoptionsB=simoptions;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptionsB);

% Drop the d2 row from Policy1B to align with Policy1A (which has no d2 since side A's n_d = n_d_withoutsemiz = n_d2 only)
% Actually Policy1A has shape [1=d2, n_a, n_z, N_j]; Policy1B has [2=d2+d3, n_a, n_semiz, N_j]. Drop d3 to compare.
Policy1Bshort=Policy1B(1,:,:,:);

fprintf('Cross test 2 (noa1+semiz): semiz as z, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V1A(:)-V1B(:))),max(abs(Policy1A(:)-Policy1Bshort(:))),max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

output=struct();

end
