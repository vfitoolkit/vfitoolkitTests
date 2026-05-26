function output=CoreFHorzExpAsset_CrossTests2_d1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 2 noa1+semiz d1: semiz as plain Markov should equal ExpAsset noa1+d1 with z.
% n_d=[n_d1, n_d2, n_d3]; d_grid=[d1_grid; d2_grid; d3_grid].
% PENDING TOOLKIT SUPPORT: side B (semiz) errors at VFI call because ExpAsset+SemiExo+noa1 not implemented.

n_d_withoutsemiz=n_d(1:2); % keep d1, d2; drop d3
d_grid_withoutsemiz=d_grid(1:n_d(1)+n_d(2));

Params.uempbenefit=0;
Params.searcheffortcost=0;

n_z=2;
z_grid=[0.6;1.4];
pi_z=[1-Params.probfindjob, Params.probfindjob;...
      Params.problosejob, 1-Params.problosejob];
Params.z1=z_grid(1);
Params.z2=z_grid(2);

SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzExpAssetSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);

vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;
vfoptions.n_semiz=n_z;
vfoptions.semiz_grid=z_grid;
vfoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.n_semiz=n_z;
simoptions.semiz_grid=z_grid;
simoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
n_semiz=vfoptions.n_semiz;

ReturnFn_semiz=@(d1,d2,d3,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_noa1_semiz(d1,d2,d3,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z=@(d1,d2,a,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_noe_noa1_nosemiz(d1,d2,a,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

%% Side A: plain markov z (no semiz)
jequaloneDist=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist(1,ceil(n_semiz/2))=1;

vfoptionsA=struct();
simoptionsA=struct();
vfoptionsA.experienceasset=vfoptions.experienceasset;
simoptionsA.experienceasset=simoptions.experienceasset;
vfoptionsA.aprimeFn=vfoptions.aprimeFn;
simoptionsA.aprimeFn=simoptions.aprimeFn;
simoptionsA.d_grid=d_grid_withoutsemiz;
simoptionsA.a_grid=simoptions.a_grid;
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(n_d_withoutsemiz,n_a,n_z,N_j,d_grid_withoutsemiz,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1A,n_d_withoutsemiz,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Side B: semiz that mimics the markov
vfoptionsB=vfoptions;
simoptionsB=simoptions;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptionsB);

% Drop d3 from Policy1B to align (Policy1A has d1,d2; Policy1B has d1,d2,d3)
Policy1Bshort=Policy1B(1:2,:,:,:);

fprintf('Cross test 2 (noa1+d1+semiz): semiz as z, this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f \n',max(abs(V1A(:)-V1B(:))),max(abs(Policy1A(:)-Policy1Bshort(:))),max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

output=struct();

end
