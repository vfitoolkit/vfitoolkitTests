function output=CoreFHorzExpAsset_CrossTests2_d1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% For models without semiz we still need
n_d_withoutsemiz=n_d(1:2); % d3 is semi-exo
d_grid_withoutsemiz=d_grid(1:sum(n_d(1:2)));

Params.uempbenefit=0; % Need this to make return fns the same
Params.searcheffortcost=0; % No actual need to set this to zero, as d2 does
% nothing in these CrossTests2 examples anyway, and d2=0 is one of the
% choices you can make in which case the return fns are the same.

% WARNING: THE z and semiz IN HERE ARE PRETTY HARDCODED!
n_z=2;
z_grid=[0.6;1.4];
pi_z=[1-Params.probfindjob, Params.probfindjob;...
    Params.problosejob, 1-Params.problosejob];
Params.z1=z_grid(1);
Params.z2=z_grid(2);

% First, do a test in which the semiz is just a duplicate of z
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzExpAssetSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);
vfoptions.n_semiz=n_z;
vfoptions.semiz_grid=z_grid;
vfoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.n_semiz=n_z;
simoptions.semiz_grid=z_grid;
simoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.d_grid=d_grid;
% For convenience
n_semiz=vfoptionsbaseline.n_semiz;

vfoptions.divideandconquer=1;
% Experience asset
vfoptions.experienceasset=1;
simoptions.experienceasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

ReturnFn_semiz=@(d1,d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_noe_semiz(d1,d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_semize=@(d1,d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_e_semiz(d1,d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

ReturnFn_z=@(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)...
    ReturnFn_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_ze=@(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)...
    ReturnFn_d1_z_e_nosemiz(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);


%% Solving for model with one markov
jequaloneDist1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist1(1,1,ceil(n_semiz/2))=1; % no assets

vfoptions1A.divideandconquer=vfoptions.divideandconquer;
vfoptions1A.experienceasset=vfoptions.experienceasset;
simoptions1A.experienceasset=simoptions.experienceasset;
vfoptions1A.aprimeFn=vfoptions.aprimeFn;
simoptions1A.aprimeFn=simoptions.aprimeFn;
simoptions1A.d_grid=simoptions.d_grid;
simoptions1A.a_grid=simoptions.a_grid;
% First, just use z (without semiz)
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(n_d_withoutsemiz,n_a,n_z,N_j,d_grid_withoutsemiz,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions1A);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1A,n_d_withoutsemiz,n_a,n_z,N_j,pi_z,Params,simoptions1A);

vfoptions1B.divideandconquer=vfoptions.divideandconquer;
vfoptions1B.experienceasset=vfoptions.experienceasset;
simoptions1B.experienceasset=simoptions.experienceasset;
vfoptions1B.aprimeFn=vfoptions.aprimeFn;
simoptions1B.aprimeFn=simoptions.aprimeFn;
simoptions1B.d_grid=simoptions.d_grid;
simoptions1B.a_grid=simoptions.a_grid;
% Second, use semiz (without z)
vfoptions1B.n_semiz=vfoptions.n_semiz;
vfoptions1B.semiz_grid=vfoptions.semiz_grid;
vfoptions1B.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions1B.n_semiz=simoptions.n_semiz;
simoptions1B.semiz_grid=simoptions.semiz_grid;
simoptions1B.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptions1B.d_grid=simoptions.d_grid;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,[],vfoptions1B);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptions1B);

Policy1Bshort=[Policy1B(1:2,:,:,:); Policy1B(4,:,:,:)]; % remove the d2 policy  (as it is not relevant, and is not in Policy1A)

% size(Policy1A)
% size(Policy1B)

fprintf('Cross test 2: semiz as z, this should be zero: %.3e \n',max(abs(V1A(:)-V1B(:))))
fprintf('Cross test 2: semiz as z, this should be zero: %.3e \n',max(abs(Policy1A(:)-Policy1Bshort(:))))
fprintf('Cross test 2: semiz as z, this should be zero: %.3e \n',max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

% squeeze(abs(sum(StationaryDist1A,1)-sum(StationaryDist1B,1))) % Directly check shocks without Policy.

%% Solving for model with one markov and one e
jequaloneDist1=zeros([n_a,n_semiz,vfoptions.n_e],'gpuArray');
jequaloneDist1(1,1,ceil(n_semiz/2),ceil(vfoptions.n_e/2))=1; % no assets

vfoptions2A.divideandconquer=vfoptions.divideandconquer;
vfoptions2A.experienceasset=vfoptions.experienceasset;
simoptions2A.experienceasset=simoptions.experienceasset;
vfoptions2A.aprimeFn=vfoptions.aprimeFn;
simoptions2A.aprimeFn=simoptions.aprimeFn;
simoptions2A.d_grid=simoptions.d_grid;
simoptions2A.a_grid=simoptions.a_grid;
% First, just use z and e (without semiz)
vfoptions2A.n_e=vfoptions.n_e;
vfoptions2A.e_grid=vfoptions.e_grid;
vfoptions2A.pi_e=vfoptions.pi_e;
simoptions2A.n_e=simoptions.n_e;
simoptions2A.e_grid=simoptions.e_grid;
simoptions2A.pi_e=simoptions.pi_e;
[V2A,Policy2A]=ValueFnIter_Case1_FHorz(n_d_withoutsemiz,n_a,n_z,N_j,d_grid_withoutsemiz,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions2A);
StationaryDist2A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy2A,n_d_withoutsemiz,n_a,n_z,N_j,pi_z,Params,simoptions2A);

vfoptions2B.divideandconquer=vfoptions.divideandconquer;
vfoptions2B.experienceasset=vfoptions.experienceasset;
simoptions2B.experienceasset=simoptions.experienceasset;
vfoptions2B.aprimeFn=vfoptions.aprimeFn;
simoptions2B.aprimeFn=simoptions.aprimeFn;
simoptions2B.d_grid=simoptions.d_grid;
simoptions2B.a_grid=simoptions.a_grid;
% Second, use semiz and e (without z)
vfoptions2B.n_semiz=vfoptions.n_semiz;
vfoptions2B.semiz_grid=vfoptions.semiz_grid;
vfoptions2B.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions2B.n_semiz=simoptions.n_semiz;
simoptions2B.semiz_grid=simoptions.semiz_grid;
simoptions2B.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptions2B.d_grid=simoptions.d_grid;
vfoptions2B.n_e=vfoptions.n_e;
vfoptions2B.e_grid=vfoptions.e_grid;
vfoptions2B.pi_e=vfoptions.pi_e;
simoptions2B.n_e=simoptions.n_e;
simoptions2B.e_grid=simoptions.e_grid;
simoptions2B.pi_e=simoptions.pi_e;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semize,Params,DiscountFactorParamNames,[],vfoptions2B);
StationaryDist2B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy2B,n_d,n_a,0,N_j,[],Params,simoptions2B);

Policy2Bshort=[Policy2B(1:2,:,:,:,:); Policy2B(4,:,:,:,:)]; % remove the d2 policy  (as it is not relevant, and is not in Policy2A)

% size(Policy2A)
% size(Policy2B)

fprintf('Cross test 2: semiz as z (with e), this should be zero: %.3e \n',max(abs(V2A(:)-V2B(:))))
fprintf('Cross test 2: semiz as z (with e), this should be zero: %.3e \n',max(abs(Policy2A(:)-Policy2Bshort(:))))
fprintf('Cross test 2: semiz as z (with e), this should be zero: %.3e \n',max(abs(StationaryDist2A(:)-StationaryDist2B(:))))

% squeeze(abs(sum(StationaryDist2A(:,1,:,:),1)-sum(StationaryDist2B(:,1,:,:),1))) % Directly check shocks without Policy.
% squeeze(abs(sum(StationaryDist2A(:,2,:,:),1)-sum(StationaryDist2B(:,2,:,:),1))) % Directly check shocks without Policy.


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end