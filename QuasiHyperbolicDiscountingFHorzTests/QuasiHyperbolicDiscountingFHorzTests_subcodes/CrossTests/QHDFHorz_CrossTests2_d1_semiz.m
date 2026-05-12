function output=CoreFHorz_CrossTests2_d1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% n_d=n_d_semiz;
% d_grid=d_grid_semiz;

% For models without semiz we still need
n_d1=n_d(1);
d1_grid=d_grid(1:n_d(1));

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
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);
vfoptions.n_semiz=n_z;
vfoptions.semiz_grid=z_grid;
vfoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.n_semiz=n_z;
simoptions.semiz_grid=z_grid;
simoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.d_grid=d_grid;
% For convenience
n_semiz=vfoptionsbaseline.n_semiz;

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

ReturnFn_semiz=@(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_semize=@(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);

ReturnFn_z=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)...
    ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_ze=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)...
    ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

% Setup some FnsToEvaluate
% FnsToEvaluate_z.assets=@(aprime,a,z) a;
% FnsToEvaluate_z.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
% FnsToEvaluate_e.assets=@(aprime,a,e) a;
% FnsToEvaluate_e.earnings=@(aprime,a,e,w,kappa_j) w*kappa_j*e;

%% Solving for model with one markov
jequaloneDist1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2))=1; % no assets

% First, just use z (without semiz)
vfoptions1A.divideandconquer=0;
simoptions1A=struct();
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(n_d1,n_a,n_z,N_j,d1_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions1A);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1A,n_d1,n_a,n_z,N_j,pi_z,Params,simoptions1A);

% Second, use semiz (without z)
vfoptions1B.divideandconquer=0;
vfoptions1B.n_semiz=vfoptions.n_semiz;
vfoptions1B.semiz_grid=vfoptions.semiz_grid;
vfoptions1B.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions1B.n_semiz=simoptions.n_semiz;
simoptions1B.semiz_grid=simoptions.semiz_grid;
simoptions1B.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptions1B.d_grid=simoptions.d_grid;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,[],vfoptions1B);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptions1B);

Policy1Bshort=[Policy1B(1,:,:,:); Policy1B(3,:,:,:)]; % remove the d2 policy  (as it is not relevant, and is not in Policy1A)

% size(Policy1A)
% size(Policy1B)

fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(V1A(:)-V1B(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(Policy1A(:)-Policy1Bshort(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

% squeeze(abs(sum(StationaryDist1A,1)-sum(StationaryDist1B,1))) % Directly check shocks without Policy.

%% Solving for model with one markov and one e
jequaloneDist1=zeros([n_a,n_semiz,vfoptions.n_e],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2),ceil(vfoptions.n_e/2))=1; % no assets

% First, just use z and e (without semiz)
vfoptions2A.n_e=vfoptions.n_e;
vfoptions2A.e_grid=vfoptions.e_grid;
vfoptions2A.pi_e=vfoptions.pi_e;
simoptions2A.n_e=simoptions.n_e;
simoptions2A.e_grid=simoptions.e_grid;
simoptions2A.pi_e=simoptions.pi_e;
[V2A,Policy2A]=ValueFnIter_Case1_FHorz(n_d1,n_a,n_z,N_j,d1_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions2A);
StationaryDist2A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy2A,n_d1,n_a,n_z,N_j,pi_z,Params,simoptions2A);

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

Policy2Bshort=[Policy2B(1,:,:,:,:); Policy2B(3,:,:,:,:)]; % remove the d2 policy  (as it is not relevant, and is not in Policy2A)

% size(Policy2A)
% size(Policy2B)

fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(V2A(:)-V2B(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(Policy2A(:)-Policy2Bshort(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(StationaryDist2A(:)-StationaryDist2B(:))))

% squeeze(abs(sum(StationaryDist2A(:,1,:,:),1)-sum(StationaryDist2B(:,1,:,:),1))) % Directly check shocks without Policy.
% squeeze(abs(sum(StationaryDist2A(:,2,:,:),1)-sum(StationaryDist2B(:,2,:,:),1))) % Directly check shocks without Policy.


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end