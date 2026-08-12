function output=EZFHorz_CrossTests2_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Epstein-Zin mirror of QHDFHorz_CrossTests2_nod1_semiz.m, run for the three EZ cases:
% a semi-exogenous state that is really just a markov must give the same answer as the markov.
% For EZ this checks the certainty-equivalent is taken correctly over the semiz transitions.
% TEST-FIRST: EZ with semi-exogenous shocks is NOT yet implemented in the toolkit;
% the semiz-side solves in this file are EXPECTED TO ERROR until the EZ SemiExo solvers exist.

% n_d=n_d2_semiz;
% d_grid=d2_grid_semiz;

Params.uempbenefit=0; % Need this to make return fns the same
Params.searcheffortcost=0; % Makes the (1-searcheffortcost*d2) factor in the EZ semiz return fns
% equal to 1, so the semiz return fns coincide with the nosemiz ones. (Also, d2 does nothing in
% these CrossTests2 examples anyway, and d2=0 is one of the choices you can make.)

% WARNING: THE z and semiz IN HERE ARE PRETTY HARDCODED!
n_z=2;
z_grid=[0.6;1.4];
pi_z=[1-Params.probfindjob, Params.probfindjob;...
    Params.problosejob, 1-Params.problosejob];
Params.z1=z_grid(1);
Params.z2=z_grid(2);

% First, do a test in which the semiz is just a duplicate of z
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);
vfoptions_semiz.n_semiz=n_z;
vfoptions_semiz.semiz_grid=z_grid;
vfoptions_semiz.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions_semiz.n_semiz=n_z;
simoptions_semiz.semiz_grid=z_grid;
simoptions_semiz.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions_semiz.d_grid=d_grid;
% For convenience
n_semiz=vfoptionsbaseline.n_semiz;

ReturnFn_semiz_cons=@(d2,aprime,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_cons_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_semize_cons=@(d2,aprime,a,semiz,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_cons_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension)...
    EZReturnFn_cons_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_ze_cons=@(aprime,a,z,e,r,w,kappa_j,agej,Jr,pension)...
    EZReturnFn_cons_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,agej,Jr,pension);
ReturnFn_semiz_posU=@(d2,aprime,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_positiveUtils_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_semize_posU=@(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_positiveUtils_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_posU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension)...
    EZReturnFn_positiveUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_ze_posU=@(aprime,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension)...
    EZReturnFn_positiveUtils_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_semiz_negU=@(d2,aprime,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_negativeUtils_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_semize_negU=@(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_negativeUtils_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_negU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension)...
    EZReturnFn_negativeUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_ze_negU=@(aprime,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension)...
    EZReturnFn_negativeUtils_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension);


%% Epstein-Zin case 1: Consumption-units (EZutils=0)
%% Solving for model with one markov
jequaloneDist1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2))=1; % no assets

% First, just use z (without semiz)
vfoptions1A=struct();
vfoptions1A.divideandconquer=0;
vfoptions1A.exoticpreferences='EpsteinZin';
vfoptions1A.EZutils=0;
vfoptions1A.EZriskaversion='ezgamma';
vfoptions1A.EZeis='ezphi';
simoptions1A=struct();
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_z_cons,Params,DiscountFactorParamNames,[],vfoptions1A);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1A,0,n_a,n_z,N_j,pi_z,Params,simoptions1A);

% Second, use semiz (without z)
vfoptions1B=vfoptions_semiz;
vfoptions1B.divideandconquer=0;
vfoptions1B.exoticpreferences='EpsteinZin';
vfoptions1B.EZutils=0;
vfoptions1B.EZriskaversion='ezgamma';
vfoptions1B.EZeis='ezphi';
simoptions1B=simoptions_semiz;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz_cons,Params,DiscountFactorParamNames,[],vfoptions1B);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptions1B);

Policy1Bshort=Policy1B(2,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy1A)

fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(V1A(:)-V1B(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(Policy1A(:)-Policy1Bshort(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

%% Solving for model with one markov and one e
jequaloneDist2=zeros([n_a,n_semiz,vfoptionsbaseline.n_e],'gpuArray');
jequaloneDist2(1,ceil(n_semiz/2),ceil(vfoptionsbaseline.n_e/2))=1; % no assets

% First, just use z and e (without semiz)
vfoptions2A=struct();
vfoptions2A.exoticpreferences='EpsteinZin';
vfoptions2A.EZutils=0;
vfoptions2A.EZriskaversion='ezgamma';
vfoptions2A.EZeis='ezphi';
vfoptions2A.n_e=vfoptionsbaseline.n_e;
vfoptions2A.e_grid=vfoptionsbaseline.e_grid;
vfoptions2A.pi_e=vfoptionsbaseline.pi_e;
simoptions2A=struct();
simoptions2A.n_e=simoptionsbaseline.n_e;
simoptions2A.e_grid=simoptionsbaseline.e_grid;
simoptions2A.pi_e=simoptionsbaseline.pi_e;
[V2A,Policy2A]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_ze_cons,Params,DiscountFactorParamNames,[],vfoptions2A);
StationaryDist2A=StationaryDist_FHorz_Case1(jequaloneDist2,AgeWeightParamNames,Policy2A,0,n_a,n_z,N_j,pi_z,Params,simoptions2A);

% Second, use semiz and e (without z)
vfoptions2B=vfoptions_semiz;
vfoptions2B.exoticpreferences='EpsteinZin';
vfoptions2B.EZutils=0;
vfoptions2B.EZriskaversion='ezgamma';
vfoptions2B.EZeis='ezphi';
vfoptions2B.n_e=vfoptionsbaseline.n_e;
vfoptions2B.e_grid=vfoptionsbaseline.e_grid;
vfoptions2B.pi_e=vfoptionsbaseline.pi_e;
simoptions2B=simoptions_semiz;
simoptions2B.n_e=simoptionsbaseline.n_e;
simoptions2B.e_grid=simoptionsbaseline.e_grid;
simoptions2B.pi_e=simoptionsbaseline.pi_e;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semize_cons,Params,DiscountFactorParamNames,[],vfoptions2B);
StationaryDist2B=StationaryDist_FHorz_Case1(jequaloneDist2,AgeWeightParamNames,Policy2B,n_d,n_a,0,N_j,[],Params,simoptions2B);

Policy2Bshort=Policy2B(2,:,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy2A)

fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(V2A(:)-V2B(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(Policy2A(:)-Policy2Bshort(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(StationaryDist2A(:)-StationaryDist2B(:))))

clear V1A V1B V2A V2B Policy1A Policy1B Policy2A Policy2B Policy1Bshort Policy2Bshort StationaryDist1A StationaryDist1B StationaryDist2A StationaryDist2B

%% Epstein-Zin case 2: Utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%% Solving for model with one markov
jequaloneDist1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2))=1; % no assets

% First, just use z (without semiz)
vfoptions1A=struct();
vfoptions1A.divideandconquer=0;
vfoptions1A.exoticpreferences='EpsteinZin';
vfoptions1A.EZutils=1;
vfoptions1A.EZpositiveutility=1;
vfoptions1A.EZriskaversion='ezrisk';
simoptions1A=struct();
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_z_posU,Params,DiscountFactorParamNames,[],vfoptions1A);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1A,0,n_a,n_z,N_j,pi_z,Params,simoptions1A);

% Second, use semiz (without z)
vfoptions1B=vfoptions_semiz;
vfoptions1B.divideandconquer=0;
vfoptions1B.exoticpreferences='EpsteinZin';
vfoptions1B.EZutils=1;
vfoptions1B.EZpositiveutility=1;
vfoptions1B.EZriskaversion='ezrisk';
simoptions1B=simoptions_semiz;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz_posU,Params,DiscountFactorParamNames,[],vfoptions1B);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptions1B);

Policy1Bshort=Policy1B(2,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy1A)

fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(V1A(:)-V1B(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(Policy1A(:)-Policy1Bshort(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

%% Solving for model with one markov and one e
jequaloneDist2=zeros([n_a,n_semiz,vfoptionsbaseline.n_e],'gpuArray');
jequaloneDist2(1,ceil(n_semiz/2),ceil(vfoptionsbaseline.n_e/2))=1; % no assets

% First, just use z and e (without semiz)
vfoptions2A=struct();
vfoptions2A.exoticpreferences='EpsteinZin';
vfoptions2A.EZutils=1;
vfoptions2A.EZpositiveutility=1;
vfoptions2A.EZriskaversion='ezrisk';
vfoptions2A.n_e=vfoptionsbaseline.n_e;
vfoptions2A.e_grid=vfoptionsbaseline.e_grid;
vfoptions2A.pi_e=vfoptionsbaseline.pi_e;
simoptions2A=struct();
simoptions2A.n_e=simoptionsbaseline.n_e;
simoptions2A.e_grid=simoptionsbaseline.e_grid;
simoptions2A.pi_e=simoptionsbaseline.pi_e;
[V2A,Policy2A]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_ze_posU,Params,DiscountFactorParamNames,[],vfoptions2A);
StationaryDist2A=StationaryDist_FHorz_Case1(jequaloneDist2,AgeWeightParamNames,Policy2A,0,n_a,n_z,N_j,pi_z,Params,simoptions2A);

% Second, use semiz and e (without z)
vfoptions2B=vfoptions_semiz;
vfoptions2B.exoticpreferences='EpsteinZin';
vfoptions2B.EZutils=1;
vfoptions2B.EZpositiveutility=1;
vfoptions2B.EZriskaversion='ezrisk';
vfoptions2B.n_e=vfoptionsbaseline.n_e;
vfoptions2B.e_grid=vfoptionsbaseline.e_grid;
vfoptions2B.pi_e=vfoptionsbaseline.pi_e;
simoptions2B=simoptions_semiz;
simoptions2B.n_e=simoptionsbaseline.n_e;
simoptions2B.e_grid=simoptionsbaseline.e_grid;
simoptions2B.pi_e=simoptionsbaseline.pi_e;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semize_posU,Params,DiscountFactorParamNames,[],vfoptions2B);
StationaryDist2B=StationaryDist_FHorz_Case1(jequaloneDist2,AgeWeightParamNames,Policy2B,n_d,n_a,0,N_j,[],Params,simoptions2B);

Policy2Bshort=Policy2B(2,:,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy2A)

fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(V2A(:)-V2B(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(Policy2A(:)-Policy2Bshort(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(StationaryDist2A(:)-StationaryDist2B(:))))

clear V1A V1B V2A V2B Policy1A Policy1B Policy2A Policy2B Policy1Bshort Policy2Bshort StationaryDist1A StationaryDist1B StationaryDist2A StationaryDist2B

%% Epstein-Zin case 3: Utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
%% Solving for model with one markov
jequaloneDist1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2))=1; % no assets

% First, just use z (without semiz)
vfoptions1A=struct();
vfoptions1A.divideandconquer=0;
vfoptions1A.exoticpreferences='EpsteinZin';
vfoptions1A.EZutils=1;
vfoptions1A.EZpositiveutility=0;
vfoptions1A.EZriskaversion='ezrisk';
simoptions1A=struct();
[V1A,Policy1A]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_z_negU,Params,DiscountFactorParamNames,[],vfoptions1A);
StationaryDist1A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1A,0,n_a,n_z,N_j,pi_z,Params,simoptions1A);

% Second, use semiz (without z)
vfoptions1B=vfoptions_semiz;
vfoptions1B.divideandconquer=0;
vfoptions1B.exoticpreferences='EpsteinZin';
vfoptions1B.EZutils=1;
vfoptions1B.EZpositiveutility=0;
vfoptions1B.EZriskaversion='ezrisk';
simoptions1B=simoptions_semiz;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semiz_negU,Params,DiscountFactorParamNames,[],vfoptions1B);
StationaryDist1B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,Policy1B,n_d,n_a,0,N_j,[],Params,simoptions1B);

Policy1Bshort=Policy1B(2,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy1A)

fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(V1A(:)-V1B(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(Policy1A(:)-Policy1Bshort(:))))
fprintf('Cross test: semiz as z, this should be zero: %2.8f \n',max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

%% Solving for model with one markov and one e
jequaloneDist2=zeros([n_a,n_semiz,vfoptionsbaseline.n_e],'gpuArray');
jequaloneDist2(1,ceil(n_semiz/2),ceil(vfoptionsbaseline.n_e/2))=1; % no assets

% First, just use z and e (without semiz)
vfoptions2A=struct();
vfoptions2A.exoticpreferences='EpsteinZin';
vfoptions2A.EZutils=1;
vfoptions2A.EZpositiveutility=0;
vfoptions2A.EZriskaversion='ezrisk';
vfoptions2A.n_e=vfoptionsbaseline.n_e;
vfoptions2A.e_grid=vfoptionsbaseline.e_grid;
vfoptions2A.pi_e=vfoptionsbaseline.pi_e;
simoptions2A=struct();
simoptions2A.n_e=simoptionsbaseline.n_e;
simoptions2A.e_grid=simoptionsbaseline.e_grid;
simoptions2A.pi_e=simoptionsbaseline.pi_e;
[V2A,Policy2A]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_ze_negU,Params,DiscountFactorParamNames,[],vfoptions2A);
StationaryDist2A=StationaryDist_FHorz_Case1(jequaloneDist2,AgeWeightParamNames,Policy2A,0,n_a,n_z,N_j,pi_z,Params,simoptions2A);

% Second, use semiz and e (without z)
vfoptions2B=vfoptions_semiz;
vfoptions2B.exoticpreferences='EpsteinZin';
vfoptions2B.EZutils=1;
vfoptions2B.EZpositiveutility=0;
vfoptions2B.EZriskaversion='ezrisk';
vfoptions2B.n_e=vfoptionsbaseline.n_e;
vfoptions2B.e_grid=vfoptionsbaseline.e_grid;
vfoptions2B.pi_e=vfoptionsbaseline.pi_e;
simoptions2B=simoptions_semiz;
simoptions2B.n_e=simoptionsbaseline.n_e;
simoptions2B.e_grid=simoptionsbaseline.e_grid;
simoptions2B.pi_e=simoptionsbaseline.pi_e;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_semize_negU,Params,DiscountFactorParamNames,[],vfoptions2B);
StationaryDist2B=StationaryDist_FHorz_Case1(jequaloneDist2,AgeWeightParamNames,Policy2B,n_d,n_a,0,N_j,[],Params,simoptions2B);

Policy2Bshort=Policy2B(2,:,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy2A)

fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(V2A(:)-V2B(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(Policy2A(:)-Policy2Bshort(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %2.8f \n',max(abs(StationaryDist2A(:)-StationaryDist2B(:))))

clear V1A V1B V2A V2B Policy1A Policy1B Policy2A Policy2B Policy1Bshort Policy2Bshort StationaryDist1A StationaryDist1B StationaryDist2A StationaryDist2B

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
