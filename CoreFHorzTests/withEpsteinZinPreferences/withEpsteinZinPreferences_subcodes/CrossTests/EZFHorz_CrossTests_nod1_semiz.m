function output=EZFHorz_CrossTests_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Epstein-Zin mirror of QHDFHorz_CrossTests_nod1_semiz.m, run for the three EZ cases.
% TEST-FIRST: EZ with semi-exogenous shocks is NOT yet implemented in the toolkit;
% everything in this file is EXPECTED TO ERROR until the EZ SemiExo solvers exist.
% NOTE: the QH version starts by comparing against a model with no z/e at all; Epstein-Zin
% currently errors on that even with semiz present, so the first leg compares a single-point z
% against a single-point e instead (both with semiz).

% n_d=n_d2_semiz;
% d_grid=d2_grid_semiz;

% Setup semiz
vfoptions_semiz.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions_semiz.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions_semiz.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions_semiz.n_semiz=simoptionsbaseline.n_semiz;
simoptions_semiz.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions_semiz.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions_semiz.d_grid=d_grid;
% For convenience
n_semiz=vfoptionsbaseline.n_semiz;

% For crosstests, set up z to just be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in same place in earnings

ReturnFn_z_cons=@(d2,aprime,a,semiz,z,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_cons_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_cons=@(d2,aprime,a,semiz,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_cons_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze_cons=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_cons_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_posU=@(d2,aprime,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_positiveUtils_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_posU=@(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_positiveUtils_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze_posU=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_positiveUtils_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_negU=@(d2,aprime,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_negativeUtils_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_negU=@(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_negativeUtils_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze_negU=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    EZReturnFn_negativeUtils_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost);


%% Epstein-Zin case 1: Consumption-units (EZutils=0)
%% Single point for z (value 1, prob 1) vs single point for e (both with semiz)
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(n_semiz/2))=1; % no assets

vfoptionsA=vfoptions_semiz;
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=0;
vfoptionsA.EZriskaversion='ezgamma';
vfoptionsA.EZeis='ezphi';
simoptionsA=simoptions_semiz;
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z_cons,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

vfoptionsB=vfoptions_semiz;
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=0;
vfoptionsB.EZriskaversion='ezgamma';
vfoptionsB.EZeis='ezphi';
vfoptionsB.n_e=1;
vfoptionsB.e_grid=1;
vfoptionsB.pi_e=1;
simoptionsB=simoptions_semiz;
simoptionsB.n_e=1;
simoptionsB.e_grid=1;
simoptionsB.pi_e=1;
[V0e,Policy0e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_cons,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist0e=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0e,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(V0z(:)-V0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(Policy0z(:)-Policy0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(StationaryDist0z(:)-StationaryDist0e(:))))

clear V0z V0e Policy0z Policy0e StationaryDist0z StationaryDist0e

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid as e
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1; % no assets, midpoint shock

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z_cons,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsC=vfoptions_semiz;
vfoptionsC.exoticpreferences='EpsteinZin';
vfoptionsC.EZutils=0;
vfoptionsC.EZriskaversion='ezgamma';
vfoptionsC.EZeis='ezphi';
vfoptionsC.n_e=vfoptionsbaseline.n_e;
vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
vfoptionsC.pi_e=vfoptionsbaseline.pi_e;
simoptionsC=simoptions_semiz;
simoptionsC.n_e=simoptionsbaseline.n_e;
simoptionsC.e_grid=simoptionsbaseline.e_grid;
simoptionsC.pi_e=simoptionsbaseline.pi_e;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_cons,Params,DiscountFactorParamNames,[],vfoptionsC);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptionsC);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist2(:))))

clear V2 Policy2 StationaryDist2

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze_cons,Params,DiscountFactorParamNames,[],vfoptionsC);
jequaloneDist3=zeros([n_a,n_semiz,1,vfoptionsC.n_e],'gpuArray');
jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptionsC.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsC);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1 (with semiz)
vfoptionsD=vfoptions_semiz;
vfoptionsD.exoticpreferences='EpsteinZin';
vfoptionsD.EZutils=0;
vfoptionsD.EZriskaversion='ezgamma';
vfoptionsD.EZeis='ezphi';
vfoptionsD.n_e=1;
vfoptionsD.e_grid=1;
vfoptionsD.pi_e=1;
simoptionsD=simoptions_semiz;
simoptionsD.n_e=1;
simoptionsD.e_grid=1;
simoptionsD.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_cons,Params,DiscountFactorParamNames,[],vfoptionsD);
jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray');
jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsD);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))

clear V1 V3 V4 Policy1 Policy3 Policy4 StationaryDist1 StationaryDist3 StationaryDist4

%% Epstein-Zin case 2: Utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%% Single point for z (value 1, prob 1) vs single point for e (both with semiz)
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(n_semiz/2))=1; % no assets

vfoptionsA=vfoptions_semiz;
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=1;
vfoptionsA.EZpositiveutility=1;
vfoptionsA.EZriskaversion='ezrisk';
simoptionsA=simoptions_semiz;
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z_posU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

vfoptionsB=vfoptions_semiz;
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=1;
vfoptionsB.EZpositiveutility=1;
vfoptionsB.EZriskaversion='ezrisk';
vfoptionsB.n_e=1;
vfoptionsB.e_grid=1;
vfoptionsB.pi_e=1;
simoptionsB=simoptions_semiz;
simoptionsB.n_e=1;
simoptionsB.e_grid=1;
simoptionsB.pi_e=1;
[V0e,Policy0e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_posU,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist0e=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0e,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(V0z(:)-V0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(Policy0z(:)-Policy0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(StationaryDist0z(:)-StationaryDist0e(:))))

clear V0z V0e Policy0z Policy0e StationaryDist0z StationaryDist0e

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid as e
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1; % no assets, midpoint shock

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z_posU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsC=vfoptions_semiz;
vfoptionsC.exoticpreferences='EpsteinZin';
vfoptionsC.EZutils=1;
vfoptionsC.EZpositiveutility=1;
vfoptionsC.EZriskaversion='ezrisk';
vfoptionsC.n_e=vfoptionsbaseline.n_e;
vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
vfoptionsC.pi_e=vfoptionsbaseline.pi_e;
simoptionsC=simoptions_semiz;
simoptionsC.n_e=simoptionsbaseline.n_e;
simoptionsC.e_grid=simoptionsbaseline.e_grid;
simoptionsC.pi_e=simoptionsbaseline.pi_e;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_posU,Params,DiscountFactorParamNames,[],vfoptionsC);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptionsC);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist2(:))))

clear V2 Policy2 StationaryDist2

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze_posU,Params,DiscountFactorParamNames,[],vfoptionsC);
jequaloneDist3=zeros([n_a,n_semiz,1,vfoptionsC.n_e],'gpuArray');
jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptionsC.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsC);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1 (with semiz)
vfoptionsD=vfoptions_semiz;
vfoptionsD.exoticpreferences='EpsteinZin';
vfoptionsD.EZutils=1;
vfoptionsD.EZpositiveutility=1;
vfoptionsD.EZriskaversion='ezrisk';
vfoptionsD.n_e=1;
vfoptionsD.e_grid=1;
vfoptionsD.pi_e=1;
simoptionsD=simoptions_semiz;
simoptionsD.n_e=1;
simoptionsD.e_grid=1;
simoptionsD.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_posU,Params,DiscountFactorParamNames,[],vfoptionsD);
jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray');
jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsD);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))

clear V1 V3 V4 Policy1 Policy3 Policy4 StationaryDist1 StationaryDist3 StationaryDist4

%% Epstein-Zin case 3: Utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
%% Single point for z (value 1, prob 1) vs single point for e (both with semiz)
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(n_semiz/2))=1; % no assets

vfoptionsA=vfoptions_semiz;
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=1;
vfoptionsA.EZpositiveutility=0;
vfoptionsA.EZriskaversion='ezrisk';
simoptionsA=simoptions_semiz;
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z_negU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

vfoptionsB=vfoptions_semiz;
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=1;
vfoptionsB.EZpositiveutility=0;
vfoptionsB.EZriskaversion='ezrisk';
vfoptionsB.n_e=1;
vfoptionsB.e_grid=1;
vfoptionsB.pi_e=1;
simoptionsB=simoptions_semiz;
simoptionsB.n_e=1;
simoptionsB.e_grid=1;
simoptionsB.pi_e=1;
[V0e,Policy0e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_negU,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist0e=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0e,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(V0z(:)-V0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(Policy0z(:)-Policy0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(StationaryDist0z(:)-StationaryDist0e(:))))

clear V0z V0e Policy0z Policy0e StationaryDist0z StationaryDist0e

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid as e
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1; % no assets, midpoint shock

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z_negU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsC=vfoptions_semiz;
vfoptionsC.exoticpreferences='EpsteinZin';
vfoptionsC.EZutils=1;
vfoptionsC.EZpositiveutility=0;
vfoptionsC.EZriskaversion='ezrisk';
vfoptionsC.n_e=vfoptionsbaseline.n_e;
vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
vfoptionsC.pi_e=vfoptionsbaseline.pi_e;
simoptionsC=simoptions_semiz;
simoptionsC.n_e=simoptionsbaseline.n_e;
simoptionsC.e_grid=simoptionsbaseline.e_grid;
simoptionsC.pi_e=simoptionsbaseline.pi_e;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_negU,Params,DiscountFactorParamNames,[],vfoptionsC);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptionsC);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist2(:))))

clear V2 Policy2 StationaryDist2

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze_negU,Params,DiscountFactorParamNames,[],vfoptionsC);
jequaloneDist3=zeros([n_a,n_semiz,1,vfoptionsC.n_e],'gpuArray');
jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptionsC.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsC);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1 (with semiz)
vfoptionsD=vfoptions_semiz;
vfoptionsD.exoticpreferences='EpsteinZin';
vfoptionsD.EZutils=1;
vfoptionsD.EZpositiveutility=0;
vfoptionsD.EZriskaversion='ezrisk';
vfoptionsD.n_e=1;
vfoptionsD.e_grid=1;
vfoptionsD.pi_e=1;
simoptionsD=simoptions_semiz;
simoptionsD.n_e=1;
simoptionsD.e_grid=1;
simoptionsD.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_negU,Params,DiscountFactorParamNames,[],vfoptionsD);
jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray');
jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsD);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))

clear V1 V3 V4 Policy1 Policy3 Policy4 StationaryDist1 StationaryDist3 StationaryDist4

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
