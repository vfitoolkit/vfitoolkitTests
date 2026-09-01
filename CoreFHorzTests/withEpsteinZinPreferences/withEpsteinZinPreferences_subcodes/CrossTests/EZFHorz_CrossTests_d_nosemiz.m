function output=EZFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Epstein-Zin mirror of QHDFHorz_CrossTests_d_nosemiz.m, run for the three EZ cases.
% For EZ these cross tests are NOT redundant (unlike QH, where they mostly re-check the vNM
% machinery): they check that an iid shock enters the same single certainty-equivalent whether
% it lives in z or in e. This is exactly the nested-vs-joint CE bug class.
% NOTE: the no-shock legs print the (expected) toolkit warning that Epstein-Zin does not make
% much sense without shocks.

% For crosstests, set up z to just be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in same place in earnings

ReturnFn_none_cons=@(d,aprime,a,r,w,kappa_j,varphi,agej,Jr,pension)...
    EZReturnFn_cons_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_z_cons=@(d,aprime,a,z,r,w,kappa_j,varphi,agej,Jr,pension)...
    EZReturnFn_cons_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_e_cons=@(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension)...
    EZReturnFn_cons_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_ze_cons=@(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension)...
    EZReturnFn_cons_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension);
ReturnFn_none_posU=@(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_positiveUtils_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_z_posU=@(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_positiveUtils_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_e_posU=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_positiveUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_ze_posU=@(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_positiveUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_none_negU=@(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_negativeUtils_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_z_negU=@(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_negativeUtils_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_e_negU=@(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_negativeUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);
ReturnFn_ze_negU=@(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)...
    EZReturnFn_negativeUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension);


%% Epstein-Zin case 1: Consumption-units (EZutils=0)
%% Solving with just a single point for z with value 1 and prob 1 gives us same as no shocks
jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1; % no assets

vfoptionsA=struct();
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=0;
vfoptionsA.EZriskaversion='ezgamma';
vfoptionsA.EZeis='ezphi';
simoptionsA=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none_cons,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z_cons,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% And a single point for z vs a single point for e
vfoptionsB=struct();
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=0;
vfoptionsB.EZriskaversion='ezgamma';
vfoptionsB.EZeis='ezphi';
vfoptionsB.n_e=1;
vfoptionsB.e_grid=1;
vfoptionsB.pi_e=1;
simoptionsB=struct();
simoptionsB.n_e=1;
simoptionsB.e_grid=1;
simoptionsB.pi_e=1;
[V0e,Policy0e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_cons,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist0e=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0e,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(V0z(:)-V0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(Policy0z(:)-Policy0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(StationaryDist0z(:)-StationaryDist0e(:))))

clear V0 V0z V0e Policy0 Policy0z Policy0e StationaryDist0 StationaryDist0z StationaryDist0e

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid as e
% zeros assets, mid points for any shocks
jequaloneDist_z=zeros(n_a,n_z,'gpuArray');
jequaloneDist_z(1,ceil(n_z/2))=1; % no assets, midpoint shock

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z_cons,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsC=struct();
vfoptionsC.exoticpreferences='EpsteinZin';
vfoptionsC.EZutils=0;
vfoptionsC.EZriskaversion='ezgamma';
vfoptionsC.EZeis='ezphi';
vfoptionsC.n_e=vfoptionsbaseline.n_e;
vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
vfoptionsC.pi_e=vfoptionsbaseline.pi_e;
simoptionsC=struct();
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
jequaloneDist3=zeros(n_a,1,vfoptionsC.n_e,'gpuArray');
jequaloneDist3(1,1,ceil(vfoptionsC.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsC);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1
vfoptionsD=struct();
vfoptionsD.exoticpreferences='EpsteinZin';
vfoptionsD.EZutils=0;
vfoptionsD.EZriskaversion='ezgamma';
vfoptionsD.EZeis='ezphi';
vfoptionsD.n_e=1;
vfoptionsD.e_grid=1;
vfoptionsD.pi_e=1;
simoptionsD=struct();
simoptionsD.n_e=1;
simoptionsD.e_grid=1;
simoptionsD.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_cons,Params,DiscountFactorParamNames,[],vfoptionsD);
jequaloneDist4=zeros(n_a,n_z,1,'gpuArray');
jequaloneDist4(1,ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsD);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))

clear V1 V3 V4 Policy1 Policy3 Policy4 StationaryDist1 StationaryDist3 StationaryDist4

%% Epstein-Zin case 2: Utility-units, positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%% Solving with just a single point for z with value 1 and prob 1 gives us same as no shocks
jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1; % no assets

vfoptionsA=struct();
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=1;
vfoptionsA.EZpositiveutility=1;
vfoptionsA.EZriskaversion='ezrisk';
simoptionsA=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none_posU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z_posU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% And a single point for z vs a single point for e
vfoptionsB=struct();
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=1;
vfoptionsB.EZpositiveutility=1;
vfoptionsB.EZriskaversion='ezrisk';
vfoptionsB.n_e=1;
vfoptionsB.e_grid=1;
vfoptionsB.pi_e=1;
simoptionsB=struct();
simoptionsB.n_e=1;
simoptionsB.e_grid=1;
simoptionsB.pi_e=1;
[V0e,Policy0e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_posU,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist0e=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0e,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(V0z(:)-V0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(Policy0z(:)-Policy0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(StationaryDist0z(:)-StationaryDist0e(:))))

clear V0 V0z V0e Policy0 Policy0z Policy0e StationaryDist0 StationaryDist0z StationaryDist0e

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid as e
% zeros assets, mid points for any shocks
jequaloneDist_z=zeros(n_a,n_z,'gpuArray');
jequaloneDist_z(1,ceil(n_z/2))=1; % no assets, midpoint shock

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z_posU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsC=struct();
vfoptionsC.exoticpreferences='EpsteinZin';
vfoptionsC.EZutils=1;
vfoptionsC.EZpositiveutility=1;
vfoptionsC.EZriskaversion='ezrisk';
vfoptionsC.n_e=vfoptionsbaseline.n_e;
vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
vfoptionsC.pi_e=vfoptionsbaseline.pi_e;
simoptionsC=struct();
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
jequaloneDist3=zeros(n_a,1,vfoptionsC.n_e,'gpuArray');
jequaloneDist3(1,1,ceil(vfoptionsC.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsC);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1
vfoptionsD=struct();
vfoptionsD.exoticpreferences='EpsteinZin';
vfoptionsD.EZutils=1;
vfoptionsD.EZpositiveutility=1;
vfoptionsD.EZriskaversion='ezrisk';
vfoptionsD.n_e=1;
vfoptionsD.e_grid=1;
vfoptionsD.pi_e=1;
simoptionsD=struct();
simoptionsD.n_e=1;
simoptionsD.e_grid=1;
simoptionsD.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_posU,Params,DiscountFactorParamNames,[],vfoptionsD);
jequaloneDist4=zeros(n_a,n_z,1,'gpuArray');
jequaloneDist4(1,ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsD);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))

clear V1 V3 V4 Policy1 Policy3 Policy4 StationaryDist1 StationaryDist3 StationaryDist4

%% Epstein-Zin case 3: Utility-units, negative-valued utility fn (EZutils=1, EZpositiveutility=0)
%% Solving with just a single point for z with value 1 and prob 1 gives us same as no shocks
jequaloneDist_none=zeros(n_a,1,'gpuArray');
jequaloneDist_none(1)=1; % no assets

vfoptionsA=struct();
vfoptionsA.exoticpreferences='EpsteinZin';
vfoptionsA.EZutils=1;
vfoptionsA.EZpositiveutility=0;
vfoptionsA.EZriskaversion='ezrisk';
simoptionsA=struct();
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none_negU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z_negU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test: no shocks vs single-point z, this should be zero: %.3e \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% And a single point for z vs a single point for e
vfoptionsB=struct();
vfoptionsB.exoticpreferences='EpsteinZin';
vfoptionsB.EZutils=1;
vfoptionsB.EZpositiveutility=0;
vfoptionsB.EZriskaversion='ezrisk';
vfoptionsB.n_e=1;
vfoptionsB.e_grid=1;
vfoptionsB.pi_e=1;
simoptionsB=struct();
simoptionsB.n_e=1;
simoptionsB.e_grid=1;
simoptionsB.pi_e=1;
[V0e,Policy0e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e_negU,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist0e=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0e,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(V0z(:)-V0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(Policy0z(:)-Policy0e(:))))
fprintf('Cross test: single-point z vs single-point e, this should be zero: %.3e \n',max(abs(StationaryDist0z(:)-StationaryDist0e(:))))

clear V0 V0z V0e Policy0 Policy0z Policy0e StationaryDist0 StationaryDist0z StationaryDist0e

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid as e
% zeros assets, mid points for any shocks
jequaloneDist_z=zeros(n_a,n_z,'gpuArray');
jequaloneDist_z(1,ceil(n_z/2))=1; % no assets, midpoint shock

[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z_negU,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

vfoptionsC=struct();
vfoptionsC.exoticpreferences='EpsteinZin';
vfoptionsC.EZutils=1;
vfoptionsC.EZpositiveutility=0;
vfoptionsC.EZriskaversion='ezrisk';
vfoptionsC.n_e=vfoptionsbaseline.n_e;
vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
vfoptionsC.pi_e=vfoptionsbaseline.pi_e;
simoptionsC=struct();
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
jequaloneDist3=zeros(n_a,1,vfoptionsC.n_e,'gpuArray');
jequaloneDist3(1,1,ceil(vfoptionsC.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsC);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1
vfoptionsD=struct();
vfoptionsD.exoticpreferences='EpsteinZin';
vfoptionsD.EZutils=1;
vfoptionsD.EZpositiveutility=0;
vfoptionsD.EZriskaversion='ezrisk';
vfoptionsD.n_e=1;
vfoptionsD.e_grid=1;
vfoptionsD.pi_e=1;
simoptionsD=struct();
simoptionsD.n_e=1;
simoptionsD.e_grid=1;
simoptionsD.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze_negU,Params,DiscountFactorParamNames,[],vfoptionsD);
jequaloneDist4=zeros(n_a,n_z,1,'gpuArray');
jequaloneDist4(1,ceil(n_z/2),1)=1; % no assets, midpoint shock
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
