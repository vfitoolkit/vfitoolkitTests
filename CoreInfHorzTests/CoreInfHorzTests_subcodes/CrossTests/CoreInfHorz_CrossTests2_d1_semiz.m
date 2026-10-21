function output=CoreInfHorz_CrossTests2_d1_semiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for InfHorz using a semi-exogenous state that is really just a markov (d2 does nothing to its
% transitions), with d1. Declared as a semiz, it has to give the same answer as when declared as a markov z.
% This is the independent check of the semiz code (value fn, agent dist and AggVars): the z code paths are already validated.
% Each of these is also done with the grid interpolation layer, with an iid e (as in the FHorz version), and with both.

Params.uempbenefit=0; % Need this to make return fns the same
Params.searcheffortcost=0; % d2 does nothing to the semiz transitions, and with no cost it does nothing at all

% The model declared without semiz has only the d1 (d2 does nothing)
n_d1=n_d(1);
d1_grid=d_grid(1:n_d1);

% Keep the original z (for the 'with z' cross-test)
n_zB=n_z;
z_gridB=z_grid;
pi_zB=pi_z;

% WARNING: THE z and semiz IN HERE ARE PRETTY HARDCODED!
n_z=2;
z_grid=[0.6;1.4];
pi_z=[1-Params.probfindjob, Params.probfindjob;...
    Params.problosejob, 1-Params.problosejob];
Params.z1=z_grid(1);
Params.z2=z_grid(2);

% The semiz is just a duplicate of z
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreInfHorzSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);
vfoptions=struct();
simoptions=struct();
vfoptions.n_semiz=n_z;
vfoptions.semiz_grid=z_grid;
vfoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.n_semiz=n_z;
simoptions.semiz_grid=z_grid;
simoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.d_grid=d_grid;

ReturnFn_semiz=@(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_z=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
% with uempbenefit=0 and searcheffortcost=0, ReturnFn_semiz is ReturnFn_z with semiz in place of z

ReturnFn_semizz=@(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_zz=@(d,aprime,a,zA,zB,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,zA,zB,r,w,sigma,eta,varphi);
% with uempbenefit=0 and searcheffortcost=0, ReturnFn_semizz is ReturnFn_zz with (semiz,z) in place of (zA,zB) (both are w*d*shock1*shock2)

% FnsToEvaluate, the same functions declared each way
FnsToEvaluate_z.assets=@(d,aprime,a,z) a;
FnsToEvaluate_z.earnings=@(d,aprime,a,z,w) w*d*z;
FnsToEvaluate_semiz.assets=@(d1,d2,aprime,a,semiz) a;
FnsToEvaluate_semiz.earnings=@(d1,d2,aprime,a,semiz,w) w*d1*semiz;
FnsToEvaluate_zz.assets=@(d,aprime,a,zA,zB) a;
FnsToEvaluate_zz.earnings=@(d,aprime,a,zA,zB,w) w*d*zA*zB;
FnsToEvaluate_semizz.assets=@(d1,d2,aprime,a,semiz,z) a;
FnsToEvaluate_semizz.earnings=@(d1,d2,aprime,a,semiz,z,w) w*d1*semiz*z;
% The AggVars below are all evaluated on the agent dist of the model declared without semiz: the two agent dists only
% agree to the iteration tolerance (checked above as StationaryDist), so evaluating both on the same one makes the
% AggVars comparison an exact identity of the semiz AggVars code.

%% Solving for model with one markov
% First, just use z (without semiz). Note only d1 (as d2 does nothing).
[V1A,Policy1A]=ValueFnIter_InfHorz(n_d1,n_a,n_z,d1_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],struct());
StationaryDist1A=StationaryDist_InfHorz(Policy1A,n_d1,n_a,n_z,pi_z,struct(),Params,[]);

% Second, use semiz (without z)
[V1B,Policy1B]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist1B=StationaryDist_InfHorz(Policy1B,n_d,n_a,0,[],simoptions,Params,[]);

Policy1Bshort=Policy1B([1,3],:,:); % remove the d2 policy (as it is not relevant, and is not in Policy1A)

fprintf('Cross test (semiz as z), V should be zero:             %.3e \n',max(abs(V1A(:)-V1B(:))))
fprintf('Cross test (semiz as z), Policy should be zero:        %.3e \n',max(abs(Policy1A(:)-Policy1Bshort(:))))
fprintf('Cross test (semiz as z), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist1A(:)-StationaryDist1B(:))))

AggVars1A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1A,Policy1A,FnsToEvaluate_z,Params,[],n_d1,n_a,n_z,d1_grid,a_grid,z_grid,struct());
AggVars1B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1A,Policy1B,FnsToEvaluate_semiz,Params,[],n_d,n_a,0,d_grid,a_grid,[],simoptions);
fprintf('Cross test (semiz as z), AggVars assets mean should be zero:   %.3e \n',abs(AggVars1A.assets.Mean-AggVars1B.assets.Mean))
fprintf('Cross test (semiz as z), AggVars earnings mean should be zero: %.3e \n',abs(AggVars1A.earnings.Mean-AggVars1B.earnings.Mean))

clear V1A Policy1A StationaryDist1A V1B Policy1B Policy1Bshort StationaryDist1B AggVars1A AggVars1B

%% Solving for model with one markov, plus the original z
% The semiz (really a markov) and the z together are the same thing as a two-dimensional markov (zA,zB), with zA
% varying first (as semiz does), whose transition matrix is kron(pi_zB,pi_zA).
n_zAB=[n_z,n_zB];
z_gridAB=[z_grid;z_gridB];
pi_zAB=kron(pi_zB,pi_z);

% First, just use the two-dimensional z (without semiz). Note only d1 (as d2 does nothing).
[V2A,Policy2A]=ValueFnIter_InfHorz(n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,pi_zAB,ReturnFn_zz,Params,DiscountFactorParamNames,[],struct());
StationaryDist2A=StationaryDist_InfHorz(Policy2A,n_d1,n_a,n_zAB,pi_zAB,struct(),Params,[]);

% Second, use semiz and z
[V2B,Policy2B]=ValueFnIter_InfHorz(n_d,n_a,n_zB,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_semizz,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist2B=StationaryDist_InfHorz(Policy2B,n_d,n_a,n_zB,pi_zB,simoptions,Params,[]);

Policy2Bshort=Policy2B([1,3],:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy2A)

fprintf('Cross test (semiz as z, with z), V should be zero:             %.3e \n',max(abs(V2A(:)-V2B(:))))
fprintf('Cross test (semiz as z, with z), Policy should be zero:        %.3e \n',max(abs(Policy2A(:)-Policy2Bshort(:))))
fprintf('Cross test (semiz as z, with z), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist2A(:)-StationaryDist2B(:))))

AggVars2A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist2A,Policy2A,FnsToEvaluate_zz,Params,[],n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,struct());
AggVars2B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist2A,Policy2B,FnsToEvaluate_semizz,Params,[],n_d,n_a,n_zB,d_grid,a_grid,z_gridB,simoptions);
fprintf('Cross test (semiz as z, with z), AggVars assets mean should be zero:   %.3e \n',abs(AggVars2A.assets.Mean-AggVars2B.assets.Mean))
fprintf('Cross test (semiz as z, with z), AggVars earnings mean should be zero: %.3e \n',abs(AggVars2A.earnings.Mean-AggVars2B.earnings.Mean))

clear V2A Policy2A StationaryDist2A V2B Policy2B Policy2Bshort StationaryDist2B AggVars2A AggVars2B

%% The same two, with the grid interpolation layer
% Both are solved with the same vfoptions.maxaprimediff (the defaults differ between models with and without a decision
% variable), so the two fine-grid problems are the same problem
vfoptionsGI=struct();
vfoptionsGI.gridinterplayer=1;
vfoptionsGI.ngridinterp=5;
vfoptionsGI.maxaprimediff=ceil(n_a/5);
simoptionsGI=struct();
simoptionsGI.gridinterplayer=1;
simoptionsGI.ngridinterp=5;
vfoptionsSGI=vfoptions; % semiz and GI
vfoptionsSGI.gridinterplayer=vfoptionsGI.gridinterplayer;
vfoptionsSGI.ngridinterp=vfoptionsGI.ngridinterp;
vfoptionsSGI.maxaprimediff=vfoptionsGI.maxaprimediff;
simoptionsSGI=simoptions;
simoptionsSGI.gridinterplayer=simoptionsGI.gridinterplayer;
simoptionsSGI.ngridinterp=simoptionsGI.ngridinterp;

% One markov: first z (without semiz), second semiz
[V5A,Policy5A]=ValueFnIter_InfHorz(n_d1,n_a,n_z,d1_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptionsGI);
StationaryDist5A=StationaryDist_InfHorz(Policy5A,n_d1,n_a,n_z,pi_z,simoptionsGI,Params,[]);
[V5B,Policy5B]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,[],vfoptionsSGI);
StationaryDist5B=StationaryDist_InfHorz(Policy5B,n_d,n_a,0,[],simoptionsSGI,Params,[]);

Policy5Bshort=Policy5B([1,3,4,5],:,:); % remove the d2 policy (as it is not relevant, and is not in Policy5A)

fprintf('Cross test (semiz as z, with GI), V should be zero:             %.3e \n',max(abs(V5A(:)-V5B(:))))
fprintf('Cross test (semiz as z, with GI), Policy should be zero:        %.3e \n',max(abs(Policy5A(:)-Policy5Bshort(:))))
fprintf('Cross test (semiz as z, with GI), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist5A(:)-StationaryDist5B(:))))

AggVars5A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist5A,Policy5A,FnsToEvaluate_z,Params,[],n_d1,n_a,n_z,d1_grid,a_grid,z_grid,simoptionsGI);
AggVars5B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist5A,Policy5B,FnsToEvaluate_semiz,Params,[],n_d,n_a,0,d_grid,a_grid,[],simoptionsSGI);
fprintf('Cross test (semiz as z, with GI), AggVars assets mean should be zero:   %.3e \n',abs(AggVars5A.assets.Mean-AggVars5B.assets.Mean))
fprintf('Cross test (semiz as z, with GI), AggVars earnings mean should be zero: %.3e \n',abs(AggVars5A.earnings.Mean-AggVars5B.earnings.Mean))

V5Bfrompolicy=ValueFnFromPolicy_InfHorz(Policy5B,n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_semiz,Params,DiscountFactorParamNames,vfoptionsSGI);
fprintf('Cross test (semiz as z, with GI), ValueFnFromPolicy of the semiz solve vs the z solve, this should be zero: %.3e \n',max(abs(V5A(:)-V5Bfrompolicy(:))))

clear V5A Policy5A StationaryDist5A V5B Policy5B Policy5Bshort StationaryDist5B AggVars5A AggVars5B V5Bfrompolicy

% One markov and the original z: first the two-dimensional z (without semiz), second semiz and z
[V6A,Policy6A]=ValueFnIter_InfHorz(n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,pi_zAB,ReturnFn_zz,Params,DiscountFactorParamNames,[],vfoptionsGI);
StationaryDist6A=StationaryDist_InfHorz(Policy6A,n_d1,n_a,n_zAB,pi_zAB,simoptionsGI,Params,[]);
[V6B,Policy6B]=ValueFnIter_InfHorz(n_d,n_a,n_zB,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_semizz,Params,DiscountFactorParamNames,[],vfoptionsSGI);
StationaryDist6B=StationaryDist_InfHorz(Policy6B,n_d,n_a,n_zB,pi_zB,simoptionsSGI,Params,[]);

Policy6Bshort=Policy6B([1,3,4,5],:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy6A)

fprintf('Cross test (semiz as z, with z, with GI), V should be zero:             %.3e \n',max(abs(V6A(:)-V6B(:))))
fprintf('Cross test (semiz as z, with z, with GI), Policy should be zero:        %.3e \n',max(abs(Policy6A(:)-Policy6Bshort(:))))
fprintf('Cross test (semiz as z, with z, with GI), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist6A(:)-StationaryDist6B(:))))

AggVars6A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist6A,Policy6A,FnsToEvaluate_zz,Params,[],n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,simoptionsGI);
AggVars6B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist6A,Policy6B,FnsToEvaluate_semizz,Params,[],n_d,n_a,n_zB,d_grid,a_grid,z_gridB,simoptionsSGI);
fprintf('Cross test (semiz as z, with z, with GI), AggVars assets mean should be zero:   %.3e \n',abs(AggVars6A.assets.Mean-AggVars6B.assets.Mean))
fprintf('Cross test (semiz as z, with z, with GI), AggVars earnings mean should be zero: %.3e \n',abs(AggVars6A.earnings.Mean-AggVars6B.earnings.Mean))

V6Bfrompolicy=ValueFnFromPolicy_InfHorz(Policy6B,n_d,n_a,n_zB,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_semizz,Params,DiscountFactorParamNames,vfoptionsSGI);
fprintf('Cross test (semiz as z, with z, with GI), ValueFnFromPolicy of the semiz solve vs the z solve, this should be zero: %.3e \n',max(abs(V6A(:)-V6Bfrompolicy(:))))

clear V6A Policy6A StationaryDist6A V6B Policy6B Policy6Bshort StationaryDist6B AggVars6A AggVars6B V6Bfrompolicy

%% The same two, with an iid e as well
vfoptionsE=struct(); % just e
vfoptionsE.n_e=vfoptionsbaseline.n_e;
vfoptionsE.e_grid=vfoptionsbaseline.e_grid;
vfoptionsE.pi_e=vfoptionsbaseline.pi_e;
simoptionsE=struct();
simoptionsE.n_e=simoptionsbaseline.n_e;
simoptionsE.e_grid=simoptionsbaseline.e_grid;
simoptionsE.pi_e=simoptionsbaseline.pi_e;
vfoptionsSE=vfoptions; % semiz and e
vfoptionsSE.n_e=vfoptionsE.n_e;
vfoptionsSE.e_grid=vfoptionsE.e_grid;
vfoptionsSE.pi_e=vfoptionsE.pi_e;
simoptionsSE=simoptions;
simoptionsSE.n_e=simoptionsE.n_e;
simoptionsSE.e_grid=simoptionsE.e_grid;
simoptionsSE.pi_e=simoptionsE.pi_e;

ReturnFn_semize=@(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_ze=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
ReturnFn_semizze=@(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_e_semiz(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_zze=@(d,aprime,a,zA,zB,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,zA*zB,e,r,w,sigma,eta,varphi);
% with uempbenefit=0 and searcheffortcost=0, ReturnFn_semize is ReturnFn_ze with semiz in place of z, and ReturnFn_semizze is
% ReturnFn_zze with (semiz,z) in place of (zA,zB) (all are w*d*shocks*e)
FnsToEvaluate_ze.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate_ze.earnings=@(d,aprime,a,z,e,w) w*d*z*e;
FnsToEvaluate_semize.assets=@(d1,d2,aprime,a,semiz,e) a;
FnsToEvaluate_semize.earnings=@(d1,d2,aprime,a,semiz,e,w) w*d1*semiz*e;
FnsToEvaluate_zze.assets=@(d,aprime,a,zA,zB,e) a;
FnsToEvaluate_zze.earnings=@(d,aprime,a,zA,zB,e,w) w*d*zA*zB*e;
FnsToEvaluate_semizze.assets=@(d1,d2,aprime,a,semiz,z,e) a;
FnsToEvaluate_semizze.earnings=@(d1,d2,aprime,a,semiz,z,e,w) w*d1*semiz*z*e;

% One markov and e: first z and e (without semiz), second semiz and e
[V3A,Policy3A]=ValueFnIter_InfHorz(n_d1,n_a,n_z,d1_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptionsE);
StationaryDist3A=StationaryDist_InfHorz(Policy3A,n_d1,n_a,n_z,pi_z,simoptionsE,Params,[]);
[V3B,Policy3B]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_semize,Params,DiscountFactorParamNames,[],vfoptionsSE);
StationaryDist3B=StationaryDist_InfHorz(Policy3B,n_d,n_a,0,[],simoptionsSE,Params,[]);

Policy3Bshort=Policy3B([1,3],:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy3A)

fprintf('Cross test (semiz as z, with e), V should be zero:             %.3e \n',max(abs(V3A(:)-V3B(:))))
fprintf('Cross test (semiz as z, with e), Policy should be zero:        %.3e \n',max(abs(Policy3A(:)-Policy3Bshort(:))))
fprintf('Cross test (semiz as z, with e), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist3A(:)-StationaryDist3B(:))))

AggVars3A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist3A,Policy3A,FnsToEvaluate_ze,Params,[],n_d1,n_a,n_z,d1_grid,a_grid,z_grid,simoptionsE);
AggVars3B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist3A,Policy3B,FnsToEvaluate_semize,Params,[],n_d,n_a,0,d_grid,a_grid,[],simoptionsSE);
fprintf('Cross test (semiz as z, with e), AggVars assets mean should be zero:   %.3e \n',abs(AggVars3A.assets.Mean-AggVars3B.assets.Mean))
fprintf('Cross test (semiz as z, with e), AggVars earnings mean should be zero: %.3e \n',abs(AggVars3A.earnings.Mean-AggVars3B.earnings.Mean))

clear V3A Policy3A StationaryDist3A V3B Policy3B Policy3Bshort StationaryDist3B AggVars3A AggVars3B

% One markov, the original z and e: first the two-dimensional z and e (without semiz), second semiz, z and e
[V4A,Policy4A]=ValueFnIter_InfHorz(n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,pi_zAB,ReturnFn_zze,Params,DiscountFactorParamNames,[],vfoptionsE);
StationaryDist4A=StationaryDist_InfHorz(Policy4A,n_d1,n_a,n_zAB,pi_zAB,simoptionsE,Params,[]);
[V4B,Policy4B]=ValueFnIter_InfHorz(n_d,n_a,n_zB,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_semizze,Params,DiscountFactorParamNames,[],vfoptionsSE);
StationaryDist4B=StationaryDist_InfHorz(Policy4B,n_d,n_a,n_zB,pi_zB,simoptionsSE,Params,[]);

Policy4Bshort=Policy4B([1,3],:,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy4A)

fprintf('Cross test (semiz as z, with z and e), V should be zero:             %.3e \n',max(abs(V4A(:)-V4B(:))))
fprintf('Cross test (semiz as z, with z and e), Policy should be zero:        %.3e \n',max(abs(Policy4A(:)-Policy4Bshort(:))))
fprintf('Cross test (semiz as z, with z and e), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist4A(:)-StationaryDist4B(:))))

AggVars4A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist4A,Policy4A,FnsToEvaluate_zze,Params,[],n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,simoptionsE);
AggVars4B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist4A,Policy4B,FnsToEvaluate_semizze,Params,[],n_d,n_a,n_zB,d_grid,a_grid,z_gridB,simoptionsSE);
fprintf('Cross test (semiz as z, with z and e), AggVars assets mean should be zero:   %.3e \n',abs(AggVars4A.assets.Mean-AggVars4B.assets.Mean))
fprintf('Cross test (semiz as z, with z and e), AggVars earnings mean should be zero: %.3e \n',abs(AggVars4A.earnings.Mean-AggVars4B.earnings.Mean))

clear V4A Policy4A StationaryDist4A V4B Policy4B Policy4Bshort StationaryDist4B AggVars4A AggVars4B

%% The same two with e, with the grid interpolation layer as well
vfoptionsEGI=vfoptionsE;
vfoptionsEGI.gridinterplayer=vfoptionsGI.gridinterplayer;
vfoptionsEGI.ngridinterp=vfoptionsGI.ngridinterp;
vfoptionsEGI.maxaprimediff=vfoptionsGI.maxaprimediff;
simoptionsEGI=simoptionsE;
simoptionsEGI.gridinterplayer=simoptionsGI.gridinterplayer;
simoptionsEGI.ngridinterp=simoptionsGI.ngridinterp;
vfoptionsSEGI=vfoptionsSE;
vfoptionsSEGI.gridinterplayer=vfoptionsGI.gridinterplayer;
vfoptionsSEGI.ngridinterp=vfoptionsGI.ngridinterp;
vfoptionsSEGI.maxaprimediff=vfoptionsGI.maxaprimediff;
simoptionsSEGI=simoptionsSE;
simoptionsSEGI.gridinterplayer=simoptionsGI.gridinterplayer;
simoptionsSEGI.ngridinterp=simoptionsGI.ngridinterp;

% One markov and e: first z and e (without semiz), second semiz and e
[V7A,Policy7A]=ValueFnIter_InfHorz(n_d1,n_a,n_z,d1_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptionsEGI);
StationaryDist7A=StationaryDist_InfHorz(Policy7A,n_d1,n_a,n_z,pi_z,simoptionsEGI,Params,[]);
[V7B,Policy7B]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_semize,Params,DiscountFactorParamNames,[],vfoptionsSEGI);
StationaryDist7B=StationaryDist_InfHorz(Policy7B,n_d,n_a,0,[],simoptionsSEGI,Params,[]);

Policy7Bshort=Policy7B([1,3,4,5],:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy7A)

fprintf('Cross test (semiz as z, with e, with GI), V should be zero:             %.3e \n',max(abs(V7A(:)-V7B(:))))
fprintf('Cross test (semiz as z, with e, with GI), Policy should be zero:        %.3e \n',max(abs(Policy7A(:)-Policy7Bshort(:))))
fprintf('Cross test (semiz as z, with e, with GI), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist7A(:)-StationaryDist7B(:))))

AggVars7A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist7A,Policy7A,FnsToEvaluate_ze,Params,[],n_d1,n_a,n_z,d1_grid,a_grid,z_grid,simoptionsEGI);
AggVars7B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist7A,Policy7B,FnsToEvaluate_semize,Params,[],n_d,n_a,0,d_grid,a_grid,[],simoptionsSEGI);
fprintf('Cross test (semiz as z, with e, with GI), AggVars assets mean should be zero:   %.3e \n',abs(AggVars7A.assets.Mean-AggVars7B.assets.Mean))
fprintf('Cross test (semiz as z, with e, with GI), AggVars earnings mean should be zero: %.3e \n',abs(AggVars7A.earnings.Mean-AggVars7B.earnings.Mean))

V7Bfrompolicy=ValueFnFromPolicy_InfHorz(Policy7B,n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_semize,Params,DiscountFactorParamNames,vfoptionsSEGI);
fprintf('Cross test (semiz as z, with e, with GI), ValueFnFromPolicy of the semiz solve vs the z solve, this should be zero: %.3e \n',max(abs(V7A(:)-V7Bfrompolicy(:))))

clear V7A Policy7A StationaryDist7A V7B Policy7B Policy7Bshort StationaryDist7B AggVars7A AggVars7B V7Bfrompolicy

% One markov, the original z and e: first the two-dimensional z and e (without semiz), second semiz, z and e
[V8A,Policy8A]=ValueFnIter_InfHorz(n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,pi_zAB,ReturnFn_zze,Params,DiscountFactorParamNames,[],vfoptionsEGI);
StationaryDist8A=StationaryDist_InfHorz(Policy8A,n_d1,n_a,n_zAB,pi_zAB,simoptionsEGI,Params,[]);
[V8B,Policy8B]=ValueFnIter_InfHorz(n_d,n_a,n_zB,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_semizze,Params,DiscountFactorParamNames,[],vfoptionsSEGI);
StationaryDist8B=StationaryDist_InfHorz(Policy8B,n_d,n_a,n_zB,pi_zB,simoptionsSEGI,Params,[]);

Policy8Bshort=Policy8B([1,3,4,5],:,:,:,:); % remove the d2 policy (as it is not relevant, and is not in Policy8A)

fprintf('Cross test (semiz as z, with z and e, with GI), V should be zero:             %.3e \n',max(abs(V8A(:)-V8B(:))))
fprintf('Cross test (semiz as z, with z and e, with GI), Policy should be zero:        %.3e \n',max(abs(Policy8A(:)-Policy8Bshort(:))))
fprintf('Cross test (semiz as z, with z and e, with GI), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist8A(:)-StationaryDist8B(:))))

AggVars8A=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist8A,Policy8A,FnsToEvaluate_zze,Params,[],n_d1,n_a,n_zAB,d1_grid,a_grid,z_gridAB,simoptionsEGI);
AggVars8B=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist8A,Policy8B,FnsToEvaluate_semizze,Params,[],n_d,n_a,n_zB,d_grid,a_grid,z_gridB,simoptionsSEGI);
fprintf('Cross test (semiz as z, with z and e, with GI), AggVars assets mean should be zero:   %.3e \n',abs(AggVars8A.assets.Mean-AggVars8B.assets.Mean))
fprintf('Cross test (semiz as z, with z and e, with GI), AggVars earnings mean should be zero: %.3e \n',abs(AggVars8A.earnings.Mean-AggVars8B.earnings.Mean))

V8Bfrompolicy=ValueFnFromPolicy_InfHorz(Policy8B,n_d,n_a,n_zB,d_grid,a_grid,z_gridB,pi_zB,ReturnFn_semizze,Params,DiscountFactorParamNames,vfoptionsSEGI);
fprintf('Cross test (semiz as z, with z and e, with GI), ValueFnFromPolicy of the semiz solve vs the z solve, this should be zero: %.3e \n',max(abs(V8A(:)-V8Bfrompolicy(:))))

clear V8A Policy8A StationaryDist8A V8B Policy8B Policy8Bshort StationaryDist8B AggVars8A AggVars8B V8Bfrompolicy

%%
output=struct();

end
