function output=EstimFHorz_ShocksCustomStats(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 3: shocks that depend on estimated parameters, jequaloneDist as a function, CustomModelStats, through
% EstimateLifeCycleModel_MethodOfMoments. The grid interpolation layer is on.
%   (a) vfoptions.ExogShockFn (rho_z, sigma_z) with rho_z estimated: at the truth (skipestimation=1) the objective is zero and the
%       Jacobian equals finite differences the test takes itself with the grids rebuilt at each step (exact); started at rho_z=0.8
%       (constrain0to1) it is recovered (close)
%   (b) vfoptions.EiidShockFn (sigma_e) with sigma_e estimated, an iid shock e with three points: the same two checks (truth 0.1,
%       started at 0.15)
%   (c) jequaloneDist given as a function, under simoptions.jequaloneDist_usergrids=1 (the user's z_grid) and =0 (the internal joint
%       grid): at the truth the objective is zero and J equals the J of the matrix jequaloneDist (exact)
%   (d) TargetMoments.CustomModelStats (a by-hand mean of assets, and a check of which z_grid is handed to the custom function) under
%       both CustomModelStats_usergrids settings: at the truth the objective is zero and the current moments equal the targets (exact)
% The shock, initial-distribution and custom-statistics functions are the CalibFHorzTests ones (on the path). Each driver run is wrapped
% so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
estimoptions=struct(); estimoptions.verbose=0;

%% (a) ExogShockFn with rho_z estimated
Params_a=Params; Params_a.rho_z=0.9; Params_a.sigma_z=0.03; % (the bank's z grid is the Farmer-Toda discretisation at these values)
vfoptions_a=vfoptions; vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibFHorz_ExogShockFn(rho_z,sigma_z);
simoptions_a=simoptions; simoptions_a.ExogShockFn=vfoptions_a.ExogShockFn;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params_a,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params_a,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params_a,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TMa=struct();
t=gather(LCPFull.earnings.StdDeviation); t(Jr:end)=NaN; TMa.AgeConditionalStats.earnings.StdDeviation=t; % (earnings are zero in retirement)
TMa.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
nT=(Jr-1)+1; W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
estimoptions_a=estimoptions; estimoptions_a.skipestimation=1;
estimoptions_a.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
ran=0; objv=NaN; maxd=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'rho_z'},TMa,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions_a,simoptions_a);
    ran=1; objv=estsummary.objectivefnval;
    dc=estsummary.doublechecks.epsilon10000;
    mu=EstimFHorz_ModelMoments(dc.estimparamsvecup(:),{'rho_z'},TMa,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions_a,simoptions_a,0,[]);
    md=EstimFHorz_ModelMoments(dc.estimparamsvecdown(:),{'rho_z'},TMa,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions_a,simoptions_a,0,[]);
    Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
    maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
catch err
    fprintf('ShocksCustomStats (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats (a) ExogShockFn, rho_z, skipestimation at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J_centered vs finite differences with the grids rebuilt, this should be zero: %.3e \n',1-ran,objv,maxd)
Params_start=Params_a; Params_start.rho_z=0.8;
estimoptions_a2=estimoptions; estimoptions_a2.fminalgo=8; estimoptions_a2.constrain0to1={'rho_z'};
ran=0; perr=NaN; tt=NaN; objv=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'rho_z'},TMa,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_a2,vfoptions_a,simoptions_a);
    tt=toc(t0); ran=1; perr=abs(EstimParams.rho_z-0.9); objv=estsummary.objectivefnval;
catch err
    fprintf('ShocksCustomStats (a) from 0.8: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats (a) ExogShockFn, rho_z started at 0.8 (truth 0.9) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)

%% (b) EiidShockFn with sigma_e estimated
ReturnFn_e=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate_e.assets=@(aprime,a,z,e) a;
FnsToEvaluate_e.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
n_e=3;
Params_b=Params_a; Params_b.sigma_e=0.1;
[e_grid,pi_e]=CalibFHorz_EiidShockFn(Params_b.sigma_e);
vfoptions_e=vfoptions; vfoptions_e.n_e=n_e; vfoptions_e.e_grid=e_grid; vfoptions_e.pi_e=pi_e;
simoptions_e=simoptions; simoptions_e.n_e=n_e; simoptions_e.e_grid=e_grid; simoptions_e.pi_e=pi_e;
simoptions_e_full=simoptions_e; simoptions_e_full.whichstats=ones(1,7);
jequaloneDist_e=zeros(n_a,n_z,n_e,'gpuArray'); jequaloneDist_e(1,ceil(n_z/2),ceil(n_e/2))=1;
[~,Policy_e]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,[],vfoptions_e);
StationaryDist_e=StationaryDist_FHorz_Case1(jequaloneDist_e,AgeWeightParamNames,Policy_e,n_d,n_a,n_z,N_j,pi_z,Params_b,simoptions_e);
AllFull_e=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_e,Policy_e,FnsToEvaluate_e,Params_b,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_e_full);
LCPFull_e=LifeCycleProfiles_FHorz_Case1(StationaryDist_e,Policy_e,FnsToEvaluate_e,Params_b,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_e_full);
TMb=struct();
t=gather(LCPFull_e.earnings.StdDeviation); t(Jr:end)=NaN; TMb.AgeConditionalStats.earnings.StdDeviation=t;
TMb.AllStats.earnings.Gini=gather(AllFull_e.earnings.Gini);
vfoptions_c=vfoptions; vfoptions_c.n_e=n_e; vfoptions_c.EiidShockFn=@(sigma_e) CalibFHorz_EiidShockFn(sigma_e);
simoptions_c=simoptions; simoptions_c.n_e=n_e; simoptions_c.EiidShockFn=vfoptions_c.EiidShockFn;
ran=0; objv=NaN; maxd=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'sigma_e'},TMb,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,[],FnsToEvaluate_e,estimoptions_a,vfoptions_c,simoptions_c);
    ran=1; objv=estsummary.objectivefnval;
    dc=estsummary.doublechecks.epsilon10000;
    mu=EstimFHorz_ModelMoments(dc.estimparamsvecup(:),{'sigma_e'},TMb,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,FnsToEvaluate_e,vfoptions_c,simoptions_c,0,[]);
    md=EstimFHorz_ModelMoments(dc.estimparamsvecdown(:),{'sigma_e'},TMb,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,FnsToEvaluate_e,vfoptions_c,simoptions_c,0,[]);
    Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
    maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
catch err
    fprintf('ShocksCustomStats (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats (b) EiidShockFn, sigma_e, skipestimation at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J_centered vs finite differences with the grids rebuilt, this should be zero: %.3e \n',1-ran,objv,maxd)
Params_start=Params_b; Params_start.sigma_e=0.15;
estimoptions_b2=estimoptions; estimoptions_b2.fminalgo=8;
ran=0; perr=NaN; tt=NaN; objv=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'sigma_e'},TMb,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_start,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,[],FnsToEvaluate_e,estimoptions_b2,vfoptions_c,simoptions_c);
    tt=toc(t0); ran=1; perr=abs(EstimParams.sigma_e-0.1); objv=estsummary.objectivefnval;
catch err
    fprintf('ShocksCustomStats (b) from 0.15: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats (b) EiidShockFn, sigma_e started at 0.15 (truth 0.1) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; sigma_e error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)

%% (c) jequaloneDist as a function, both usergrids settings; (d) CustomModelStats: both at the truth, beta and pension
TMc=struct();
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TMc.AgeConditionalStats.assets.Mean=t;
TMc.AllStats.assets.Gini=gather(AllFull.assets.Gini);
nTc=sum(~isnan(t))+1; Wc=eye(nTc); Omegac=diag(0.01*(1+(1:nTc)'/nTc));
EstimParamNames={'beta','pension'};
estimoptions_s=estimoptions; estimoptions_s.skipestimation=1;
estimoptions_s.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
[~,~,estsummary_mat]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMc,Wc,Omegac,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_s,vfoptions,simoptions);
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z) CalibFHorz_jequaloneDistFn(a_grid,z_grid,n_a,n_z);
for ug=[1,0]
    simoptions_ug=simoptions; simoptions_ug.jequaloneDist_usergrids=ug;
    simoptions_ug.a_grid=a_grid; simoptions_ug.z_grid=z_grid; % (a functional jequaloneDist needs the grids in simoptions, as the agent-distribution command says)
    ran=0; objv=NaN; dJ=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMc,Wc,Omegac,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_s,vfoptions,simoptions_ug);
        ran=1; objv=estsummary.objectivefnval; dJ=CoreMaxAbs(estsummary.variousmatrices.J(:)-estsummary_mat.variousmatrices.J(:));
    catch err
        fprintf('ShocksCustomStats (c) jequaloneDist_usergrids=%i: the estimation command errored with: %s \n',ug,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats (c) jequaloneDist as a function, jequaloneDist_usergrids=%i, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J vs the matrix jequaloneDist J, this should be zero: %.3e \n',ug,1-ran,objv,dJ)
end
TMd=struct();
TMd.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TMd.CustomModelStats.meanassets=gather(AllFull.assets.Mean);
TMd.CustomModelStats.zgridcheck=0;
nTd=3; Wd=eye(nTd); Omegad=diag(0.01*(1+(1:nTd)'/nTd));
[z_gridvals_J_expected,~,~]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0); % the internal joint grid
for ug=[1,0]
    estimoptions_d=estimoptions_s;
    estimoptions_d.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibFHorz_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
    estimoptions_d.CustomModelStats_usergrids=ug;
    if ug==1
        estimoptions_d.CalibFHorz_expected_zgrid=z_grid;
    else
        estimoptions_d.CalibFHorz_expected_zgrid=z_gridvals_J_expected;
    end
    ran=0; objv=NaN; dmom=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMd,Wd,Omegad,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_d,vfoptions,simoptions);
        ran=1; objv=estsummary.objectivefnval; dmom=CoreMaxAbs(estsummary.currentmomentvec(:)-estsummary.targetmomentvec(:));
    catch err
        fprintf('ShocksCustomStats (d) CustomModelStats_usergrids=%i: the estimation command errored with: %s \n',ug,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats (d) CustomModelStats, CustomModelStats_usergrids=%i, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (by-hand mean of assets vs AllStats, the z_grid check), this should be zero: %.3e \n',ug,1-ran,objv,dmom)
end

output=struct();
output.ran=ran;

end
