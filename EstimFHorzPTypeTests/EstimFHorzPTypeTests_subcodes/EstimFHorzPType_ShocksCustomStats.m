function output=EstimFHorzPType_ShocksCustomStats(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% EstimFHorzPTypeTests part 3: shocks that depend on estimated parameters, jequaloneDist as a function, CustomModelStats, through
% EstimateLifeCycleModel_PType_MethodOfMoments. The grid interpolation layer is on.
%   (a) vfoptions.ExogShockFn common to the types, rho_z estimated: at the truth (skipestimation=1) the objective is zero and J equals
%       finite differences the test takes itself with the grids rebuilt at each step (exact); started at 0.8 (constrain0to1), recovered (close)
%   (a2) vfoptions.ExogShockFn as a structure over the types (.low with rho_z_low, .high with rho_z_high), rho_z_low estimated: the same
%        two checks (targets: the low type's earnings profile)
%   (b) vfoptions.EiidShockFn with sigma_e estimated: the same two checks (truth 0.1, started at 0.15)
%   (c) jequaloneDist as a function under both jequaloneDist_usergrids settings: at the truth the objective is zero and J equals the
%       matrix jequaloneDist's J (exact)
%   (d) CustomModelStats under both CustomModelStats_usergrids settings: at the truth the objective is zero and the current moments
%       equal the targets (exact)
% The function files are the CalibFHorzPTypeTests ones (on the path). Each driver run is wrapped so an error in the plumbing is reported as a
% failed check.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
estimoptions=struct(); estimoptions.verbose=0;
estimoptions_s=estimoptions; estimoptions_s.skipestimation=1;
estimoptions_s.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')

%% (a) ExogShockFn common to the types, rho_z estimated; (a2) per type as a structure
Params_a=Params; Params_a.rho_z=0.9; Params_a.sigma_z=0.03;
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_a,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params_a,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params_a,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMa=struct();
t=gather(LCPFull.earnings.StdDeviation); t(Jr:end)=NaN; TMa.AgeConditionalStats.earnings.StdDeviation=t; % grouped
TMa.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
nT=(Jr-1)+1; W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
TMa2=struct();
t=gather(LCPFull.earnings.low.StdDeviation); t(Jr:end)=NaN; TMa2.AgeConditionalStats.earnings.low.StdDeviation=t; % the low type
TMa2.AllStats.earnings.low.Gini=gather(AllFull.earnings.low.Gini);
for form=1:2
    if form==1
        vfoptions_a=vfoptions; vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibFHorzPType_ExogShockFn(rho_z,sigma_z);
        pname='rho_z'; Params_t=Params_a; TMt=TMa; lab='common ExogShockFn, rho_z';
    else
        vfoptions_a=vfoptions; vfoptions_a.ExogShockFn=struct();
        vfoptions_a.ExogShockFn.low=@(rho_z_low,sigma_z) CalibFHorzPType_ExogShockFn(rho_z_low,sigma_z);
        vfoptions_a.ExogShockFn.high=@(rho_z_high,sigma_z) CalibFHorzPType_ExogShockFn(rho_z_high,sigma_z);
        Params_t=Params_a; Params_t.rho_z_low=0.9; Params_t.rho_z_high=0.9; pname='rho_z_low'; TMt=TMa2; lab='per-type ExogShockFn (structure), rho_z_low';
    end
    simoptions_a=simoptions; simoptions_a.ExogShockFn=vfoptions_a.ExogShockFn;
    ran=0; objv=NaN; maxd=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments({pname},TMt,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_t,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_s,vfoptions_a,simoptions_a);
        ran=1; objv=estsummary.objectivefnval;
        dc=estsummary.doublechecks.epsilon10000;
        mu=EstimFHorzPType_ModelMoments(dc.estimparamsvecup(:),{pname},TMt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_t,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions_a,simoptions_a,0,[]);
        md=EstimFHorzPType_ModelMoments(dc.estimparamsvecdown(:),{pname},TMt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_t,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions_a,simoptions_a,0,[]);
        Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
        maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
    catch err
        fprintf('ShocksCustomStats PType (a) %s: the estimation command errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats PType (a) %s, skipestimation at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J_centered vs finite differences with the grids rebuilt, this should be zero: %.3e \n',lab,1-ran,objv,maxd)
    Params_start=Params_t; Params_start.(pname)=0.8;
    estimoptions_a2=estimoptions; estimoptions_a2.fminalgo=8; estimoptions_a2.constrain0to1={pname};
    ran=0; perr=NaN; tt=NaN; objv=NaN;
    try
        t0=tic;
        [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments({pname},TMt,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_a2,vfoptions_a,simoptions_a);
        tt=toc(t0); ran=1; perr=abs(EstimParams.(pname)-0.9); objv=estsummary.objectivefnval;
    catch err
        fprintf('ShocksCustomStats PType (a) %s from 0.8: the estimation command errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats PType (a) %s started at 0.8 (truth 0.9) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; parameter error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
end

%% (b) EiidShockFn with sigma_e estimated
ReturnFn_e=@(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension);
FnsToEvaluate_e.assets=@(aprime,a,z,e) a;
FnsToEvaluate_e.earnings=@(aprime,a,z,e,w,kappa_j_pt) w*kappa_j_pt*z*e;
n_e=3;
Params_b=Params_a; Params_b.sigma_e=0.1;
[e_grid,pi_e]=CalibFHorzPType_EiidShockFn(Params_b.sigma_e);
vfoptions_e=vfoptions; vfoptions_e.n_e=n_e; vfoptions_e.e_grid=e_grid; vfoptions_e.pi_e=pi_e;
simoptions_e=simoptions; simoptions_e.n_e=n_e; simoptions_e.e_grid=e_grid; simoptions_e.pi_e=pi_e;
simoptions_e_full=simoptions_e; simoptions_e_full.whichstats=ones(1,7);
jequaloneDist_e=zeros(n_a,n_z,n_e,'gpuArray'); jequaloneDist_e(1,ceil(n_z/2),ceil(n_e/2))=1;
[~,Policy_e]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,vfoptions_e);
StationaryDist_e=StationaryDist_Case1_FHorz_PType(jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,Policy_e,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_b,simoptions_e);
AllFull_e=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_e,Policy_e,FnsToEvaluate_e,Params_b,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_e_full);
LCPFull_e=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_e,Policy_e,FnsToEvaluate_e,Params_b,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_e_full);
TMb=struct();
t=gather(LCPFull_e.earnings.StdDeviation); t(Jr:end)=NaN; TMb.AgeConditionalStats.earnings.StdDeviation=t;
TMb.AllStats.earnings.Gini=gather(AllFull_e.earnings.Gini);
vfoptions_c=vfoptions; vfoptions_c.n_e=n_e; vfoptions_c.EiidShockFn=@(sigma_e) CalibFHorzPType_EiidShockFn(sigma_e);
simoptions_c=simoptions; simoptions_c.n_e=n_e; simoptions_c.EiidShockFn=vfoptions_c.EiidShockFn;
ran=0; objv=NaN; maxd=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments({'sigma_e'},TMb,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_e,estimoptions_s,vfoptions_c,simoptions_c);
    ran=1; objv=estsummary.objectivefnval;
    dc=estsummary.doublechecks.epsilon10000;
    mu=EstimFHorzPType_ModelMoments(dc.estimparamsvecup(:),{'sigma_e'},TMb,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate_e,vfoptions_c,simoptions_c,0,[]);
    md=EstimFHorzPType_ModelMoments(dc.estimparamsvecdown(:),{'sigma_e'},TMb,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_b,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate_e,vfoptions_c,simoptions_c,0,[]);
    Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
    maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
catch err
    fprintf('ShocksCustomStats PType (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats PType (b) EiidShockFn, sigma_e, skipestimation at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J_centered vs finite differences with the grids rebuilt, this should be zero: %.3e \n',1-ran,objv,maxd)
Params_start=Params_b; Params_start.sigma_e=0.15;
estimoptions_b2=estimoptions; estimoptions_b2.fminalgo=8;
ran=0; perr=NaN; tt=NaN; objv=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments({'sigma_e'},TMb,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_start,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_e,estimoptions_b2,vfoptions_c,simoptions_c);
    tt=toc(t0); ran=1; perr=abs(EstimParams.sigma_e-0.1); objv=estsummary.objectivefnval;
catch err
    fprintf('ShocksCustomStats PType (b) from 0.15: the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats PType (b) EiidShockFn, sigma_e started at 0.15 (truth 0.1) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; sigma_e error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)

%% (c) jequaloneDist as a function; (d) CustomModelStats: at the truth with beta and pension
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMc=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TMc.AgeConditionalStats.assets.Mean=t;
TMc.AllStats.assets.Gini=gather(AllFull.assets.Gini);
nTc=sum(~isnan(t))+1; Wc=eye(nTc); Omegac=diag(0.01*(1+(1:nTc)'/nTc));
EstimParamNames={'beta','pension'};
[~,~,estsummary_mat]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMc,Wc,Omegac,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_s,vfoptions,simoptions);
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z) CalibFHorzPType_jequaloneDistFn(a_grid,z_grid,n_a,n_z);
for ug=[1,0]
    simoptions_ug=simoptions; simoptions_ug.jequaloneDist_usergrids=ug;
    simoptions_ug.a_grid=a_grid; simoptions_ug.z_grid=z_grid; % (a functional jequaloneDist needs the grids in simoptions, as the agent-distribution command says)
    ran=0; objv=NaN; dJ=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMc,Wc,Omegac,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_s,vfoptions,simoptions_ug);
        ran=1; objv=estsummary.objectivefnval; dJ=CoreMaxAbs(estsummary.variousmatrices.J(:)-estsummary_mat.variousmatrices.J(:));
    catch err
        fprintf('ShocksCustomStats PType (c) jequaloneDist_usergrids=%i: the estimation command errored with: %s \n',ug,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats PType (c) jequaloneDist as a function, jequaloneDist_usergrids=%i, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J vs the matrix jequaloneDist J, this should be zero: %.3e \n',ug,1-ran,objv,dJ)
end
TMd=struct();
TMd.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TMd.CustomModelStats.meanassets=gather(AllFull.assets.Mean);
TMd.CustomModelStats.zgridcheck=0;
nTd=3; Wd=eye(nTd); Omegad=diag(0.01*(1+(1:nTd)'/nTd));
[z_gridvals_J_expected,~,~]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0); % the internal form (usergrids=0)
[~,~,vfoptions_ug1]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,1); % the user's own grid in the PType form (usergrids=1)
for ug=[1,0]
    estimoptions_d=estimoptions_s;
    estimoptions_d.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibFHorzPType_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
    estimoptions_d.CustomModelStats_usergrids=ug;
    if ug==1
        estimoptions_d.CalibFHorz_expected_zgrid=vfoptions_ug1.user_z_grid;
    else
        estimoptions_d.CalibFHorz_expected_zgrid=z_gridvals_J_expected;
    end
    ran=0; objv=NaN; dmom=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMd,Wd,Omegad,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_d,vfoptions,simoptions);
        ran=1; objv=estsummary.objectivefnval; dmom=CoreMaxAbs(estsummary.currentmomentvec(:)-estsummary.targetmomentvec(:));
    catch err
        fprintf('ShocksCustomStats PType (d) CustomModelStats_usergrids=%i: the estimation command errored with: %s \n',ug,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats PType (d) CustomModelStats, CustomModelStats_usergrids=%i, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (by-hand population mean of assets vs grouped AllStats, the z_grid check), this should be zero: %.3e \n',ug,1-ran,objv,dmom)
end

output=struct();
output.ran=ran;

end
