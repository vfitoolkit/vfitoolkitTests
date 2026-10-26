function output=CalibFHorzPType_ShocksCustomStats(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 5: input forms of CalibrateLifeCycleModel_PType that change what the objective rebuilds each evaluation.
%
%   (a) vfoptions.ExogShockFn common to the types, with rho_z calibrated (the grids are rebuilt every evaluation): at the truth the
%       objective is zero and rho_z does not move (exact); started at 0.8 it is recovered (close). Targets: the per-type age-conditional
%       standard deviation of earnings while working and the grouped Gini of earnings.
%   (a2) vfoptions.ExogShockFn as a structure over the types (.low with rho_z_low, .high with rho_z_high), rho_z_low calibrated: the
%       calibrated-shocks detection has to look inside the structure (it errored before 2026-10-07). At the truth zero; from 0.8 recovered.
%   (b) jequaloneDist given as a function, with simoptions.jequaloneDist_usergrids=1 and =0: at the truth the objective is zero.
%   (c) CustomModelStats (PType signature, with Names_i): a by-hand population mean of assets (target: the grouped AllStats mean at the
%       truth) and a check that the z_grid handed to the custom function is the user's grid (CustomModelStats_usergrids=1) or the
%       internal form (=0); target 0. At the truth the objective is zero in both settings. With permanent types the user's grid
%       arrives as a structure over the types (each the z_grid as input), which is what ExogShockSetup_FHorz_PType returns as user_z_grid.
% The runs of (b) and (c) are wrapped so that an error in the plumbing is reported as a failed check rather than ending the bank.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
Jr=Params.Jr;

vfoptions=struct();
simoptions=struct();
% The grid interpolation layer makes the moments vary smoothly with the parameters; the targets are generated with the same options
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);

%% Targets at the truth from the driver's grids
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);

caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;

%% (a) ExogShockFn common to the types, rho_z calibrated
TMa=struct();
t=gather(LCPFull.earnings.low.StdDeviation); t(Jr:end)=NaN; TMa.AgeConditionalStats.earnings.low.StdDeviation=t;
t=gather(LCPFull.earnings.high.StdDeviation); t(Jr:end)=NaN; TMa.AgeConditionalStats.earnings.high.StdDeviation=t;
TMa.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
Params_a=Params;
Params_a.rho_z=0.9;
Params_a.sigma_z=0.03;
vfoptions_a=vfoptions;
vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibFHorzPType_ExogShockFn(rho_z,sigma_z);
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'rho_z'},TMa,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a,simoptions);
tt=toc(t0);
fprintf('ShocksCustomStats PType (a) common ExogShockFn, rho_z started at the truth (%.1f s): objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e; current vs target moments, this should be zero: %.3e \n',tt,calibsummary.objvalue,abs(CalibParams.rho_z-0.9),max(abs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:))))
Params_a.rho_z=0.8;
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'rho_z'},TMa,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a,simoptions);
tt=toc(t0);
fprintf('ShocksCustomStats PType (a) common ExogShockFn, rho_z started at 0.8 (%.1f s): recovered %.6f (truth 0.9), objective %.3e; error, this should be close to zero: %.3e \n',tt,CalibParams.rho_z,calibsummary.objvalue,abs(CalibParams.rho_z-0.9))

%% (a2) ExogShockFn as a structure over the types, rho_z_low calibrated
Params_a2=Params;
Params_a2.rho_z_low=0.9; Params_a2.rho_z_high=0.9; Params_a2.sigma_z=0.03;
vfoptions_a2=vfoptions;
vfoptions_a2.ExogShockFn=struct();
vfoptions_a2.ExogShockFn.low=@(rho_z_low,sigma_z) CalibFHorzPType_ExogShockFn(rho_z_low,sigma_z);
vfoptions_a2.ExogShockFn.high=@(rho_z_high,sigma_z) CalibFHorzPType_ExogShockFn(rho_z_high,sigma_z);
TMa2=struct();
t=gather(LCPFull.earnings.low.StdDeviation); t(Jr:end)=NaN; TMa2.AgeConditionalStats.earnings.low.StdDeviation=t;
for ss=1:2
    if ss==2
        Params_a2.rho_z_low=0.8;
    end
    ran=0; objv=NaN; moved=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'rho_z_low'},TMa2,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a2,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a2,simoptions);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.rho_z_low-0.9);
    catch err
        fprintf('ShocksCustomStats PType (a2) per-type ExogShockFn: the calibration errored with: %s \n',err.message)
    end
    if ss==1
        fprintf('ShocksCustomStats PType (a2) per-type ExogShockFn (structure), rho_z_low started at the truth (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z_low moved, this should be zero: %.3e \n',tt,1-ran,objv,moved)
    else
        fprintf('ShocksCustomStats PType (a2) per-type ExogShockFn (structure), rho_z_low started at 0.8 (%.1f s), objective %.3e: did not run, this should be zero: %.3e; error, this should be close to zero: %.3e \n',tt,objv,1-ran,moved)
    end
end

%% (b) jequaloneDist as a function
TMb=struct();
t=gather(LCPFull.assets.low.Mean); t(1)=NaN; TMb.AgeConditionalStats.assets.low.Mean=t;
t=gather(LCPFull.assets.high.Mean); t(1)=NaN; TMb.AgeConditionalStats.assets.high.Mean=t;
TMb.AllStats.assets.Gini=gather(AllFull.assets.Gini);
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z) CalibFHorzPType_jequaloneDistFn(a_grid,z_grid,n_a,n_z);
for ug=[1,0]
    simoptions_b=simoptions;
    simoptions_b.a_grid=a_grid;
    simoptions_b.z_grid=z_grid;
    simoptions_b.jequaloneDist_usergrids=ug;
    ran=0; objv=NaN; moved=NaN;
    try
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'beta'},TMb,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions,simoptions_b);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta);
    catch err
        fprintf('ShocksCustomStats PType (b) jequaloneDist as a function, jequaloneDist_usergrids=%i: the calibration errored with: %s \n',ug,err.message)
    end
    fprintf('ShocksCustomStats PType (b) jequaloneDist as a function, jequaloneDist_usergrids=%i, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',ug,1-ran,objv,moved)
end

%% (c) CustomModelStats under both usergrids settings
TMc=struct();
TMc.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TMc.CustomModelStats.meanassets=gather(AllFull.assets.Mean);
TMc.CustomModelStats.zgridcheck=0;
[z_gridvals_J_expected,~,~]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0); % the internal form (usergrids=0)
[~,~,vfoptions_ug1]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,1); % the user's own grid in the PType form (usergrids=1): a structure over the types, each the z_grid as input
for ug=[1,0]
    caliboptions_c=caliboptions;
    caliboptions_c.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibFHorzPType_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
    caliboptions_c.CustomModelStats_usergrids=ug;
    if ug==1
        caliboptions_c.CalibFHorz_expected_zgrid=vfoptions_ug1.user_z_grid;
    else
        caliboptions_c.CalibFHorz_expected_zgrid=z_gridvals_J_expected;
    end
    ran=0; objv=NaN; moved=NaN; dmom=NaN;
    try
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'beta'},TMc,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_c,vfoptions,simoptions);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta); dmom=max(abs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:)));
    catch err
        fprintf('ShocksCustomStats PType (c) CustomModelStats, CustomModelStats_usergrids=%i: the calibration errored with: %s \n',ug,err.message)
    end
    fprintf('ShocksCustomStats PType (c) CustomModelStats, CustomModelStats_usergrids=%i, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (by-hand population mean of assets vs grouped AllStats, z_grid check), this should be zero: %.3e; beta moved, this should be zero: %.3e \n',ug,1-ran,objv,dmom,moved)
end

output=struct();
output.TMa=TMa;

end
