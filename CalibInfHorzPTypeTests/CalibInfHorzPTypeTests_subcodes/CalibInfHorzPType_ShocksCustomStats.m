function output=CalibInfHorzPType_ShocksCustomStats(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn)
% CalibInfHorzPTypeTests part 5: input forms of CalibrateInfHorzAgentModel_PType that change what the objective rebuilds each evaluation.
%
%   (a) vfoptions.ExogShockFn common to the types, with rho_z calibrated (the grids are rebuilt every evaluation): at the truth the
%       objective is zero and rho_z does not move (exact); started at 0.8 it is recovered (close). Targets: the grouped autocorrelation
%       of assets, the low type's autocorrelation of earnings, the grouped Gini of earnings.
%   (a2) vfoptions.ExogShockFn as a structure over the types (.low with rho_z_low, .high with rho_z_high), rho_z_low calibrated: the
%       calibrated-shocks detection has to look inside the structure (the unextended driver reads the structure as one function). At the
%       truth zero; from 0.8 recovered. Targets: the low type's autocorrelations of assets and earnings.
%   (b) CustomModelStats (PType signature, with Names_i): a by-hand population mean of assets (target: the grouped AllStats mean at
%       the truth) and a check that the z_grid handed to the custom function is the user's grid (CustomModelStats_usergrids=1) or the
%       internal form (=0); target 0. At the truth the objective is zero in both settings. The internal form is what
%       ExogShockSetup_InfHorz_PType returns; the user's grid in the PType form is its user_z_grid output with KeepOriginalGrid=1.
% Each run is wrapped so an error in the plumbing is reported as a failed check.

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;

vfoptions=struct();
simoptions=struct();
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);

%% Targets at the truth from the driver's grids
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_full);
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);

caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;
output=struct();

%% (a) ExogShockFn common to the types, rho_z calibrated
TMa=struct();
TMa.AutoCorrTransProbs.assets.AutoCorrelation=gather(ACFull.assets.AutoCorrelation);
TMa.AutoCorrTransProbs.earnings.low.AutoCorrelation=gather(ACFull.earnings.low.AutoCorrelation);
TMa.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
fprintf('ShocksCustomStats PType (a) targets at the truth: grouped autocorrelation of assets %.4f, low type autocorrelation of earnings %.4f, grouped Gini of earnings %.4f \n',TMa.AutoCorrTransProbs.assets.AutoCorrelation,TMa.AutoCorrTransProbs.earnings.low.AutoCorrelation,TMa.AllStats.earnings.Gini)
Params_a=Params; Params_a.rho_z=0.9; Params_a.sigma_z=0.1;
vfoptions_a=vfoptions; vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibInfHorzPType_ExogShockFn(rho_z,sigma_z);
simoptions_a=simoptions; simoptions_a.ExogShockFn=vfoptions_a.ExogShockFn;
for ss=1:2
    Params_start=Params_a; lab='rho_z started at the truth';
    if ss==2
        Params_start.rho_z=0.8; lab='rho_z started at 0.8 (truth 0.9)';
    end
    ran=0; objv=NaN; perr=NaN; dmom=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType({'rho_z'},TMa,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a,simoptions_a);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.rho_z-0.9); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
    catch err
        fprintf('ShocksCustomStats PType (a) common ExogShockFn, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('ShocksCustomStats PType (a) common ExogShockFn, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e; current vs target moments, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr,dmom)
    else
        fprintf('ShocksCustomStats PType (a) common ExogShockFn, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
    end
end
output.TMa=TMa;

%% (a2) ExogShockFn as a structure over the types, rho_z_low calibrated
Params_a2=Params; Params_a2.rho_z_low=0.9; Params_a2.rho_z_high=0.9; Params_a2.sigma_z=0.1;
vfoptions_a2=vfoptions;
vfoptions_a2.ExogShockFn=struct();
vfoptions_a2.ExogShockFn.low=@(rho_z_low,sigma_z) CalibInfHorzPType_ExogShockFn(rho_z_low,sigma_z);
vfoptions_a2.ExogShockFn.high=@(rho_z_high,sigma_z) CalibInfHorzPType_ExogShockFn(rho_z_high,sigma_z);
simoptions_a2=simoptions; simoptions_a2.ExogShockFn=vfoptions_a2.ExogShockFn;
TMa2=struct();
TMa2.AutoCorrTransProbs.assets.low.AutoCorrelation=gather(ACFull.assets.low.AutoCorrelation);
TMa2.AutoCorrTransProbs.earnings.low.AutoCorrelation=gather(ACFull.earnings.low.AutoCorrelation);
for ss=1:2
    Params_start=Params_a2; lab='rho_z_low started at the truth';
    if ss==2
        Params_start.rho_z_low=0.8; lab='rho_z_low started at 0.8 (truth 0.9)';
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType({'rho_z_low'},TMa2,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a2,simoptions_a2);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.rho_z_low-0.9);
    catch err
        fprintf('ShocksCustomStats PType (a2) per-type ExogShockFn, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('ShocksCustomStats PType (a2) per-type ExogShockFn (structure), %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z_low moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
    else
        fprintf('ShocksCustomStats PType (a2) per-type ExogShockFn (structure), %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z_low error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
    end
end

%% (b) CustomModelStats under both usergrids settings
TMb=struct();
TMb.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TMb.CustomModelStats.meanassets=gather(AllFull.assets.Mean);
TMb.CustomModelStats.zgridcheck=0;
[z_gridvals_expected,~,~]=ExogShockSetup_InfHorz_PType(n_z,z_grid,pi_z,Names_i,Params,vfoptions,3,0); % the internal form (usergrids=0)
try
    [~,~,vfoptions_ug1]=ExogShockSetup_InfHorz_PType(n_z,z_grid,pi_z,Names_i,Params,vfoptions,3,1); % the user's own grid in the PType form (usergrids=1)
    user_expected=vfoptions_ug1.user_z_grid;
catch
    user_expected=z_grid; % (if the PType setup does not report the user's grid, the plain z_grid is what the custom function should receive)
end
for ug=[1,0]
    caliboptions_b=caliboptions;
    caliboptions_b.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibInfHorzPType_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
    caliboptions_b.CustomModelStats_usergrids=ug;
    if ug==1
        caliboptions_b.CalibInfHorz_expected_zgrid=user_expected;
    else
        caliboptions_b.CalibInfHorz_expected_zgrid=z_gridvals_expected;
    end
    ran=0; objv=NaN; moved=NaN; dmom=NaN;
    try
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType({'beta'},TMb,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_b,vfoptions,simoptions);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
    catch err
        fprintf('ShocksCustomStats PType (b) CustomModelStats, CustomModelStats_usergrids=%i: the calibration errored with: %s \n',ug,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats PType (b) CustomModelStats, CustomModelStats_usergrids=%i, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (by-hand population mean of assets vs grouped AllStats, z_grid check), this should be zero: %.3e; beta moved, this should be zero: %.3e \n',ug,1-ran,objv,dmom,moved)
end

end
