function output=CalibInfHorz_ShocksCustomStats(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn)
% CalibInfHorzTests part 5: two input forms of CalibrateInfHorzAgentModel that change what the objective has to rebuild each evaluation.
%
%   (a) vfoptions.ExogShockFn with the autocorrelation rho_z of the shock a calibrated parameter (the grids are rebuilt every
%       evaluation, caliboptions.calibrateshocks=1). At rho_z=0.9, sigma_z=0.1 the function reproduces the driver's z_grid and pi_z,
%       so targets generated with the driver's grids are the model's own: started at the truth the objective is zero and rho_z does
%       not move (exact); started at rho_z=0.8 it is recovered (close). Targets: the autocorrelation of assets and of earnings
%       (AutoCorrTransProbs, both depend on rho_z) and the Gini of earnings (AllStats).
%   (b) CustomModelStats: a by-hand mean of assets (target: the AllStats mean at the truth) and a check that the z_grid handed to the
%       custom function is the user's grid (CustomModelStats_usergrids=1) or the internal joint grid (=0); target 0. At the truth the
%       objective is zero in both settings, and nonzero if the wrong grid form arrives.
%   (c) CustomModelStats_usergrids=1 together with ExogShockFn (rho_z calibrated): the grid handed to the custom function must be the
%       user's own grid as rebuilt from the current rho_z (vfoptions.user_z_grid), at the truth the driver's z_grid; objective zero.
% (There is no jequaloneDist in the infinite-horizon model, so the FHorz part's function-form test has no counterpart here.)
% Each run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;

vfoptions=struct();
simoptions=struct();
vfoptions.gridinterplayer=1; % (as in part 4: the moments then vary smoothly with the parameters)
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;

%% Targets at the truth from the driver's grids
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_full);
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions);

caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;
output=struct();

%% (a) ExogShockFn with rho_z calibrated
TMa=struct();
TMa.AutoCorrTransProbs.assets.AutoCorrelation=gather(ACFull.assets.AutoCorrelation);
TMa.AutoCorrTransProbs.earnings.AutoCorrelation=gather(ACFull.earnings.AutoCorrelation);
TMa.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
fprintf('ShocksCustomStats (a) targets at the truth: autocorrelation of assets %.4f and of earnings %.4f, Gini of earnings %.4f \n',TMa.AutoCorrTransProbs.assets.AutoCorrelation,TMa.AutoCorrTransProbs.earnings.AutoCorrelation,TMa.AllStats.earnings.Gini)
Params_a=Params;
Params_a.rho_z=0.9;
Params_a.sigma_z=0.1;
vfoptions_a=vfoptions;
vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibInfHorz_ExogShockFn(rho_z,sigma_z);
simoptions_a=simoptions;
simoptions_a.ExogShockFn=vfoptions_a.ExogShockFn;
for ss=1:2
    Params_start=Params_a; lab='rho_z started at the truth';
    if ss==2
        Params_start.rho_z=0.8; lab='rho_z started at 0.8 (truth 0.9)';
    end
    ran=0; objv=NaN; perr=NaN; dmom=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'rho_z'},TMa,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a,simoptions_a);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.rho_z-0.9); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
    catch err
        fprintf('ShocksCustomStats (a) ExogShockFn, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('ShocksCustomStats (a) ExogShockFn, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e; current vs target moments, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr,dmom)
    else
        fprintf('ShocksCustomStats (a) ExogShockFn, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
    end
end
output.TMa=TMa;

%% (b) CustomModelStats under both usergrids settings
TMb=struct();
TMb.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TMb.CustomModelStats.meanassets=gather(AllFull.assets.Mean);
TMb.CustomModelStats.zgridcheck=0;
[z_gridvals_expected,~,~]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0); % the internal joint grid
for ug=[1,0]
    caliboptions_b=caliboptions;
    caliboptions_b.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibInfHorz_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
    caliboptions_b.CustomModelStats_usergrids=ug;
    if ug==1
        caliboptions_b.CalibInfHorz_expected_zgrid=z_grid;
    else
        caliboptions_b.CalibInfHorz_expected_zgrid=z_gridvals_expected;
    end
    ran=0; objv=NaN; moved=NaN; dmom=NaN;
    try
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'beta'},TMb,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_b,vfoptions,simoptions);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
    catch err
        fprintf('ShocksCustomStats (b) CustomModelStats, CustomModelStats_usergrids=%i: the calibration errored with: %s \n',ug,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ShocksCustomStats (b) CustomModelStats, CustomModelStats_usergrids=%i, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (by-hand mean of assets vs AllStats, z_grid check), this should be zero: %.3e; beta moved, this should be zero: %.3e \n',ug,1-ran,objv,dmom,moved)
end

%% (c) CustomModelStats_usergrids=1 together with ExogShockFn: the custom function gets the user's grid as rebuilt from the current rho_z
caliboptions_c=caliboptions;
caliboptions_c.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibInfHorz_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
caliboptions_c.CustomModelStats_usergrids=1;
caliboptions_c.CalibInfHorz_expected_zgrid=z_grid; % at the true rho_z the rebuilt user grid is the driver's
TMc=TMb;
TMc.AutoCorrTransProbs.assets.AutoCorrelation=gather(ACFull.assets.AutoCorrelation);
ran=0; objv=NaN; moved=NaN; dmom=NaN;
try
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'rho_z'},TMc,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_c,vfoptions_a,simoptions_a);
    ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.rho_z-0.9); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
catch err
    fprintf('ShocksCustomStats (c) CustomModelStats_usergrids=1 with ExogShockFn: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ShocksCustomStats (c) CustomModelStats_usergrids=1 with ExogShockFn, rho_z at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (incl. the rebuilt user grid handed to the custom function), this should be zero: %.3e; rho_z moved, this should be zero: %.3e \n',1-ran,objv,dmom,moved)

end
