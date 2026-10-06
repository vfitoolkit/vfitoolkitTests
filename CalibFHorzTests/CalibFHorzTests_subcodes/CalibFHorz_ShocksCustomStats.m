function output=CalibFHorz_ShocksCustomStats(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 5: three input forms of CalibrateLifeCycleModel that change what the objective has to rebuild each evaluation.
%
%   (a) vfoptions.ExogShockFn with the autocorrelation rho_z of the shock a calibrated parameter (the grids are rebuilt every
%       evaluation, caliboptions.calibrateshocks=1). At rho_z=0.9, sigma_z=0.03 the function reproduces the driver's z_grid and pi_z,
%       so targets generated with the driver's grids are the model's own: started at the truth the objective is zero and rho_z does
%       not move (exact); started at rho_z=0.8 it is recovered (close). Targets: the age-conditional standard deviation of earnings
%       while working and the Gini of earnings, both of which depend on rho_z.
%   (b) jequaloneDist given as a function, with simoptions.jequaloneDist_usergrids=1 (the user's z_grid) and =0 (the internal joint
%       grid): started at the truth with the matrix-version targets, the objective is zero. Each run is wrapped so that an error in
%       the plumbing is reported as a failed check rather than ending the bank.
%   (c) CustomModelStats: a by-hand mean of assets (target: the AllStats mean at the truth) and a check that the z_grid handed to the
%       custom function is the user's grid (CustomModelStats_usergrids=1) or the internal joint grid (=0); target 0. At the truth the
%       objective is zero in both settings, and nonzero if the wrong grid form arrives.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
Jr=Params.Jr;

vfoptions=struct();
simoptions=struct();
% The grid interpolation layer makes the policy, and so every moment, vary smoothly with the parameters. Without it the moments
% are piecewise constant in beta on this coarse grid (run of 2026-10-07: lsqnonlin stalled on its step-size tolerance with the
% objective still at 0.13 in the constrained cases, and fminsearch found an exact zero away from the truth). The targets are
% generated with the same options, so the truth is still exact.
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;

%% Targets at the truth from the driver's grids
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);

caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;

%% (a) ExogShockFn with rho_z calibrated
TMa=struct();
t=gather(LCPFull.earnings.StdDeviation); t(Jr:end)=NaN; TMa.AgeConditionalStats.earnings.StdDeviation=t;
TMa.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
Params_a=Params;
Params_a.rho_z=0.9;
Params_a.sigma_z=0.03;
vfoptions_a=vfoptions;
vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibFHorz_ExogShockFn(rho_z,sigma_z);
tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel({'rho_z'},TMa,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a,simoptions);
tt=toc;
fprintf('ShocksCustomStats (a) ExogShockFn, rho_z started at the truth (%.1f s): objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e; current vs target moments, this should be zero: %.3e \n',tt,calibsummary.objvalue,abs(CalibParams.rho_z-0.9),max(abs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:))))
Params_a.rho_z=0.8;
tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel({'rho_z'},TMa,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions,vfoptions_a,simoptions);
tt=toc;
fprintf('ShocksCustomStats (a) ExogShockFn, rho_z started at 0.8 (%.1f s): recovered %.6f (truth 0.9), objective %.3e; error, this should be close to zero: %.3e \n',tt,CalibParams.rho_z,calibsummary.objvalue,abs(CalibParams.rho_z-0.9))

%% (b) jequaloneDist as a function
TMb=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; TMb.AgeConditionalStats.assets.Mean=t;
TMb.AllStats.assets.Gini=gather(AllFull.assets.Gini);
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z) CalibFHorz_jequaloneDistFn(a_grid,z_grid,n_a,n_z);
for ug=[1,0]
    simoptions_b=simoptions;
    simoptions_b.a_grid=a_grid;
    simoptions_b.z_grid=z_grid;
    simoptions_b.jequaloneDist_usergrids=ug;
    ran=0; objv=NaN; moved=NaN;
    try
        [CalibParams,calibsummary]=CalibrateLifeCycleModel({'beta'},TMb,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,[],FnsToEvaluate,caliboptions,vfoptions,simoptions_b);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta);
    catch err
        fprintf('ShocksCustomStats (b) jequaloneDist as a function, jequaloneDist_usergrids=%i: the calibration errored with: %s \n',ug,err.message)
    end
    fprintf('ShocksCustomStats (b) jequaloneDist as a function, jequaloneDist_usergrids=%i, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',ug,1-ran,objv,moved)
end

%% (c) CustomModelStats under both usergrids settings
TMc=struct();
TMc.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TMc.CustomModelStats.meanassets=gather(AllFull.assets.Mean);
TMc.CustomModelStats.zgridcheck=0;
[z_gridvals_J_expected,~,~]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0); % the internal joint grid
for ug=[1,0]
    caliboptions_c=caliboptions;
    caliboptions_c.CustomModelStats=@(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibFHorz_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
    caliboptions_c.CustomModelStats_usergrids=ug;
    if ug==1
        caliboptions_c.CalibFHorz_expected_zgrid=z_grid;
    else
        caliboptions_c.CalibFHorz_expected_zgrid=z_gridvals_J_expected;
    end
    ran=0; objv=NaN; moved=NaN; dmom=NaN;
    try
        [CalibParams,calibsummary]=CalibrateLifeCycleModel({'beta'},TMc,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions_c,vfoptions,simoptions);
        ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta); dmom=max(abs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:)));
    catch err
        fprintf('ShocksCustomStats (c) CustomModelStats, CustomModelStats_usergrids=%i: the calibration errored with: %s \n',ug,err.message)
    end
    fprintf('ShocksCustomStats (c) CustomModelStats, CustomModelStats_usergrids=%i, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments (by-hand mean of assets vs AllStats, z_grid check), this should be zero: %.3e; beta moved, this should be zero: %.3e \n',ug,1-ran,objv,dmom,moved)
end

output=struct();
output.TMa=TMa;

end
