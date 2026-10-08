function output=CalibInfHorz_RestrictedTargets(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn)
% CalibInfHorzTests part 6: round trips through CalibrateInfHorzAgentModel on the AutoCorrTransProbs and CrossSectionCovarCorr targets,
% restricted and with horizons and matrices, and caliboptions.logmoments by their three- and four-level names. (Their parsing, selectors
% and moment vector are covered in parts 1 and 2.) Grid interpolation layer as in part 4; targets generated at the truth.
%   (a) rho_z (ExogShockFn) from the autocorrelation of assets at horizons 1 and 2 among the rich (restricted, with a horizon suffix)
%       and the unrestricted autocorrelation of earnings: at the truth (exact), from 0.8 (close)
%   (b) beta from the assets-consumption correlation among the rich (a four-level cross-section name) and the unrestricted
%       assets-earnings covariance: at the truth (exact), from 0.92 (close); whichcombos=0 from the same start gives exactly the same
%       parameters, objective and moments as =1
%   (c) the unrestricted CorrelationMatrix (diagonal NaN) and the CovarianceMatrix among the poor (one entry NaN) as targets: at the
%       truth the objective is zero and the current moments equal the targets (the matrices vectorised with their omissions dropped)
%   (d) logmoments by name on the four-level cross-section name and the three-level restricted AutoCorr name, the other targets
%       unnamed: calibsummary.logmoments is the expected vector, the objective zero on the log targets
% Each run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w) (1+r)*a+w*z-aprime;

vfoptions=struct();
simoptions=struct();
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);

%% Targets at the truth
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
simoptions_ac=simoptions;
simoptions_ac.timehorizons=2;
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;
output=struct();

%% (a) rho_z from restricted autocorrelations with a horizon suffix
TMa=struct();
TMa.AutoCorrTransProbs.rich.assets.AutoCorrelation=gather(ACFull.rich.assets.AutoCorrelation);
TMa.AutoCorrTransProbs.rich.assets.AutoCorrelation_k2=gather(ACFull.rich.assets.tperiods2.AutoCorrelation);
TMa.AutoCorrTransProbs.earnings.AutoCorrelation=gather(ACFull.earnings.AutoCorrelation);
fprintf('RestrictedTargets (a) targets at the truth: autocorrelation of assets among the rich %.4f (horizon one) and %.4f (horizon two), of earnings %.4f \n',TMa.AutoCorrTransProbs.rich.assets.AutoCorrelation,TMa.AutoCorrTransProbs.rich.assets.AutoCorrelation_k2,TMa.AutoCorrTransProbs.earnings.AutoCorrelation)
Params_a=Params; Params_a.rho_z=0.9; Params_a.sigma_z=0.1;
vfoptions_a=vfoptions; vfoptions_a.ExogShockFn=@(rho_z,sigma_z) CalibInfHorz_ExogShockFn(rho_z,sigma_z);
simoptions_a=simoptions; simoptions_a.ExogShockFn=vfoptions_a.ExogShockFn;
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
        fprintf('RestrictedTargets (a) %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('RestrictedTargets (a) restricted AutoCorr targets with a horizon, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; rho_z moved, this should be zero: %.3e; current vs target moments, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr,dmom)
    else
        fprintf('RestrictedTargets (a) restricted AutoCorr targets with a horizon, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; rho_z error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
    end
end

%% (b) beta from a restricted cross-section correlation (four-level name) and an unrestricted covariance
TMb=struct();
TMb.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=gather(CSFull.rich.assets.CorrelationWith.consumption);
TMb.CrossSectionCovarCorr.assets.CovarianceWith.earnings=gather(CSFull.assets.CovarianceWith.earnings);
fprintf('RestrictedTargets (b) targets at the truth: assets-consumption correlation among the rich %.4f, assets-earnings covariance %.4f \n',TMb.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption,TMb.CrossSectionCovarCorr.assets.CovarianceWith.earnings)
res=cell(3,1);
for ss=1:3
    Params_start=Params; caliboptions_b=caliboptions; lab='beta started at the truth';
    if ss>=2
        Params_start.beta=0.92; lab='beta started at 0.92 (truth 0.95)';
    end
    if ss==3
        caliboptions_b.whichcombos=0; lab='beta started at 0.92, whichcombos=0';
    end
    ran=0; objv=NaN; perr=NaN; dmom=NaN; tt=NaN; rec=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'beta'},TMb,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_b,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; rec=CalibParams.beta; perr=abs(rec-Params.beta); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
        res{ss}=struct('beta',CalibParams.beta,'objv',objv,'cur',calibsummary.currentmomentvec(:));
    catch err
        fprintf('RestrictedTargets (b) %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('RestrictedTargets (b) cross-section targets, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e; current vs target moments, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr,dmom)
    else
        fprintf('RestrictedTargets (b) cross-section targets, %s (%.1f s), recovered beta %.6f, objective %.3e: did not run, this should be zero: %.3e; beta error, this should be close to zero: %.3e \n',lab,tt,rec,objv,1-ran,perr)
    end
end
if ~isempty(res{2}) && ~isempty(res{3})
    fprintf('RestrictedTargets (b) whichcombos=0 vs =1 from the same start: beta, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',abs(res{3}.beta-res{2}.beta),abs(res{3}.objv-res{2}.objv),CoreMaxAbs(res{3}.cur-res{2}.cur))
else
    fprintf('RestrictedTargets (b) whichcombos=0 vs =1: one of the two did not run, this should be zero: %.3e \n',1)
end

%% (c) matrix targets at the truth
TMc=struct();
t=gather(CSFull.CorrelationMatrix); t(logical(eye(3)))=NaN; TMc.CrossSectionCovarCorr.CorrelationMatrix=t;
t=gather(CSFull.poor.CovarianceMatrix); t(1,3)=NaN; t(3,1)=NaN; TMc.CrossSectionCovarCorr.poor.CovarianceMatrix=t;
nTc=6+7; % the six off-diagonal correlations and the seven kept covariance entries
ran=0; objv=NaN; dmom=NaN; nbad=NaN; moved=NaN;
try
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'beta'},TMc,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions,vfoptions,simoptions);
    ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:)); nbad=abs(length(calibsummary.targetmomentvec)-nTc);
catch err
    fprintf('RestrictedTargets (c) matrix targets: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RestrictedTargets (c) CorrelationMatrix (diagonal omitted) and the CovarianceMatrix among the poor (one pair omitted), at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current vs target moments, this should be zero: %.3e; number of targeted entries wrong, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',1-ran,objv,dmom,nbad,moved)

%% (d) logmoments by name on a four-level and a three-level name
TMd=struct();
TMd.AutoCorrTransProbs.rich.assets.AutoCorrelation=log(gather(ACFull.rich.assets.AutoCorrelation)); % logged
TMd.AutoCorrTransProbs.earnings.AutoCorrelation=gather(ACFull.earnings.AutoCorrelation); % not logged
TMd.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=log(gather(CSFull.rich.assets.CorrelationWith.consumption)); % logged
TMd.CrossSectionCovarCorr.assets.CovarianceWith.earnings=gather(CSFull.assets.CovarianceWith.earnings); % not logged
expected_d=[1;0;1;0];
caliboptions_d=caliboptions;
caliboptions_d.logmoments=struct();
caliboptions_d.logmoments.AutoCorrTransProbs.rich.assets.AutoCorrelation=1;
caliboptions_d.logmoments.CrossSectionCovarCorr.rich.assets.CorrelationWith.consumption=1;
ran=0; objv=NaN; dl=NaN; dmom=NaN; moved=NaN;
try
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'beta'},TMd,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_d,vfoptions,simoptions);
    ran=1; objv=calibsummary.objvalue; moved=abs(CalibParams.beta-Params.beta); dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:));
    if length(calibsummary.logmoments)==4
        dl=CoreMaxAbs(calibsummary.logmoments(:)-expected_d);
    else
        dl=abs(length(calibsummary.logmoments)-4);
    end
catch err
    fprintf('RestrictedTargets (d) logmoments by name: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RestrictedTargets (d) logmoments by name on a three-level AutoCorr and a four-level cross-section name, at the truth: did not run, this should be zero: %.3e; logmoments vector vs expected, this should be zero: %.3e; objective, this should be zero: %.3e; current vs (log) target moments, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',1-ran,dl,objv,dmom,moved)

output.TMa=TMa;
output.TMb=TMb;

end
