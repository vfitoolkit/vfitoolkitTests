function output=CalibInfHorz_RoundTrips(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn)
% CalibInfHorzTests part 4: round trips through CalibrateInfHorzAgentModel.
%
% Targets are generated from the model at the true parameters, so the truth is known exactly. A calibration started AT the truth
% must return the truth with an objective of zero (exact: the optimisers cannot improve on a zero residual and return their start).
% One started AWAY from it must get back within the tolerances ("close to zero": only as exact as caliboptions.toleranceparams).
% The two parameters are beta and the wage w; the targets are the mean of assets (beta) and the mean of earnings (w times the mean
% of z, so w is identified exactly), plus the Gini of assets.
%   (a) beta and w, started at the truth, lsqnonlin (fminalgo=8): objective and parameter moves are zero
%   (b) the same with fminsearch (fminalgo=1)
%   (c) started away (beta 0.92, w 1.30), lsqnonlin: recovers (close)
%   (d) started away, fminsearch: recovers (close)
%   (e) started away with caliboptions.whichcombos=0: exactly the same parameters and objective as (c)
%   (f) constrainpositive on w, softplus; (g) the same with the log method; (h) constrain0to1 on beta; (i) constrainAtoB on beta
%       in [0.9,0.99]: each started away recovers (close)
%   (j) omitcalibparam: a vector parameter theta=[log w; discount rate] mapped to (w, beta) by ParametrizeParamsFn, theta(1)
%       calibrated and theta(2) fixed by the mask; started at log(1.3): theta(1) recovered (close), theta(2) unchanged (exact)
%   (k) ParametrizeParamsFn: both entries of theta calibrated from a wrong start (close)
% Each run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;

vfoptions=struct();
simoptions=struct();
% The grid interpolation layer makes the policy, and so every moment, vary smoothly with the parameters (without it the moments are
% piecewise constant in beta on this coarse grid and the optimisers stall; the FHorz bank found this on 2026-10-07). The targets are
% generated with the same options, so the truth is still exact.
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;

%% Targets at the truth
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
TM.AllStats.assets.Mean=gather(AllFull.assets.Mean);
TM.AllStats.earnings.Mean=gather(AllFull.earnings.Mean);
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
fprintf('RoundTrips: targets at the truth, mean assets %.4f, mean earnings %.4f, Gini of assets %.4f \n',TM.AllStats.assets.Mean,TM.AllStats.earnings.Mean,TM.AllStats.assets.Gini)

CalibParamNames={'beta','w'};
truth=[Params.beta;Params.w];
caliboptions=struct();
caliboptions.verbose=0;
output=struct();
output.TargetMoments=TM;

%% (a),(b) started at the truth
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    ran=0; objv=NaN; pmove=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; pmove=CoreMaxAbs([CalibParams.beta;CalibParams.w]-truth);
    catch err
        fprintf('RoundTrips (%s): the calibration errored with: %s \n',char('a'+(alg==1)),[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips (%s) started at the truth, fminalgo=%i (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; max parameter move, this should be zero: %.3e \n',char('a'+(alg==1)),alg,tt,1-ran,objv,pmove)
end

%% (c),(d) started away
Params_away=Params; Params_away.beta=0.92; Params_away.w=1.30; % (not too far below the truth: the further beta falls, the closer the economy gets to everyone at the borrowing constraint, where the moments are flat in beta)
CalibParams_c=[]; calibsummary_c=[];
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    ran=0; objv=NaN; perr=NaN; tt=NaN; rec=[NaN;NaN];
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; rec=[CalibParams.beta;CalibParams.w]; perr=CoreMaxAbs(rec-truth);
        if alg==8
            CalibParams_c=CalibParams; calibsummary_c=calibsummary;
        end
    catch err
        fprintf('RoundTrips (%s): the calibration errored with: %s \n',char('c'+(alg==1)),[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips (%s) started at beta=0.92, w=1.30, fminalgo=%i (%.1f s): recovered beta %.6f and w %.6f (truth %.2f, %.2f), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',char('c'+(alg==1)),alg,tt,rec(1),rec(2),truth(1),truth(2),objv,1-ran,perr)
end

%% (e) whichcombos=0 from the same wrong start, lsqnonlin: identical to (c)
caliboptions_e=caliboptions; caliboptions_e.fminalgo=8; caliboptions_e.whichcombos=0;
ran=0; dp=NaN; dobj=NaN; dmom=NaN;
try
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_e,vfoptions,simoptions);
    ran=1;
    if ~isempty(CalibParams_c)
        dp=CoreMaxAbs([CalibParams.beta;CalibParams.w]-[CalibParams_c.beta;CalibParams_c.w]); dobj=abs(calibsummary.objvalue-calibsummary_c.objvalue);
        dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary_c.currentmomentvec(:));
    end
catch err
    fprintf('RoundTrips (e): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips (e) whichcombos=0 vs =1 from the same start: did not run, this should be zero: %.3e; parameters, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dp,dobj,dmom)

%% (f)-(i) constraints, started away, lsqnonlin
labels={'(f) constrainpositive w, softplus','(g) constrainpositive w, log','(h) constrain0to1 beta','(i) constrainAtoB beta in [0.9,0.99]'};
for cc=1:4
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=8;
    Params_start=Params_away;
    if cc==1
        caliboptions_t.constrainpositive={'w'};
        caliboptions_t.constrainpositivemethod='softplus';
    elseif cc==2
        caliboptions_t.constrainpositive={'w'};
        caliboptions_t.constrainpositivemethod='log';
    elseif cc==3
        caliboptions_t.constrain0to1={'beta'};
    elseif cc==4
        caliboptions_t.constrainAtoB={'beta'};
        caliboptions_t.constrainAtoBlimits.beta=[0.9,0.99];
        Params_start.beta=0.93; % (the start must satisfy the constraint, and not sit on its edge)
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN; rec=[NaN;NaN];
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; rec=[CalibParams.beta;CalibParams.w]; perr=CoreMaxAbs(rec-truth);
    catch err
        fprintf('RoundTrips %s: the calibration errored with: %s \n',labels{cc},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips %s, started away (%.1f s): recovered beta %.6f and w %.6f, objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',labels{cc},tt,rec(1),rec(2),objv,1-ran,perr)
end

%% (j) omitcalibparam on a vector parameter: theta=[log w; discount rate], mapped to (w, beta) by ParametrizeParamsFn; theta(2) fixed
ParametrizeParamsFn=@(Parameters) CalibInfHorz_ParametrizeParamsFn(Parameters);
theta_truth=[log(Params.w); 1/Params.beta-1];
Params_th=Params; Params_th.theta=theta_truth; Params_th=CalibInfHorz_ParametrizeParamsFn(Params_th); % (w and beta unchanged at the truth)
fprintf('RoundTrips (j),(k): ParametrizeParamsFn at the true theta reproduces w and beta, this should be zero: %.3e \n',CoreMaxAbs([Params_th.w-Params.w;Params_th.beta-Params.beta]))
Params_start=Params_th;
Params_start.theta(1)=log(1.3);
caliboptions_j=caliboptions; caliboptions_j.fminalgo=8;
caliboptions_j.omitcalibparam.theta=theta_truth; % the fixed values, NaN where calibrated
caliboptions_j.omitcalibparam.theta(1)=NaN;
ran=0; objv=NaN; e1=NaN; e2=NaN; tt=NaN;
try
    t0=tic;
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'theta'},TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,ParametrizeParamsFn,FnsToEvaluate,caliboptions_j,vfoptions,simoptions);
    tt=toc(t0); ran=1; objv=calibsummary.objvalue;
    th=CalibParams.theta(:);
    e1=abs(th(1)-theta_truth(1)); e2=abs(th(2)-theta_truth(2))+abs(length(th)-2);
catch err
    fprintf('RoundTrips (j): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips (j) omitcalibparam, theta(1)=log w from log(1.3), theta(2) fixed by the mask (%.1f s), objective %.3e: did not run, this should be zero: %.3e; theta(1) error, this should be close to zero: %.3e; theta(2) moved (or the wrong length), this should be zero: %.3e \n',tt,objv,1-ran,e1,e2)

%% (k) ParametrizeParamsFn: both entries of theta calibrated
Params_start=Params_th;
Params_start.theta=[log(1.3); 0.087]; % (beta 0.92)
caliboptions_k=caliboptions; caliboptions_k.fminalgo=8;
ran=0; objv=NaN; perr=NaN; tt=NaN; th=[NaN;NaN];
try
    t0=tic;
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel({'theta'},TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,ParametrizeParamsFn,FnsToEvaluate,caliboptions_k,vfoptions,simoptions);
    tt=toc(t0); ran=1; objv=calibsummary.objvalue; th=CalibParams.theta(:); perr=CoreMaxAbs(th-theta_truth);
catch err
    fprintf('RoundTrips (k): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips (k) ParametrizeParamsFn, theta from (log 1.3, 0.087) (%.1f s): recovered %.6f and %.6f (truth %.6f, %.6f), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,th(1),th(2),theta_truth(1),theta_truth(2),objv,1-ran,perr)

output.CalibParams_c=CalibParams_c;

end
