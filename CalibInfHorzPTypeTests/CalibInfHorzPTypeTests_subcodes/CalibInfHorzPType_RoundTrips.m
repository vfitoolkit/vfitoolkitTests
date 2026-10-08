function output=CalibInfHorzPType_RoundTrips(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn)
% CalibInfHorzPTypeTests part 4: round trips through CalibrateInfHorzAgentModel_PType.
%
% Targets are generated from the model at the true parameters, so the truth is known exactly: the grouped mean and Gini of assets
% (beta) and the per-type means of earnings (w_pt per type: earnings are w_pt*z, so each type's wage is identified by its own mean).
% A calibration started AT the truth returns the truth with an objective of zero (exact); one started AWAY gets back within the
% tolerances ("close to zero"). Grid interpolation layer as in the single-type bank.
%   (a),(b) beta (common) and w_pt (per type, given as the vector [1;1.4], N_i as its first dimension), started at the truth,
%           lsqnonlin and fminsearch: zero; the result is a numeric vector of N_i
%   (c),(d) started away (beta 0.92, w_pt [1.2;1.2]): recover (close)
%   (e) caliboptions.whichcombos=0 from the same start: exactly the same parameters, objective and moments as (c)
%   (f) w_pt as a structure over the types: started away, recovered (close); the result is a structure
%   (g) w_pt as a row [1,N_i] (N_i as the second dimension): recovered (close); the result is a numeric vector
%   (h) constrainpositive on the per-type w_pt (softplus); (i) constrain0to1 on beta: recovered (close)
%   (j) omitcalibparam on a common vector parameter theta=[log w_low; log w_high; discount rate] mapped by ParametrizeParamsFn to
%       w_pt (a structure) and beta; theta(1:2) calibrated from log(1.2), theta(3) fixed by the mask: recovered (close), fixed (exact)
%   (k) ParametrizeParamsFn with all three entries of theta calibrated from a wrong start (close)
%   (l) the type weights calibrated: ptypemass=[sharelow;1-sharelow] through ParametrizeParamsFn, sharelow from 0.5 (truth 0.6),
%       identified by the grouped mean of earnings (the types differ by wage): recovered (close)
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
caliboptions=struct();
caliboptions.verbose=0;
output=struct();

%% Targets at the truth
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
TM.AllStats.assets.Mean=gather(AllFull.assets.Mean);
TM.AllStats.earnings.low.Mean=gather(AllFull.earnings.low.Mean);
TM.AllStats.earnings.high.Mean=gather(AllFull.earnings.high.Mean);
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
fprintf('RoundTrips PType: targets at the truth, grouped mean assets %.4f, mean earnings low %.4f and high %.4f, grouped Gini of assets %.4f \n',TM.AllStats.assets.Mean,TM.AllStats.earnings.low.Mean,TM.AllStats.earnings.high.Mean,TM.AllStats.assets.Gini)
output.TargetMoments=TM;

CalibParamNames={'beta','w_pt'};
truth=[Params.beta;Params.w_pt(:)];

%% (a),(b) started at the truth
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    ran=0; objv=NaN; pmove=NaN; tt=NaN; formbad=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType(CalibParamNames,TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue;
        formbad=~(isnumeric(CalibParams.w_pt) && numel(CalibParams.w_pt)==N_i);
        if formbad==0
            pmove=CoreMaxAbs([CalibParams.beta;CalibParams.w_pt(:)]-truth);
        end
    catch err
        fprintf('RoundTrips PType (%s): the calibration errored with: %s \n',char('a'+(alg==1)),[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips PType (%s) beta and the per-type w_pt (vector, N_i first), started at the truth, fminalgo=%i (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; result not a numeric vector of N_i, this should be zero: %.3e; max parameter move, this should be zero: %.3e \n',char('a'+(alg==1)),alg,tt,1-ran,objv,formbad,pmove)
end

%% (c),(d) started away
Params_away=Params; Params_away.beta=0.92; Params_away.w_pt=[1.2;1.2];
CalibParams_c=[]; calibsummary_c=[];
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    ran=0; objv=NaN; perr=NaN; tt=NaN; rec=nan(3,1);
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType(CalibParamNames,TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; rec=[CalibParams.beta;CalibParams.w_pt(:)]; perr=CoreMaxAbs(rec-truth);
        if alg==8
            CalibParams_c=CalibParams; calibsummary_c=calibsummary;
        end
    catch err
        fprintf('RoundTrips PType (%s): the calibration errored with: %s \n',char('c'+(alg==1)),[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips PType (%s) started at beta=0.92, w_pt=[1.2;1.2], fminalgo=%i (%.1f s): recovered beta %.6f, w_pt %.6f and %.6f (truth %.2f, %.2f, %.2f), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',char('c'+(alg==1)),alg,tt,rec(1),rec(2),rec(3),truth(1),truth(2),truth(3),objv,1-ran,perr)
end

%% (e) whichcombos=0 from the same wrong start, lsqnonlin: identical to (c)
caliboptions_e=caliboptions; caliboptions_e.fminalgo=8; caliboptions_e.whichcombos=0;
ran=0; dp=NaN; dobj=NaN; dmom=NaN;
try
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType(CalibParamNames,TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_e,vfoptions,simoptions);
    ran=1;
    if ~isempty(CalibParams_c)
        dp=CoreMaxAbs([CalibParams.beta;CalibParams.w_pt(:)]-[CalibParams_c.beta;CalibParams_c.w_pt(:)]); dobj=abs(calibsummary.objvalue-calibsummary_c.objvalue);
        dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary_c.currentmomentvec(:));
    end
catch err
    fprintf('RoundTrips PType (e): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips PType (e) whichcombos=0 vs =1 from the same start: did not run, this should be zero: %.3e; parameters, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',1-ran,dp,dobj,dmom)

%% (f) w_pt as a structure over the types; (g) as a row [1,N_i]
for ff=1:2
    Params_start=Params; Params_start.beta=0.92;
    if ff==1
        Params_start.w_pt=struct(); Params_start.w_pt.low=1.2; Params_start.w_pt.high=1.2; lab='(f) w_pt as a structure over the types';
    else
        Params_start.w_pt=[1.2,1.2]; lab='(g) w_pt as a row [1,N_i]';
    end
    caliboptions_f=caliboptions; caliboptions_f.fminalgo=8;
    ran=0; objv=NaN; perr=NaN; tt=NaN; formbad=NaN; rec=nan(3,1);
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType(CalibParamNames,TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_f,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue;
        if ff==1
            formbad=~isstruct(CalibParams.w_pt);
            if formbad==0
                rec=[CalibParams.beta;CalibParams.w_pt.low;CalibParams.w_pt.high];
            end
        else
            formbad=~(isnumeric(CalibParams.w_pt) && numel(CalibParams.w_pt)==N_i);
            if formbad==0
                rec=[CalibParams.beta;CalibParams.w_pt(:)];
            end
        end
        perr=CoreMaxAbs(rec-truth);
    catch err
        fprintf('RoundTrips PType %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips PType %s, started away (%.1f s): recovered beta %.6f, w_pt %.6f and %.6f, objective %.3e: did not run, this should be zero: %.3e; result in the wrong form, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',lab,tt,rec(1),rec(2),rec(3),objv,1-ran,formbad,perr)
end

%% (h),(i) constraints, started away (the structure form)
Params_start=Params; Params_start.beta=0.92; Params_start.w_pt=struct(); Params_start.w_pt.low=1.2; Params_start.w_pt.high=1.2;
labels={'(h) constrainpositive on the per-type w_pt, softplus','(i) constrain0to1 on beta'};
for cc=1:2
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=8;
    if cc==1
        caliboptions_t.constrainpositive={'w_pt'}; caliboptions_t.constrainpositivemethod='softplus';
    else
        caliboptions_t.constrain0to1={'beta'};
    end
    ran=0; objv=NaN; perr=NaN; tt=NaN; rec=nan(3,1);
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType(CalibParamNames,TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
        tt=toc(t0); ran=1; objv=calibsummary.objvalue; rec=[CalibParams.beta;CalibParams.w_pt.low;CalibParams.w_pt.high]; perr=CoreMaxAbs(rec-truth);
    catch err
        fprintf('RoundTrips PType %s: the calibration errored with: %s \n',labels{cc},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('RoundTrips PType %s, started away (%.1f s): recovered beta %.6f, w_pt %.6f and %.6f, objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',labels{cc},tt,rec(1),rec(2),rec(3),objv,1-ran,perr)
end

%% (j) omitcalibparam on a common vector parameter theta, mapped by ParametrizeParamsFn to w_pt (structure) and beta; theta(3) fixed
ParametrizeParamsFn=@(Parameters) CalibInfHorzPType_ParametrizeParamsFn(Parameters);
theta_truth=[log(Params.w_pt(1)); log(Params.w_pt(2)); 1/Params.beta-1];
Params_th=Params; Params_th.theta=theta_truth; Params_th=CalibInfHorzPType_ParametrizeParamsFn(Params_th);
fprintf('RoundTrips PType (j),(k): ParametrizeParamsFn at the true theta reproduces w_pt and beta, this should be zero: %.3e \n',CoreMaxAbs([Params_th.w_pt.low-Params.w_pt(1);Params_th.w_pt.high-Params.w_pt(2);Params_th.beta-Params.beta]))
Params_start=Params_th;
Params_start.theta(1:2)=log(1.2);
caliboptions_j=caliboptions; caliboptions_j.fminalgo=8;
caliboptions_j.omitcalibparam.theta=theta_truth; % the fixed values, NaN where calibrated
caliboptions_j.omitcalibparam.theta(1:2)=NaN;
ran=0; objv=NaN; e12=NaN; e3=NaN; tt=NaN;
try
    t0=tic;
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType({'theta'},TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,caliboptions_j,vfoptions,simoptions);
    tt=toc(t0); ran=1; objv=calibsummary.objvalue;
    th=CalibParams.theta(:);
    e12=CoreMaxAbs(th(1:2)-theta_truth(1:2)); e3=abs(th(3)-theta_truth(3))+abs(length(th)-3);
catch err
    fprintf('RoundTrips PType (j): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips PType (j) omitcalibparam on theta, theta(1:2) (the log wages) from log(1.2), theta(3) fixed by the mask (%.1f s), objective %.3e: did not run, this should be zero: %.3e; theta(1:2) error, this should be close to zero: %.3e; theta(3) moved (or the wrong length), this should be zero: %.3e \n',tt,objv,1-ran,e12,e3)

%% (k) ParametrizeParamsFn: all three entries of theta calibrated
Params_start=Params_th;
Params_start.theta=[log(1.2); log(1.2); 0.087];
caliboptions_k=caliboptions; caliboptions_k.fminalgo=8;
ran=0; objv=NaN; perr=NaN; tt=NaN; th=nan(3,1);
try
    t0=tic;
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType({'theta'},TM,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,caliboptions_k,vfoptions,simoptions);
    tt=toc(t0); ran=1; objv=calibsummary.objvalue; th=CalibParams.theta(:); perr=CoreMaxAbs(th-theta_truth);
catch err
    fprintf('RoundTrips PType (k): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips PType (k) ParametrizeParamsFn, theta from (log 1.2, log 1.2, 0.087) (%.1f s): recovered %.6f, %.6f, %.6f (truth %.6f, %.6f, %.6f), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,th(1),th(2),th(3),theta_truth(1),theta_truth(2),theta_truth(3),objv,1-ran,perr)

%% (l) the type weights calibrated: ptypemass=[sharelow;1-sharelow] through ParametrizeParamsFn
ParametrizeWeightsFn=@(Parameters) CalibInfHorzPType_ParametrizeWeightsFn(Parameters);
TMl=struct();
TMl.AllStats.earnings.Mean=gather(AllFull.earnings.Mean); % the grouped mean of earnings: 0.6*1*E[z]+0.4*1.4*E[z], so it pins the share
TMl.AllStats.assets.Mean=gather(AllFull.assets.Mean);
Params_start=Params; Params_start.sharelow=0.5;
caliboptions_l=caliboptions; caliboptions_l.fminalgo=8;
ran=0; objv=NaN; perr=NaN; tt=NaN; rec=NaN;
try
    t0=tic;
    [CalibParams,calibsummary]=CalibrateInfHorzAgentModel_PType({'sharelow'},TMl,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,PTypeDistParamNames,ParametrizeWeightsFn,FnsToEvaluate,caliboptions_l,vfoptions,simoptions);
    tt=toc(t0); ran=1; objv=calibsummary.objvalue; rec=CalibParams.sharelow; perr=abs(rec-0.6);
catch err
    fprintf('RoundTrips PType (l): the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('RoundTrips PType (l) the type weights calibrated through ParametrizeParamsFn, sharelow from 0.5 (truth 0.6) (%.1f s): recovered %.6f, objective %.3e: did not run, this should be zero: %.3e; error, this should be close to zero: %.3e \n',tt,rec,objv,1-ran,perr)

output.CalibParams_c=CalibParams_c;

end
