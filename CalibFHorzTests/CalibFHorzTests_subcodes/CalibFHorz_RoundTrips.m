function output=CalibFHorz_RoundTrips(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% CalibFHorzTests part 4: round trips through CalibrateLifeCycleModel.
%
% Targets are generated from the model at the true parameters, so the truth is known exactly. A calibration started AT the truth
% must return the truth with an objective of zero (exact: the optimisers cannot improve on a zero residual and return their start).
% One started AWAY from it must get back within the tolerances ("close to zero": only as exact as caliboptions.toleranceparams).
%   (a) beta and pension, started at the truth, lsqnonlin (fminalgo=8): objective and parameter moves are zero
%   (b) the same with fminsearch (fminalgo=1)
%   (c) started away (beta 0.90, pension 0.30), lsqnonlin: recovers (close)
%   (d) started away, fminsearch: recovers (close)
%   (e) started away with caliboptions.whichcombos=0: exactly the same parameters and objective as (c)
%   (f) constrainpositive on pension, softplus; (g) the same with the log method; (h) constrain0to1 on beta; (i) constrainAtoB on
%       beta in [0.9,0.99]: each started away recovers (close)
%   (j) omitcalibparam: kappa_j calibrated at ages 1-6 only (the other ages fixed by the mask), targets the age-conditional mean of
%       earnings; started at 0.8 times the truth at those ages: ages 1-6 recovered (close), the other ages unchanged (exact)
%   (k) ParametrizeParamsFn: kappa_j = kappa0 + kappaslope*(agej-1) while working; kappa0 and kappaslope calibrated from a wrong
%       start (close)

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

%% Targets at the truth: the age-conditional mean of assets (age 1 omitted, everyone has none) and the Gini of assets
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; TM.AgeConditionalStats.assets.Mean=t;
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);

CalibParamNames={'beta','pension'};
truth=[Params.beta;Params.pension];
caliboptions=struct();
caliboptions.verbose=0;

%% (a),(b) started at the truth
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    tic;
    [CalibParams,calibsummary]=CalibrateLifeCycleModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
    tt=toc;
    fprintf('RoundTrips (%s) started at the truth, fminalgo=%i (%.1f s): objective, this should be zero: %.3e; max parameter move, this should be zero: %.3e \n',char('a'+(alg==1)),alg,tt,calibsummary.objvalue,max(abs([CalibParams.beta;CalibParams.pension]-truth)))
end

%% (c),(d) started away
Params_away=Params; Params_away.beta=0.90; Params_away.pension=0.30;
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    tic;
    [CalibParams,calibsummary]=CalibrateLifeCycleModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
    tt=toc;
    fprintf('RoundTrips (%s) started at beta=0.90, pension=0.30, fminalgo=%i (%.1f s): recovered beta %.6f and pension %.6f (truth %.2f, %.2f), objective %.3e; max parameter error, this should be close to zero: %.3e \n',char('c'+(alg==1)),alg,tt,CalibParams.beta,CalibParams.pension,truth(1),truth(2),calibsummary.objvalue,max(abs([CalibParams.beta;CalibParams.pension]-truth)))
    if alg==8
        CalibParams_c=CalibParams; calibsummary_c=calibsummary;
    end
end

%% (e) whichcombos=0 from the same wrong start, lsqnonlin: identical to (c)
caliboptions_e=caliboptions; caliboptions_e.fminalgo=8; caliboptions_e.whichcombos=0;
[CalibParams,calibsummary]=CalibrateLifeCycleModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions_e,vfoptions,simoptions);
fprintf('RoundTrips (e) whichcombos=0 vs =1 from the same start: parameters, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',max(abs([CalibParams.beta;CalibParams.pension]-[CalibParams_c.beta;CalibParams_c.pension])),abs(calibsummary.objvalue-calibsummary_c.objvalue),max(abs(calibsummary.currentmomentvec(:)-calibsummary_c.currentmomentvec(:))))

%% (f)-(i) constraints, started away, lsqnonlin
labels={'(f) constrainpositive pension, softplus','(g) constrainpositive pension, log','(h) constrain0to1 beta','(i) constrainAtoB beta in [0.9,0.99]'};
for cc=1:4
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=8;
    Params_start=Params_away;
    if cc==1
        caliboptions_t.constrainpositive={'pension'};
        caliboptions_t.constrainpositivemethod='softplus';
    elseif cc==2
        caliboptions_t.constrainpositive={'pension'};
        caliboptions_t.constrainpositivemethod='log';
    elseif cc==3
        caliboptions_t.constrain0to1={'beta'};
    elseif cc==4
        caliboptions_t.constrainAtoB={'beta'};
        caliboptions_t.constrainAtoBlimits.beta=[0.9,0.99];
        Params_start.beta=0.92; % (the start must satisfy the constraint)
    end
    tic;
    [CalibParams,calibsummary]=CalibrateLifeCycleModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
    tt=toc;
    fprintf('RoundTrips %s, started away (%.1f s): recovered beta %.6f and pension %.6f, objective %.3e; max parameter error, this should be close to zero: %.3e \n',labels{cc},tt,CalibParams.beta,CalibParams.pension,calibsummary.objvalue,max(abs([CalibParams.beta;CalibParams.pension]-truth)))
end

%% (j) omitcalibparam: kappa_j at ages 1-6 calibrated, the rest fixed; targets the age-conditional mean of earnings
TMj=struct();
t=gather(LCPFull.earnings.Mean); t(Jr:end)=NaN; TMj.AgeConditionalStats.earnings.Mean=t;
Params_start=Params;
Params_start.kappa_j(1:6)=0.8*Params.kappa_j(1:6);
caliboptions_j=caliboptions; caliboptions_j.fminalgo=8;
caliboptions_j.omitcalibparam.kappa_j=Params.kappa_j; % the fixed values, NaN where calibrated
caliboptions_j.omitcalibparam.kappa_j(1:6)=NaN;
tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel({'kappa_j'},TMj,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,caliboptions_j,vfoptions,simoptions);
tt=toc;
kj=reshape(CalibParams.kappa_j,[1,N_j]);
fprintf('RoundTrips (j) omitcalibparam, kappa_j at ages 1-6 from 0.8x the truth, the other %i ages fixed (%.1f s), objective %.3e: calibrated ages max error, this should be close to zero: %.3e; fixed ages moved, this should be zero: %.3e \n',N_j-6,tt,calibsummary.objvalue,max(abs(kj(1:6)-Params.kappa_j(1:6))),max(abs(kj(7:end)-Params.kappa_j(7:end))))

%% (k) ParametrizeParamsFn: kappa_j a straight line in age from kappa0 and kappaslope
Params_pf=Params;
Params_pf.kappa0=0.5;
Params_pf.kappaslope=0.05;
Params_pf=CalibFHorz_ParametrizeParamsFn(Params_pf); % the truth has kappa_j from the line
[~,Policy_pf]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_pf,DiscountFactorParamNames,[],vfoptions);
StationaryDist_pf=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_pf,n_d,n_a,n_z,N_j,pi_z,Params_pf,simoptions);
LCP_pf=LifeCycleProfiles_FHorz_Case1(StationaryDist_pf,Policy_pf,FnsToEvaluate,Params_pf,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TMk=struct();
t=gather(LCP_pf.earnings.Mean); t(Jr:end)=NaN; TMk.AgeConditionalStats.earnings.Mean=t;
Params_start=Params_pf;
Params_start.kappa0=0.4;
Params_start.kappaslope=0.08;
caliboptions_k=caliboptions; caliboptions_k.fminalgo=8;
ParametrizeParamsFn=@(Parameters) CalibFHorz_ParametrizeParamsFn(Parameters);
tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel({'kappa0','kappaslope'},TMk,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,caliboptions_k,vfoptions,simoptions);
tt=toc;
fprintf('RoundTrips (k) ParametrizeParamsFn, kappa0 and kappaslope from (0.4, 0.08) (%.1f s): recovered %.6f and %.6f (truth 0.5, 0.05), objective %.3e; max parameter error, this should be close to zero: %.3e \n',tt,CalibParams.kappa0,CalibParams.kappaslope,calibsummary.objvalue,max(abs([CalibParams.kappa0;CalibParams.kappaslope]-[0.5;0.05])))

output=struct();
output.TargetMoments=TM;
output.CalibParams_c=CalibParams_c;

end
