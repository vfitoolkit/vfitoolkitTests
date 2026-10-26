function output=CalibFHorzPType_RoundTrips(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 4: round trips through CalibrateLifeCycleModel_PType.
%
% Targets are generated from the model at the true parameters. A calibration started AT the truth returns the truth with an
% objective of zero (exact); one started AWAY gets back within the tolerances ("close to zero"). The grid interpolation layer is on,
% so the moments vary smoothly with the parameters.
%   (a),(b) beta and pension (both common), started at the truth, lsqnonlin and fminsearch: zero
%   (c),(d) started away (beta 0.90, pension 0.30): recover (close)
%   (e) caliboptions.whichcombos=0 from the same start: exactly the same parameters, objective and moments as (c)
%   (f) pension PER TYPE as a structure (truth low 0.5, high 0.4), with beta: started away, recovered (close); the result is a structure
%   (g) pension per type as a matrix with N_i as its second dimension: the same, result returned as a matrix
%   (h) constrainpositive on the per-type pension (softplus); (i) constrain0to1 on beta: recovered (close)
%   (j) omitcalibparam on the per-type kappa_j_pt, calibrated at ages 1-6 for both types (the truth has the two types equal, as one
%       omit mask serves every type): calibrated ages recovered (close), the fixed ages unchanged (exact)
%   (k) ParametrizeParamsFn: kappa_j_pt from a per-type intercept (structure) and a common slope, both calibrated from a wrong start

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
Jr=Params.Jr;

vfoptions=struct();
simoptions=struct();
% The grid interpolation layer makes the policy, and so every moment, vary smoothly with the parameters (on the coarse grid alone they
% are piecewise constant and lsqnonlin stalls). The targets are generated with the same options, so the truth is still exact.
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
caliboptions=struct();
caliboptions.verbose=0;

%% Targets at the truth (common parameters): the per-type age-conditional mean of assets (age 1 omitted) and the grouped Gini of assets
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
t=gather(LCPFull.assets.low.Mean); t(1)=NaN; TM.AgeConditionalStats.assets.low.Mean=t;
t=gather(LCPFull.assets.high.Mean); t(1)=NaN; TM.AgeConditionalStats.assets.high.Mean=t;
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);

CalibParamNames={'beta','pension'};
truth=[Params.beta;Params.pension];

%% (a),(b) started at the truth
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    t0=tic;
    [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
    tt=toc(t0);
    fprintf('RoundTrips PType (%s) started at the truth, fminalgo=%i (%.1f s): objective, this should be zero: %.3e; max parameter move, this should be zero: %.3e \n',char('a'+(alg==1)),alg,tt,calibsummary.objvalue,max(abs([CalibParams.beta;CalibParams.pension]-truth)))
end

%% (c),(d) started away
Params_away=Params; Params_away.beta=0.90; Params_away.pension=0.30;
for alg=[8,1]
    caliboptions_t=caliboptions; caliboptions_t.fminalgo=alg;
    t0=tic;
    [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_t,vfoptions,simoptions);
    tt=toc(t0);
    fprintf('RoundTrips PType (%s) started at beta=0.90, pension=0.30, fminalgo=%i (%.1f s): recovered beta %.6f and pension %.6f (truth %.2f, %.2f), objective %.3e; max parameter error, this should be close to zero: %.3e \n',char('c'+(alg==1)),alg,tt,CalibParams.beta,CalibParams.pension,truth(1),truth(2),calibsummary.objvalue,max(abs([CalibParams.beta;CalibParams.pension]-truth)))
    if alg==8
        CalibParams_c=CalibParams; calibsummary_c=calibsummary;
    end
end

%% (e) whichcombos=0 from the same wrong start: identical to (c)
caliboptions_e=caliboptions; caliboptions_e.fminalgo=8; caliboptions_e.whichcombos=0;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_e,vfoptions,simoptions);
fprintf('RoundTrips PType (e) whichcombos=0 vs =1 from the same start: parameters, this should be zero: %.3e; objective, this should be zero: %.3e; current moments, this should be zero: %.3e \n',max(abs([CalibParams.beta;CalibParams.pension]-[CalibParams_c.beta;CalibParams_c.pension])),abs(calibsummary.objvalue-calibsummary_c.objvalue),max(abs(calibsummary.currentmomentvec(:)-calibsummary_c.currentmomentvec(:))))

%% (f),(g) a per-type pension: truth low 0.5, high 0.4; targets regenerated
Params_pt=Params;
Params_pt.pension=struct(); Params_pt.pension.low=0.5; Params_pt.pension.high=0.4;
[~,Policy_pt]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_pt,DiscountFactorParamNames,vfoptions);
StationaryDist_pt=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_pt,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_pt,simoptions);
AllFull_pt=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_pt,Policy_pt,FnsToEvaluate,Params_pt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull_pt=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_pt,Policy_pt,FnsToEvaluate,Params_pt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMpt=struct();
t=gather(LCPFull_pt.assets.low.Mean); t(1)=NaN; TMpt.AgeConditionalStats.assets.low.Mean=t;
t=gather(LCPFull_pt.assets.high.Mean); t(1)=NaN; TMpt.AgeConditionalStats.assets.high.Mean=t;
TMpt.AllStats.assets.Gini=gather(AllFull_pt.assets.Gini);
truth_pt=[0.95;0.5;0.4];
Params_start=Params_pt; Params_start.beta=0.90; Params_start.pension.low=0.30; Params_start.pension.high=0.30;
caliboptions_f=caliboptions; caliboptions_f.fminalgo=8;
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TMpt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_f,vfoptions,simoptions);
tt=toc(t0);
if isstruct(CalibParams.pension)
    rec=[CalibParams.beta;CalibParams.pension.low;CalibParams.pension.high]; formbad=0;
else
    rec=[CalibParams.beta;NaN;NaN]; formbad=1;
end
fprintf('RoundTrips PType (f) pension per type as a structure (truth low 0.5, high 0.4), with beta, started away (%.1f s): recovered beta %.6f, pension low %.6f, high %.6f, objective %.3e; result not a structure, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),rec(3),calibsummary.objvalue,formbad,max(abs(rec-truth_pt)))
Params_start=Params; Params_start.pension=[0.30,0.30]; Params_start.beta=0.90; % the matrix form: N_i as the second dimension
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TMpt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_f,vfoptions,simoptions);
tt=toc(t0);
if isnumeric(CalibParams.pension) && numel(CalibParams.pension)==N_i
    rec=[CalibParams.beta;reshape(CalibParams.pension,[N_i,1])]; formbad=0;
else
    rec=[CalibParams.beta;NaN;NaN]; formbad=1;
end
fprintf('RoundTrips PType (g) pension per type as a matrix [1,N_i], started away (%.1f s): recovered beta %.6f, pension %.6f, %.6f, objective %.3e; result not a numeric vector of N_i, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),rec(3),calibsummary.objvalue,formbad,max(abs(rec-truth_pt)))

%% (h),(i) constraints, started away
Params_start=Params_pt; Params_start.beta=0.90; Params_start.pension.low=0.30; Params_start.pension.high=0.30;
caliboptions_h=caliboptions; caliboptions_h.fminalgo=8; caliboptions_h.constrainpositive={'pension'}; caliboptions_h.constrainpositivemethod='softplus';
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TMpt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_h,vfoptions,simoptions);
tt=toc(t0);
rec=[CalibParams.beta;CalibParams.pension.low;CalibParams.pension.high];
fprintf('RoundTrips PType (h) constrainpositive on the per-type pension, softplus, started away (%.1f s): recovered beta %.6f, pension low %.6f, high %.6f, objective %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),rec(3),calibsummary.objvalue,max(abs(rec-truth_pt)))
caliboptions_i=caliboptions; caliboptions_i.fminalgo=8; caliboptions_i.constrain0to1={'beta'};
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TMpt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_i,vfoptions,simoptions);
tt=toc(t0);
rec=[CalibParams.beta;CalibParams.pension.low;CalibParams.pension.high];
fprintf('RoundTrips PType (i) constrain0to1 on beta, per-type pension, started away (%.1f s): recovered beta %.6f, pension low %.6f, high %.6f, objective %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),rec(3),calibsummary.objvalue,max(abs(rec-truth_pt)))

%% (j) omitcalibparam on kappa_j_pt: the truth has both types with the same profile (one omit mask serves every type)
Params_om=Params;
Params_om.kappa_j_pt=[Params.kappa_j; Params.kappa_j];
[~,Policy_om]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_om,DiscountFactorParamNames,vfoptions);
StationaryDist_om=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_om,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_om,simoptions);
LCPFull_om=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_om,Policy_om,FnsToEvaluate,Params_om,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMj=struct();
t=gather(LCPFull_om.earnings.low.Mean); t(Jr:end)=NaN; TMj.AgeConditionalStats.earnings.low.Mean=t;
t=gather(LCPFull_om.earnings.high.Mean); t(Jr:end)=NaN; TMj.AgeConditionalStats.earnings.high.Mean=t;
Params_start=Params_om;
Params_start.kappa_j_pt(:,1:6)=0.8*Params_om.kappa_j_pt(:,1:6);
caliboptions_j=caliboptions; caliboptions_j.fminalgo=8;
caliboptions_j.omitcalibparam.kappa_j_pt=Params.kappa_j; % the fixed values (the same for both types), NaN where calibrated
caliboptions_j.omitcalibparam.kappa_j_pt(1:6)=NaN;
ran=0; objv=NaN; errcal=NaN; errfix=NaN;
try
    t0=tic;
    [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'kappa_j_pt'},TMj,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_j,vfoptions,simoptions);
    tt=toc(t0);
    ran=1; objv=calibsummary.objvalue;
    if isstruct(CalibParams.kappa_j_pt)
        kj=[reshape(CalibParams.kappa_j_pt.low,[1,N_j]); reshape(CalibParams.kappa_j_pt.high,[1,N_j])];
    else
        kj=reshape(CalibParams.kappa_j_pt,[N_i,N_j]);
    end
    errcal=max(max(abs(kj(:,1:6)-Params_om.kappa_j_pt(:,1:6)))); errfix=max(max(abs(kj(:,7:end)-Params_om.kappa_j_pt(:,7:end))));
catch err
    tt=NaN;
    fprintf('RoundTrips PType (j) omitcalibparam: the calibration errored with: %s \n',err.message)
end
fprintf('RoundTrips PType (j) omitcalibparam on kappa_j_pt, ages 1-6 of both types from 0.8x the truth, the other ages fixed (%.1f s), objective %.3e: did not run, this should be zero: %.3e; calibrated ages max error, this should be close to zero: %.3e; fixed ages moved, this should be zero: %.3e \n',tt,objv,1-ran,errcal,errfix)

%% (k) ParametrizeParamsFn: kappa_j_pt from a per-type intercept (structure) and a common slope
Params_pf=Params;
Params_pf.kappa0=struct(); Params_pf.kappa0.low=0.5; Params_pf.kappa0.high=0.75;
Params_pf.kappaslope=0.05;
Params_pf=CalibFHorzPType_ParametrizeParamsFn(Params_pf); % the truth has kappa_j_pt from the lines
[~,Policy_pf]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_pf,DiscountFactorParamNames,vfoptions);
StationaryDist_pf=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_pf,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_pf,simoptions);
LCP_pf=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_pf,Policy_pf,FnsToEvaluate,Params_pf,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMk=struct();
t=gather(LCP_pf.earnings.low.Mean); t(Jr:end)=NaN; TMk.AgeConditionalStats.earnings.low.Mean=t;
t=gather(LCP_pf.earnings.high.Mean); t(Jr:end)=NaN; TMk.AgeConditionalStats.earnings.high.Mean=t;
Params_start=Params_pf;
Params_start.kappa0.low=0.4; Params_start.kappa0.high=0.6;
Params_start.kappaslope=0.08;
caliboptions_k=caliboptions; caliboptions_k.fminalgo=8;
ParametrizeParamsFn=@(Parameters) CalibFHorzPType_ParametrizeParamsFn(Parameters);
t0=tic;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'kappa0','kappaslope'},TMk,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,caliboptions_k,vfoptions,simoptions);
tt=toc(t0);
rec=[CalibParams.kappa0.low;CalibParams.kappa0.high;CalibParams.kappaslope];
fprintf('RoundTrips PType (k) ParametrizeParamsFn, per-type kappa0 (truth 0.5, 0.75) and common kappaslope (truth 0.05) from (0.4, 0.6, 0.08) (%.1f s): recovered %.6f, %.6f, %.6f, objective %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),rec(3),calibsummary.objvalue,max(abs(rec-[0.5;0.75;0.05])))

output=struct();
output.TargetMoments=TM;
output.CalibParams_c=CalibParams_c;

end
