function output=EstimFHorz_Runtime(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 8: the runtime of EstimateLifeCycleModel_MethodOfMoments with estimoptions.skipestimation=1 (one model solve and the
% Jacobian block: 1+2*nParams evaluations at each of four step sizes) under estimoptions.whichcombos=1 and =0, for a sparse target set
% (one age-conditional statistic of one function under one restriction at five ages, plus one AllStats mean; four FnsToEvaluate and two
% restrictions in simoptions). Informational: the moments were shown identical in part 4. (b): estimoptions.skipjacobian=1 (one evaluation,
% no Jacobian) vs skipestimation alone, and ignored without skipestimation; timed.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+w*kappa_j*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
TM=struct();
t=nan(1,N_j); t(6:10)=1; TM.AgeConditionalStats.poor.earnings.Mean=t; % (values do not matter for the timing; the poor have mass at every age, the rich may not, and a NaN moment breaks rank(J))
TM.AllStats.assets.Mean=1;
nT=6; W=eye(nT); Omega=eye(nT);
estimoptions=struct(); estimoptions.verbose=0; estimoptions.skipestimation=1;
times=zeros(2,1);
for wc=[1,0]
    estimoptions_t=estimoptions; estimoptions_t.whichcombos=wc;
    EstimateLifeCycleModel_MethodOfMoments({'beta','pension'},TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions); % warm-up
    t0=tic;
    EstimateLifeCycleModel_MethodOfMoments({'beta','pension'},TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions);
    times(2-wc)=toc(t0);
end
fprintf('Runtime: skipestimation (21 objective evaluations), sparse targets, whichcombos=1 %.2f s, whichcombos=0 %.2f s, ratio %.2f (informational) \n',times(1),times(2),times(1)/times(2))

output=struct();
output.times=times;

%% (b) estimoptions.skipjacobian=1 (with skipestimation=1): one evaluation at the input parameters, no Jacobian
% The same objective and moment vectors as skipestimation alone (the Jacobian block's centre-point evaluation), the standard deviations and
% confidence intervals NaN, no J, Sigma, sensitivity matrix or local identification, the warning field present; and with skipestimation=0
% skipjacobian is ignored (an estimation from the truth, lsqnonlin, reports J). Timed against the skipestimation run above (informational).
estimoptions_t=estimoptions; estimoptions_t.whichcombos=1;
[EP1,CI1,es1]=EstimateLifeCycleModel_MethodOfMoments({'beta','pension'},TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions); % skipestimation=1, the Jacobian
ran=0; d=nan(1,7); tsj=NaN;
try
    estimoptions_t.skipjacobian=1;
    t0=tic;
    [EP2,CI2,es2]=EstimateLifeCycleModel_MethodOfMoments({'beta','pension'},TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions);
    tsj=toc(t0);
    ran=1;
    d(1)=CoreMaxAbs([es2.objectivefnval-es1.objectivefnval; es2.currentmomentvec(:)-es1.currentmomentvec(:); es2.targetmomentvec(:)-es1.targetmomentvec(:)]); % the same objective and moment vectors
    d(2)=CoreMaxAbs([EP2.beta;EP2.pension]-[EP1.beta;EP1.pension]); % the same (input) parameters
    d(3)=double(~(all(isnan(es2.EstimParamsStdDev.beta(:))) && all(isnan(es2.EstimParamsStdDev.pension(:))))); % standard deviations NaN
    d(4)=double(~(all(isnan(CI2.beta(:))) && all(isnan(CI2.pension(:))) && isequal(size(CI2.beta),size(CI1.beta)))); % confidence intervals NaN, same shape
    d(5)=double(isfield(es2.variousmatrices,'J')+isfield(es2.variousmatrices,'Sigma')+isfield(es2,'sensitivitymatrix')+isfield(es2,'localidentification')+isfield(es2,'doublechecks')); % nothing built on the Jacobian
    d(6)=double(~isfield(es2,'warningskipjacobian')); % the warning
catch err
    fprintf('Runtime (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
try
    estimoptions_u=estimoptions_t; estimoptions_u.skipestimation=0; estimoptions_u.fminalgo=8; % skipjacobian=1 is ignored without skipestimation=1
    TM_u=TM; cm=es1.currentmomentvec; TM_u.AllStats.assets.Mean=cm(1); TM_u.AgeConditionalStats.poor.earnings.Mean(6:10)=cm(2:6); % the targets at the truth (the vector is AllStats first, then the ages), so lsqnonlin from the truth stops at once
    [~,~,es3]=EstimateLifeCycleModel_MethodOfMoments({'beta','pension'},TM_u,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_u,vfoptions,simoptions);
    d(7)=double(~isfield(es3.variousmatrices,'J') + isfield(es3,'warningskipjacobian')); % J reported, no skipjacobian warning
catch err
    fprintf('Runtime (b): the estimation from the truth with skipjacobian=1 (ignored) errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Runtime (b) skipjacobian=1 with skipestimation=1: did not run, this should be zero: %.3e; objective and moment vectors vs skipestimation alone, this should be zero: %.3e; EstimParams, this should be zero: %.3e; std devs not NaN, this should be zero: %.3e; conf ints not NaN, this should be zero: %.3e; Jacobian-based fields present, this should be zero: %.3e; warning absent, this should be zero: %.3e; skipjacobian not ignored with skipestimation=0, this should be zero: %.3e \n',1-ran,d(1),d(2),d(3),d(4),d(5),d(6),d(7))
fprintf('Runtime (b) skipjacobian=1 %.2f s vs skipestimation with the Jacobian (whichcombos=1) %.2f s, ratio %.3f (informational; one evaluation vs 1+2*nParams at each of four steps) \n',tsj,times(1),tsj/times(1))
output.skipjacobiantime=tsj;
output.skipjacobianchecks=d;

end
