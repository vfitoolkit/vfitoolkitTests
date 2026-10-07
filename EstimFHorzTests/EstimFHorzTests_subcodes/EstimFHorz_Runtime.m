function output=EstimFHorz_Runtime(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 8: the runtime of EstimateLifeCycleModel_MethodOfMoments with estimoptions.skipestimation=1 (one model solve and the
% Jacobian block: 1+2*nParams evaluations at each of four step sizes) under estimoptions.whichcombos=1 and =0, for a sparse target set
% (one age-conditional statistic of one function under one restriction at five ages, plus one AllStats mean; four FnsToEvaluate and two
% restrictions in simoptions). Informational: the moments were shown identical in part 4.

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

end
