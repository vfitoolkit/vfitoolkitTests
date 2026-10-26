function output=EstimFHorzPTypeMLE_PointEstimates(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(3): EstimateLifeCycleModel_PType_MaximumLikelihood point estimates, on one panel (N=1000: 700 of the first type, 300 of the
% second) simulated by hand from the model at the truth (EstimFHorzPTypeMLE_SimulatePTypePanel, fixed seed). Estimates beta, the s.d.
% of the consumption measurement error, and pshare_high (the mass of the second type, via EstimFHorzPTypeMLE_Parametrize).
% NOT identities (sampling error, optimizer tolerance): "close to zero" checks, with the distance in standard errors informational.
%   (a) CMA-ES started at the truth: recovers the truth; the log-likelihood at the estimate is at least that at the truth
%   (b) CMA-ES started away (beta=0.92, sig_c=0.2, pshare_high=0.5), constrained: the same estimate as (a)
%   (c) the posterior type probabilities at the estimate pick the true type

Params.sig_c=0.1; Params.sig_y=0.05; Params.sig_a=0.2; % the truth (measurement error)
Params.pshare_high=Params.(PTypeDistParamNames{1})(2);
truth=[Params.beta; Params.sig_c; Params.pshare_high];
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j_pt*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt,agej,Jr) (agej<Jr)*w*kappa_j_pt*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';

[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
[Data,truetype]=EstimFHorzPTypeMLE_SimulatePTypePanel(Policy,n_a,n_z,N_j,Names_i,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,1000,20261018);
N=size(Data.cons,1);
LogLtruth=gather(LogLikelihood_FHorz_PType(Data,MeasurementError,jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions));

EstimParamNames={'beta','sig_c','pshare_high'};
ParametrizeParamsFn=@(Params) EstimFHorzPTypeMLE_Parametrize(Params);
estimoptions.verbose=0;
estimoptions.constrainpositive={'sig_c'};
estimoptions.constrain0to1={'beta','pshare_high'};
estimoptions.inopts.TolX=1e-6; % CMA-ES tolerances (on the unconstrained parameters, and on minus the log-likelihood)
estimoptions.inopts.TolFun=1e-2;
estimoptions.inopts.MaxFunEvals=800; % each evaluation solves both types, so keep the run bounded
estimoptions.inopts.DispFinal='off';
estimoptions.inopts.DispModulo=Inf;
estimoptions.inopts.SaveVariables='off';
estimoptions.inopts.LogModulo=0;
estimoptions.inopts.Seed=20261019; % so CMA-ES draws the same candidates every run

%% (a) CMA-ES from the truth
tic;
[EP_a,~,es_a]=EstimateLifeCycleModel_PType_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions,vfoptions,simoptions);
t_a=toc;
est_a=[EP_a.beta; EP_a.sig_c; EP_a.pshare_high];
sd_a=[es_a.EstimParamsStdDev.beta; es_a.EstimParamsStdDev.sig_c; es_a.EstimParamsStdDev.pshare_high];
fprintf('PTypePointEstimates (a) CMA-ES from the truth (%.1f s, %i evaluations): beta %.5f (truth %.5f), sig_c %.5f (truth %.5f), pshare_high %.4f (truth %.4f); distance in standard errors %.2f, %.2f, %.2f (informational) \n',t_a,es_a.cmaesoutputs.counteval,est_a(1),truth(1),est_a(2),truth(2),est_a(3),truth(3),(est_a-truth)./sd_a)
fprintf('PTypePointEstimates (a) max parameter error, this should be close to zero: %.3e; LogL(truth)-LogL(estimate) if positive, this should be close to zero: %.3e \n',max(abs(est_a-truth)),max(0,LogLtruth-es_a.LogL))
fprintf('PTypePointEstimates (a) type masses at the estimate vs pshare_high, this should be zero: %.3e \n',CoreMaxAbs(es_a.ptypemasses-[1-est_a(3); est_a(3)]))

%% (b) CMA-ES from away
Params_b=Params; Params_b.beta=0.92; Params_b.sig_c=0.2; Params_b.pshare_high=0.5;
Params_b=EstimFHorzPTypeMLE_Parametrize(Params_b);
tic;
[EP_b,~,es_b]=EstimateLifeCycleModel_PType_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions,vfoptions,simoptions);
t_b=toc;
est_b=[EP_b.beta; EP_b.sig_c; EP_b.pshare_high];
fprintf('PTypePointEstimates (b) CMA-ES from beta=0.92, sig_c=0.2, pshare_high=0.5 (%.1f s, %i evaluations): beta %.5f, sig_c %.5f, pshare_high %.4f, LogL %.4f (informational) \n',t_b,es_b.cmaesoutputs.counteval,est_b(1),est_b(2),est_b(3),es_b.LogL)
fprintf('PTypePointEstimates (b) max parameter error vs the truth, this should be close to zero: %.3e; vs the estimate of (a), this should be close to zero: %.3e \n',max(abs(est_b-truth)),max(abs(est_b-est_a)))

%% (c) the posterior picks the true type
posttrue=es_a.PTypePosterior(sub2ind([N,length(Names_i)],(1:N)',truetype));
fprintf('PTypePointEstimates (c) one minus the mean posterior probability of the true type at the estimate, this should be close to zero: %.3e (share classified correctly %.4f, informational) \n',1-mean(posttrue),mean(posttrue>0.5))

output.est_a=est_a; output.est_b=est_b;
end
