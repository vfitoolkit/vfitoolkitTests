function output=EstimFHorzMLE_PointEstimates(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(5): EstimateLifeCycleModel_MaximumLikelihood point estimates, on one panel (N=1000) simulated by hand from the model at
% the truth (EstimFHorzMLE_SimulatePanel, fixed seed). Estimates beta and the s.d. of the consumption measurement error.
% These are NOT identities: the estimate differs from the truth by sampling error (a finite panel), and the optimizers stop at a
% tolerance. So the checks are "close to zero", and the informational lines give the distance in standard errors.
%   (a) CMA-ES (the default) started at the truth: recovers the truth; the log-likelihood at the estimate is at least that at the truth
%   (b) CMA-ES started away (beta=0.92, sig_c=0.2), with beta constrained to 0-1 and sig_c positive: recovers the truth, and the same
%       estimate as (a)
%   (c) fminsearch started away: the same estimate as (a)

Params.sig_c=0.1; Params.sig_y=0.05; Params.sig_a=0.2; % the truth (measurement error)
truth=[Params.beta; Params.sig_c];
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
N=1000;
Data=EstimFHorzMLE_SimulatePanel(Policy,n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,N,20261016);
LogLtruth=gather(LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions));

EstimParamNames={'beta','sig_c'};
estimoptions.verbose=0;
estimoptions.inopts.TolX=1e-6; % CMA-ES tolerances (on the unconstrained parameters, and on minus the log-likelihood)
estimoptions.inopts.TolFun=1e-3;
estimoptions.inopts.MaxFunEvals=2000;
estimoptions.inopts.DispFinal='off';
estimoptions.inopts.DispModulo=Inf;
estimoptions.inopts.SaveVariables='off';
estimoptions.inopts.LogModulo=0;
estimoptions.inopts.Seed=20261017; % so CMA-ES draws the same candidates every run

%% (a) CMA-ES from the truth
tic;
[EP_a,~,es_a]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);
t_a=toc;
est_a=[EP_a.beta; EP_a.sig_c];
sd_a=[es_a.EstimParamsStdDev.beta; es_a.EstimParamsStdDev.sig_c];
fprintf('PointEstimates (a) CMA-ES from the truth (%.1f s, %i evaluations): beta %.5f (truth %.5f), sig_c %.5f (truth %.5f); distance in standard errors %.2f, %.2f (informational) \n',t_a,es_a.cmaesoutputs.counteval,est_a(1),truth(1),est_a(2),truth(2),(est_a(1)-truth(1))/sd_a(1),(est_a(2)-truth(2))/sd_a(2))
fprintf('PointEstimates (a) max parameter error, this should be close to zero: %.3e; LogL(truth)-LogL(estimate) if positive, this should be close to zero: %.3e \n',max(abs(est_a-truth)),max(0,LogLtruth-es_a.LogL))
fprintf('PointEstimates (a) LogL at the estimate %.4f, at the truth %.4f; gradient at the estimate %s (informational) \n',es_a.LogL,LogLtruth,mat2str(es_a.gradient',4))

%% (b) CMA-ES from away, with constraints
Params_b=Params; Params_b.beta=0.92; Params_b.sig_c=0.2;
estimoptions_b=estimoptions;
estimoptions_b.constrain0to1={'beta'};
estimoptions_b.constrainpositive={'sig_c'};
tic;
[EP_b,~,es_b]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_b,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_b,vfoptions,simoptions);
t_b=toc;
est_b=[EP_b.beta; EP_b.sig_c];
fprintf('PointEstimates (b) CMA-ES from beta=0.92, sig_c=0.2, constrained (%.1f s, %i evaluations): beta %.5f, sig_c %.5f, LogL %.4f (informational) \n',t_b,es_b.cmaesoutputs.counteval,est_b(1),est_b(2),es_b.LogL)
fprintf('PointEstimates (b) max parameter error vs the truth, this should be close to zero: %.3e; vs the estimate of (a), this should be close to zero: %.3e \n',max(abs(est_b-truth)),max(abs(est_b-est_a)))

%% (c) fminsearch from away
Params_c=Params; Params_c.beta=0.92; Params_c.sig_c=0.2;
estimoptions_c=estimoptions;
estimoptions_c.fminalgo=1;
estimoptions_c.constrainpositive={'sig_c'};
tic;
[EP_c,~,es_c]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_c,vfoptions,simoptions);
t_c=toc;
est_c=[EP_c.beta; EP_c.sig_c];
fprintf('PointEstimates (c) fminsearch from beta=0.92, sig_c=0.2 (%.1f s): beta %.5f, sig_c %.5f, LogL %.4f (informational) \n',t_c,est_c(1),est_c(2),es_c.LogL)
fprintf('PointEstimates (c) max parameter error vs the estimate of (a), this should be close to zero: %.3e \n',max(abs(est_c-est_a)))

output.est_a=est_a; output.est_b=est_b; output.est_c=est_c;
end
