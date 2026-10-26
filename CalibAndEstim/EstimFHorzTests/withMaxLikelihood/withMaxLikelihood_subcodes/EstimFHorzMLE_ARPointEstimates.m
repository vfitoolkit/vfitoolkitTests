function output=EstimFHorzMLE_ARPointEstimates(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(8): EstimateLifeCycleModel_MaximumLikelihood with AR(1) measurement error: estimates rho_c, sigAR_c (the AR(1) component of
% the consumption measurement error) and sig_c (its iid component), on one panel (N=1000) simulated by hand with a CONTINUOUS AR(1)
% component (EstimFHorzMLE_SimulatePanel), while the likelihood uses its Rouwenhorst discretization (nARpoints=7). So the estimates
% differ from the truth by sampling error and by the discretization: "close to zero" checks, the distance in standard errors informational.
%   (a) CMA-ES started away (rho_c=0.3, sigAR_c=0.1, sig_c=0.2): recovers the truth
%   (c) the same with nARpoints=15 and 21: the 7-point grid is coarser than the iid s.d., which biases the estimates; the distance from
%       the truth in standard errors should shrink as the grid gets finer (informational), and the estimates stay close (checks)
%   (b) the same with the AR(1) component left out of the model: the log-likelihood at its estimate is lower (informational), and the
%       iid s.d. it estimates absorbs the AR(1) variance (informational)

Params.sig_c=0.05; Params.sig_y=0.05; Params.sig_a=0.2;
Params.rho_c=0.7; Params.sigAR_c=0.08; % the truth
truth=[Params.rho_c; Params.sigAR_c; Params.sig_c];
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c'; MeasurementError.cons.rho='rho_c'; MeasurementError.cons.sigmaAR='sigAR_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
Data=EstimFHorzMLE_SimulatePanel(Policy,n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,1000,20261021);

estimoptions.verbose=0;
estimoptions.nARpoints=7;
estimoptions.constrain0to1={'rho_c'};
estimoptions.constrainpositive={'sigAR_c','sig_c'};
estimoptions.inopts.TolX=1e-6;
estimoptions.inopts.TolFun=1e-2;
estimoptions.inopts.MaxFunEvals=1000;
estimoptions.inopts.DispFinal='off';
estimoptions.inopts.DispModulo=Inf;
estimoptions.inopts.SaveVariables='off';
estimoptions.inopts.LogModulo=0;
estimoptions.inopts.Seed=20261022; % so CMA-ES draws the same candidates every run

%% (a) CMA-ES from away
Params_a=Params; Params_a.rho_c=0.3; Params_a.sigAR_c=0.1; Params_a.sig_c=0.2;
tic;
[EP_a,~,es_a]=EstimateLifeCycleModel_MaximumLikelihood({'rho_c','sigAR_c','sig_c'},Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);
t_a=toc;
est_a=[EP_a.rho_c; EP_a.sigAR_c; EP_a.sig_c];
sd_a=[es_a.EstimParamsStdDev.rho_c; es_a.EstimParamsStdDev.sigAR_c; es_a.EstimParamsStdDev.sig_c];
fprintf('ARPointEstimates (a) CMA-ES from rho_c=0.3, sigAR_c=0.1, sig_c=0.2 (%.1f s, %i evaluations): rho_c %.4f (truth %.4f), sigAR_c %.4f (truth %.4f), sig_c %.4f (truth %.4f); distance in standard errors %.2f, %.2f, %.2f (informational) \n',t_a,es_a.cmaesoutputs.counteval,est_a(1),truth(1),est_a(2),truth(2),est_a(3),truth(3),(est_a-truth)./sd_a)
fprintf('ARPointEstimates (a) max parameter error, this should be close to zero: %.3e \n',max(abs(est_a-truth)))

%% (b) without the AR(1) component
ME_b=MeasurementError; ME_b.cons=rmfield(ME_b.cons,{'rho','sigmaAR'});
estimoptions_b=estimoptions; estimoptions_b.constrain0to1={}; estimoptions_b.constrainpositive={'sig_c'};
[EP_b,~,es_b]=EstimateLifeCycleModel_MaximumLikelihood({'sig_c'},Data,ME_b,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_b,vfoptions,simoptions);
fprintf('ARPointEstimates (b) without the AR(1) component: sig_c %.4f vs sqrt(sig_c^2+sigAR_c^2/(1-rho_c^2)) %.4f at the truth; LogL %.3f vs %.3f with it; AIC %.3f vs %.3f with it (informational) \n',EP_b.sig_c,sqrt(Params.sig_c^2+Params.sigAR_c^2/(1-Params.rho_c^2)),es_b.LogL,es_a.LogL,es_b.AIC,es_a.AIC)

%% (c) finer AR(1) grids
for nn=[15,21]
    estimoptions_c=estimoptions; estimoptions_c.nARpoints=nn;
    tic;
    [EP_c,~,es_c]=EstimateLifeCycleModel_MaximumLikelihood({'rho_c','sigAR_c','sig_c'},Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_a,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_c,vfoptions,simoptions);
    t_c=toc;
    est_c=[EP_c.rho_c; EP_c.sigAR_c; EP_c.sig_c];
    sd_c=[es_c.EstimParamsStdDev.rho_c; es_c.EstimParamsStdDev.sigAR_c; es_c.EstimParamsStdDev.sig_c];
    fprintf('ARPointEstimates (c) nARpoints=%i (%.1f s, %i evaluations): rho_c %.4f, sigAR_c %.4f, sig_c %.4f; distance in standard errors %.2f, %.2f, %.2f; iid s.d. over AR grid spacing %.2f (informational) \n',nn,t_c,es_c.cmaesoutputs.counteval,est_c(1),est_c(2),est_c(3),(est_c-truth)./sd_c,es_c.measurementerror.cons.sigmaovergridspacing)
    fprintf('ARPointEstimates (c) max parameter error with the finer grid, this should be close to zero: %.3e \n',max(abs(est_c-truth)))
end
fprintf('ARPointEstimates (c) with nARpoints=7 the iid s.d. over AR grid spacing was %.2f (informational) \n',es_a.measurementerror.cons.sigmaovergridspacing)

output.est_a=est_a;
end
