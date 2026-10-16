function output=EstimFHorzMLE_StdErrors(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(4): EstimateLifeCycleModel_MaximumLikelihood at the truth with estimoptions.skipestimation=1 (no optimizer), so that every
% reported number can be recomputed by hand from direct LogLikelihood_FHorz calls. Identities throughout:
%   (a) LogL and LogLi at the estimate vs a direct solve + LogLikelihood_FHorz
%   (b) the scores (and the individual log-likelihoods behind them) vs direct calls at the reported up/down parameter vectors
%   (c) OPG=Scores'*Scores, Sigma_OPG, the standard errors and the confidence intervals from it; the gradient; local identification
%   (d) the Hessian vs direct calls at the reported points (both parameters stepped up / down), Sigma_Hessian and the sandwich
%   (e) AIC, BIC and the number of observations; the measurement error noise share
%   (f) estimoptions.stderrors='sandwich' and 'Hessian' put those standard errors in EstimParamsStdDev and the confidence intervals
%   (g) omitestimparam: an age-dependent measurement error s.d. with only two ages estimated: standard errors NaN exactly at the fixed ages
%   (h) refinecheck_ngridinterp: the log-likelihood at the estimate with ngridinterp=7 vs a direct solve with ngridinterp=7
%   (i) refusals: fminalgo=8, a grid interpolation layer mismatch, an unknown stderrors
% The data: one panel (N=500) simulated by hand from the model at the truth (EstimFHorzMLE_SimulatePanel, fixed seed).

Params.sig_c=0.1; Params.sig_y=0.05; Params.sig_a=0.2; % the truth (measurement error)
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j,agej,Jr) (agej<Jr)*w*kappa_j*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';

[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
N=500;
Data=EstimFHorzMLE_SimulatePanel(Policy,n_a,n_z,N_j,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,N,20261015);

EstimParamNames={'beta','sig_c'};
estimoptions.verbose=0;
estimoptions.skipestimation=1;
estimoptions.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
estimoptions.constrainpositive={'sig_c'};
[EstimParams,EstimParamsConfInts,estsummary]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);

% A direct solve + log-likelihood at a parameter vector [beta; sig_c]
K=2;

%% (a)
[LogLhand,LogLihand]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf('StdErrors (a) LogL at the estimate vs direct, this should be zero: %.3e; LogLi, this should be zero: %.3e \n',abs(estsummary.LogL-gather(LogLhand)),CoreMaxAbs(estsummary.LogLi-gather(LogLihand)))
fprintf('StdErrors (a) EstimParams vs the truth (skipestimation), this should be zero: %.3e \n',CoreMaxAbs([EstimParams.beta-Params.beta; EstimParams.sig_c-Params.sig_c]))

%% (b) the scores
dc=estsummary.doublechecks.epsilon10000; % the step (eedefault=3, set above)
LogLi_up=zeros(N,K); LogLi_down=zeros(N,K);
for pp=1:K
    for updown=1:2
        if updown==1
            vec=dc.estimparamsvec; vec(pp)=dc.estimparamsvecup(pp);
        else
            vec=dc.estimparamsvec; vec(pp)=dc.estimparamsvecdown(pp);
        end
        P=Params; P.beta=vec(1); P.sig_c=vec(2);
        [~,Pol]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,P,DiscountFactorParamNames,[],vfoptions);
        [~,temp]=LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Pol,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
        if updown==1
            LogLi_up(:,pp)=gather(temp);
        else
            LogLi_down(:,pp)=gather(temp);
        end
    end
end
Shand=(LogLi_up-LogLi_down)./((dc.estimparamsvecup-dc.estimparamsvecdown)');
fprintf('StdErrors (b) LogLi at the up steps vs direct, this should be zero: %.3e; at the down steps, this should be zero: %.3e \n',CoreMaxAbs(dc.LogLi_up-LogLi_up),CoreMaxAbs(dc.LogLi_down-LogLi_down))
fprintf('StdErrors (b) Scores vs centered differences of direct calls, this should be zero: %.3e (relative to the largest score %.3e) \n',CoreMaxAbs(estsummary.variousmatrices.Scores-Shand),max(abs(Shand(:))))

%% (c) OPG and the standard errors from it
S=estsummary.variousmatrices.Scores;
fprintf('StdErrors (c) OPG vs Scores''*Scores, this should be zero: %.3e; Sigma_OPG vs inv(OPG), this should be zero: %.3e \n',CoreMaxAbs(estsummary.variousmatrices.OPG-S'*S),CoreMaxAbs(estsummary.variousmatrices.Sigma_OPG-inv(S'*S)))
sd=sqrt(diag(inv(S'*S)));
fprintf('StdErrors (c) EstimParamsStdDev (default OPG) vs sqrt(diag(inv(OPG))), this should be zero: %.3e \n',CoreMaxAbs([estsummary.EstimParamsStdDev.beta;estsummary.EstimParamsStdDev.sig_c]-sd))
ci=[[Params.beta;Params.sig_c]-1.645*sd, [Params.beta;Params.sig_c]+1.645*sd];
fprintf('StdErrors (c) 90-percent confidence intervals vs estimate -/+ 1.645 sd, this should be zero: %.3e \n',CoreMaxAbs([EstimParamsConfInts.beta;EstimParamsConfInts.sig_c]-ci))
fprintf('StdErrors (c) gradient vs sum of the scores, this should be zero: %.3e; rank(OPG) is full, this should be one: %i \n',CoreMaxAbs(estsummary.gradient-sum(S,1)'),double(estsummary.localidentification.yesidentified))
fprintf('StdErrors (c) standard errors (OPG): beta %.5f, sig_c %.5f (informational) \n',sd(1),sd(2))

%% (d) the Hessian
hs=estsummary.doublechecks.hessian;
L0=sum(LogLihand);
Lup=sum(LogLi_up,1)'; Ldown=sum(LogLi_down,1)';
x=dc.estimparamsvec; xup=dc.estimparamsvecup; xdown=dc.estimparamsvecdown;
Hhand=zeros(K,K);
for jj=1:K
    Hhand(jj,jj)=(Lup(jj)-2*L0+Ldown(jj))/((xup(jj)-x(jj))*(x(jj)-xdown(jj)));
end
vec=x; vec(1)=xup(1); vec(2)=xup(2);
P=Params; P.beta=vec(1); P.sig_c=vec(2);
[~,Pol]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,P,DiscountFactorParamNames,[],vfoptions);
Lupup=gather(LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Pol,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions));
vec=x; vec(1)=xdown(1); vec(2)=xdown(2);
P=Params; P.beta=vec(1); P.sig_c=vec(2);
[~,Pol]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,P,DiscountFactorParamNames,[],vfoptions);
Ldowndown=gather(LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Pol,FnsToEvaluate,P,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions));
Hhand(1,2)=(Lupup-Lup(1)-Lup(2)+2*L0-Ldown(1)-Ldown(2)+Ldowndown)/(2*(xup(1)-x(1))*(xup(2)-x(2)));
Hhand(2,1)=Hhand(1,2);
fprintf('StdErrors (d) the points: LogL_upup vs direct, this should be zero: %.3e; LogL_downdown vs direct, this should be zero: %.3e \n',abs(hs.LogL_upup(1,2)-Lupup),abs(hs.LogL_downdown(1,2)-Ldowndown))
fprintf('StdErrors (d) Hessian vs second differences of direct calls, this should be zero: %.3e (relative to the largest entry %.3e) \n',CoreMaxAbs(estsummary.variousmatrices.Hessian-Hhand),max(abs(Hhand(:))))
SH=inv(-Hhand);
fprintf('StdErrors (d) Sigma_Hessian vs inv(-Hessian), this should be zero: %.3e; Sigma_sandwich vs inv(-H)*OPG*inv(-H), this should be zero: %.3e \n', ...
    CoreMaxAbs(estsummary.variousmatrices.Sigma_Hessian-inv(-estsummary.variousmatrices.Hessian)),CoreMaxAbs(estsummary.variousmatrices.Sigma_sandwich-SH*(S'*S)*SH))
fprintf('StdErrors (d) standard errors: Hessian beta %.5f sig_c %.5f; sandwich beta %.5f sig_c %.5f (informational) \n',sqrt(SH(1,1)),sqrt(SH(2,2)),sqrt(estsummary.variousmatrices.Sigma_sandwich(1,1)),sqrt(estsummary.variousmatrices.Sigma_sandwich(2,2)))

%% (e) information criteria, measurement error diagnostics
nobs=sum(~isnan(Data.cons(:)))+sum(~isnan(Data.earnings(:)))+sum(~isnan(Data.assets(:)));
fprintf('StdErrors (e) nobservations vs count, this should be zero: %i; AIC vs its formula, this should be zero: %.3e; BIC vs its formula, this should be zero: %.3e \n',estsummary.nobservations-nobs,abs(estsummary.AIC-(2*K-2*L0)),abs(estsummary.BIC-(K*log(N)-2*L0)))
tempd=Data.cons(~isnan(Data.cons));
fprintf('StdErrors (e) noise share of cons vs sig_c^2/var(cons data), this should be zero: %.3e; sig_y reported, this should be zero: %.3e \n',abs(estsummary.measurementerror.cons.noiseshare-Params.sig_c^2/var(tempd)),abs(estsummary.measurementerror.earnings.sigma-Params.sig_y))

%% (f) stderrors='sandwich' and 'Hessian'
stderrs={'sandwich','Hessian'};
for ss=1:2
    estimoptions_f=estimoptions; estimoptions_f.stderrors=stderrs{ss};
    [~,CI_f,es_f]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_f,vfoptions,simoptions);
    sd_f=[es_f.(['EstimParamsStdDev_',stderrs{ss}]).beta; es_f.(['EstimParamsStdDev_',stderrs{ss}]).sig_c];
    fprintf('StdErrors (f) stderrors=%s: EstimParamsStdDev vs EstimParamsStdDev_%s, this should be zero: %.3e; vs this run''s Sigma_%s, this should be zero: %.3e; confidence intervals, this should be zero: %.3e \n',stderrs{ss},stderrs{ss}, ...
        CoreMaxAbs([es_f.EstimParamsStdDev.beta;es_f.EstimParamsStdDev.sig_c]-sd_f),stderrs{ss},CoreMaxAbs(sd_f-sqrt(diag(es_f.variousmatrices.(['Sigma_',stderrs{ss}])))), ...
        CoreMaxAbs([CI_f.beta;CI_f.sig_c]-[[Params.beta;Params.sig_c]-1.645*sd_f, [Params.beta;Params.sig_c]+1.645*sd_f]))
end

%% (g) omitestimparam: an age-dependent s.d. of the consumption measurement error, only ages 3 and 12 estimated
Params_g=Params; Params_g.sig_c_j=Params.sig_c*ones(1,N_j);
ME_g=MeasurementError; ME_g.cons.sigma='sig_c_j';
estimoptions_g=estimoptions; estimoptions_g.constrainpositive={};
omitmask=Params_g.sig_c_j; omitmask([3,12])=NaN;
estimoptions_g.omitestimparam.sig_c_j=omitmask;
[EP_g,CI_g,es_g]=EstimateLifeCycleModel_MaximumLikelihood({'beta','sig_c_j'},Data,ME_g,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_g,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_g,vfoptions,simoptions);
sdg=es_g.EstimParamsStdDev.sig_c_j;
fprintf('StdErrors (g) omitestimparam: std devs NaN exactly at the fixed ages (pattern mismatches), this should be zero: %i; size of the std devs vs the parameter, this should be zero: %i \n',sum(isnan(sdg(:))'~=~isnan(omitmask)),numel(sdg)-N_j)
tempsd=sqrt(diag(es_g.variousmatrices.Sigma_OPG));
fprintf('StdErrors (g) omitestimparam: the two estimated std devs vs sqrt(diag(Sigma)), this should be zero: %.3e; scores columns minus the estimated entries, this should be zero: %i; confidence interval rows minus the ages, this should be zero: %i \n',CoreMaxAbs(reshape(sdg([3,12]),[],1)-tempsd(2:3)),size(es_g.variousmatrices.Scores,2)-3,size(CI_g.sig_c_j,1)-N_j)
fprintf('StdErrors (g) EstimParams.sig_c_j vs the input (skipestimation), this should be zero: %.3e \n',CoreMaxAbs(EP_g.sig_c_j(:)-Params_g.sig_c_j(:)))

%% (h) refinecheck_ngridinterp
estimoptions_h=estimoptions; estimoptions_h.refinecheck_ngridinterp=7;
[~,~,es_h]=EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions);
vfo7=vfoptions; vfo7.ngridinterp=7; sio7=simoptions; sio7.ngridinterp=7;
[~,Pol7]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfo7);
L7=gather(LogLikelihood_FHorz(Data,MeasurementError,jequaloneDist,Pol7,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,sio7));
fprintf('StdErrors (h) refinecheck: LogL with the bank''s ngridinterp vs the estimate''s LogL, this should be zero: %.3e; with the coarser ngridinterp vs direct, this should be zero: %.3e \n',abs(es_h.refinecheck.LogL(1)-es_h.LogL),abs(es_h.refinecheck.LogL(2)-L7))
fprintf('StdErrors (h) LogL with ngridinterp=%i: %.4f, with 7: %.4f (informational) \n',vfoptions.ngridinterp,es_h.refinecheck.LogL(1),es_h.refinecheck.LogL(2))

%% (i) refusals
errcases={'fminalgo=8','gridinterplayer mismatch','ngridinterp mismatch','unknown stderrors'};
errmsgs={'estimoptions.fminalgo must be 4','gridinterplayer must be the same','ngridinterp must both be set','estimoptions.stderrors must be'};
for ee=1:length(errcases)
    esto=estimoptions; vfo=vfoptions; sio=simoptions;
    switch ee
        case 1
            esto.fminalgo=8;
        case 2
            sio.gridinterplayer=0;
        case 3
            sio.ngridinterp=vfoptions.ngridinterp+1;
        case 4
            esto.stderrors='bootstrap';
    end
    errok=0;
    try
        EstimateLifeCycleModel_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,[],FnsToEvaluate,esto,vfo,sio);
    catch err
        errok=contains(err.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',err.message);
        end
    end
    fprintf('StdErrors (i) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.estsummary=estsummary;
end
