function output=EstimFHorzPTypeMLE_StdErrors(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist,vfoptions,simoptions)
% doPart(2): EstimateLifeCycleModel_PType_MaximumLikelihood at the truth with estimoptions.skipestimation=1, every reported number
% recomputed from direct solves + LogLikelihood_FHorz_PType (EstimFHorzPTypeMLE_LogLiAt). Estimated: beta (common to the types), the
% s.d. of the consumption measurement error by type (a structure over the types, so two blocks), and pshare_high, the mass of the
% second type, which ParametrizeParamsFn (EstimFHorzPTypeMLE_Parametrize) turns into the type masses. Identities throughout:
%   (a) LogL, LogLi, PTypePosterior and the type masses at the estimate vs direct
%   (b) the scores vs direct calls at the reported up/down parameter vectors
%   (c) OPG, Sigma, the standard errors and confidence intervals in the per-type layout (a structure over the types)
%   (d) the Hessian: the points with both parameters stepped vs direct, and the entries from the reported values
%   (e) AIC and BIC with four parameters
%   (f) the per-type parameter given as a matrix (N_i as its first dimension): the same standard errors, as a matrix
%   (g) refusals
% The data: one panel (N=400) simulated by hand from the model at the truth (EstimFHorzPTypeMLE_SimulatePTypePanel, fixed seed).

Params.sig_c.(Names_i{1})=0.1; Params.sig_c.(Names_i{2})=0.1; Params.sig_y=0.05; Params.sig_a=0.2; % the truth (measurement error)
Params.pshare_high=Params.(PTypeDistParamNames{1})(2);
FnsToEvaluate.cons=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+(agej<Jr)*w*kappa_j_pt*z+(agej>=Jr)*pension-aprime;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt,agej,Jr) (agej<Jr)*w*kappa_j_pt*z;
FnsToEvaluate.assets=@(aprime,a,z) a;
MeasurementError.cons.distribution='normal'; MeasurementError.cons.sigma='sig_c';
MeasurementError.earnings.distribution='normal'; MeasurementError.earnings.sigma='sig_y';
MeasurementError.assets.distribution='normal'; MeasurementError.assets.sigma='sig_a';

[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
Data=EstimFHorzPTypeMLE_SimulatePTypePanel(Policy,n_a,n_z,N_j,Names_i,a_grid,z_grid,pi_z,jequaloneDist,Params,vfoptions.ngridinterp,400,20261017);
N=size(Data.cons,1);

EstimParamNames={'beta','sig_c','pshare_high'};
estimoptions.verbose=0;
estimoptions.skipestimation=1;
estimoptions.constrainpositive={'sig_c'};
estimoptions.constrain0to1={'pshare_high'};
ParametrizeParamsFn=@(Params) EstimFHorzPTypeMLE_Parametrize(Params);
[EstimParams,EstimParamsConfInts,estsummary]=EstimateLifeCycleModel_PType_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions,vfoptions,simoptions);
K=4;
x=[Params.beta; Params.sig_c.(Names_i{1}); Params.sig_c.(Names_i{2}); Params.pshare_high];
LogLiAt=@(vec) EstimFHorzPTypeMLE_LogLiAt(vec,Params,Data,MeasurementError,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,vfoptions,simoptions);

dc=estsummary.doublechecks.epsilon10000; % the default step (eedefault=3)
xe=dc.estimparamsvec; % the point the estimation evaluated at: the truth after the round trip through the constraint transforms

%% (a)
[LogLi0,Post0]=LogLiAt(xe);
fprintf('PTypeStdErrors (a) LogL vs direct, this should be zero: %.3e; LogLi, this should be zero: %.3e; PTypePosterior, this should be zero: %.3e \n',abs(estsummary.LogL-sum(LogLi0)),CoreMaxAbs(estsummary.LogLi-LogLi0),CoreMaxAbs(estsummary.PTypePosterior-gather(Post0)))
fprintf('PTypeStdErrors (a) ptypemasses vs the truth, this should be zero: %.3e; the parameter vector of the estimation vs the truth in block order, this should be zero: %.3e \n',CoreMaxAbs(estsummary.ptypemasses-Params.(PTypeDistParamNames{1})(:)),CoreMaxAbs(estsummary.doublechecks.epsilon10000.estimparamsvec-x))

%% (b) the scores
LogLi_up=zeros(N,K); LogLi_down=zeros(N,K);
for pp=1:K
    vec=xe; vec(pp)=dc.estimparamsvecup(pp);
    LogLi_up(:,pp)=LogLiAt(vec);
    vec=xe; vec(pp)=dc.estimparamsvecdown(pp);
    LogLi_down(:,pp)=LogLiAt(vec);
end
Shand=(LogLi_up-LogLi_down)./((dc.estimparamsvecup-dc.estimparamsvecdown)');
fprintf('PTypeStdErrors (b) LogLi at the up steps vs direct, this should be zero: %.3e; at the down steps, this should be zero: %.3e \n',CoreMaxAbs(dc.LogLi_up-LogLi_up),CoreMaxAbs(dc.LogLi_down-LogLi_down))
fprintf('PTypeStdErrors (b) Scores vs centered differences of direct calls, this should be zero: %.3e (relative to the largest score %.3e) \n',CoreMaxAbs(estsummary.variousmatrices.Scores-Shand),max(abs(Shand(:))))

%% (c) OPG and the per-type layout
S=estsummary.variousmatrices.Scores;
sd=sqrt(diag(inv(S'*S)));
fprintf('PTypeStdErrors (c) OPG vs Scores''*Scores, this should be zero: %.3e; Sigma_OPG vs inv(OPG), this should be zero: %.3e \n',CoreMaxAbs(estsummary.variousmatrices.OPG-S'*S),CoreMaxAbs(estsummary.variousmatrices.Sigma_OPG-inv(S'*S)))
sdrep=[estsummary.EstimParamsStdDev.beta; estsummary.EstimParamsStdDev.sig_c.(Names_i{1}); estsummary.EstimParamsStdDev.sig_c.(Names_i{2}); estsummary.EstimParamsStdDev.pshare_high];
fprintf('PTypeStdErrors (c) EstimParamsStdDev (sig_c a structure over the types) vs sqrt(diag(inv(OPG))), this should be zero: %.3e; vs EstimParamsStdDev_OPG, this should be zero: %.3e \n',CoreMaxAbs(sdrep-sd),CoreMaxAbs(sdrep-[estsummary.EstimParamsStdDev_OPG.beta; estsummary.EstimParamsStdDev_OPG.sig_c.(Names_i{1}); estsummary.EstimParamsStdDev_OPG.sig_c.(Names_i{2}); estsummary.EstimParamsStdDev_OPG.pshare_high]))
cirep=[EstimParamsConfInts.beta; EstimParamsConfInts.sig_c.(Names_i{1}); EstimParamsConfInts.sig_c.(Names_i{2}); EstimParamsConfInts.pshare_high];
fprintf('PTypeStdErrors (c) confidence intervals vs estimate -/+ z sd, this should be zero: %.3e; EstimParams.sig_c by type vs the truth, this should be zero: %.3e \n',CoreMaxAbs(cirep-[x-1.645*sd, x+1.645*sd]),CoreMaxAbs([EstimParams.sig_c.(Names_i{1});EstimParams.sig_c.(Names_i{2})]-x(2:3)))
fprintf('PTypeStdErrors (c) rank(OPG) is full, this should be one: %i \n',double(estsummary.localidentification.yesidentified))
fprintf('PTypeStdErrors (c) standard errors (OPG): beta %.5f, sig_c by type %.5f and %.5f, pshare_high %.5f (informational) \n',sd(1),sd(2),sd(3),sd(4))

%% (d) the Hessian
hs=estsummary.doublechecks.hessian;
vec=xe; vec(1)=dc.estimparamsvecup(1); vec(4)=dc.estimparamsvecup(4);
Lupup=sum(LogLiAt(vec));
vec=xe; vec(1)=dc.estimparamsvecdown(1); vec(4)=dc.estimparamsvecdown(4);
Ldowndown=sum(LogLiAt(vec));
fprintf('PTypeStdErrors (d) the points with beta and pshare_high both stepped: up vs direct, this should be zero: %.3e; down vs direct, this should be zero: %.3e \n',abs(hs.LogL_upup(1,4)-Lupup),abs(hs.LogL_downdown(1,4)-Ldowndown))
L0=sum(LogLi0); Lup=sum(LogLi_up,1)'; Ldown=sum(LogLi_down,1)';
H14=(Lupup-Lup(1)-Lup(4)+2*L0-Ldown(1)-Ldown(4)+Ldowndown)/(2*hs.hup(1)*hs.hup(4));
H22=(Lup(2)-2*L0+Ldown(2))/(hs.hup(2)*hs.hdown(2));
fprintf('PTypeStdErrors (d) Hessian off-diagonal entry for beta and pshare_high vs by hand, this should be zero: %.3e; diagonal entry for sig_c of the first type, this should be zero: %.3e (relative to %.3e) \n',abs(estsummary.variousmatrices.Hessian(1,4)-H14),abs(estsummary.variousmatrices.Hessian(2,2)-H22),abs(H22))
SH=inv(-estsummary.variousmatrices.Hessian);
fprintf('PTypeStdErrors (d) Sigma_sandwich vs inv(-H)*OPG*inv(-H), this should be zero: %.3e \n',CoreMaxAbs(estsummary.variousmatrices.Sigma_sandwich-SH*(S'*S)*SH))

%% (e) information criteria
fprintf('PTypeStdErrors (e) AIC vs its formula, this should be zero: %.3e; BIC vs its formula, this should be zero: %.3e \n',abs(estsummary.AIC-(2*K-2*L0)),abs(estsummary.BIC-(K*log(N)-2*L0)))

%% (f) the per-type parameter as a matrix
Params_f=Params;
Params_f.sig_c=[Params.sig_c.(Names_i{1}); Params.sig_c.(Names_i{2})];
[EP_f,CI_f,es_f]=EstimateLifeCycleModel_PType_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_f,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,estimoptions,vfoptions,simoptions);
fprintf('PTypeStdErrors (f) matrix form: LogL vs the structure form, this should be zero: %.3e; std devs of sig_c (a matrix) vs the structure form, this should be zero: %.3e; EstimParams.sig_c size vs a column over the types, this should be zero: %i \n', ...
    abs(es_f.LogL-estsummary.LogL),CoreMaxAbs(es_f.EstimParamsStdDev.sig_c(:)-sdrep(2:3)),sum(abs(size(EP_f.sig_c)-[2,1])))
fprintf('PTypeStdErrors (f) matrix form: confidence intervals by type vs the structure form, this should be zero: %.3e \n',CoreMaxAbs([CI_f.sig_c.(Names_i{1});CI_f.sig_c.(Names_i{2})]-cirep(2:3,:)))

%% (g) refusals
errcases={'fminalgo=8','unknown stderrors'};
errmsgs={'estimoptions.fminalgo must be 4','estimoptions.stderrors must be'};
for ee=1:length(errcases)
    esto=estimoptions;
    switch ee
        case 1
            esto.fminalgo=8;
        case 2
            esto.stderrors='bootstrap';
    end
    errok=0;
    try
        EstimateLifeCycleModel_PType_MaximumLikelihood(EstimParamNames,Data,MeasurementError,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,esto,vfoptions,simoptions);
    catch err
        errok=contains(err.message,errmsgs{ee});
        if ~errok
            fprintf('   (the error was: %s) \n',err.message);
        end
    end
    fprintf('PTypeStdErrors (g) %s errors with the intended message, this should be zero: %i \n',errcases{ee},1-errok);
end

output.estsummary=estsummary;
end
