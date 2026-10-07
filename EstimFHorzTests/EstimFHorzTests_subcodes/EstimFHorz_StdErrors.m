function output=EstimFHorz_StdErrors(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 1: the standard-error machinery of EstimateLifeCycleModel_MethodOfMoments at given parameters
% (estimoptions.skipestimation=1 at the truth, so no optimiser runs). Targets generated from the model at the truth; W the identity,
% Omega (CoVarMatrixDataMoments) an arbitrary positive-definite diagonal.
%   (a) ran; EstimParams are the input parameters (skipestimation), objectivefnval zero at the truth (exact)
%   (b) the Jacobian: at each of the four step sizes the reported J_centered, J_up and J_down equal finite differences the test takes
%       itself from direct objective calls at the reported up/down parameter vectors (exact); J is the default step's J_centered
%   (c) Sigma equals inv(J'WJ) J'W Omega W J inv(J'WJ) from the reported J, W, Omega (exact up to roundoff)
%   (d) the standard deviations are sqrt(diag(Sigma)), the confidence intervals are estimate +/- z*sd with z=1.645 (90, the main
%       output), 1.96 (95) and 1 (68) (exact)
%   (e) the sensitivity matrix is -inv(J'WJ) J'W; rank(J) is the number of parameters and the estimate is locally identified
%   (f) efficientW=1 with W=inv(Omega): Sigma is inv(J'WJ), Omega is reported empty
%   (g) must error: a weighting matrix of the wrong size; bootstrapStdErrors=1; iterateGMM=2
%   (h) logmoments with log targets: scalar 1, a per-target-entry vector, a per-target-name vector and the by-name structure each give
%       a zero objective at the truth and a J equal to finite differences of the logged moments (exact)
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
Jr=Params.Jr;

vfoptions=struct(); simoptions=struct();
% The grid interpolation layer makes the moments vary smoothly with the parameters, so finite differences at a relative step of
% 1.5e-4 are meaningful (on the coarse grid alone the moments are piecewise constant and J would be zero)
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);

%% Targets at the truth
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TM.AgeConditionalStats.assets.Mean=t; % (no assets at age 1; an age with no assets at all is omitted too, as its log is -Inf: run 1 of 2026-10-07 had age 2 at zero)
nomit=sum(isnan(t)); % the omitted ages (the target vector keeps them as NaN entries)
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
TM.AllStats.earnings.Mean=gather(AllFull.earnings.Mean);
nT=(N_j-nomit)+2;
W=eye(nT);
Omega=diag(0.01*(1+(1:nT)'/nT)); % any positive-definite covariance of the data moments

EstimParamNames={'beta','pension'};
nP=2;
truth=[Params.beta;Params.pension];
estimoptions=struct();
estimoptions.verbose=0;
estimoptions.skipestimation=1;

%% (a) the run
ran=0; objv=NaN; pmove=NaN;
try
    [EstimParams,EstimParamsConfInts,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);
    ran=1; objv=estsummary.objectivefnval; pmove=CoreMaxAbs([EstimParams.beta;EstimParams.pension]-truth);
catch err
    fprintf('StdErrors (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('StdErrors (a) skipestimation at the truth, %i target entries of which %i ages omitted: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; EstimParams vs the input parameters, this should be zero: %.3e \n',nT,nomit,1-ran,objv,pmove)
if ran==0
    output=struct(); output.ran=0;
    return
end

%% (b) the Jacobian vs the test's own finite differences at the reported up/down parameter vectors
epsnames={'epsilon1','epsilon100','epsilon10000','epsilon1000000'};
maxd=0; bad=0;
for ee=1:4
    dc=estsummary.doublechecks.(epsnames{ee});
    x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
    m0=EstimFHorz_ModelMoments(x0,EstimParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
    Ju=zeros(nT,nP); Jd=zeros(nT,nP); Jc=zeros(nT,nP);
    for pp=1:nP
        xp=x0; xp(pp)=xu(pp); mu=EstimFHorz_ModelMoments(xp,EstimParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
        xp=x0; xp(pp)=xd(pp); md=EstimFHorz_ModelMoments(xp,EstimParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
        Ju(:,pp)=(mu-m0)/(xu(pp)-x0(pp)); Jd(:,pp)=(m0-md)/(x0(pp)-xd(pp)); Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
    end
    if isequal(size(dc.J_centered),[nT,nP])
        maxd=max([maxd; abs(dc.J_centered(:)-Jc(:)); abs(dc.J_up(:)-Ju(:)); abs(dc.J_down(:)-Jd(:)); abs(dc.J(:)-Jc(:))]); % (no constraints, so J is J_centered)
    else
        bad=bad+1;
    end
    if ee==3
        bad=bad+~isequal(estsummary.variousmatrices.J,dc.J); % the default step is epsilon10000 (eedefault=3)
    end
end
fprintf('StdErrors (b) Jacobian at the four step sizes: J_centered, J_up, J_down and J vs the test''s finite differences at the reported up/down parameter vectors, this should be zero: %.3e; wrong sizes or the default J not the epsilon10000 one, this should be zero: %.3e \n',maxd,bad)

%% (c) Sigma from the reported J, W, Omega
J=estsummary.variousmatrices.J;
SigmaHand=((J'*W*J)\(J'*W))*Omega*((J'*W*J)\(J'*W))';
Sigma=estsummary.variousmatrices.Sigma;
fprintf('StdErrors (c) Sigma vs inv(J''WJ) J''W Omega W J inv(J''WJ) from the reported J, W, Omega, relative, this should be zero: %.3e; W and Omega reported as input, this should be zero: %.3e \n',CoreMaxAbs(Sigma(:)-SigmaHand(:))/CoreMaxAbs(SigmaHand(:)),CoreMaxAbs(estsummary.variousmatrices.W(:)-W(:))+CoreMaxAbs(estsummary.variousmatrices.Omega(:)-Omega(:)))

%% (d) standard deviations and confidence intervals
sd=sqrt(diag(Sigma));
d=[abs(estsummary.EstimParamsStdDev.beta-sd(1)); abs(estsummary.EstimParamsStdDev.pension-sd(2))];
d=[d; abs(EstimParamsConfInts.beta(:)-(EstimParams.beta+[-1;1]*1.645*sd(1))); abs(EstimParamsConfInts.pension(:)-(EstimParams.pension+[-1;1]*1.645*sd(2)))];
d=[d; abs(estsummary.confidenceintervals.confint95.EstimParamsConfInts.beta(:)-(EstimParams.beta+[-1;1]*1.96*sd(1))); abs(estsummary.confidenceintervals.confint68.EstimParamsConfInts.pension(:)-(EstimParams.pension+[-1;1]*1*sd(2)))];
fprintf('StdErrors (d) standard deviations vs sqrt(diag(Sigma)), 90/95/68 confidence intervals vs estimate +/- z*sd, this should be zero: %.3e \n',max(d))

%% (e) sensitivity and identification
Lambda=-((J'*W*J)\(J'*W));
fprintf('StdErrors (e) sensitivity matrix vs -inv(J''WJ) J''W, relative, this should be zero: %.3e; rank(J) is %i (hand %i), identified flag true, this should be zero: %.3e \n',CoreMaxAbs(estsummary.sensitivitymatrix(:)-Lambda(:))/CoreMaxAbs(Lambda(:)),estsummary.localidentification.rankJ,nP,abs(estsummary.localidentification.rankJ-nP)+~estsummary.localidentification.yesidentified)

%% (f) efficient weighting
ran=0; d=NaN; omegaempty=NaN;
try
    estimoptions_f=estimoptions; estimoptions_f.efficientW=1;
    Winv=inv(Omega);
    [~,~,estsummary_f]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,Winv,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_f,vfoptions,simoptions);
    ran=1;
    Jf=estsummary_f.variousmatrices.J;
    SigmaHand_f=inv(Jf'*Winv*Jf);
    d=CoreMaxAbs(estsummary_f.variousmatrices.Sigma(:)-SigmaHand_f(:))/CoreMaxAbs(SigmaHand_f(:));
    omegaempty=~isempty(estsummary_f.variousmatrices.Omega);
catch err
    fprintf('StdErrors (f): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('StdErrors (f) efficientW=1 with W=inv(Omega): did not run, this should be zero: %.3e; Sigma vs inv(J''WJ), relative, this should be zero: %.3e; Omega not reported empty, this should be zero: %.3e \n',1-ran,d,omegaempty)

%% (g) must error
errored=0;
try
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,eye(nT+1),Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);
catch
    errored=errored+1;
end
try
    estimoptions_g=estimoptions; estimoptions_g.bootstrapStdErrors=1;
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_g,vfoptions,simoptions);
catch
    errored=errored+1;
end
try
    estimoptions_g=estimoptions; estimoptions_g.iterateGMM=2; estimoptions_g.skipestimation=0; estimoptions_g.fminalgo=8; % (iterated GMM is checked after the estimation step, so estimate; at the truth lsqnonlin stops at once)
    EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_g,vfoptions,simoptions);
catch
    errored=errored+1;
end
fprintf('StdErrors (g) a wrong-size weighting matrix, bootstrapStdErrors=1 and iterateGMM=2 must all error, this should be zero: %.3e \n',3-errored)

%% (h) logmoments, four forms, with log targets
TMlog=TM;
TMlog.AgeConditionalStats.assets.Mean=log(TM.AgeConditionalStats.assets.Mean); % (NaN stays NaN)
TMlog.AllStats.assets.Gini=log(TM.AllStats.assets.Gini);
TMlog.AllStats.earnings.Mean=log(TM.AllStats.earnings.Mean);
TMhalf=TM; TMhalf.AgeConditionalStats.assets.Mean=log(TM.AgeConditionalStats.assets.Mean); % only the age profile logged
% the target vector order is AllStats (Gini, earnings Mean) then AgeConditionalStats (assets Mean, N_j entries incl. the NaN)
forms={1, [1;1;ones(N_j,1)], [0;0;1], struct()};
forms{4}.AgeConditionalStats.assets.Mean=1;
labels={'scalar 1 (all logged)','a vector per target entry (all logged)','a vector per target name [0,0,1] (the age profile logged)','by name (the age profile logged)'};
logvecs={ones(N_j+2,1), ones(N_j+2,1), [0;0;ones(N_j,1)], [0;0;ones(N_j,1)]}; % one entry per target entry (N_j+2 of them), the NaN entries included
for ff=1:4
    estimoptions_h=estimoptions; estimoptions_h.logmoments=forms{ff};
    if ff<=2
        TMh=TMlog;
    else
        TMh=TMhalf;
    end
    ran=0; objv=NaN; maxd=NaN;
    try
        [~,~,estsummary_h]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TMh,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions);
        ran=1; objv=estsummary_h.objectivefnval;
        dc=estsummary_h.doublechecks.epsilon10000;
        x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
        Jc=zeros(nT,nP);
        for pp=1:nP
            xp=x0; xp(pp)=xu(pp); mu=EstimFHorz_ModelMoments(xp,EstimParamNames,TMh,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,logvecs{ff},[]);
            xp=x0; xp(pp)=xd(pp); md=EstimFHorz_ModelMoments(xp,EstimParamNames,TMh,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,logvecs{ff},[]);
            Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
        end
        if isequal(size(dc.J_centered),[nT,nP])
            maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
        else
            maxd=1;
        end
    catch err
        fprintf('StdErrors (h) logmoments %s: the estimation command errored with: %s \n',labels{ff},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('StdErrors (h) logmoments %s, log targets at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; J_centered vs finite differences of the logged moments, this should be zero: %.3e \n',labels{ff},1-ran,objv,maxd)
end

output=struct();
output.ran=1;
output.estsummary=estsummary;

end
