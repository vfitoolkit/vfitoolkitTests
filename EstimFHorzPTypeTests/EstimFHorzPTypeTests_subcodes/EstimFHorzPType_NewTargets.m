function output=EstimFHorzPType_NewTargets(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% EstimFHorzPTypeTests part 4: the AutoCorrTransProbs, CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets (including one
% under a conditional restriction) through EstimateLifeCycleModel_MethodOfMoments, alongside AllStats and AgeConditionalStats targets.
%   (a) at the truth (skipestimation=1): the objective is zero, the current moments equal the targets, and J equals the test's finite
%       differences at the reported up/down vectors (exact)
%   (b) estimoptions.whichcombos=0 (every statistic of every targeted function computed): the same J, Sigma and current moments as
%       whichcombos=1 (exact)
%   (c) started away (beta 0.90, pension 0.30), lsqnonlin, on the unrestricted targets: recovered (close)
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+w*kappa_j_pt*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);

%% Targets at the truth from full runs of the five commands
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AXFull=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
TM=struct();
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TM.AgeConditionalStats.assets.Mean=t;
t=gather(ACFull.earnings.AutoCorrelation); t(Jr-1:end)=NaN; TM.AutoCorrTransProbs.earnings.AutoCorrelation=t;              % grouped, the working ages
t=gather(ACFull.earnings.low.AutoCorrelation); t(Jr-1:end)=NaN; TM.AutoCorrTransProbs.earnings.low.AutoCorrelation=t;      % the low type
t=gather(ACFull.rich.assets.high.AutoCorrelation); t(1:5)=NaN; TM.AutoCorrTransProbs.rich.assets.high.AutoCorrelation=t;  % the high type under a restriction
TM.CrossSectionCovarCorr.assets.CorrelationWith.earnings=gather(CSFull.assets.CorrelationWith.earnings);                   % grouped pair
TM.CrossSectionCovarCorr.assets.low.CorrelationWith.earnings=gather(CSFull.assets.low.CorrelationWith.earnings);           % the low type's pair
TM.CrossSectionCovarCorr.rich.assets.CovarianceWith.consumption=gather(CSFull.rich.assets.CovarianceWith.consumption);     % grouped under a restriction
t=gather(AXFull.assets.CovarianceWith.consumption); t(1)=NaN; TM.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption=t;         % grouped
t=gather(AXFull.assets.high.CovarianceWith.consumption); t(1)=NaN; TM.AgeConditionalCrossSectionCovarCorr.assets.high.CovarianceWith.consumption=t; % the high type
% the targets as a vector, NaN dropped, in the kinds' order (AllStats, AgeConditionalStats, AutoCorr, CrossSec, AgeCrossSec), each in its field order
tvec=[TM.AllStats.assets.Gini; TM.AgeConditionalStats.assets.Mean(:); TM.AutoCorrTransProbs.earnings.AutoCorrelation(:); TM.AutoCorrTransProbs.earnings.low.AutoCorrelation(:); TM.AutoCorrTransProbs.rich.assets.high.AutoCorrelation(:); ...
    TM.CrossSectionCovarCorr.assets.CorrelationWith.earnings; TM.CrossSectionCovarCorr.assets.low.CorrelationWith.earnings; TM.CrossSectionCovarCorr.rich.assets.CovarianceWith.consumption; ...
    TM.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption(:); TM.AgeConditionalCrossSectionCovarCorr.assets.high.CovarianceWith.consumption(:)];
tvec=tvec(~isnan(tvec));
nT=length(tvec); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'}; nP=2;
truth=[Params.beta;Params.pension];
estimoptions=struct(); estimoptions.verbose=0;

%% (a) at the truth, (b) whichcombos=0
res=cell(2,1);
for wc=[1,0]
    estimoptions_t=estimoptions; estimoptions_t.skipestimation=1; estimoptions_t.whichcombos=wc;
    ran=0; objv=NaN; dmom=NaN; maxd=NaN;
    try
        [~,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions);
        ran=1; objv=estsummary.objectivefnval;
        if length(estsummary.currentmomentvec)==nT
            dmom=CoreMaxAbs(estsummary.currentmomentvec(:)-tvec);
        else
            dmom=abs(length(estsummary.currentmomentvec)-nT);
        end
        res{2-wc}=estsummary;
        if wc==1
            dc=estsummary.doublechecks.epsilon10000;
            x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
            Jc=zeros(nT,nP);
            for pp=1:nP
                xp=x0; xp(pp)=xu(pp); mu=EstimFHorzPType_ModelMoments(xp,EstimParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
                xp=x0; xp(pp)=xd(pp); md=EstimFHorzPType_ModelMoments(xp,EstimParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
                Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
            end
            if isequal(size(dc.J_centered),[nT,nP])
                maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
            else
                maxd=1;
            end
        end
    catch err
        fprintf('NewTargets PType whichcombos=%i: the estimation command errored with: %s \n',wc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if wc==1
        fprintf('NewTargets PType (a) AutoCorr, CrossSection and age-conditional CrossSection targets (grouped, per type, one restricted), %i target entries, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current moments vs targets, this should be zero: %.3e; J_centered vs finite differences, this should be zero: %.3e \n',nT,1-ran,objv,dmom,maxd)
    else
        fprintf('NewTargets PType (b) the same with estimoptions.whichcombos=0: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current moments vs targets, this should be zero: %.3e \n',1-ran,objv,dmom)
    end
end
if ~isempty(res{1}) && ~isempty(res{2})
    fprintf('NewTargets PType (b) whichcombos=0 vs =1: J, this should be zero: %.3e; Sigma, this should be zero: %.3e; current moments, this should be zero: %.3e \n',CoreMaxAbs(res{1}.variousmatrices.J(:)-res{2}.variousmatrices.J(:)),CoreMaxAbs(res{1}.variousmatrices.Sigma(:)-res{2}.variousmatrices.Sigma(:)),CoreMaxAbs(res{1}.currentmomentvec(:)-res{2}.currentmomentvec(:)))
else
    fprintf('NewTargets PType (b) whichcombos=0 vs =1: one of the two did not run, this should be zero: %.3e \n',1)
end

%% (c) started away (without the restricted targets: at a wrong start the rich may have no mass at some ages, and a NaN moment stops lsqnonlin before it begins, run 1 of 2026-10-07)
TMc=TM;
TMc.AutoCorrTransProbs=rmfield(TMc.AutoCorrTransProbs,'rich');
TMc.CrossSectionCovarCorr=rmfield(TMc.CrossSectionCovarCorr,'rich');
nTc=nT-sum(~isnan(TM.AutoCorrTransProbs.rich.assets.high.AutoCorrelation))-1;
Wc=eye(nTc); Omegac=diag(0.01*(1+(1:nTc)'/nTc));
Params_away=Params; Params_away.beta=0.90; Params_away.pension=0.30;
estimoptions_c=estimoptions; estimoptions_c.fminalgo=8;
ran=0; perr=NaN; tt=NaN; objv=NaN; rec=[NaN;NaN];
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMc,Wc,Omegac,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_c,vfoptions,simoptions);
    tt=toc(t0); ran=1; rec=[EstimParams.beta;EstimParams.pension]; perr=CoreMaxAbs(rec-truth); objv=estsummary.objectivefnval;
catch err
    fprintf('NewTargets PType (c): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('NewTargets PType (c) started at beta=0.90, pension=0.30 (%.1f s), recovered %.6f, %.6f, objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),objv,1-ran,perr)

output=struct();
output.TM=TM;

end
