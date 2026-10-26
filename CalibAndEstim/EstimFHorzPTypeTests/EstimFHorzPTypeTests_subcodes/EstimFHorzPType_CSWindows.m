function output=EstimFHorzPType_CSWindows(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% EstimFHorzPTypeTests part 11 (permanent types: grouped targets, and per type assets_2periodAverage.low.CorrelationWith and earnings_2periodAverage.high.StdDeviation): CrossSectionCovarCorr and AgeConditionalCrossSectionCovarCorr targets of the windows of simoptions.ntimeperiods
% (added 2026-10-14: fnname_2periodAverage.CorrelationWith.(fn2_2periodAverage), ..., one under a conditional restriction, and an
% age-conditional window target by start age alongside a single-period one) through EstimateLifeCycleModel_PType_MethodOfMoments, with the
% AllStats and AgeConditionalStats targets of part 4 for identification.
%   (a) at the truth (skipestimation=1): the objective is zero, the current moments equal the targets, and J equals the test's finite
%       differences at the reported up/down vectors (exact)
%   (b) estimoptions.whichcombos=0: the same J, Sigma and current moments as whichcombos=1 (exact)
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

%% Targets at the truth from full runs
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
simoptions_w=simoptions; simoptions_w.ntimeperiods=2; simoptions_w.pi_z=pi_z;
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
simoptions_w.ntimeperiods=[1,2];
AXFull=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
TM=struct();
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TM.AgeConditionalStats.assets.Mean=t;
TM.CrossSectionCovarCorr.assets_2periodAverage.CorrelationWith.earnings_2periodAverage=gather(CSFull.assets_2periodAverage.CorrelationWith.earnings_2periodAverage);
TM.CrossSectionCovarCorr.assets_2periodAverage.low.CorrelationWith.earnings_2periodAverage=gather(CSFull.assets_2periodAverage.low.CorrelationWith.earnings_2periodAverage); % the low type's window pair
TM.CrossSectionCovarCorr.rich.assets_2periodAverage.CovarianceWith.consumption_2periodAverage=gather(CSFull.rich.assets_2periodAverage.CovarianceWith.consumption_2periodAverage); % under a restriction
t=gather(AXFull.assets.CovarianceWith.consumption); t(1)=NaN; TM.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption=t; % single period
t=gather(AXFull.earnings_2periodAverage.StdDeviation); t(Jr-1:end)=NaN; TM.AgeConditionalCrossSectionCovarCorr.earnings_2periodAverage.StdDeviation=t; % the 2-period windows of the working ages, by start age
t=gather(AXFull.earnings_2periodAverage.high.StdDeviation); t(Jr-1:end)=NaN; TM.AgeConditionalCrossSectionCovarCorr.earnings_2periodAverage.high.StdDeviation=t; % the high type
% the targets as a vector, NaN dropped, in the kinds' order (AllStats, AgeConditionalStats, CrossSec, AgeCrossSec)
tvec=[TM.AllStats.assets.Gini; TM.AgeConditionalStats.assets.Mean(:); ...
    TM.CrossSectionCovarCorr.assets_2periodAverage.CorrelationWith.earnings_2periodAverage; TM.CrossSectionCovarCorr.assets_2periodAverage.low.CorrelationWith.earnings_2periodAverage; TM.CrossSectionCovarCorr.rich.assets_2periodAverage.CovarianceWith.consumption_2periodAverage; ...
    TM.AgeConditionalCrossSectionCovarCorr.assets.CovarianceWith.consumption(:); TM.AgeConditionalCrossSectionCovarCorr.earnings_2periodAverage.StdDeviation(:); TM.AgeConditionalCrossSectionCovarCorr.earnings_2periodAverage.high.StdDeviation(:)];
tvec=tvec(~isnan(tvec));
nT=length(tvec); W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'}; nP=2;
truth=[Params.beta;Params.pension];
estimoptions=struct(); estimoptions.verbose=0;

%% (a) at the truth, (b) whichcombos=0
res=cell(2,1);
for wc=[1,0]
    estimoptions_t=estimoptions; estimoptions_t.skipestimation=1; estimoptions_t.whichcombos=wc;
    estimoptions_t.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
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
        fprintf('CSWindows PType whichcombos=%i: the estimation command errored with: %s \n',wc,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if wc==1
        fprintf('CSWindows PType (a) CrossSection and age-conditional CrossSection window targets (one restricted), %i target entries, at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current moments vs targets, this should be zero: %.3e; J_centered vs finite differences, this should be zero: %.3e \n',nT,1-ran,objv,dmom,maxd)
    else
        fprintf('CSWindows PType (b) the same with estimoptions.whichcombos=0: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; current moments vs targets, this should be zero: %.3e \n',1-ran,objv,dmom)
    end
end
if ~isempty(res{1}) && ~isempty(res{2})
    fprintf('CSWindows PType (b) whichcombos=0 vs =1: J, this should be zero: %.3e; Sigma, this should be zero: %.3e; current moments, this should be zero: %.3e \n',CoreMaxAbs(res{1}.variousmatrices.J(:)-res{2}.variousmatrices.J(:)),CoreMaxAbs(res{1}.variousmatrices.Sigma(:)-res{2}.variousmatrices.Sigma(:)),CoreMaxAbs(res{1}.currentmomentvec(:)-res{2}.currentmomentvec(:)))
else
    fprintf('CSWindows PType (b) whichcombos=0 vs =1: one of the two did not run, this should be zero: %.3e \n',1)
end

%% (c) started away, without the restricted target (at a wrong start the rich may have no mass, and a NaN moment stops lsqnonlin)
TMc=TM;
TMc.CrossSectionCovarCorr=rmfield(TMc.CrossSectionCovarCorr,'rich');
nTc=nT-1;
Wc=eye(nTc); Omegac=diag(0.01*(1+(1:nTc)'/nTc));
Params_away=Params; Params_away.beta=0.90; Params_away.pension=0.30;
estimoptions_c=estimoptions; estimoptions_c.fminalgo=8;
ran=0; perr=NaN; tt=NaN; objv=NaN; rec=[NaN;NaN];
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TMc,Wc,Omegac,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_c,vfoptions,simoptions);
    tt=toc(t0); ran=1; rec=[EstimParams.beta;EstimParams.pension]; perr=CoreMaxAbs(rec-truth); objv=estsummary.objectivefnval;
catch err
    fprintf('CSWindows PType (c): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('CSWindows PType (c) started at beta=0.90, pension=0.30 (%.1f s), recovered %.6f, %.6f, objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),objv,1-ran,perr)

output=struct();
output.TM=TM;

end
