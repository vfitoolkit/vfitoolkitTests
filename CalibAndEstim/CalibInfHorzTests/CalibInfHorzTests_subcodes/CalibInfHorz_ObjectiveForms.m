function output=CalibInfHorz_ObjectiveForms(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,ReturnFn)
% CalibInfHorzTests part 3: the forms of the objective, and the logmoments/weights plumbing of CalibrateInfHorzAgentModel.
%
% First, the objective called directly (as CalibrateInfHorzAgentModel calls it) with TWO calibrated parameters, so that the division
% of the scalar objective by the number of parameters shows, and with targets perturbed away from the model so the objective is not
% zero. The model moment vector m comes from vectoroutput=1 (validated against a hand-built vector in part 2); t is the target
% vector; w the weights. Checks (all exact, against the formulas on m and t):
%   (a) sum_squared:            sum(w.*(m-t).^2)/2
%   (b) sum_logratiosquared:    sum(w.*log(m./t).^2)/2 (on a set of positive targets)
%   (c) MethodOfMoments:        (m-t)'*W*(m-t)/2 with a weight matrix W
%   (d) the lsqnonlin residual vector (vectoroutput=2, weights sqrt(w) as the driver sets them): r(1:end-1)=sqrt(w).*(m-t), r(end)=0
%       (the penalty, zero without constraints), so sum(r(1:end-1).^2) is 2 times the sum_squared objective
%   (e) the same residual for sum_logratiosquared
%   (f) logmoments all ones: the moment vector is log(m) and sum_squared is sum(w.*(log(m)-t).^2)/2 with t the log targets
% Then, through CalibrateInfHorzAgentModel started at the truth with model-generated targets (the optimiser stops where it starts):
%   (g) the objective is zero, calibsummary.currentmomentvec equals calibsummary.targetmomentvec, beta is unchanged
%   (h) logmoments=1 with log targets: zero objective, calibsummary.logmoments all ones, current moments equal the log targets
%   (i) logmoments per target name [1,0,1,0]: calibsummary.logmoments equals the expansion over the entries of each target
%   (j) logmoments by name, a restricted (three-level) AllStats name and a CrossSectionCovarCorr name: the expected vector
%   (k) a non-uniform weights vector of the right length: zero objective; one of the wrong length must error
%   (l) caliboptions.whichcombos=0: the same current moments as =1, exactly
%   (m) MethodOfMoments with a weight matrix W under lsqnonlin, perturbed targets: calibsummary.objvalue equals (m-t)'*W*(m-t) at the
%       returned point (the driver hands lsqnonlin chol(W)); the old elementwise sqrt(W) would not
% The four targets (one AllStats, one restricted AllStats, one AutoCorrTransProbs, one CrossSectionCovarCorr) are all scalars here
% (InfHorz), so the per-target-name and per-entry forms of logmoments coincide in length; (i) still exercises the name form.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w) (1+r)*a+w*z-aprime;
nParams=2;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% Grids as CalibrateInfHorzAgentModel sets them up, one solve, the full stats (to build targets from)
[z_gridvals,pi_zg,vfoptions_g]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0);
simoptions_g=simoptions;
simoptions_g.e_gridvals=vfoptions_g.e_gridvals;
simoptions_g.pi_e=vfoptions_g.pi_e;
vfoptions_g.alreadygridvals=1;
simoptions_g.alreadygridvals=1;
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_g);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_zg,simoptions_g,Params);
simoptions_full=simoptions_g;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_full);
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,simoptions_g);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_g);

% A target set of strictly positive moments (so logs and log-ratios are finite): four target names, one of each kind plus a restricted one
TMpos=struct();
TMpos.AllStats.assets.Mean=gather(AllFull.assets.Mean);
TMpos.AllStats.rich.earnings.Mean=gather(AllFull.rich.earnings.Mean);
TMpos.AutoCorrTransProbs.assets.AutoCorrelation=gather(ACFull.assets.AutoCorrelation);
TMpos.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=gather(CSFull.earnings.CorrelationWith.consumption);
fprintf('ObjectiveForms: the four targets at the truth are %.4f, %.4f, %.4f, %.4f (all must be positive for the log forms); not positive, this should be zero: %.3e \n',TMpos.AllStats.assets.Mean,TMpos.AllStats.rich.earnings.Mean,TMpos.AutoCorrTransProbs.assets.AutoCorrelation,TMpos.CrossSectionCovarCorr.earnings.CorrelationWith.consumption,sum([TMpos.AllStats.assets.Mean,TMpos.AllStats.rich.earnings.Mean,TMpos.AutoCorrTransProbs.assets.AutoCorrelation,TMpos.CrossSectionCovarCorr.earnings.CorrelationWith.consumption]<=0))
% The same names with the targets perturbed away from the model
TMpert=TMpos;
TMpert.AllStats.assets.Mean=1.1*TMpert.AllStats.assets.Mean+0.01;
TMpert.AllStats.rich.earnings.Mean=0.9*TMpert.AllStats.rich.earnings.Mean+0.02;
TMpert.AutoCorrTransProbs.assets.AutoCorrelation=0.95*TMpert.AutoCorrTransProbs.assets.AutoCorrelation+0.01;
TMpert.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=0.9*TMpert.CrossSectionCovarCorr.earnings.CorrelationWith.consumption+0.01;

%% The objective called directly: two parameters
CalibParamNames={'beta','w'};
calibparamsvec0=[Params.beta;Params.w];
calibparamsvecindex=[0;1;2];
calibomitparams_counter=[0;0];
calibomitparamsmatrix={}; % (a cell since 2026-10-08, one omitted-values vector per masked parameter; empty here)
caliboptions=struct();
caliboptions.verbose=0;
caliboptions.constrainpositive={};
caliboptions.constrainpositivemethod='softplus';
caliboptions.constrain0to1={};
caliboptions.constrainAtoB={};
caliboptions.logmoments=0;
caliboptions.metric='sum_squared';
caliboptions.simulatemoments=0;
caliboptions.vectoroutput=1;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=0;
caliboptions.whichcombos=1;
[calibparamsvec0,caliboptions]=ParameterConstraints_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,caliboptions,1);
output=struct();
ran=0;
try
    [targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorr,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSec,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TMpert,FnsToEvaluate,0,simoptions_g);
    ran=1;
catch err
    fprintf('ObjectiveForms: SetupTargetMoments_InfHorz errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('ObjectiveForms: SetupTargetMoments_InfHorz with simoptions, did not run, this should be zero: %.3e \n',1-ran)
if ran==1
    caliboptions.selectors=selectors;
    caliboptions.autocorrtimehorizons=autocorrtimehorizons;
    actual=~isnan(targetmomentvec);
    nT=sum(actual);
    tvec=targetmomentvec(actual);
    w=linspace(0.5,1.5,nT)'; % non-uniform weights
    argcell={CalibParamNames,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zg,ReturnFn,[],Params,DiscountFactorParamNames,[],FnsToEvaluate,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorr,FnsToEvaluate_CrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix};
    ranm=0; m=[];
    try
        % the model moment vector (vectoroutput=1; weights unused)
        caliboptions.weights=w;
        m=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions,vfoptions_g,simoptions_g));
        m=m(:);
        ranm=1;
    catch err
        fprintf('ObjectiveForms: the objective (vectoroutput=1) errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ObjectiveForms: %i targeted entries; the objective with vectoroutput=1 did not run, this should be zero: %.3e \n',nT,1-ranm)
    if ranm==1 && length(m)==nT
        fprintf('ObjectiveForms: the perturbed targets differ from the model by at least %.3e (so the objective is not zero); zero, this should be zero: %.3e \n',min(abs(m-tvec)),min(abs(m-tvec))==0)
        % (a) sum_squared
        caliboptions_a=caliboptions; caliboptions_a.vectoroutput=0; caliboptions_a.metric='sum_squared';
        Obj=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_a,vfoptions_g,simoptions_g));
        handObj=sum(w.*(m-tvec).^2)/nParams;
        fprintf('ObjectiveForms (a) sum_squared vs sum(w.*(m-t).^2)/nParams, relative, this should be zero: %.3e \n',abs(Obj-handObj)/abs(handObj))
        % (b) sum_logratiosquared
        caliboptions_b=caliboptions_a; caliboptions_b.metric='sum_logratiosquared';
        Obj=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_b,vfoptions_g,simoptions_g));
        handObj=sum(w.*(log(m./tvec)).^2)/nParams;
        fprintf('ObjectiveForms (b) sum_logratiosquared vs sum(w.*log(m./t).^2)/nParams, relative, this should be zero: %.3e \n',abs(Obj-handObj)/abs(handObj))
        % (c) MethodOfMoments with a weight matrix
        W=diag(w)+0.1*(ones(nT)-eye(nT))/nT; % symmetric, positive definite (diagonally dominant)
        caliboptions_c=caliboptions_a; caliboptions_c.metric='MethodOfMoments'; caliboptions_c.weights=W;
        Obj=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_c,vfoptions_g,simoptions_g));
        handObj=(m-tvec)'*W*(m-tvec)/nParams;
        fprintf('ObjectiveForms (c) MethodOfMoments vs (m-t)''*W*(m-t)/nParams, relative, this should be zero: %.3e \n',abs(Obj-handObj)/abs(handObj))
        % (d) the residual vector for lsqnonlin (the driver takes the square roots of the weights before creating the handle)
        caliboptions_d=caliboptions; caliboptions_d.vectoroutput=2; caliboptions_d.metric='sum_squared'; caliboptions_d.weights=sqrt(w);
        r=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_d,vfoptions_g,simoptions_g));
        r=r(:);
        if length(r)==nT+1
            fprintf('ObjectiveForms (d) residual vector: length is nT+1, this should be zero: %.3e; r(1:nT) vs sqrt(w).*(m-t), this should be zero: %.3e; r(end) (the penalty, no constraints), this should be zero: %.3e; sum(r(1:nT).^2) vs nParams*sum_squared objective, relative, this should be zero: %.3e \n',0,CoreMaxAbs(r(1:nT)-sqrt(w).*(m-tvec)),abs(r(end)),abs(sum(r(1:nT).^2)-sum(w.*(m-tvec).^2))/sum(w.*(m-tvec).^2))
        else
            fprintf('ObjectiveForms (d) residual vector: length is nT+1, this should be zero: %.3e \n',abs(length(r)-(nT+1)))
        end
        % (e) the residual for sum_logratiosquared
        caliboptions_e=caliboptions_d; caliboptions_e.metric='sum_logratiosquared';
        r=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell{:},caliboptions_e,vfoptions_g,simoptions_g));
        r=r(:);
        if length(r)==nT+1
            fprintf('ObjectiveForms (e) sum_logratiosquared residual r(1:nT) vs sqrt(w).*log(m./t), this should be zero: %.3e \n',CoreMaxAbs(r(1:nT)-sqrt(w).*log(m./tvec)))
        else
            fprintf('ObjectiveForms (e) sum_logratiosquared residual: length is nT+1, this should be zero: %.3e \n',abs(length(r)-(nT+1)))
        end
        % (f) logmoments all ones (the targets are then log targets; the driver would have expanded the scalar 1 to this vector)
        caliboptions_f=caliboptions_a; caliboptions_f.logmoments=ones(length(targetmomentvec),1);
        argcell_f=argcell; argcell_f{19}=log(targetmomentvec); % (the 19th entry is targetmomentvec)
        Obj=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell_f{:},caliboptions_f,vfoptions_g,simoptions_g));
        handObj=sum(w.*(log(m)-log(tvec)).^2)/nParams;
        caliboptions_f1=caliboptions_f; caliboptions_f1.vectoroutput=1;
        mlog=gather(CalibrateInfHorzAgentModel_objectivefn(calibparamsvec0,argcell_f{:},caliboptions_f1,vfoptions_g,simoptions_g));
        fprintf('ObjectiveForms (f) logmoments all ones: moment vector vs log(m), this should be zero: %.3e; sum_squared vs sum(w.*(log(m)-log(t)).^2)/nParams, relative, this should be zero: %.3e \n',CoreMaxAbs(mlog(:)-log(m)),abs(Obj-handObj)/abs(handObj))
        output.m=m;
        output.tvec=tvec;
    elseif ranm==1
        fprintf('ObjectiveForms: the moment vector has %i entries, not nT=%i, this should be zero: %.3e \n',length(m),nT,abs(length(m)-nT))
    end
end

%% Through CalibrateInfHorzAgentModel, started at the truth with the model as the targets (lsqnonlin stops where it starts)
CalibParamNames1={'beta'};
caliboptions_g=struct();
caliboptions_g.verbose=0;
caliboptions_g.fminalgo=8;
nAll=4; % entries of all four targets
labs={'(g) at the truth','(h) logmoments=1, log targets','(i) logmoments per target name [1,0,1,0]','(j) logmoments by name (the restricted AllStats name and the CrossSectionCovarCorr name)','(k) non-uniform weights of the right length','(l) whichcombos=0'};
expected={zeros(nAll,1),ones(nAll,1),[1;0;1;0],[0;1;0;1],zeros(nAll,1),zeros(nAll,1)};
TMs=cell(6,1); opts=cell(6,1);
TMs{1}=TMpos; opts{1}=caliboptions_g;
TMlog=TMpos; % (h)
TMlog.AllStats.assets.Mean=log(TMlog.AllStats.assets.Mean);
TMlog.AllStats.rich.earnings.Mean=log(TMlog.AllStats.rich.earnings.Mean);
TMlog.AutoCorrTransProbs.assets.AutoCorrelation=log(TMlog.AutoCorrTransProbs.assets.AutoCorrelation);
TMlog.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=log(TMlog.CrossSectionCovarCorr.earnings.CorrelationWith.consumption);
TMs{2}=TMlog; opts{2}=caliboptions_g; opts{2}.logmoments=1;
TMi=TMpos; % (i): the first and third targets logged
TMi.AllStats.assets.Mean=log(TMi.AllStats.assets.Mean);
TMi.AutoCorrTransProbs.assets.AutoCorrelation=log(TMi.AutoCorrTransProbs.assets.AutoCorrelation);
TMs{3}=TMi; opts{3}=caliboptions_g; opts{3}.logmoments=[1,0,1,0];
TMj=TMpos; % (j): the second and fourth targets logged, named
TMj.AllStats.rich.earnings.Mean=log(TMj.AllStats.rich.earnings.Mean);
TMj.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=log(TMj.CrossSectionCovarCorr.earnings.CorrelationWith.consumption);
TMs{4}=TMj; opts{4}=caliboptions_g; opts{4}.logmoments=struct(); opts{4}.logmoments.AllStats.rich.earnings.Mean=1; opts{4}.logmoments.CrossSectionCovarCorr.earnings.CorrelationWith.consumption=1;
TMs{5}=TMpos; opts{5}=caliboptions_g; opts{5}.weights=linspace(0.5,1.5,nAll);
TMs{6}=TMpos; opts{6}=caliboptions_g; opts{6}.whichcombos=0;
cur=cell(6,1);
for cc=1:6
    ranc=0; objv=NaN; dmom=NaN; dl=NaN; bmove=NaN;
    try
        [CalibParams,calibsummary]=CalibrateInfHorzAgentModel(CalibParamNames1,TMs{cc},n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,opts{cc},vfoptions,simoptions);
        ranc=1; objv=calibsummary.objvalue; bmove=abs(CalibParams.beta-Params.beta);
        cur{cc}=calibsummary.currentmomentvec(:);
        dmom=CoreMaxAbs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:))+abs(length(calibsummary.currentmomentvec)-nAll);
        if length(calibsummary.logmoments)==nAll
            dl=CoreMaxAbs(calibsummary.logmoments(:)-expected{cc});
        else
            dl=abs(length(calibsummary.logmoments)-nAll);
        end
    catch err
        fprintf('ObjectiveForms %s: CalibrateInfHorzAgentModel errored with: %s \n',labs{cc},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('ObjectiveForms %s: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; calibsummary current vs target moments (%i entries), this should be zero: %.3e; calibsummary.logmoments vs expected, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',labs{cc},1-ranc,objv,nAll,dmom,dl,bmove)
end
% (k) weights of the wrong length must error
errored=0;
caliboptions_bad=caliboptions_g; caliboptions_bad.weights=ones(nAll+1,1);
try
    CalibrateInfHorzAgentModel(CalibParamNames1,TMpos,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_bad,vfoptions,simoptions);
catch
    errored=1;
end
fprintf('ObjectiveForms (k) weights of the wrong length must error, this should be zero: %.3e \n',1-errored)
% (l) whichcombos=0 vs =1: the same current moments
if ~isempty(cur{1}) && ~isempty(cur{6}) && length(cur{1})==length(cur{6})
    fprintf('ObjectiveForms (l) whichcombos=0 vs =1 through CalibrateInfHorzAgentModel: current moments, this should be zero: %.3e \n',CoreMaxAbs(cur{6}-cur{1}))
else
    fprintf('ObjectiveForms (l) whichcombos=0 vs =1: one of the two did not run or the lengths differ, this should be zero: %.3e \n',1)
end
% (m) MethodOfMoments with a weight matrix through the driver under lsqnonlin: the driver hands lsqnonlin R=chol(W) as the weights, so
% the residual is R*(m-t) and lsqnonlin's resnorm (returned as calibsummary.objvalue) is (m-t)'*R'*R*(m-t)=(m-t)'*W*(m-t). Perturbed
% targets so m~=t at the returned point. Under the old elementwise sqrt(W) the resnorm is (m-t)'*S'*S*(m-t) with S=sqrt(W), which
% differs from (m-t)'*W*(m-t) as W is not diagonal: the inverted check confirms this W tells the two apart.
Wm=diag(linspace(0.5,1.5,nAll))+0.1*(ones(nAll)-eye(nAll))/nAll; % as the W of (c): symmetric, positive definite (diagonally dominant)
caliboptions_m=caliboptions_g; caliboptions_m.metric='MethodOfMoments'; caliboptions_m.weights=Wm;
ranm=0;
try
    [~,calibsummary_m]=CalibrateInfHorzAgentModel(CalibParamNames1,TMpert,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],FnsToEvaluate,caliboptions_m,vfoptions,simoptions);
    ranm=1;
catch err
    fprintf('ObjectiveForms (m) MethodOfMoments through CalibrateInfHorzAgentModel errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
if ranm==1
    cm=calibsummary_m.currentmomentvec(:); tm=calibsummary_m.targetmomentvec(:);
    if length(cm)==nAll && length(tm)==nAll
        dm=cm-tm;
        handObj=dm'*Wm*dm;
        Sm=sqrt(Wm);
        oldObj=dm'*(Sm'*Sm)*dm;
        fprintf('ObjectiveForms (m) MethodOfMoments with a weight matrix through CalibrateInfHorzAgentModel (lsqnonlin, perturbed targets): did not run, this should be zero: %.3e; calibsummary.objvalue vs (m-t)''*W*(m-t) at the returned point, relative, this should be zero: %.3e; the old elementwise sqrt(W) resnorm vs (m-t)''*W*(m-t), relative, should NOT be zero: %.3e \n',0,abs(calibsummary_m.objvalue-handObj)/abs(handObj),abs(oldObj-handObj)/abs(handObj))
    else
        fprintf('ObjectiveForms (m) MethodOfMoments through CalibrateInfHorzAgentModel: calibsummary has %i current and %i target entries, not %i, this should be zero: %.3e \n',length(cm),length(tm),nAll,1)
    end
else
    fprintf('ObjectiveForms (m) MethodOfMoments with a weight matrix through CalibrateInfHorzAgentModel: did not run, this should be zero: %.3e \n',1)
end

end
