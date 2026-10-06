function output=CalibFHorzPType_ObjectiveForms(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 3: the forms of the objective, and the logmoments/weights plumbing of CalibrateLifeCycleModel_PType.
%
% First the objective called directly (as the driver calls it) with TWO common parameters and targets perturbed away from the model.
% m is the model moment vector from vectoroutput=1 (validated in part 2), t the targets, w the weights:
%   (a) sum_squared: sum(w.*(m-t).^2)/2   (b) sum_logratiosquared: sum(w.*log(m./t).^2)/2 (positive targets)
%   (c) MethodOfMoments: (m-t)'*W*(m-t)/2   (d) the lsqnonlin residual: r(1:nT)=sqrt(w).*(m-t), r(end)=0
%   (e) the same residual for sum_logratiosquared   (f) logmoments all ones: the moment vector is log(m)
% Then through CalibrateLifeCycleModel_PType started at the truth with model-generated targets (the optimiser stops where it starts):
%   (g) zero objective, calibsummary current moments equal the targets, beta unchanged
%   (h) logmoments=1 with log targets   (i) logmoments per target name [1,0,1,0]
%   (j) logmoments by name, including a four-level name (restriction.fn.type.stat) and a per-type age-conditional name
%   (k) a non-uniform weights vector of the right length; one of the wrong length must error
%   (l) caliboptions.whichcombos=0 gives the same current moments as =1, exactly

N_i=length(Names_i);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.consumption=@(aprime,a,z,r,w,kappa_j_pt,agej,Jr,pension) (1+r)*a+w*kappa_j_pt*z*(agej<Jr)+pension*(agej>=Jr)-aprime;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
Jr=Params.Jr;
nParams=2;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);

%% Grids as the driver sets them up, one solve, the full stats (to build targets from)
[z_gridvals_J,pi_z_J,vfoptions_g]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0);
simoptions_g=simoptions;
simoptions_g.e_gridvals_J=vfoptions_g.e_gridvals_J;
simoptions_g.pi_e_J=vfoptions_g.pi_e_J;
vfoptions_g.alreadygridvals=1;
simoptions_g.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,vfoptions_g);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z_J,Params,simoptions_g);
simoptions_full=simoptions_g;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,simoptions_full);

% A target set of strictly positive moments: four target names in every layout kind
TMpos=struct();
TMpos.AllStats.assets.Mean=gather(AllFull.assets.Mean);                            % grouped
TMpos.AllStats.rich.earnings.low.Mean=gather(AllFull.rich.earnings.low.Mean);      % restriction.fn.type.stat (four levels)
t=gather(LCPFull.earnings.Mean); t(Jr:end)=NaN; TMpos.AgeConditionalStats.earnings.Mean=t;    % grouped profile, NaN in retirement
TMpos.AgeConditionalStats.consumption.high.Mean=gather(LCPFull.consumption.high.Mean);         % per-type profile
TMpert=TMpos;
TMpert.AllStats.assets.Mean=1.1*TMpert.AllStats.assets.Mean+0.01;
TMpert.AllStats.rich.earnings.low.Mean=0.9*TMpert.AllStats.rich.earnings.low.Mean+0.02;
TMpert.AgeConditionalStats.earnings.Mean=1.05*TMpert.AgeConditionalStats.earnings.Mean+0.01;
TMpert.AgeConditionalStats.consumption.high.Mean=0.95*TMpert.AgeConditionalStats.consumption.high.Mean-0.01;

%% The objective called directly: two common parameters
CalibParamNames={'beta','pension'};
calibparamsvec0=[Params.beta;Params.pension];
nCalibParams=2;
nCalibParamsFinder=[1,0;2,0];
calibparamsvecindex=[0;1;2];
calibparamssizes=[1,1;1,1];
calibomitparams_counter=[0;0];
calibomitparamsmatrix=zeros(N_j,1);
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
caliboptions.calibsemiexo=0;
caliboptions.whichcombos=1;
[calibparamsvec0,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,nCalibParamsFinder,caliboptions,1);
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors]=SetupTargetMoments_FHorz(TMpert,FnsToEvaluate,1,N_j,simoptions_g,Names_i);
caliboptions.selectors=selectors;
actual=~isnan(targetmomentvec);
nT=sum(actual);
tvec=targetmomentvec(actual);
w=linspace(0.5,1.5,nT)';
argcell={CalibParamNames,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix};

caliboptions.weights=w;
m=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions,vfoptions_g,simoptions_g));
m=m(:);
fprintf('ObjectiveForms PType: %i targeted entries; the perturbed targets differ from the model by at least %.3e (so the objective is not zero) \n',nT,min(abs(m-tvec)))

caliboptions_a=caliboptions; caliboptions_a.vectoroutput=0; caliboptions_a.metric='sum_squared';
Obj=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_a,vfoptions_g,simoptions_g));
handObj=sum(w.*(m-tvec).^2)/nParams;
fprintf('ObjectiveForms PType (a) sum_squared vs sum(w.*(m-t).^2)/nParams, relative, this should be zero: %.3e \n',abs(Obj-handObj)/abs(handObj))
caliboptions_b=caliboptions_a; caliboptions_b.metric='sum_logratiosquared';
Obj=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_b,vfoptions_g,simoptions_g));
handObj=sum(w.*(log(m./tvec)).^2)/nParams;
fprintf('ObjectiveForms PType (b) sum_logratiosquared vs sum(w.*log(m./t).^2)/nParams, relative, this should be zero: %.3e \n',abs(Obj-handObj)/abs(handObj))
W=diag(w)+0.1*(ones(nT)-eye(nT))/nT;
caliboptions_c=caliboptions_a; caliboptions_c.metric='MethodOfMoments'; caliboptions_c.weights=W;
Obj=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_c,vfoptions_g,simoptions_g));
handObj=(m-tvec)'*W*(m-tvec)/nParams;
fprintf('ObjectiveForms PType (c) MethodOfMoments vs (m-t)''*W*(m-t)/nParams, relative, this should be zero: %.3e \n',abs(Obj-handObj)/abs(handObj))
caliboptions_d=caliboptions; caliboptions_d.vectoroutput=2; caliboptions_d.metric='sum_squared'; caliboptions_d.weights=sqrt(w);
r=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_d,vfoptions_g,simoptions_g));
r=r(:);
if length(r)==nT+1
    dres=max(abs(r(1:nT)-sqrt(w).*(m-tvec))); rend=abs(r(end)); drel=abs(sum(r(1:nT).^2)-sum(w.*(m-tvec).^2))/sum(w.*(m-tvec).^2);
else
    dres=1; rend=1; drel=1;
end
fprintf('ObjectiveForms PType (d) residual vector: length is nT+1, this should be zero: %.3e; r(1:nT) vs sqrt(w).*(m-t), this should be zero: %.3e; r(end) (the penalty, no constraints), this should be zero: %.3e; sum(r(1:nT).^2) vs nParams*sum_squared objective, relative, this should be zero: %.3e \n',abs(length(r)-(nT+1)),dres,rend,drel)
caliboptions_e=caliboptions_d; caliboptions_e.metric='sum_logratiosquared';
r=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell{:},caliboptions_e,vfoptions_g,simoptions_g));
r=r(:);
if length(r)==nT+1
    dres=max(abs(r(1:nT)-sqrt(w).*log(m./tvec)));
else
    dres=1;
end
fprintf('ObjectiveForms PType (e) sum_logratiosquared residual r(1:nT) vs sqrt(w).*log(m./t), this should be zero: %.3e \n',dres)
caliboptions_f=caliboptions_a; caliboptions_f.logmoments=ones(length(targetmomentvec),1);
argcell_f=argcell; argcell_f{22}=log(targetmomentvec); % (the 22nd entry is targetmomentvec)
Obj=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell_f{:},caliboptions_f,vfoptions_g,simoptions_g));
handObj=sum(w.*(log(m)-log(tvec)).^2)/nParams;
caliboptions_f1=caliboptions_f; caliboptions_f1.vectoroutput=1;
mlog=gather(CalibrateLifeCycleModel_PType_objectivefn(calibparamsvec0,argcell_f{:},caliboptions_f1,vfoptions_g,simoptions_g));
fprintf('ObjectiveForms PType (f) logmoments all ones: moment vector vs log(m), this should be zero: %.3e; sum_squared vs sum(w.*(log(m)-log(t)).^2)/nParams, relative, this should be zero: %.3e \n',max(abs(mlog(:)-log(m))),abs(Obj-handObj)/abs(handObj))

%% Through CalibrateLifeCycleModel_PType, started at the truth with the model as the targets
CalibParamNames1={'beta'};
caliboptions_g=struct();
caliboptions_g.verbose=0;
caliboptions_g.fminalgo=8;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMpos,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_g,vfoptions,simoptions);
fprintf('ObjectiveForms PType (g) at the truth, targets kept %i of %i entries: objective, this should be zero: %.3e; calibsummary current vs target moments, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',length(calibsummary.targetmomentvec),nT,calibsummary.objvalue,max(abs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:))),abs(CalibParams.beta-Params.beta))
% (h) logmoments=1 with log targets
TMlog=TMpos;
TMlog.AllStats.assets.Mean=log(TMlog.AllStats.assets.Mean);
TMlog.AllStats.rich.earnings.low.Mean=log(TMlog.AllStats.rich.earnings.low.Mean);
TMlog.AgeConditionalStats.earnings.Mean=log(TMlog.AgeConditionalStats.earnings.Mean);
TMlog.AgeConditionalStats.consumption.high.Mean=log(TMlog.AgeConditionalStats.consumption.high.Mean);
caliboptions_h=caliboptions_g; caliboptions_h.logmoments=1;
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMlog,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_h,vfoptions,simoptions);
nAll=2+2*N_j;
if length(calibsummary.logmoments)==nAll
    dl=max(abs(calibsummary.logmoments(:)-ones(nAll,1)));
else
    dl=abs(length(calibsummary.logmoments)-nAll);
end
fprintf('ObjectiveForms PType (h) logmoments=1, log targets, %i target entries: objective, this should be zero: %.3e; logmoments vector vs all ones, this should be zero: %.3e; current vs target (log) moments, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',nAll,calibsummary.objvalue,dl,max(abs(calibsummary.currentmomentvec(:)-calibsummary.targetmomentvec(:))),abs(CalibParams.beta-Params.beta))
% (i) logmoments per target name [1,0,1,0]
TMi=TMpos;
TMi.AllStats.assets.Mean=log(TMi.AllStats.assets.Mean);
TMi.AgeConditionalStats.earnings.Mean=log(TMi.AgeConditionalStats.earnings.Mean);
caliboptions_i=caliboptions_g; caliboptions_i.logmoments=[1,0,1,0];
expected_i=[1;0;ones(N_j,1);zeros(N_j,1)];
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMi,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_i,vfoptions,simoptions);
if length(calibsummary.logmoments)==nAll
    dl=max(abs(calibsummary.logmoments(:)-expected_i));
else
    dl=abs(length(calibsummary.logmoments)-nAll);
end
fprintf('ObjectiveForms PType (i) logmoments per target name [1,0,1,0]: logmoments vector vs its expansion over the entries, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',dl,calibsummary.objvalue,abs(CalibParams.beta-Params.beta))
% (j) logmoments by name: the four-level AllStats name and the per-type age-conditional name; the other two unmentioned
TMj=TMpos;
TMj.AllStats.rich.earnings.low.Mean=log(TMj.AllStats.rich.earnings.low.Mean);
TMj.AgeConditionalStats.consumption.high.Mean=log(TMj.AgeConditionalStats.consumption.high.Mean);
caliboptions_j=caliboptions_g;
caliboptions_j.logmoments=struct();
caliboptions_j.logmoments.AllStats.rich.earnings.low.Mean=1;
caliboptions_j.logmoments.AgeConditionalStats.consumption.high.Mean=1;
expected_j=[0;1;zeros(N_j,1);ones(N_j,1)];
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMj,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_j,vfoptions,simoptions);
if length(calibsummary.logmoments)==nAll
    dl=max(abs(calibsummary.logmoments(:)-expected_j));
else
    dl=abs(length(calibsummary.logmoments)-nAll);
end
fprintf('ObjectiveForms PType (j) logmoments by name (a four-level restricted name and a per-type age-conditional name): logmoments vector vs expected, this should be zero: %.3e; objective, this should be zero: %.3e; beta moved, this should be zero: %.3e \n',dl,calibsummary.objvalue,abs(CalibParams.beta-Params.beta))
% (k) weights
caliboptions_k=caliboptions_g; caliboptions_k.weights=linspace(0.5,1.5,nT);
[CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMpos,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_k,vfoptions,simoptions);
errored=0;
caliboptions_bad=caliboptions_g; caliboptions_bad.weights=ones(nT+1,1);
try
    CalibrateLifeCycleModel_PType(CalibParamNames1,TMpos,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_bad,vfoptions,simoptions);
catch
    errored=1;
end
fprintf('ObjectiveForms PType (k) non-uniform weights of the right length: objective, this should be zero: %.3e; beta moved, this should be zero: %.3e; weights of the wrong length must error, this should be zero: %.3e \n',calibsummary.objvalue,abs(CalibParams.beta-Params.beta),1-errored)
% (l) whichcombos=0 gives the same moments
caliboptions_l=caliboptions_g; caliboptions_l.whichcombos=0;
[~,calibsummary_l]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMpos,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_l,vfoptions,simoptions);
[~,calibsummary_1]=CalibrateLifeCycleModel_PType(CalibParamNames1,TMpos,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions_g,vfoptions,simoptions);
fprintf('ObjectiveForms PType (l) whichcombos=0 vs =1 through CalibrateLifeCycleModel_PType: current moments, this should be zero: %.3e; objectives, this should be zero: %.3e \n',max(abs(calibsummary_l.currentmomentvec(:)-calibsummary_1.currentmomentvec(:))),abs(calibsummary_l.objvalue-calibsummary_1.objvalue))

output=struct();
output.m=m;
output.tvec=tvec;

end
