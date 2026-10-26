function m=EstimFHorz_ModelMoments(paramvec,EstimParamNames,TargetMoments,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,logmoments,omitestimparam,ParametrizeParamsFn,forcecalibrateshocks)
% EstimFHorzTests: the vector of model moments at the targeted entries for a parameter vector, by the same route as the estimation
% command's Jacobian handle: SetupTargetMoments_FHorz, the grids as the drivers set them up, then CalibrateLifeCycleModel_objectivefn
% with vectoroutput=1 and no constraints (the Jacobian is taken on the model parameters). paramvec is the stacked vector of the
% estimated parameters in EstimParamNames order (an omitted parameter contributes only its NaN entries of omitestimparam).
% logmoments is a vector with one entry per target entry (NaN entries included), or 0. omitestimparam is a structure as
% estimoptions.omitestimparam, or []. ParametrizeParamsFn (optional) is applied inside the objective, as the drivers do. When a shock
% function (vfoptions.ExogShockFn / EiidShockFn) takes an estimated parameter the grids are rebuilt inside the objective at every call
% (calibrateshocks=1), as the drivers do; when jequaloneDist is a function the user's grids are kept for it (jequaloneDist_usergrids).
if nargin<22
    ParametrizeParamsFn=[];
end
vfoptions_m=vfoptions; simoptions_m=simoptions;
calibrateshocks=0;
if isfield(vfoptions_m,'ExogShockFn') && ~isempty(intersect(getAnonymousFnInputNames(vfoptions_m.ExogShockFn),EstimParamNames))
    calibrateshocks=1;
end
if isfield(vfoptions_m,'EiidShockFn') && ~isempty(intersect(getAnonymousFnInputNames(vfoptions_m.EiidShockFn),EstimParamNames))
    calibrateshocks=1;
end
if nargin>=23 && ~isempty(forcecalibrateshocks)
    % Caller states whether the shocks depend on the estimated parameters. Needed when the dependence runs through parameters
    % DERIVED by ParametrizeParamsFn (part 9): this helper is the test's independent reference, so it must not share the drivers'
    % detection machinery (ParseParametrizeParamsFn) -- the test knows the truth and forces it.
    calibrateshocks=forcecalibrateshocks;
end
if ~isfield(simoptions_m,'jequaloneDist_usergrids')
    simoptions_m.jequaloneDist_usergrids=1;
end
if calibrateshocks==0
    KeepOriginalGrid=(simoptions_m.jequaloneDist_usergrids==1);
    [z_gridvals_J,pi_z_J,vfoptions_m]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions_m,3,KeepOriginalGrid);
    simoptions_m.e_gridvals_J=vfoptions_m.e_gridvals_J;
    simoptions_m.pi_e_J=vfoptions_m.pi_e_J;
    if simoptions_m.jequaloneDist_usergrids==1
        simoptions_m.user_z_grid=vfoptions_m.user_z_grid;
        simoptions_m.user_pi_z=vfoptions_m.user_pi_z;
    end
else
    z_gridvals_J=z_grid; pi_z_J=pi_z; % passed through as the drivers do: the objective rebuilds the grids that depend on the current parameters, and with only an EiidShockFn these ARE the z grids (empty placeholders broke the setup, 2026-10-07)
end
vfoptions_m.alreadygridvals=1;
simoptions_m.alreadygridvals=1;
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TargetMoments,FnsToEvaluate,0,N_j,simoptions_m);
nP=length(EstimParamNames);
calibparamsvecindex=zeros(nP+1,1);
calibomitparams_counter=zeros(nP,1);
calibomitparamsmatrix=zeros(N_j,1);
for pp=1:nP
    if ~isempty(omitestimparam) && isfield(omitestimparam,EstimParamNames{pp})
        mask=omitestimparam.(EstimParamNames{pp}); mask=mask(:);
        calibomitparams_counter(pp)=1;
        calibomitparamsmatrix(:,sum(calibomitparams_counter))=mask;
        calibparamsvecindex(pp+1)=calibparamsvecindex(pp)+sum(isnan(mask));
    else
        calibparamsvecindex(pp+1)=calibparamsvecindex(pp)+numel(Params.(EstimParamNames{pp}));
    end
end
caliboptions=struct();
caliboptions.verbose=0;
caliboptions.constrainpositive={};
caliboptions.constrainpositivemethod='softplus';
caliboptions.constrain0to1={};
caliboptions.constrainAtoB={};
caliboptions.metric='sum_squared';
caliboptions.simulatemoments=0;
caliboptions.vectoroutput=1;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=calibrateshocks;
caliboptions.whichcombos=1;
caliboptions.selectors=selectors;
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
if isscalar(logmoments)
    caliboptions.logmoments=logmoments*ones(length(targetmomentvec),1);
else
    caliboptions.logmoments=reshape(logmoments,[],1);
end
[~,caliboptions]=ParameterConstraints_TransformParamsToUnconstrained(paramvec(:),calibparamsvecindex,EstimParamNames,caliboptions,1); % (no constraints: this just turns the empty name lists into zero vectors)
m=gather(CalibrateLifeCycleModel_objectivefn(paramvec(:),EstimParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,[],Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,ParametrizeParamsFn,FnsToEvaluate,[],usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix,caliboptions,vfoptions_m,simoptions_m));
m=m(:);

end
