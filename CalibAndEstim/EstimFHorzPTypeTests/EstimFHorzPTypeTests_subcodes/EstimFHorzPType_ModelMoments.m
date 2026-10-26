function m=EstimFHorzPType_ModelMoments(paramvec,EstimParamNames,TargetMoments,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,logmoments,omitestimparam,ParametrizeParamsFn)
% EstimFHorzPTypeTests: the vector of model moments at the targeted entries for a parameter vector, by the same route as the PType
% estimation command's Jacobian handle: SetupTargetMoments_FHorz (useptype=1), the grids as the driver sets them up, then
% CalibrateLifeCycleModel_PType_objectivefn with vectoroutput=1 and no constraints. paramvec stacks the estimated parameters in
% EstimParamNames order, a per-type parameter (a structure over the types, or a matrix with N_i as one dimension) contributing one block
% per type in Names_i order (an omitted parameter contributes only its NaN entries of omitestimparam, the same mask for every type).
% logmoments is a vector with one entry per target entry (NaN entries included), or 0. omitestimparam is a structure as
% estimoptions.omitestimparam, or []. ParametrizeParamsFn (optional) is applied inside the objective. When a shock function
% (vfoptions.ExogShockFn / EiidShockFn, possibly a structure over the types) takes an estimated parameter the grids are rebuilt inside
% the objective at every call (calibrateshocks=1), as the driver does; when jequaloneDist is a function the user's grids are kept for it.
if nargin<24
    ParametrizeParamsFn=[];
end
N_i=length(Names_i);
vfoptions_m=vfoptions; simoptions_m=simoptions;
calibrateshocks=0;
for fnm={'ExogShockFn','EiidShockFn'}
    if isfield(vfoptions_m,fnm{1})
        F=vfoptions_m.(fnm{1});
        if isstruct(F)
            fns=fieldnames(F);
            for ii=1:length(fns)
                if ~isempty(intersect(getAnonymousFnInputNames(F.(fns{ii})),EstimParamNames))
                    calibrateshocks=1;
                end
            end
        elseif ~isempty(intersect(getAnonymousFnInputNames(F),EstimParamNames))
            calibrateshocks=1;
        end
    end
end
if ~isfield(simoptions_m,'jequaloneDist_usergrids')
    simoptions_m.jequaloneDist_usergrids=1;
end
if calibrateshocks==0
    KeepOriginalGrid=(simoptions_m.jequaloneDist_usergrids==1);
    [z_gridvals_J,pi_z_J,vfoptions_m]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions_m,3,KeepOriginalGrid);
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
if ~isfield(simoptions_m,'warnjequaloneptypeasdim')
    simoptions_m.warnjequaloneptypeasdim=0;
end
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TargetMoments,FnsToEvaluate,1,N_j,simoptions_m,Names_i);
% the parameters as the PType driver lays them out: a matrix with N_i as a dimension becomes a structure over the types
nCalibParams=0; nCalibParamsFinder=zeros(0,2);
for pp=1:length(EstimParamNames)
    P=Params.(EstimParamNames{pp});
    if isstruct(P)
        for ii=1:N_i
            if isfield(P,Names_i{ii})
                nCalibParams=nCalibParams+1; nCalibParamsFinder(nCalibParams,:)=[pp,ii];
            end
        end
    elseif any(size(P)==N_i)
        if size(P,1)==N_i
            P=P';
        end
        Params=rmfield(Params,EstimParamNames{pp});
        for ii=1:N_i
            nCalibParams=nCalibParams+1; nCalibParamsFinder(nCalibParams,:)=[pp,ii];
            Params.(EstimParamNames{pp}).(Names_i{ii})=P(:,ii);
        end
    else
        nCalibParams=nCalibParams+1; nCalibParamsFinder(nCalibParams,:)=[pp,0];
    end
end
calibparamsvecindex=zeros(nCalibParams+1,1);
calibparamssizes=zeros(nCalibParams,2);
calibomitparams_counter=zeros(nCalibParams,1);
calibomitparamsmatrix=zeros(N_j,1);
for pp=1:nCalibParams
    if nCalibParamsFinder(pp,2)==0
        cur=Params.(EstimParamNames{nCalibParamsFinder(pp,1)});
    else
        cur=Params.(EstimParamNames{nCalibParamsFinder(pp,1)}).(Names_i{nCalibParamsFinder(pp,2)});
    end
    calibparamssizes(pp,:)=size(cur);
    if ~isempty(omitestimparam) && isfield(omitestimparam,EstimParamNames{nCalibParamsFinder(pp,1)})
        mask=omitestimparam.(EstimParamNames{nCalibParamsFinder(pp,1)}); mask=mask(:);
        calibomitparams_counter(pp)=1;
        calibomitparamsmatrix(:,sum(calibomitparams_counter))=mask;
        calibparamsvecindex(pp+1)=calibparamsvecindex(pp)+sum(isnan(mask));
    else
        calibparamsvecindex(pp+1)=calibparamsvecindex(pp)+numel(cur);
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
caliboptions.calibsemiexo=0;
caliboptions.whichcombos=1;
caliboptions.selectors=selectors;
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
if isscalar(logmoments)
    caliboptions.logmoments=logmoments*ones(length(targetmomentvec),1);
else
    caliboptions.logmoments=reshape(logmoments,[],1);
end
[~,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(paramvec(:),calibparamsvecindex,EstimParamNames,nCalibParamsFinder,caliboptions,1); % (no constraints: the empty name lists become zero vectors)
m=gather(CalibrateLifeCycleModel_PType_objectivefn(paramvec(:),EstimParamNames,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,ParametrizeParamsFn,FnsToEvaluate,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,nCalibParams,nCalibParamsFinder,calibparamsvecindex,calibparamssizes,calibomitparams_counter,calibomitparamsmatrix,caliboptions,vfoptions_m,simoptions_m));
m=m(:);

end
