function output=EstimFHorz_Sensitivity(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 5: the sensitivity to calibrated parameters (Jorgensen 2023) and local identification, at the truth with
% estimoptions.skipestimation=1.
%   (a) estimoptions.CalibParamsNames={'r','w'}: the reported Jcalib_centered, Jcalib_up and Jcalib_down equal finite differences the
%       test takes itself at the reported up/down values of r and w (exact); sensitivitytocalibrationmatrix equals the sensitivity
%       matrix times Jcalib_centered (exact)
%   (b) an unidentified pair: w and pension estimated from the age profile of mean earnings while working, which pension does not
%       enter: rank(J) is 1 of 2 and the identification flag is false; the pension column of J is exactly zero
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
Jr=Params.Jr;
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
estimoptions=struct(); estimoptions.verbose=0; estimoptions.skipestimation=1;

%% (a) sensitivity to the calibrated r and w
TM=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TM.AgeConditionalStats.assets.Mean=t;
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
nT=sum(~isnan(t))+1; W=eye(nT); Omega=diag(0.01*(1+(1:nT)'/nT));
EstimParamNames={'beta','pension'};
CalibParamsNames={'r','w'};
estimoptions_a=estimoptions; estimoptions_a.CalibParamsNames=CalibParamsNames;
ran=0; maxd=NaN; dsens=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions_a,vfoptions,simoptions);
    ran=1;
    dc=estsummary.doublechecks;
    x0=[Params.beta;Params.pension];
    m0=EstimFHorz_ModelMoments(x0,EstimParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
    Jcu=zeros(nT,2); Jcd=zeros(nT,2); Jcc=zeros(nT,2);
    for pp=1:2
        Pu=Params; Pu.(CalibParamsNames{pp})=dc.calibparamsvecup(pp);
        Pd=Params; Pd.(CalibParamsNames{pp})=dc.calibparamsvecdown(pp);
        mu=EstimFHorz_ModelMoments(x0,EstimParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Pu,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
        md=EstimFHorz_ModelMoments(x0,EstimParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Pd,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
        c0=Params.(CalibParamsNames{pp}); cu=dc.calibparamsvecup(pp); cd=dc.calibparamsvecdown(pp);
        Jcu(:,pp)=(mu-m0)/(cu-c0); Jcd(:,pp)=(m0-md)/(c0-cd); Jcc(:,pp)=(mu-md)/(cu-cd);
    end
    if isequal(size(dc.Jcalib),[nT,2])
        maxd=max([abs(dc.Jcalib(:)-Jcc(:)); abs(dc.Jcalib_up(:)-Jcu(:)); abs(dc.Jcalib_down(:)-Jcd(:))]);
        hand=estsummary.sensitivitymatrix*dc.Jcalib;
        dsens=CoreMaxAbs(estsummary.sensitivitytocalibrationmatrix(:)-hand(:))/CoreMaxAbs(hand(:));
    else
        maxd=1; dsens=1;
    end
catch err
    fprintf('Sensitivity (a): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Sensitivity (a) CalibParamsNames r and w at the truth: did not run, this should be zero: %.3e; Jcalib centered/up/down vs finite differences at the reported steps, this should be zero: %.3e; sensitivity-to-calibration matrix vs Lambda*Jcalib, relative, this should be zero: %.3e \n',1-ran,maxd,dsens)

%% (b) an unidentified pair: w and pension from the earnings profile while working
TMb=struct();
t=gather(LCPFull.earnings.Mean); t(Jr:end)=NaN; TMb.AgeConditionalStats.earnings.Mean=t;
nTb=Jr-1; Wb=eye(nTb); Omegab=diag(0.01*(1+(1:nTb)'/nTb));
ran=0; rk=NaN; flag=NaN; colpension=NaN;
try
    [~,~,estsummary]=EstimateLifeCycleModel_MethodOfMoments({'w','pension'},TMb,Wb,Omegab,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);
    ran=1; rk=estsummary.localidentification.rankJ; flag=estsummary.localidentification.yesidentified; colpension=CoreMaxAbs(estsummary.variousmatrices.J(:,2));
catch err
    fprintf('Sensitivity (b): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Sensitivity (b) w and pension from the working-age earnings profile, rank(J) is %i (hand 1): did not run, this should be zero: %.3e; identified flag (should be false), this should be zero: %.3e; rank minus one, this should be zero: %.3e; the pension column of J, this should be zero: %.3e \n',rk,1-ran,double(flag),abs(rk-1),colpension)

output=struct();
output.ran=ran;

end
