function output=CalibSGE_FHorz_JointNested(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 3: joint vs nested, and the inner general eqm solve of the nested mode (CalibrateOLGModel).
%
% (a) CalibrateOLGModel_Nested_objectivefn (vectoroutput=1) at the truth, its inner GE solve started from equilibrium prices 10% off:
%     the GE prices it reports against the equilibrium, and its GE conditions (both near zero)
% (b) its moments against hand-built ones (VFI, agent dist, AllStats, LifeCycleProfiles) at its own GE prices (exact)
% (c) the joint objective at (truth, the nested objective's GE prices) gives the same moments and the same GE conditions (exact)
% (d) the GE-price constraints reach the inner solver (they used to arrive as 0/1 vectors where names are expected): constrainpositive
%     on Tr with 'log' and with 'softplus', inner tolerance 1e-4, both find the equilibrium (near zero) but along different paths
%     (should NOT be zero); constrain0to1 on r also finds it
% (e) CalibrateOLGModel, nested, returns the GE prices in CalibParams (they were missing), equal to calibsummary's
% (f) CalibrateOLGModel, joint, with heteroagentoptions.constrainpositive on Tr by 'log' while caliboptions uses 'softplus': runs (the
%     per-parameter softplus flags of the GE prices were not concatenated, so it indexed past their end) and returns its start

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

%% The general eqm at the truth
heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
peq=[Params.r;Params.Tr;Params.tau_c];

%% Targets (their values do not matter here, only which moments they name)
TM=struct();
TM.AllStats.K.Mean=1;
TM.AllStats.H.Mean=1;
TM.AllStats.H.StdDeviation=1;
TM.AgeConditionalStats.K.Mean=ones(1,N_j);
TM.AgeConditionalStats.H.Median=ones(1,N_j);

%% Both objectives, set up as CalibrateOLGModel sets them up
[z_gridvals_J,pi_z_J,vfoptionsJ]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,1);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals_J=vfoptionsJ.e_gridvals_J;
simoptionsJ.pi_e_J=vfoptionsJ.pi_e_J;
simoptionsJ.jequaloneDist_usergrids=1;
simoptionsJ.user_z_grid=vfoptionsJ.user_z_grid;
simoptionsJ.user_pi_z=vfoptionsJ.user_pi_z;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;
CalibParamNames={'beta','varphi'};
calibparamsvec0=[Params.beta; Params.varphi];
caliboptionsN=struct();
caliboptionsN.verbose=0;
caliboptionsN.constrainpositive={};
caliboptionsN.constrainpositivemethod='softplus';
caliboptionsN.constrain0to1={};
caliboptionsN.constrainAtoB={};
caliboptionsN.relativeGEweight=10;
caliboptionsN.useCustomModelStats=0;
caliboptionsN.calibrateshocks=0;
caliboptionsN.whichcombos=1;
caliboptionsN.nCalibParamsNonGE=2;
caliboptionsN.metric='sum_squared';
[calibparamsvec0,caliboptionsN]=ParameterConstraints_TransformParamsToUnconstrained(calibparamsvec0,[0;1;2],CalibParamNames,caliboptionsN,1);
heteroagentoptionsJ=struct();
heteroagentoptionsJ.outputGEstruct=2;
heteroagentoptionsJ.multiGEweights=ones(1,3);
heteroagentoptionsJ.constrainpositive={};
heteroagentoptionsJ.constrainpositivemethod='softplus';
heteroagentoptionsJ.constrain0to1={};
heteroagentoptionsJ.constrainAtoB={};
heteroagentoptionsJ.useCustomModelStats=0;
heteroagentoptionsJ.useintermediateEqns=0;
[~,heteroagentoptionsJ]=ParameterConstraints_TransformParamsToUnconstrained(peq,0:1:3,GEPriceParamNames,heteroagentoptionsJ,1);
caliboptions=caliboptionsN;
caliboptions.constrainpositive=[caliboptionsN.constrainpositive; heteroagentoptionsJ.constrainpositive];
caliboptions.constrainAtoB=[caliboptionsN.constrainAtoB; heteroagentoptionsJ.constrainAtoB];
caliboptions.constrain0to1=[caliboptionsN.constrain0to1; heteroagentoptionsJ.constrain0to1];
caliboptions.constrainpositivesoftplus=[caliboptionsN.constrainpositivesoftplus; heteroagentoptionsJ.constrainpositivesoftplus];
caliboptions.constrainAtoBlimits=zeros(5,2);
GEeqnNames=fieldnames(GeneralEqmEqns);
GeneralEqmEqnsCell=cell(1,3);
for gg=1:3
    GeneralEqmEqnParamNames(gg).Names=getAnonymousFnInputNames(GeneralEqmEqns.(GEeqnNames{gg}));
    GeneralEqmEqnsCell{gg}=GeneralEqmEqns.(GEeqnNames{gg});
end
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
caliboptions.selectors=selectors;
caliboptionsN.selectors=selectors;
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptionsN.logmoments=caliboptions.logmoments;
caliboptions.weights=ones(length(targetmomentvec),1);
caliboptionsN.weights=caliboptions.weights;
caliboptions.vectoroutput=1;
caliboptionsN.vectoroutput=1;
nmom=length(targetmomentvec);
ParamsOff=Params; % the inner GE solve starts 10% off the equilibrium
ParamsOff.r=1.1*Params.r; ParamsOff.Tr=0.9*Params.Tr; ParamsOff.tau_c=1.1*Params.tau_c;
jointargs={[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,5,[0;1;2;3;4;5],zeros(5,1),zeros(N_j,1)};
nestedargs={CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,ParamsOff,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqns,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,2,[0;1;2],zeros(2,1),zeros(N_j,1)};
GEinner.heteroagentoptions=heteroagentoptions;
GEinner.heteroagentoptions.outputGEstruct=2;
GEinner.z_grid=z_grid;
GEinner.pi_z=pi_z;
GEinner.vfoptions=vfoptions;
GEinner.simoptions=simoptions;

%% (a) the nested objective's inner GE solve
[mn,GEn]=CalibrateOLGModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinner,vfoptionsJ,simoptionsJ);
mn=gather(mn);
fprintf('JointNested (a): nested inner GE prices vs the equilibrium, should be near zero: %.3e; its GE conditions, should be near zero: %.3e \n',max(abs(GEn.prices-peq)),max(abs(GEn.conditions)))

%% (b) its moments against hand-built ones at its own GE prices
ParamsN=Params;
ParamsN.r=GEn.prices(1); ParamsN.Tr=GEn.prices(2); ParamsN.tau_c=GEn.prices(3);
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,ParamsN,DiscountFactorParamNames,[],vfoptionsJ);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,ParamsN,simoptionsJ);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,ParamsN,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,ParamsN,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
hand=gather([AllFull.K.Mean; AllFull.H.Mean; AllFull.H.StdDeviation; LCPFull.K.Mean(:); LCPFull.H.Median(:)]);
if length(mn)==length(hand)
    d=max(abs(mn-hand));
else
    d=1;
end
fprintf('JointNested (b): nested moments vs hand-built at its own GE prices, this should be zero: %.3e \n',d)

%% (c) joint at the nested objective's prices
calibparamsvecJ=[calibparamsvec0; GEn.prices]; % (no constraints, so the GE prices enter untransformed)
[mj,GEj]=CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},caliboptions,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);
mj=gather(mj);
if length(mj)==nmom+3 && length(mn)==nmom
    d1=max(abs(mj(1:nmom)-mn)); d2=max(abs(mj(nmom+1:end)-GEn.conditions));
else
    d1=1; d2=1;
end
fprintf('JointNested (c): joint vs nested at the same GE prices, moments, this should be zero: %.3e; GE conditions, this should be zero: %.3e; joint GEout.prices vs nested, this should be zero: %.3e \n',d1,d2,max(abs(GEj.prices-GEn.prices)))

%% (d) constraints on the GE prices reach the inner solver
GEinnerL=GEinner;
GEinnerL.heteroagentoptions.toleranceGEprices=1e-4;
GEinnerL.heteroagentoptions.toleranceGEcondns=1e-4;
GEinnerL.heteroagentoptions.constrainpositive={'Tr'};
GEinnerS=GEinnerL;
GEinnerL.heteroagentoptions.constrainpositivemethod='log';
GEinnerS.heteroagentoptions.constrainpositivemethod='softplus';
GEinnerB=GEinnerS;
GEinnerB.heteroagentoptions.constrainpositive={};
GEinnerB.heteroagentoptions.constrain0to1={'r'};
try
    [~,GEl]=CalibrateOLGModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerL,vfoptionsJ,simoptionsJ);
    [~,GEs]=CalibrateOLGModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerS,vfoptionsJ,simoptionsJ);
    [~,GEb]=CalibrateOLGModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerB,vfoptionsJ,simoptionsJ);
    dl=max(abs(GEl.prices-peq)); ds=max(abs(GEs.prices-peq)); db=max(abs(GEb.prices-peq)); dls=max(abs(GEl.prices-GEs.prices));
    failed=0;
catch err
    fprintf('JointNested (d) errored: %s \n',err.message)
    dl=1; ds=1; db=1; dls=0; failed=1;
end
fprintf('JointNested (d): inner GE with Tr constrainpositive log vs equilibrium, should be near zero: %.3e; softplus, should be near zero: %.3e; r constrain0to1, should be near zero: %.3e; did not run, this should be zero: %.3e \n',dl,ds,db,failed)
fprintf('JointNested (d): log vs softplus inner GE prices (the constraint reached the inner solver), this should NOT be zero: %.3e \n',dls)

%% (e) CalibrateOLGModel, nested: the GE prices in CalibParams
% The run stops after fminsearch's initial simplex (tolerances far above any change), so it returns its start, the truth.
TMt=struct();
TMt.AllStats.K.Mean=gather(AllFull.K.Mean); % (at the nested prices, which are the equilibrium to within its tolerance)
TMt.AllStats.H.Mean=gather(AllFull.H.Mean);
caliboptionsC=struct();
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
caliboptionsC.jointoptimization=0;
heteroagentoptionsC=heteroagentoptions;
try
    [CalibParamsN,calibsummaryN]=CalibrateOLGModel(CalibParamNames,TMt,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsC,caliboptionsC,vfoptions,simoptions);
    missing=sum(~isfield(CalibParamsN,GEPriceParamNames));
    if missing==0
        dp=max(abs([CalibParamsN.r;CalibParamsN.Tr;CalibParamsN.tau_c]-peq));
    else
        dp=1;
    end
    dg=max(abs(calibsummaryN.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('JointNested (e) errored: %s \n',err.message)
    missing=3; dp=1; dg=1; failed=1;
end
fprintf('JointNested (e): nested CalibrateOLGModel, GE prices missing from CalibParams, this should be zero: %.3e; did not run, this should be zero: %.3e \n',missing,failed)
fprintf('JointNested (e): nested CalibParams GE prices vs the equilibrium, should be near zero: %.3e; calibsummary.GeneralEqmConditions, should be near zero: %.3e \n',dp,dg)

%% (f) joint with a 'log' constrainpositive GE price while the calibrated parameters use 'softplus'
caliboptionsC.jointoptimization=1;
heteroagentoptionsF=struct();
heteroagentoptionsF.constrainpositive={'Tr'};
heteroagentoptionsF.constrainpositivemethod='log';
try
    [CalibParamsF,~]=CalibrateOLGModel(CalibParamNames,TMt,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsF,caliboptionsC,vfoptions,simoptions);
    dF=max(abs([CalibParamsF.beta-Params.beta, CalibParamsF.varphi-Params.varphi, CalibParamsF.r-Params.r, CalibParamsF.Tr-Params.Tr, CalibParamsF.tau_c-Params.tau_c]));
    failed=0;
catch err
    fprintf('JointNested (f) errored: %s \n',err.message)
    dF=1; failed=1;
end
fprintf('JointNested (f): joint, Tr constrainpositive by log, calibrated params by softplus: did not run, this should be zero: %.3e; returned its start (the log round trip of Tr), this should be zero: %.3e \n',failed,dF)

output=struct();
output.GEn=GEn;

end
