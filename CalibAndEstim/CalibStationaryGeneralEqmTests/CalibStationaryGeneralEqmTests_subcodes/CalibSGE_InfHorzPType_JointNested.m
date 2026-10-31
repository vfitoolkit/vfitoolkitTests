function output=CalibSGE_InfHorzPType_JointNested(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 19: joint vs nested, and the inner general eqm solve of the nested mode (CalibrateBIHAModel_PType).
% Part 9 with two permanent types (see there and part 14). The nested PType objective could never run before 2026-10-29 (its inner solve
% used undefined Params and z_grid), so every nested line here is new ground.
% (a) the nested objective's inner GE solve, started 10% off: prices vs the equilibrium, GE conditions (near zero; loose, as in part 9:
%     fminsearch over three InfHorz GE prices started off the equilibrium stalls on the discrete hours grid)
% (b) its moments vs hand-built at its own GE prices (exact)
% (c) joint at the nested objective's prices: same moments, same GE conditions (exact)
% (d) GE-price constraints reach the inner solver: 'log' vs 'softplus' on Tr (should NOT be zero), each near the equilibrium (loose)
% (e) CalibrateBIHAModel_PType, nested, returns the GE prices in CalibParams
% (f) CalibrateBIHAModel_PType, joint, with a 'log' constrainpositive GE price while the calibrated parameters use 'softplus'
% (g) the type masses calibrated directly (CalibParamNames={'ptypemass'}): CalibParams.ptypemass comes back as the vector over the types

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
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
p_eqm=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
peq=[Params.r;Params.Tr;Params.tau_c];

%% Targets (their values do not matter here, only which moments they name): grouped and per type
TM=struct();
TM.AllStats.K.Mean=1;
TM.AllStats.H.typeA.Mean=1;
TM.AllStats.H.StdDeviation=1;
TM.CrossSectionCovarCorr.H.typeB.Mean=1;

%% Both objectives, set up as CalibrateBIHAModel_PType sets them up
[z_gridvals,pi_zJ,vfoptionsJ]=ExogShockSetup_InfHorz_PType(n_z,z_grid,pi_z,Names_i,Params,vfoptions,3,0);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals=vfoptionsJ.e_gridvals;
simoptionsJ.pi_e=vfoptionsJ.pi_e;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;
simoptionsJ.warnjequaloneptypeasdim=0;
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
[calibparamsvec0,caliboptionsN]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,[0;1;2],CalibParamNames,[1,0;2,0],caliboptionsN,1);
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
[targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorrStats,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSecStats,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,1,simoptions,Names_i);
caliboptions.selectors=selectors;
caliboptionsN.selectors=selectors;
caliboptions.autocorrtimehorizons=autocorrtimehorizons;
caliboptionsN.autocorrtimehorizons=autocorrtimehorizons;
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptionsN.logmoments=caliboptions.logmoments;
caliboptions.weights=ones(length(targetmomentvec),1);
caliboptionsN.weights=caliboptions.weights;
caliboptions.vectoroutput=1;
caliboptionsN.vectoroutput=1;
nmom=length(targetmomentvec);
ParamsOff=Params; % the inner GE solve starts 10% off the equilibrium
ParamsOff.r=1.1*Params.r; ParamsOff.Tr=0.9*Params.Tr; ParamsOff.tau_c=1.1*Params.tau_c;
jointargs={[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats,5,[1,0;2,0;3,0;4,0;5,0],[0;1;2;3;4;5],ones(5,2),zeros(5,1),{}};
nestedargs={CalibParamNames,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,ParamsOff,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats,2,[1,0;2,0],[0;1;2],ones(2,2),zeros(2,1),{}};
GEinner.heteroagentoptions=heteroagentoptions;
GEinner.heteroagentoptions.outputGEstruct=2;
GEinner.z_grid=z_grid;
GEinner.pi_z=pi_z;
GEinner.vfoptions=vfoptions;
GEinner.simoptions=simoptions;

%% (a) the nested objective's inner GE solve
try
    [mn,GEn]=CalibrateBIHAModel_PType_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinner,vfoptionsJ,simoptionsJ);
    mn=gather(mn);
    failed=0;
catch err
    fprintf('JointNested InfHorz PType (a) the nested objective errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    mn=nan(nmom,1); GEn.prices=nan(3,1); GEn.conditions=nan(3,1); failed=1;
end
fprintf('JointNested InfHorz PType (a): nested objective did not run, this should be zero: %.3e; inner GE prices vs the equilibrium, should be near zero (loose: InfHorz fminsearch over three GE prices, started off the equilibrium, stalls on the discrete hours grid): %.3e; its GE conditions, should be near zero (loose: same reason): %.3e \n',failed,max(abs(GEn.prices-peq)),max(abs(GEn.conditions)))

%% (b) its moments against hand-built ones at its own GE prices
ParamsN=Params;
ParamsN.r=GEn.prices(1); ParamsN.Tr=GEn.prices(2); ParamsN.tau_c=GEn.prices(3);
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,ParamsN,DiscountFactorParamNames,vfoptionsJ);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_zJ,ParamsN,simoptionsJ);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,ParamsN,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simoptions_full);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,ParamsN,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simoptionsJ);
hand=gather([AllFull.K.Mean; AllFull.H.typeA.Mean; AllFull.H.StdDeviation; CSFull.H.typeB.Mean]);
if length(mn)==length(hand)
    d=max(abs(mn-hand));
else
    d=1;
end
fprintf('JointNested InfHorz PType (b): nested moments vs hand-built at its own GE prices, this should be zero: %.3e \n',d)

%% (c) joint at the nested objective's prices
calibparamsvecJ=[calibparamsvec0; GEn.prices]; % (no constraints, so the GE prices enter untransformed)
[mj,GEj]=CalibrateBIHAModel_PType_Joint_objectivefn(calibparamsvecJ,jointargs{:},caliboptions,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);
mj=gather(mj);
if length(mj)==nmom+3 && length(mn)==nmom
    d1=max(abs(mj(1:nmom)-mn)); d2=max(abs(mj(nmom+1:end)-GEn.conditions));
else
    d1=1; d2=1;
end
fprintf('JointNested InfHorz PType (c): joint vs nested at the same GE prices, moments, this should be zero: %.3e; GE conditions, this should be zero: %.3e; joint GEout.prices vs nested, this should be zero: %.3e \n',d1,d2,max(abs(GEj.prices-GEn.prices)))

%% (d) constraints on the GE prices reach the inner solver
GEinnerL=GEinner;
GEinnerL.heteroagentoptions.toleranceGEprices=1e-4;
GEinnerL.heteroagentoptions.toleranceGEcondns=1e-4;
GEinnerL.heteroagentoptions.constrainpositive={'Tr'};
GEinnerS=GEinnerL;
GEinnerL.heteroagentoptions.constrainpositivemethod='log';
GEinnerS.heteroagentoptions.constrainpositivemethod='softplus';
try
    [~,GEl]=CalibrateBIHAModel_PType_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerL,vfoptionsJ,simoptionsJ);
    [~,GEs]=CalibrateBIHAModel_PType_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerS,vfoptionsJ,simoptionsJ);
    dl=max(abs(GEl.prices-peq)); ds=max(abs(GEs.prices-peq)); dls=max(abs(GEl.prices-GEs.prices));
    failed=0;
catch err
    fprintf('JointNested InfHorz PType (d) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dl=1; ds=1; dls=0; failed=1;
end
fprintf('JointNested InfHorz PType (d): inner GE with Tr constrainpositive log vs equilibrium, should be near zero (loose: a coarse inner tolerance, and the InfHorz three-price fminsearch stall): %.3e; softplus, should be near zero (loose: same reason): %.3e; did not run, this should be zero: %.3e \n',dl,ds,failed)
fprintf('JointNested InfHorz PType (d): log vs softplus inner GE prices (the constraint reached the inner solver), this should NOT be zero: %.3e \n',dls)

%% (e) CalibrateBIHAModel_PType, nested: the GE prices in CalibParams
TMt=struct();
TMt.AllStats.K.Mean=gather(AllFull.K.Mean); % (at the nested prices, which are the equilibrium to within its tolerance)
TMt.AllStats.H.typeA.Mean=gather(AllFull.H.typeA.Mean);
caliboptionsC=struct();
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
caliboptionsC.jointoptimization=0;
try
    [CalibParamsN,calibsummaryN]=CalibrateBIHAModel_PType(CalibParamNames,TMt,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsC,vfoptions,simoptions);
    missing=sum(~isfield(CalibParamsN,GEPriceParamNames));
    if missing==0
        dp=max(abs([CalibParamsN.r;CalibParamsN.Tr;CalibParamsN.tau_c]-peq));
    else
        dp=1;
    end
    dg=max(abs(calibsummaryN.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('JointNested InfHorz PType (e) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    missing=3; dp=1; dg=1; failed=1;
end
fprintf('JointNested InfHorz PType (e): nested CalibrateBIHAModel_PType, GE prices missing from CalibParams, this should be zero: %.3e; did not run, this should be zero: %.3e \n',missing,failed)
fprintf('JointNested InfHorz PType (e): nested CalibParams GE prices vs the equilibrium, should be near zero: %.3e; calibsummary.GeneralEqmConditions, should be near zero: %.3e \n',dp,dg)

%% (f) joint with a 'log' constrainpositive GE price while the calibrated parameters use 'softplus'
caliboptionsC.jointoptimization=1;
heteroagentoptionsF=struct();
heteroagentoptionsF.constrainpositive={'Tr'};
heteroagentoptionsF.constrainpositivemethod='log';
try
    [CalibParamsF,~]=CalibrateBIHAModel_PType(CalibParamNames,TMt,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsF,caliboptionsC,vfoptions,simoptions);
    dF=max(abs([CalibParamsF.beta-Params.beta, CalibParamsF.varphi-Params.varphi, CalibParamsF.r-Params.r, CalibParamsF.Tr-Params.Tr, CalibParamsF.tau_c-Params.tau_c]));
    failed=0;
catch err
    fprintf('JointNested InfHorz PType (f) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dF=1; failed=1;
end
fprintf('JointNested InfHorz PType (f): joint, Tr constrainpositive by log, calibrated params by softplus: did not run, this should be zero: %.3e; returned its start (the log round trip of Tr), this should be zero: %.3e \n',failed,dF)

%% (g) the type masses calibrated directly
TMg=TMt;
TMg.AllStats.H.typeB.Mean=gather(AllFull.H.typeB.Mean);
try
    [CalibParamsG,~]=CalibrateBIHAModel_PType({'ptypemass'},TMg,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsC,vfoptions,simoptions);
    if isnumeric(CalibParamsG.ptypemass) && numel(CalibParamsG.ptypemass)==length(Names_i)
        formbad=0; moved=max(abs(CalibParamsG.ptypemass(:)-Params.ptypemass(:)));
    else
        formbad=1; moved=1;
    end
    failed=0;
catch err
    fprintf('JointNested InfHorz PType (g) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    formbad=1; moved=1; failed=1;
end
fprintf('JointNested InfHorz PType (g) the type masses calibrated directly: did not run, this should be zero: %.3e; CalibParams.ptypemass not a vector over the types, this should be zero: %.3e; returned its start, this should be zero: %.3e \n',failed,formbad,moved)

output=struct();
output.GEn=GEn;

end
