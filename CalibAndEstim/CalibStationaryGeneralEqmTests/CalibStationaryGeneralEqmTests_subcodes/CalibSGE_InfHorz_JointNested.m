function output=CalibSGE_InfHorz_JointNested(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 9: joint vs nested, and the inner general eqm solve of the nested mode (CalibrateBIHAModel).
% The InfHorz counterpart of part 3 (see there). The nested objective could never run before (undefined Params and z_grid, and the
% cell of GE eqns where the structure is needed), so every nested line here is new ground.
% (a) the nested objective's inner GE solve, started 10% off: prices vs the equilibrium, GE conditions (near zero)
%     [loose since 2026-10-14: with 101 hours points the truth equilibrium reaches 5.5e-5, but fminsearch over three InfHorz GE prices
%     started 10% off stops near 2.5e-3 (that is the GE solver, CoreStationaryGeneralEqm's territory, not the calibration); likewise (d)]
% (b) its moments vs hand-built at its own GE prices (exact)
% (c) joint at the nested objective's prices: same moments, same GE conditions (exact)
% (d) GE-price constraints reach the inner solver: 'log' vs 'softplus' on Tr (should NOT be zero), each near the equilibrium
% (e) CalibrateBIHAModel, nested, returns the GE prices in CalibParams
% (f) CalibrateBIHAModel, joint, with a 'log' constrainpositive GE price while the calibrated parameters use 'softplus'
% (g) runtime of one nested and one joint evaluation (informational)

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
FnsToEvaluate.earnings=@(d,aprime,a,z,r,alpha,delta,A) ((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d;
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

%% The general eqm at the truth
heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;
p_eqm=HeteroAgentStationaryEqm_InfHorz(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
peq=[Params.r;Params.Tr;Params.tau_c];

%% Targets (their values do not matter here, only which moments they name)
TM=struct();
TM.AllStats.K.Mean=1;
TM.AllStats.H.Mean=1;
TM.AllStats.H.StdDeviation=1;
TM.CrossSectionCovarCorr.H.Mean=1;
TM.CrossSectionCovarCorr.K.CorrelationWith.earnings=0.5;

%% Both objectives, set up as CalibrateBIHAModel sets them up
[z_gridvals,pi_zJ,vfoptionsJ]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals=vfoptionsJ.e_gridvals;
simoptionsJ.pi_e=vfoptionsJ.pi_e;
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
[targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorrStats,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSecStats,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,0,simoptions);
caliboptions.whichcombos=1; % (CalibrateBIHAModel's default; it hands the objective the selectors and horizons from SetupTargetMoments_InfHorz)
caliboptions.selectors=selectors;
caliboptions.autocorrtimehorizons=autocorrtimehorizons;
caliboptionsN.whichcombos=1;
caliboptionsN.selectors=selectors;
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
statargs={usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats};
jointargs=[{[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames},statargs,{5,[0;1;2;3;4;5],zeros(5,1),zeros(1,1)}];
nestedargs=[{CalibParamNames,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,ParamsOff,DiscountFactorParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqns},statargs,{2,[0;1;2],zeros(2,1),zeros(1,1)}];
GEinner.heteroagentoptions=heteroagentoptions;
GEinner.heteroagentoptions.outputGEstruct=2;
GEinner.z_grid=z_grid;
GEinner.pi_z=pi_z;
GEinner.vfoptions=vfoptions;
GEinner.simoptions=simoptions;

%% (a) the nested objective's inner GE solve
try
    tt=tic;
    [mn,GEn]=CalibrateBIHAModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinner,vfoptionsJ,simoptionsJ);
    timeN=toc(tt);
    mn=gather(mn);
    failed=0;
catch err
    fprintf('JointNested InfHorz (a) the nested objective errored: %s \n',err.message)
    mn=nan(nmom,1); GEn.prices=nan(3,1); GEn.conditions=nan(3,1); failed=1; timeN=NaN;
end
fprintf('JointNested InfHorz (a): nested objective did not run, this should be zero: %.3e; inner GE prices vs the equilibrium, should be near zero (loose: InfHorz fminsearch over three GE prices, started off the equilibrium, stalls on the discrete hours grid): %.3e; its GE conditions, should be near zero (loose: same reason): %.3e \n',failed,max(abs(GEn.prices-peq)),max(abs(GEn.conditions)))

%% (b) its moments against hand-built ones at its own GE prices
ParamsN=Params;
ParamsN.r=GEn.prices(1); ParamsN.Tr=GEn.prices(2); ParamsN.tau_c=GEn.prices(3);
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,ParamsN,DiscountFactorParamNames,[],vfoptionsJ);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_zJ,simoptionsJ,ParamsN);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,ParamsN,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_full);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,ParamsN,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptionsJ);
hand=gather([AllFull.K.Mean; AllFull.H.Mean; AllFull.H.StdDeviation; CSFull.H.Mean; CSFull.K.CorrelationWith.earnings]);
if length(mn)==length(hand)
    d=max(abs(mn-hand));
else
    d=1;
end
fprintf('JointNested InfHorz (b): nested moments vs hand-built at its own GE prices (incl. the two-name CrossSectionCovarCorr target), this should be zero: %.3e \n',d)

%% (c) joint at the nested objective's prices
calibparamsvecJ=[calibparamsvec0; GEn.prices]; % (no constraints, so the GE prices enter untransformed)
tt=tic;
[mj,GEj]=CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},caliboptions,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);
timeJ=toc(tt);
mj=gather(mj);
if length(mj)==nmom+3 && length(mn)==nmom
    d1=max(abs(mj(1:nmom)-mn)); d2=max(abs(mj(nmom+1:end)-GEn.conditions));
else
    d1=1; d2=1;
end
fprintf('JointNested InfHorz (c): joint vs nested at the same GE prices, moments, this should be zero: %.3e; GE conditions, this should be zero: %.3e; joint GEout.prices vs nested, this should be zero: %.3e \n',d1,d2,max(abs(GEj.prices-GEn.prices)))

%% (d) constraints on the GE prices reach the inner solver
GEinnerL=GEinner;
GEinnerL.heteroagentoptions.toleranceGEprices=1e-4;
GEinnerL.heteroagentoptions.toleranceGEcondns=1e-4;
GEinnerL.heteroagentoptions.constrainpositive={'Tr'};
GEinnerS=GEinnerL;
GEinnerL.heteroagentoptions.constrainpositivemethod='log';
GEinnerS.heteroagentoptions.constrainpositivemethod='softplus';
try
    [~,GEl]=CalibrateBIHAModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerL,vfoptionsJ,simoptionsJ);
    [~,GEs]=CalibrateBIHAModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},caliboptionsN,GEinnerS,vfoptionsJ,simoptionsJ);
    dl=max(abs(GEl.prices-peq)); ds=max(abs(GEs.prices-peq)); dls=max(abs(GEl.prices-GEs.prices));
    failed=0;
catch err
    fprintf('JointNested InfHorz (d) errored: %s \n',err.message)
    dl=1; ds=1; dls=0; failed=1;
end
fprintf('JointNested InfHorz (d): inner GE with Tr constrainpositive log vs equilibrium, should be near zero (loose: a coarse inner tolerance, and the InfHorz three-price fminsearch stall): %.3e; softplus, should be near zero (loose: same reason): %.3e; did not run, this should be zero: %.3e \n',dl,ds,failed)
fprintf('JointNested InfHorz (d): log vs softplus inner GE prices (the constraint reached the inner solver), this should NOT be zero: %.3e \n',dls)

%% (e) CalibrateBIHAModel, nested: the GE prices in CalibParams
TMt=struct();
TMt.AllStats.K.Mean=gather(AllFull.K.Mean); % (at the nested prices, which are the equilibrium to within its tolerance)
TMt.AllStats.H.Mean=gather(AllFull.H.Mean);
caliboptionsC=struct();
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
caliboptionsC.jointoptimization=0;
try
    [CalibParamsN,calibsummaryN]=CalibrateBIHAModel(CalibParamNames,TMt,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsC,vfoptions,simoptions);
    missing=sum(~isfield(CalibParamsN,GEPriceParamNames));
    if missing==0
        dp=max(abs([CalibParamsN.r;CalibParamsN.Tr;CalibParamsN.tau_c]-peq));
    else
        dp=1;
    end
    dg=max(abs(calibsummaryN.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('JointNested InfHorz (e) errored: %s \n',err.message)
    missing=3; dp=1; dg=1; failed=1;
end
fprintf('JointNested InfHorz (e): nested CalibrateBIHAModel, GE prices missing from CalibParams, this should be zero: %.3e; did not run, this should be zero: %.3e \n',missing,failed)
fprintf('JointNested InfHorz (e): nested CalibParams GE prices vs the equilibrium, should be near zero: %.3e; calibsummary.GeneralEqmConditions, should be near zero: %.3e \n',dp,dg)

%% (f) joint with a 'log' constrainpositive GE price while the calibrated parameters use 'softplus'
caliboptionsC.jointoptimization=1;
heteroagentoptionsF=struct();
heteroagentoptionsF.constrainpositive={'Tr'};
heteroagentoptionsF.constrainpositivemethod='log';
try
    [CalibParamsF,~]=CalibrateBIHAModel(CalibParamNames,TMt,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsF,caliboptionsC,vfoptions,simoptions);
    dF=max(abs([CalibParamsF.beta-Params.beta, CalibParamsF.varphi-Params.varphi, CalibParamsF.r-Params.r, CalibParamsF.Tr-Params.Tr, CalibParamsF.tau_c-Params.tau_c]));
    failed=0;
catch err
    fprintf('JointNested InfHorz (f) errored: %s \n',err.message)
    dF=1; failed=1;
end
fprintf('JointNested InfHorz (f): joint, Tr constrainpositive by log, calibrated params by softplus: did not run, this should be zero: %.3e; returned its start (the log round trip of Tr), this should be zero: %.3e \n',failed,dF)

%% (g) runtime
fprintf('JointNested InfHorz (g) runtime (informational): one nested objective evaluation (a GE solve from 10%% off) %.2f s, one joint evaluation %.2f s \n',timeN,timeJ)

output=struct();
output.GEn=GEn;

end
