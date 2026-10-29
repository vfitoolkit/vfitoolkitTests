function output=CalibSGE_InfHorz_MomentVector(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 7: the vector that CalibrateBIHAModel_Joint_objectivefn builds (vectoroutput=1: the model moments
% at the targeted entries, then the GE conditions), against one built by hand.
%
% The general eqm is solved at the truth. Full AllStats, AutoCorrTransProbs and CrossSectionCovarCorr runs there give the targets:
% AllStats unrestricted and under a conditional restriction, an autocorrelation, and cross-sectional targets named by two names
% (H.Mean, which the objective used to read with three names) and by three (K.CorrelationWith.earnings). The joint objective is then
% called, set up exactly as CalibrateBIHAModel sets it up, at (truth, equilibrium prices).
% Checks (all exact):
%   (a) the moment block equals the hand-built vector
%   (b) the moment block equals the targets (the targets are the model)
%   (c) the GE block equals HeteroAgentStationaryEqm_InfHorz with maxiter=0 (the GE conditions evaluated at the prices in Params)
%   (d) the second output of the objective: the GE prices are the equilibrium prices it was handed, the GE conditions the GE block

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
FnsToEvaluate.earnings=@(d,aprime,a,z,r,alpha,delta,A) ((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d;
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

simoptions.conditionalrestrictions.rich=@(d,aprime,a,z) (a>2);

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
fprintf('MomentVector InfHorz: general eqm at the truth r=%.6f, Tr=%.6f, tau_c=%.6f \n',Params.r,Params.Tr,Params.tau_c)

%% Grids as CalibrateBIHAModel sets them up, used by the full runs too
[z_gridvals,pi_zJ,vfoptionsJ]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals=vfoptionsJ.e_gridvals;
simoptionsJ.pi_e=vfoptionsJ.pi_e;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;

%% Full stats at the equilibrium
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsJ);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_zJ,simoptionsJ,Params);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_full);
simoptions_unr=simoptionsJ;
simoptions_unr=rmfield(simoptions_unr,'conditionalrestrictions'); % (the AutoCorr and CrossSec targets are unrestricted; the objective drops the restrictions for them)
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,simoptions_unr);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_unr);

%% Targets from the full runs
TM=struct();
TM.AllStats.K.Mean=gather(AllFull.K.Mean);
TM.AllStats.H.Mean=gather(AllFull.H.Mean);
TM.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
TM.AllStats.K.Variance=gather(AllFull.K.Variance);
TM.AllStats.rich.H.Mean=gather(AllFull.rich.H.Mean);
TM.AutoCorrTransProbs.earnings.AutoCorrelation=gather(ACFull.earnings.AutoCorrelation);
TM.CrossSectionCovarCorr.H.Mean=gather(CSFull.H.Mean);
TM.CrossSectionCovarCorr.K.CorrelationWith.earnings=gather(CSFull.K.CorrelationWith.earnings);
hand=[TM.AllStats.K.Mean; TM.AllStats.K.Variance; TM.AllStats.H.Mean; TM.AllStats.earnings.Gini; TM.AllStats.rich.H.Mean; TM.AutoCorrTransProbs.earnings.AutoCorrelation; TM.CrossSectionCovarCorr.H.Mean; TM.CrossSectionCovarCorr.K.CorrelationWith.earnings];
% (the order of TargetMoments: AllStats fields K (Mean, Variance), H, earnings, rich; then AutoCorrTransProbs; then CrossSectionCovarCorr)

%% The joint objective, set up as CalibrateBIHAModel sets it up
CalibParamNames={'beta','varphi'};
calibparamsvec0=[Params.beta; Params.varphi];
caliboptions=struct();
caliboptions.verbose=0;
caliboptions.constrainpositive={};
caliboptions.constrainpositivemethod='softplus';
caliboptions.constrain0to1={};
caliboptions.constrainAtoB={};
caliboptions.metric='sum_squared';
caliboptions.relativeGEweight=10;
caliboptions.vectoroutput=1;
caliboptions.useCustomModelStats=0;
caliboptions.calibrateshocks=0;
caliboptions.nCalibParamsNonGE=2;
[calibparamsvec0,caliboptions]=ParameterConstraints_TransformParamsToUnconstrained(calibparamsvec0,[0;1;2],CalibParamNames,caliboptions,1);
heteroagentoptionsJ=struct();
heteroagentoptionsJ.outputGEstruct=2;
heteroagentoptionsJ.multiGEweights=ones(1,3);
heteroagentoptionsJ.constrainpositive={};
heteroagentoptionsJ.constrainpositivemethod='softplus';
heteroagentoptionsJ.constrain0to1={};
heteroagentoptionsJ.constrainAtoB={};
heteroagentoptionsJ.useCustomModelStats=0;
heteroagentoptionsJ.useintermediateEqns=0;
GEparamsvec0=[Params.r; Params.Tr; Params.tau_c];
[GEparamsvec0,heteroagentoptionsJ]=ParameterConstraints_TransformParamsToUnconstrained(GEparamsvec0,0:1:3,GEPriceParamNames,heteroagentoptionsJ,1);
calibparamsvecJ=[calibparamsvec0; GEparamsvec0];
caliboptions.constrainpositive=[caliboptions.constrainpositive; heteroagentoptionsJ.constrainpositive];
caliboptions.constrainAtoB=[caliboptions.constrainAtoB; heteroagentoptionsJ.constrainAtoB];
caliboptions.constrain0to1=[caliboptions.constrain0to1; heteroagentoptionsJ.constrain0to1];
caliboptions.constrainpositivesoftplus=[caliboptions.constrainpositivesoftplus; heteroagentoptionsJ.constrainpositivesoftplus];
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
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptions.weights=ones(length(targetmomentvec),1);
if length(targetmomentvec)==length(hand)
    d=max(abs(targetmomentvec-hand));
else
    d=1;
end
fprintf('MomentVector InfHorz: SetupTargetMoments targetmomentvec vs the targets in TargetMoments order, this should be zero: %.3e \n',d)

[m,GEout]=CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats,5,[0;1;2;3;4;5],zeros(5,1),zeros(1,1),caliboptions,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);
m=gather(m);
nmom=length(hand);
if length(m)==nmom+3
    mm=m(1:nmom); mg=m(nmom+1:end);
    d=max(abs(mm-hand)); dt=max(abs(mm-targetmomentvec));
else
    mg=nan(3,1); d=1; dt=1;
end
fprintf('MomentVector InfHorz (a): moment block vs hand-built, this should be zero: %.3e; length mismatch, this should be zero: %.3e \n',d,abs(length(m)-(nmom+3)))
fprintf('MomentVector InfHorz (b): moment block vs the targets (targets are the model), this should be zero: %.3e \n',dt)
caliboptions0=caliboptions;
caliboptions0.whichcombos=0; % every statistic of every targeted function: the same moments
m0=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats,5,[0;1;2;3;4;5],zeros(5,1),{},caliboptions0,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
if length(m0)==length(m)
    d=max(abs(m0-m));
else
    d=1;
end
fprintf('MomentVector InfHorz (b): whichcombos=0 vs =1, moments and GE conditions, this should be zero: %.3e \n',d)

%% (c) the GE block against HeteroAgentStationaryEqm with maxiter=0 at the same prices
heteroagentoptions0=heteroagentoptions;
heteroagentoptions0.maxiter=0;
[~,GEcondns0]=HeteroAgentStationaryEqm_InfHorz(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions0,simoptions,vfoptions);
g0=gather([GEcondns0.CapitalMarket; GEcondns0.GovBudget; GEcondns0.ConsTax]);
fprintf('MomentVector InfHorz (c): GE block vs HeteroAgentStationaryEqm maxiter=0 at the same prices, this should be zero: %.3e (the GE conditions themselves, max abs %.3e) \n',max(abs(mg-g0)),max(abs(g0)))

%% (d) the second output
fprintf('MomentVector InfHorz (d): GEout.prices vs the prices handed in, this should be zero: %.3e; GEout.conditions vs the GE block, this should be zero: %.3e \n',max(abs(GEout.prices-[Params.r;Params.Tr;Params.tau_c])),max(abs(GEout.conditions-mg)))

output=struct();
output.TargetMoments=TM;
output.m=m;

end
