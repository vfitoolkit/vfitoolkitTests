function output=CalibSGE_FHorz_MomentVector(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 1: the vector that CalibrateOLGModel_Joint_objectivefn builds (vectoroutput=1: the model moments
% at the targeted entries, then the GE conditions), against one built by hand.
%
% The general eqm is solved at the truth. Full AllStats, LifeCycleProfiles, AutoCorrTransProbs and CrossSectionCovarCorr runs there
% give the targets (NaN where an entry is omitted): every kind of target CalibrateOLGModel now takes, with a conditional restriction.
% The joint objective is then called, set up exactly as CalibrateOLGModel sets it up, at (truth, equilibrium prices).
% Checks (all exact):
%   (a) the moment block equals the hand-built vector (walk TargetMoments, read the same field of the full run, vectorise
%       column-major, drop the NaN targets), under caliboptions.whichcombos=1 and =0
%   (b) the moment block equals the targets (the targets are the model)
%   (c) the GE block equals HeteroAgentStationaryEqm_Case1_FHorz with maxiter=0 (the GE conditions evaluated at the prices in Params)
%   (d) the second output of the objective: the GE prices are the equilibrium prices it was handed, the GE conditions the GE block

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
FnsToEvaluate.earnings=@(d,aprime,a,z,r,kappa_j,alpha,delta,A) ((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d;
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
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
fprintf('MomentVector: general eqm at the truth r=%.6f, Tr=%.6f, tau_c=%.6f \n',Params.r,Params.Tr,Params.tau_c)

%% Grids as CalibrateOLGModel sets them up (jequaloneDist_usergrids defaults to 1, so the user's grids are kept), used by the full runs too
[z_gridvals_J,pi_z_J,vfoptionsJ]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,1);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals_J=vfoptionsJ.e_gridvals_J;
simoptionsJ.pi_e_J=vfoptionsJ.pi_e_J;
simoptionsJ.jequaloneDist_usergrids=1;
simoptionsJ.user_z_grid=vfoptionsJ.user_z_grid;
simoptionsJ.user_pi_z=vfoptionsJ.user_pi_z;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;

%% Full stats at the equilibrium
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsJ);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptionsJ);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
simoptions_ac=simoptionsJ;
simoptions_ac.timehorizons=[1,2]; % horizons 1 and 2 (exactly the horizons listed are computed; horizon 2 comes out with the suffix _horizon2)
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_ac);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptionsJ);

%% Targets from the full runs, NaN where omitted
TM=struct();
TM.AllStats.K.Mean=gather(AllFull.K.Mean);
TM.AllStats.H.Mean=gather(AllFull.H.Mean);
TM.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
TM.AllStats.earnings.QuantileCutoffs=gather(AllFull.earnings.QuantileCutoffs);
TM.AllStats.K.Variance=gather(AllFull.K.Variance);
TM.AllStats.rich.H.Mean=gather(AllFull.rich.H.Mean);
TM.AllStats.rich.earnings.Median=gather(AllFull.rich.earnings.Median);
t=gather(LCPFull.K.Mean); t(1)=NaN; TM.AgeConditionalStats.K.Mean=t; % (everyone has zero assets at age 1)
TM.AgeConditionalStats.H.Mean=gather(LCPFull.H.Mean);
t=gather(LCPFull.earnings.StdDeviation); t(1:2)=NaN; TM.AgeConditionalStats.earnings.StdDeviation=t;
t=gather(LCPFull.rich.C.Mean); t(1:5)=NaN; TM.AgeConditionalStats.rich.C.Mean=t;
t=gather(ACFull.earnings.AutoCorrelation); t(1)=NaN; TM.AutoCorrTransProbs.earnings.AutoCorrelation=t; % horizon 1, by start age
TM.AutoCorrTransProbs.H.AutoCovariance_horizon2=gather(ACFull.H.AutoCovariance_horizon2); % horizon 2
TM.CrossSectionCovarCorr.K.CorrelationWith.earnings=gather(CSFull.K.CorrelationWith.earnings);
TM.CrossSectionCovarCorr.H.Mean=gather(CSFull.H.Mean);

% The hand-built vector, walking TargetMoments in field order (the full-run structs have the same nesting as the targets)
hand=[]; tvec=[];
kinds={'AllStats','AgeConditionalStats','AutoCorrTransProbs','CrossSectionCovarCorr'};
for kk=1:length(kinds)
    if kk==1
        Full=AllFull;
    elseif kk==2
        Full=LCPFull;
    elseif kk==3
        Full=ACFull;
    else
        Full=CSFull;
    end
    a1vec=fieldnames(TM.(kinds{kk}));
    for a1=1:length(a1vec)
        a2vec=fieldnames(TM.(kinds{kk}).(a1vec{a1}));
        for a2=1:length(a2vec)
            tval=TM.(kinds{kk}).(a1vec{a1}).(a2vec{a2});
            if isstruct(tval)
                a3vec=fieldnames(tval);
                for a3=1:length(a3vec)
                    v=gather(Full.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)];
                    v=tval.(a3vec{a3}); tvec=[tvec; v(:)];
                end
            else
                v=gather(Full.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)];
                tvec=[tvec; tval(:)];
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('MomentVector: %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

%% The joint objective, set up as CalibrateOLGModel sets it up
CalibParamNames={'beta','varphi'};
nCalibParamsNonGE=2;
calibparamsvec0=[Params.beta; Params.varphi];
calibparamsvecindex=[0;1;2];
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
caliboptions.nCalibParamsNonGE=nCalibParamsNonGE;
[calibparamsvec0,caliboptions]=ParameterConstraints_TransformParamsToUnconstrained(calibparamsvec0,calibparamsvecindex,CalibParamNames,caliboptions,1);
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
calibparamsvec0=[calibparamsvec0; GEparamsvec0];
CalibParamNamesJ=[CalibParamNames,GEPriceParamNames];
caliboptions.constrainpositive=[caliboptions.constrainpositive; heteroagentoptionsJ.constrainpositive];
caliboptions.constrainAtoB=[caliboptions.constrainAtoB; heteroagentoptionsJ.constrainAtoB];
caliboptions.constrain0to1=[caliboptions.constrain0to1; heteroagentoptionsJ.constrain0to1];
caliboptions.constrainpositivesoftplus=[caliboptions.constrainpositivesoftplus; heteroagentoptionsJ.constrainpositivesoftplus];
caliboptions.constrainAtoBlimits=zeros(5,2);
calibparamsvecindex=[0;1;2;3;4;5];
calibomitparams_counter=zeros(5,1);
calibomitparamsmatrix=zeros(N_j,1);
GEeqnNames=fieldnames(GeneralEqmEqns);
GeneralEqmEqnsCell=cell(1,3);
for gg=1:3
    GeneralEqmEqnParamNames(gg).Names=getAnonymousFnInputNames(GeneralEqmEqns.(GEeqnNames{gg}));
    GeneralEqmEqnsCell{gg}=GeneralEqmEqns.(GEeqnNames{gg});
end
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TM,FnsToEvaluate,0,N_j,simoptions);
caliboptions.selectors=selectors;
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);

if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('MomentVector: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)

nmom=sum(~isnan(targetmomentvec));
handA=hand(actual);
tA=targetmomentvec(~isnan(targetmomentvec));
for wc=1:-1:0
    caliboptions_wc=caliboptions;
    caliboptions_wc.whichcombos=wc;
    [m,GEout]=CalibrateOLGModel_Joint_objectivefn(calibparamsvec0,CalibParamNamesJ,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,5,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix,caliboptions_wc,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);
    m=gather(m);
    if length(m)==nmom+3
        mm=m(1:nmom); mg=m(nmom+1:end);
        d=max(abs(mm-handA)); dt=max(abs(mm-tA));
    else
        mm=m; mg=nan(3,1); d=1; dt=1;
    end
    fprintf('MomentVector (a) whichcombos=%i: moment block vs hand-built, this should be zero: %.3e; NaN entries + length mismatch, this should be zero: %.3e \n',wc,d,sum(isnan(mm))+abs(length(m)-(nmom+3)))
    fprintf('MomentVector (b) whichcombos=%i: moment block vs the targets (targets are the model), this should be zero: %.3e \n',wc,dt)
    if wc==1
        m1=m; GEout1=GEout; mg1=mg;
    end
end

%% (c) the GE block against HeteroAgentStationaryEqm with maxiter=0 at the same prices
heteroagentoptions0=heteroagentoptions;
heteroagentoptions0.maxiter=0;
[~,GEcondns0]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions0,simoptions,vfoptions);
g0=[GEcondns0.CapitalMarket; GEcondns0.GovBudget; GEcondns0.ConsTax];
fprintf('MomentVector (c): GE block vs HeteroAgentStationaryEqm maxiter=0 at the same prices, this should be zero: %.3e (the GE conditions themselves, max abs %.3e) \n',max(abs(mg1-gather(g0))),max(abs(gather(g0))))

%% (d) the second output
fprintf('MomentVector (d): GEout.prices vs the prices handed in, this should be zero: %.3e; GEout.conditions vs the GE block, this should be zero: %.3e \n',max(abs(GEout1.prices-[Params.r;Params.Tr;Params.tau_c])),max(abs(GEout1.conditions-mg1)))

output=struct();
output.TargetMoments=TM;
output.m1=m1;
output.hand=handA;

end
