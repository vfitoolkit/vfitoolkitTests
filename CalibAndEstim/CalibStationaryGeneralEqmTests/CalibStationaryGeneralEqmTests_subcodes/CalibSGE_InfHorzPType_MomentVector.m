function output=CalibSGE_InfHorzPType_MomentVector(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 17: the vector that CalibrateBIHAModel_PType_Joint_objectivefn builds (vectoroutput=1: the model
% moments at the targeted entries, then the GE conditions), against one built by hand. Part 7 with two permanent types.
%
% The general eqm is solved at the truth (HeteroAgentStationaryEqm_InfHorz_PType). Full PType AllStats, AutoCorrTransProbs and
% CrossSectionCovarCorr runs there give the targets, grouped, per type, restricted and restricted per type. The joint objective is then
% called, set up exactly as CalibrateBIHAModel_PType sets it up, at (truth, equilibrium prices). Checks (all exact):
%   (a) the moment block equals the hand-built vector (walk TargetMoments two to four names deep, read the same field of the full run),
%       under caliboptions.whichcombos=1 and =0
%   (b) the moment block equals the targets (the targets are the model)
%   (c) the GE block equals HeteroAgentStationaryEqm_InfHorz_PType with maxiter=0 (the GE conditions at the prices in Params)
%   (d) the second output: the GE prices are the equilibrium prices it was handed, the GE conditions the GE block
% Before 2026-10-29 the BIHA PType joint objective had no whichcombos and returned no second output.

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
p_eqm=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
fprintf('MomentVector InfHorz PType: general eqm at the truth r=%.6f, Tr=%.6f, tau_c=%.6f \n',Params.r,Params.Tr,Params.tau_c)

%% Grids as CalibrateBIHAModel_PType sets them up, used by the full runs too
[z_gridvals,pi_zJ,vfoptionsJ]=ExogShockSetup_InfHorz_PType(n_z,z_grid,pi_z,Names_i,Params,vfoptions,3,0);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals=vfoptionsJ.e_gridvals;
simoptionsJ.pi_e=vfoptionsJ.pi_e;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;
simoptionsJ.warnjequaloneptypeasdim=0;

%% Full stats at the equilibrium
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,vfoptionsJ);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_zJ,Params,simoptionsJ);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simoptions_full);
ACFull=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zJ,simoptionsJ);
CSFull=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,simoptionsJ);

%% Targets from the full runs: grouped, per type, restricted, restricted per type
TM=struct();
TM.AllStats.K.Mean=gather(AllFull.K.Mean);
TM.AllStats.K.typeA.Mean=gather(AllFull.K.typeA.Mean);
TM.AllStats.H.typeB.Mean=gather(AllFull.H.typeB.Mean);
TM.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
TM.AllStats.earnings.typeB.QuantileCutoffs=gather(AllFull.earnings.typeB.QuantileCutoffs);
TM.AllStats.rich.H.Mean=gather(AllFull.rich.H.Mean);
TM.AllStats.rich.earnings.typeA.Median=gather(AllFull.rich.earnings.typeA.Median);
TM.AutoCorrTransProbs.earnings.AutoCorrelation=gather(ACFull.earnings.AutoCorrelation);
TM.AutoCorrTransProbs.earnings.typeA.AutoCorrelation=gather(ACFull.earnings.typeA.AutoCorrelation);
TM.CrossSectionCovarCorr.H.typeB.Mean=gather(CSFull.H.typeB.Mean);
TM.CrossSectionCovarCorr.K.CorrelationWith.earnings=gather(CSFull.K.CorrelationWith.earnings);
TM.CrossSectionCovarCorr.K.typeA.CorrelationWith.earnings=gather(CSFull.K.typeA.CorrelationWith.earnings);

% The hand-built vector, walking TargetMoments in field order, two to four names deep (the full-run structs have the same nesting)
hand=[]; tvec=[];
kinds={'AllStats','AutoCorrTransProbs','CrossSectionCovarCorr'};
for kk=1:length(kinds)
    if kk==1
        Full=AllFull;
    elseif kk==2
        Full=ACFull;
    else
        Full=CSFull;
    end
    a1vec=fieldnames(TM.(kinds{kk}));
    for a1=1:length(a1vec)
        a2vec=fieldnames(TM.(kinds{kk}).(a1vec{a1}));
        for a2=1:length(a2vec)
            t2=TM.(kinds{kk}).(a1vec{a1}).(a2vec{a2});
            if isstruct(t2)
                a3vec=fieldnames(t2);
                for a3=1:length(a3vec)
                    t3=t2.(a3vec{a3});
                    if isstruct(t3)
                        a4vec=fieldnames(t3);
                        for a4=1:length(a4vec)
                            v=gather(Full.(a1vec{a1}).(a2vec{a2}).(a3vec{a3}).(a4vec{a4})); hand=[hand; v(:)];
                            v=t3.(a4vec{a4}); tvec=[tvec; v(:)];
                        end
                    else
                        v=gather(Full.(a1vec{a1}).(a2vec{a2}).(a3vec{a3})); hand=[hand; v(:)];
                        tvec=[tvec; t3(:)];
                    end
                end
            else
                v=gather(Full.(a1vec{a1}).(a2vec{a2})); hand=[hand; v(:)];
                tvec=[tvec; t2(:)];
            end
        end
    end
end
actual=~isnan(tvec);
fprintf('MomentVector InfHorz PType: %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %.3e \n',length(tvec),sum(actual),sum(isnan(hand(actual))))

%% The joint objective, set up as CalibrateBIHAModel_PType sets it up (two common calibrated parameters, then the three GE prices)
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
[calibparamsvec0,caliboptions]=ParameterConstraints_PType_TransformParamsToUnconstrained(calibparamsvec0,[0;1;2],CalibParamNames,[1,0;2,0],caliboptions,1);
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
[targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorrStats,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSecStats,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TM,FnsToEvaluate,1,simoptions,Names_i);
caliboptions.selectors=selectors;
caliboptions.autocorrtimehorizons=autocorrtimehorizons;
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptions.weights=ones(sum(~isnan(targetmomentvec)),1);
if length(targetmomentvec)==length(tvec)
    d=abs(targetmomentvec-tvec); d=d(~isnan(d));
    patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
else
    d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
end
fprintf('MomentVector InfHorz PType: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %.3e \n',max([d;0]),patternbad)

nmom=sum(~isnan(targetmomentvec));
handA=hand(actual);
tA=targetmomentvec(~isnan(targetmomentvec));
jointargs={[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,Names_i,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats,5,[1,0;2,0;3,0;4,0;5,0],[0;1;2;3;4;5],ones(5,2),zeros(5,1),{}};
for wc=1:-1:0
    caliboptions_wc=caliboptions;
    caliboptions_wc.whichcombos=wc;
    [m,GEout]=CalibrateBIHAModel_PType_Joint_objectivefn(calibparamsvecJ,jointargs{:},caliboptions_wc,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);
    m=gather(m);
    if length(m)==nmom+3
        mm=m(1:nmom); mg=m(nmom+1:end);
        d=max(abs(mm-handA)); dt=max(abs(mm-tA));
    else
        mm=m; mg=nan(3,1); d=1; dt=1;
    end
    fprintf('MomentVector InfHorz PType (a) whichcombos=%i: moment block vs hand-built, this should be zero: %.3e; NaN entries + length mismatch, this should be zero: %.3e \n',wc,d,sum(isnan(mm))+abs(length(m)-(nmom+3)))
    fprintf('MomentVector InfHorz PType (b) whichcombos=%i: moment block vs the targets (targets are the model), this should be zero: %.3e \n',wc,dt)
    if wc==1
        m1=m; GEout1=GEout; mg1=mg;
    end
end

%% (c) the GE block against HeteroAgentStationaryEqm_InfHorz_PType with maxiter=0 at the same prices
heteroagentoptions0=heteroagentoptions;
heteroagentoptions0.maxiter=0;
[~,GEcondns0]=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions0,simoptions,vfoptions);
g0=gather([GEcondns0.CapitalMarket; GEcondns0.GovBudget; GEcondns0.ConsTax]);
fprintf('MomentVector InfHorz PType (c): GE block vs HeteroAgentStationaryEqm_InfHorz_PType maxiter=0 at the same prices, this should be zero: %.3e (the GE conditions themselves, max abs %.3e) \n',max(abs(mg1-g0)),max(abs(g0)))

%% (d) the second output
fprintf('MomentVector InfHorz PType (d): GEout.prices vs the prices handed in, this should be zero: %.3e; GEout.conditions vs the GE block, this should be zero: %.3e \n',max(abs(GEout1.prices-[Params.r;Params.Tr;Params.tau_c])),max(abs(GEout1.conditions-mg1)))

output=struct();
output.TargetMoments=TM;
output.m1=m1;

end
