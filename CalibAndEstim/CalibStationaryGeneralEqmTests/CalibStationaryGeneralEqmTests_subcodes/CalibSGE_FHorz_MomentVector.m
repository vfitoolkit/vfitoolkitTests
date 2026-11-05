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
%   (e) window targets (simoptions.ntimeperiods) of AgeConditionalStats, AutoCorrTransProbs, CrossSectionCovarCorr and
%       AgeConditionalCrossSectionCovarCorr, with a restriction: the joint moment block equals the hand walk through full runs with
%       the windows and the targets, whichcombos=1 and =0, with the default agegroupings (one call per command) and with
%       agegroupings=1:5:N_j (the windows in second LifeCycleProfiles and age-conditional cross-section calls); the nested moment
%       vector equals the joint one at the nested GE prices (added 2026-11-05: the OLG objectives passed the windows only to
%       AutoCorrTransProbs)

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

%% (e) window targets (simoptions.ntimeperiods) of every kind, joint and nested
% Full runs with the windows at the equilibrium; window targets of AgeConditionalStats, AutoCorrTransProbs, CrossSectionCovarCorr and
% AgeConditionalCrossSectionCovarCorr (one under the restriction), plus a single-period target of each age-conditional kind. Scenario 1
% has the default agegroupings (one call per command); scenario 2 agegroupings=1:5:N_j, the single-period age-conditional targets by
% age group, so the windows (by start age) are a second LifeCycleProfiles and a second age-conditional cross-section call.
simoptions_W=simoptionsJ; simoptions_W.whichstats=ones(1,7); simoptions_W.ntimeperiods=[1,2,3]; simoptions_W.pi_z=pi_z_J;
LCPW=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_W);
CSW=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_W);
simoptions_W.ntimeperiods=[1,2];
ACXW=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_W);
simoptions_W.timehorizons=[1,2];
ACW=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,simoptions_W);
simoptions_G=simoptionsJ; simoptions_G.whichstats=ones(1,7); simoptions_G.agegroupings=1:5:N_j;
LCPG=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_G); % the single period by age groups
ACXG=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_G);
kindsW={'AllStats','AgeConditionalStats','AutoCorrTransProbs','CrossSectionCovarCorr','AgeConditionalCrossSectionCovarCorr'};
% the nested objective: its inner GE solve (from the equilibrium) lands within tolerance of the equilibrium, not on it, so its moments
% are checked against the joint objective at the nested GE prices (exact), the joint objective against the full runs as above
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
caliboptionsN.vectoroutput=1;
[calibparamsvecN,caliboptionsN]=ParameterConstraints_TransformParamsToUnconstrained([Params.beta; Params.varphi],[0;1;2],CalibParamNames,caliboptionsN,1);
GEinner=struct();
GEinner.heteroagentoptions=heteroagentoptions;
GEinner.heteroagentoptions.outputGEstruct=2;
GEinner.z_grid=z_grid;
GEinner.pi_z=pi_z;
GEinner.vfoptions=vfoptions;
GEinner.simoptions=simoptions;
for scen=1:2
    LCPFullW=LCPW; ACXFullW=ACXW;
    simoptions_s=simoptions; simoptionsJ_s=simoptionsJ;
    if scen==1
        scenstr='one call each (default agegroupings)';
    else
        scenstr='second calls for the windows (agegroupings 1:5:N_j)';
        LCPFullW.H=LCPG.H; % the single-period age-conditional fields by age groups
        ACXFullW.earnings=ACXG.earnings;
        simoptions_s.agegroupings=1:5:N_j;
        simoptionsJ_s.agegroupings=1:5:N_j;
    end
    % the targets, the values from the full runs, NaN where omitted
    TMW=struct();
    TMW.AllStats.K.Mean=gather(AllFull.K.Mean);
    TMW.AgeConditionalStats.H.Mean=gather(LCPFullW.H.Mean);
    t=gather(LCPFullW.earnings_2periodAverage.StdDeviation); t([1,N_j])=NaN; TMW.AgeConditionalStats.earnings_2periodAverage.StdDeviation=t;
    t=gather(LCPFullW.rich.C_3periodAverage.Mean); t([1:5,N_j-1:N_j])=NaN; TMW.AgeConditionalStats.rich.C_3periodAverage.Mean=t;
    t=gather(ACW.earnings_2periodAverage.AutoCorrelation_horizon2); t(1)=NaN; TMW.AutoCorrTransProbs.earnings_2periodAverage.AutoCorrelation_horizon2=t;
    TMW.CrossSectionCovarCorr.H.Mean=gather(CSW.H.Mean);
    TMW.CrossSectionCovarCorr.earnings_2periodAverage.CorrelationWith.K_2periodAverage=gather(CSW.earnings_2periodAverage.CorrelationWith.K_2periodAverage);
    TMW.CrossSectionCovarCorr.rich.H_3periodAverage.Mean=gather(CSW.rich.H_3periodAverage.Mean);
    t=gather(ACXFullW.earnings.Mean); t(2)=NaN; TMW.AgeConditionalCrossSectionCovarCorr.earnings.Mean=t;
    t=gather(ACXFullW.K_2periodAverage.CovarianceWith.earnings_2periodAverage); t([1,11:end])=NaN; TMW.AgeConditionalCrossSectionCovarCorr.K_2periodAverage.CovarianceWith.earnings_2periodAverage=t;
    % the hand walk: every leaf of the targets in depth-first fieldnames order (as SetupTargetMoments reads them), the same path in the full runs
    hand=[]; tvec=[];
    for kk=1:length(kindsW)
        if kk==1
            Full=AllFull;
        elseif kk==2
            Full=LCPFullW;
        elseif kk==3
            Full=ACW;
        elseif kk==4
            Full=CSW;
        else
            Full=ACXFullW;
        end
        paths={{}};
        while ~isempty(paths)
            p=paths{1}; paths(1)=[];
            node=TMW.(kindsW{kk}); v=Full;
            for ll=1:length(p)
                node=node.(p{ll}); v=v.(p{ll});
            end
            if isstruct(node)
                fn=fieldnames(node);
                children=cell(1,length(fn));
                for ll=1:length(fn)
                    children{ll}=[p,fn(ll)];
                end
                paths=[children,paths];
            else
                v=gather(v); hand=[hand; v(:)]; tvec=[tvec; node(:)];
            end
        end
    end
    actual=~isnan(tvec);
    fprintf('MomentVector (e) %s: %i target entries, %i of them targeted; the model is NaN at a targeted entry, this should be zero: %i \n',scenstr,length(tvec),sum(actual),sum(isnan(hand(actual))))

    [targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMW,FnsToEvaluate,0,N_j,simoptions_s);
    if length(targetmomentvec)==length(tvec)
        d=abs(targetmomentvec-tvec); d=d(~isnan(d));
        patternbad=sum(isnan(targetmomentvec)~=isnan(tvec));
    else
        d=1; patternbad=abs(length(targetmomentvec)-length(tvec));
    end
    if scen==1
        splitbad=~(selectors.ACStats.split==0 && selectors.AgeCrossSec.split==0);
    else
        splitbad=~(selectors.ACStats.split==1 && isfield(selectors,'ACStatsW') && selectors.AgeCrossSec.split==1 && isfield(selectors,'AgeCrossSecW'));
    end
    fprintf('MomentVector (e) %s: SetupTargetMoments targetmomentvec vs the hand walk, this should be zero: %.3e; length + NaN-pattern mismatches, this should be zero: %i; split as expected, this should be zero: %i \n',scenstr,max([d;0]),patternbad,splitbad)
    caliboptionsW=caliboptions;
    caliboptionsW.selectors=selectors;
    caliboptionsW.logmoments=zeros(length(targetmomentvec),1);
    caliboptionsW.weights=ones(sum(~isnan(targetmomentvec)),1);
    nmomW=sum(~isnan(targetmomentvec));
    handA=hand(actual);
    tA=targetmomentvec(~isnan(targetmomentvec));
    for wc=1:-1:0
        caliboptions_wc=caliboptionsW;
        caliboptions_wc.whichcombos=wc;
        m=CalibrateOLGModel_Joint_objectivefn(calibparamsvec0,CalibParamNamesJ,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,5,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix,caliboptions_wc,heteroagentoptionsJ,vfoptionsJ,simoptionsJ_s);
        m=gather(m);
        if length(m)==nmomW+3
            mm=m(1:nmomW);
            d=max(abs(mm-handA)); dt=max(abs(mm-tA));
        else
            mm=m; d=1; dt=1;
        end
        fprintf('MomentVector (e) %s, joint, whichcombos=%i: moment block vs hand-built, this should be zero: %.3e; vs the targets, this should be zero: %.3e; NaN entries + length mismatch, this should be zero: %i \n',scenstr,wc,d,dt,sum(isnan(mm))+abs(length(m)-(nmomW+3)))
        caliboptionsN_wc=caliboptionsN;
        caliboptionsN_wc.selectors=selectors;
        caliboptionsN_wc.logmoments=caliboptionsW.logmoments;
        caliboptionsN_wc.weights=caliboptionsW.weights;
        caliboptionsN_wc.whichcombos=wc;
        [mn,GEn]=CalibrateOLGModel_Nested_objectivefn(calibparamsvecN,CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqns,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,2,[0;1;2],zeros(2,1),zeros(N_j,1),caliboptionsN_wc,GEinner,vfoptionsJ,simoptionsJ_s);
        mn=gather(mn);
        pN=GEn.prices; % (no GE price is constrained here, so the transformed prices are the prices)
        mj=CalibrateOLGModel_Joint_objectivefn([calibparamsvec0(1:2); pN],CalibParamNamesJ,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,5,calibparamsvecindex,calibomitparams_counter,calibomitparamsmatrix,caliboptions_wc,heteroagentoptionsJ,vfoptionsJ,simoptionsJ_s);
        mj=gather(mj);
        if length(mn)==nmomW && length(mj)==nmomW+3
            dn=max(abs(mn-mj(1:nmomW))); dnt=max(abs(mn-tA));
        else
            dn=1; dnt=1;
        end
        fprintf('MomentVector (e) %s, nested, whichcombos=%i: moment vector vs the joint one at the nested GE prices, this should be zero: %.3e; NaN entries + length mismatch, this should be zero: %i; (vs the targets, should be near zero: %.3e; inner GE prices vs the equilibrium, should be near zero: %.3e) \n',scenstr,wc,dn,sum(isnan(mn))+abs(length(mn)-nmomW),dnt,max(abs(GEn.prices-[Params.r;Params.Tr;Params.tau_c])))
    end
end

output=struct();
output.TargetMoments=TM;
output.m1=m1;
output.hand=handA;

end
