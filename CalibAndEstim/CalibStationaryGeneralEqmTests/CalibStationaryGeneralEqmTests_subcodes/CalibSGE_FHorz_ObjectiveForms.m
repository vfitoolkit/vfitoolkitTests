function output=CalibSGE_FHorz_ObjectiveForms(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 2: the forms of the OLG calibration objective, and the option plumbing of CalibrateOLGModel.
%
% (a) Direct calls of CalibrateOLGModel_Joint_objectivefn at (truth, equilibrium prices) with targets 5% off the model, weights w,
%     relativeGEweight=7 and multiGEweights=[1,2,0.5]. With m, g the objective's own moment and GE vectors (vectoroutput=1):
%       sum_squared, scalar:          sum(w.*(m-t).^2)/2 + 7*sum(mw.*g.^2)       (2 = number of calibrated, non-GE, parameters)
%       sum_squared, vector:          [sqrt(w).*(m-t)/sqrt(2); sqrt(7)*sqrt(mw).*g; 0]   (CalibrateOLGModel hands the objective sqrt(w))
%       sum_logratiosquared, scalar:  sum(w.*log(m./t).^2)/2 + 7*sum(mw.*g.^2)
%       sum_logratiosquared, vector:  [sqrt(w).*log(m./t)/sqrt(2); sqrt(7)*sqrt(mw).*g; 0]
%       MethodOfMoments, scalar:      (m-t)'*W*(m-t)/2 + 7*sum(mw.*g.^2)         (W a full symmetric positive definite matrix)
%       MethodOfMoments, vector:      sum of squares equal to the scalar, with the objective handed chol(W,'upper')
%     all exact; each vector form's sum of squares (the penalty entry aside) is its scalar form (it used to weight the GE conditions by
%     relativeGEweight^2 and multiGEweights^2, and not to divide the targets part).
% (a2) Through CalibrateOLGModel with lsqnonlin, sum_squared with weights and MethodOfMoments with a full W: the reported objective value
%     equals the scalar formula at the solution (exact). Pins the residual weighting and the square root of W (it was the elementwise
%     sqrt, right only for a diagonal W) end to end.
% (b) CalibrateOLGModel_Nested_objectivefn at the truth (inner GE solve): its scalar objective equals the targets part of the joint
%     scalar objective at the nested objective's own GE prices (the two now divide by the same number of parameters). Exact.
% (c) Through CalibrateOLGModel, each run stopped after fminsearch's initial simplex (so it returns its start, the truth): multiGEweights
%     as a column gives exactly what the row gives (it used to make a matrix); caliboptions.logmoments as a scalar, one entry per
%     target, one entry per target name, and by name (incl. AgeConditionalStats, which read the AllStats names) becomes the expected
%     vector, and calibsummary.currentmomentvec equals the (log) targets.
% (d) One evaluation with verbose=1 in the scalar and in the vector form: exactly one 'breakdown' line each (the vector form printed one
%     per residual).

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
FnsToEvaluate.earnings=@(d,aprime,a,z,r,kappa_j,alpha,delta,A) ((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d;
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

%% Targets at the truth (computed through the objective itself, vectoroutput=1, below), then perturbed by 5%
[z_gridvals_J,pi_z_J,vfoptionsJ]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,1);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals_J=vfoptionsJ.e_gridvals_J;
simoptionsJ.pi_e_J=vfoptionsJ.pi_e_J;
simoptionsJ.jequaloneDist_usergrids=1;
simoptionsJ.user_z_grid=vfoptionsJ.user_z_grid;
simoptionsJ.user_pi_z=vfoptionsJ.user_pi_z;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsJ);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z_J,Params,simoptionsJ);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,simoptions_full);
TM=struct();
TM.AllStats.K.Mean=gather(AllFull.K.Mean);
TM.AllStats.H.Mean=gather(AllFull.H.Mean);
TM.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
t=gather(LCPFull.K.Mean); t(1)=NaN; TM.AgeConditionalStats.K.Mean=t; % (everyone has zero assets at age 1)
TM.AgeConditionalStats.H.Mean=gather(LCPFull.H.Mean);
TMp=TM; % 5% off
TMp.AllStats.K.Mean=1.05*TM.AllStats.K.Mean;
TMp.AllStats.H.Mean=1.05*TM.AllStats.H.Mean;
TMp.AllStats.earnings.Gini=1.05*TM.AllStats.earnings.Gini;
TMp.AgeConditionalStats.K.Mean=1.05*TM.AgeConditionalStats.K.Mean;
TMp.AgeConditionalStats.H.Mean=1.05*TM.AgeConditionalStats.H.Mean;

%% The joint objective, set up as CalibrateOLGModel sets it up
CalibParamNames={'beta','varphi'};
calibparamsvec0=[Params.beta; Params.varphi];
caliboptionsN=struct(); % (the nested objective's: the constraint vectors for the two calibrated parameters only)
caliboptionsN.verbose=0;
caliboptionsN.constrainpositive={};
caliboptionsN.constrainpositivemethod='softplus';
caliboptionsN.constrain0to1={};
caliboptionsN.constrainAtoB={};
caliboptionsN.relativeGEweight=7;
caliboptionsN.useCustomModelStats=0;
caliboptionsN.calibrateshocks=0;
caliboptionsN.whichcombos=1;
caliboptionsN.nCalibParamsNonGE=2;
[calibparamsvec0,caliboptionsN]=ParameterConstraints_TransformParamsToUnconstrained(calibparamsvec0,[0;1;2],CalibParamNames,caliboptionsN,1);
heteroagentoptionsJ=struct();
heteroagentoptionsJ.outputGEstruct=2;
mw=[1,2,0.5];
heteroagentoptionsJ.multiGEweights=mw;
heteroagentoptionsJ.constrainpositive={};
heteroagentoptionsJ.constrainpositivemethod='softplus';
heteroagentoptionsJ.constrain0to1={};
heteroagentoptionsJ.constrainAtoB={};
heteroagentoptionsJ.useCustomModelStats=0;
heteroagentoptionsJ.useintermediateEqns=0;
GEparamsvec0=[Params.r; Params.Tr; Params.tau_c];
[GEparamsvec0,heteroagentoptionsJ]=ParameterConstraints_TransformParamsToUnconstrained(GEparamsvec0,0:1:3,GEPriceParamNames,heteroagentoptionsJ,1);
calibparamsvecJ=[calibparamsvec0; GEparamsvec0];
CalibParamNamesJ=[CalibParamNames,GEPriceParamNames];
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
[targetmomentvec,usingallstats,usinglcp,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,acsmomentnames,acscummomentsizes,ACStats_whichstats,FnsToEvaluate_ACStats,cmsmomentnames,cmscummomentsizes,selectors,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec]=SetupTargetMoments_FHorz(TMp,FnsToEvaluate,0,N_j,simoptions);
caliboptions.selectors=selectors;
caliboptionsN.selectors=selectors;
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptionsN.logmoments=caliboptions.logmoments;
actual=~isnan(targetmomentvec);
t=targetmomentvec(actual);
nmom=sum(actual);
w=1+(1:nmom)'/nmom; % arbitrary positive weights
Wfull=eye(nmom)+0.3*ones(nmom,nmom)/nmom; % a full symmetric positive definite weighting matrix

% One argument list for every joint call (everything but calibparamsvec and caliboptions), one for the nested
jointargs={CalibParamNamesJ,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,5,[0;1;2;3;4;5],zeros(5,1),zeros(N_j,1)};
GEinner.heteroagentoptions=heteroagentoptions;
GEinner.heteroagentoptions.outputGEstruct=2;
GEinner.z_grid=z_grid;
GEinner.pi_z=pi_z;
GEinner.vfoptions=vfoptions;
GEinner.simoptions=simoptions;
nestedargs={CalibParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqns,usingallstats,usinglcp,usingcustomstats,targetmomentvec,allstatmomentnames,acsmomentnames,cmsmomentnames,allstatcummomentsizes,acscummomentsizes,cmscummomentsizes,AllStats_whichstats,ACStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_ACStats,usingautocorr,autocorrmomentnames,autocorrcummomentsizes,FnsToEvaluate_AutoCorr,autocorrtimehorizons,usingcrosssec,crosssecmomentnames,crossseccummomentsizes,FnsToEvaluate_CrossSec,usingagecrosssec,agecrosssecmomentnames,agecrossseccummomentsizes,FnsToEvaluate_AgeCrossSec,2,[0;1;2],zeros(2,1),zeros(N_j,1)};

%% (a) the joint objective forms
co=caliboptions; co.vectoroutput=1; co.metric='sum_squared'; co.weights=w;
v1=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
m=v1(1:nmom); g=v1(nmom+1:end);
GEpart=7*sum(mw(:).*g.^2);
co=caliboptions; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w;
O=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms (a) sum_squared scalar vs formula (targets part /2 calibrated params, + 7*weighted GE), this should be zero: %.3e (value %.6e) \n',abs(O-(sum(w.*(m-t).^2)/2+GEpart)),O)
co=caliboptions; co.vectoroutput=2; co.metric='sum_squared'; co.weights=sqrt(w);
O=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
f=[sqrt(w).*(m-t)/sqrt(2); sqrt(7)*sqrt(mw(:)).*g; 0]; % (the square roots of the scalar form's weights and divisor)
if length(O)==length(f)
    d=max(abs(O(:)-f));
else
    d=1;
end
fprintf('ObjectiveForms (a) sum_squared residual vector vs formula, this should be zero: %.3e \n',d)
fprintf('ObjectiveForms (a) sum_squared residual vector: its sum of squares vs the scalar objective, this should be zero: %.3e \n',abs(sum(O(1:end-1).^2)-(sum(w.*(m-t).^2)/2+GEpart)))
co=caliboptions; co.vectoroutput=0; co.metric='sum_logratiosquared'; co.weights=w;
O=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms (a) sum_logratiosquared scalar vs formula, this should be zero: %.3e \n',abs(O-(sum(w.*log(m./t).^2)/2+GEpart)))
co=caliboptions; co.vectoroutput=2; co.metric='sum_logratiosquared'; co.weights=sqrt(w);
O=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
f=[sqrt(w).*log(m./t)/sqrt(2); sqrt(7)*sqrt(mw(:)).*g; 0];
if length(O)==length(f)
    d=max(abs(O(:)-f));
else
    d=1;
end
fprintf('ObjectiveForms (a) sum_logratiosquared residual vector vs formula, this should be zero: %.3e \n',d)
fprintf('ObjectiveForms (a) sum_logratiosquared residual vector: its sum of squares vs the scalar objective, this should be zero: %.3e \n',abs(sum(O(1:end-1).^2)-(sum(w.*log(m./t).^2)/2+GEpart)))
co=caliboptions; co.vectoroutput=0; co.metric='MethodOfMoments'; co.weights=Wfull;
O=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms (a) MethodOfMoments scalar (full weighting matrix) vs formula, this should be zero: %.3e \n',abs(O-((m-t)'*Wfull*(m-t)/2+GEpart)))
co=caliboptions; co.vectoroutput=2; co.metric='MethodOfMoments'; co.weights=chol(Wfull,'upper'); % (what CalibrateOLGModel now hands the objective for MethodOfMoments under lsqnonlin)
O=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms (a) MethodOfMoments residual vector (upper Cholesky factor of W): its sum of squares vs the scalar objective, this should be zero: %.3e \n',abs(sum(O(1:end-1).^2)-((m-t)'*Wfull*(m-t)/2+GEpart)))

%% (a2) through CalibrateOLGModel with lsqnonlin: the objective value it reports (lsqnonlin's sum of squares) against the scalar formula
% at the solution, from calibsummary's moments and GE conditions (default relativeGEweight=10, multiGEweights ones). This pins both the
% weighting of the GE conditions in the residual vector and the square root of a full weighting matrix, end to end.
for mm=1:2
    caliboptionsR=struct();
    caliboptionsR.verbose=0;
    caliboptionsR.fminalgo=8;
    caliboptionsR.toleranceparams=1e-2; % (the identity holds wherever lsqnonlin stops)
    caliboptionsR.toleranceobjective=1e-4;
    if mm==1
        caliboptionsR.metric='sum_squared'; caliboptionsR.weights=w; lab='sum_squared with weights';
    else
        caliboptionsR.metric='MethodOfMoments'; caliboptionsR.weights=Wfull; lab='MethodOfMoments with a full weighting matrix';
    end
    try
        [~,sR]=CalibrateOLGModel(CalibParamNames,TMp,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,struct(),caliboptionsR,vfoptions,simoptions);
        dR=sR.currentmomentvec-sR.targetmomentvec;
        if mm==1
            formula=sum(w.*dR.^2)/2+10*sum(sR.GeneralEqmConditions.^2);
        else
            formula=dR'*Wfull*dR/2+10*sum(sR.GeneralEqmConditions.^2);
        end
        d=abs(sR.objvalue-formula);
        failed=0;
    catch err
        fprintf('ObjectiveForms (a2) %s errored: %s \n',lab,err.message)
        d=1; formula=NaN; failed=1;
    end
    fprintf('ObjectiveForms (a2) lsqnonlin, %s: reported objective vs the scalar formula at the solution, this should be zero: %.3e; did not run, this should be zero: %.3e (formula %.6e) \n',lab,d,failed,formula)
end

%% (b) nested vs the targets part of joint at the nested objective's GE prices
co=caliboptionsN; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w;
[On,GEn]=CalibrateOLGModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},co,GEinner,vfoptionsJ,simoptionsJ);
On=gather(On);
calibparamsvecJn=[calibparamsvec0; GEn.prices]; % (no constraints, so the GE prices enter untransformed)
co=caliboptions; co.vectoroutput=1; co.metric='sum_squared'; co.weights=w;
vn=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJn,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
co=caliboptions; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w;
Ojn=gather(CalibrateOLGModel_Joint_objectivefn(calibparamsvecJn,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
Ojn_targets=Ojn-7*sum(mw(:).*vn(nmom+1:end).^2);
fprintf('ObjectiveForms (b) nested scalar objective vs the targets part of joint at the same GE prices, this should be zero: %.3e (nested %.6e; nested GE prices vs equilibrium, should be near zero: %.3e) \n',abs(On-Ojn_targets),On,max(abs(GEn.prices-[Params.r;Params.Tr;Params.tau_c])))

%% (c) through CalibrateOLGModel: multiGEweights as a column, and the four forms of caliboptions.logmoments
% Each run stops after fminsearch's initial simplex (tolerances far above any change), so it returns its start, the truth.
caliboptionsC=struct();
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
heteroagentoptionsC=struct();
heteroagentoptionsC.multiGEweights=mw;
[CalibParamsRow,calibsummaryRow]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsC,caliboptionsC,vfoptions,simoptions);
heteroagentoptionsC.multiGEweights=mw'; % a column
try
    [CalibParamsCol,calibsummaryCol]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsC,caliboptionsC,vfoptions,simoptions);
    d=max(abs([CalibParamsCol.beta-CalibParamsRow.beta, CalibParamsCol.varphi-CalibParamsRow.varphi, CalibParamsCol.r-CalibParamsRow.r, CalibParamsCol.Tr-CalibParamsRow.Tr, CalibParamsCol.tau_c-CalibParamsRow.tau_c, calibsummaryCol.objvalue-calibsummaryRow.objvalue]));
    failed=0;
catch err
    fprintf('ObjectiveForms (c) multiGEweights as a column errored: %s \n',err.message)
    d=1; failed=1;
end
fprintf('ObjectiveForms (c) multiGEweights column vs row through CalibrateOLGModel, this should be zero: %.3e; did not run, this should be zero: %.3e \n',d,failed)
fprintf('ObjectiveForms (c) the run returned its start (the truth), this should be zero: %.3e; currentmomentvec vs targets, this should be zero: %.3e; GeneralEqmConditions at the truth, should be near zero: %.3e \n',max(abs([CalibParamsRow.beta-Params.beta,CalibParamsRow.varphi-Params.varphi])),max(abs(calibsummaryRow.currentmomentvec-calibsummaryRow.targetmomentvec)),max(abs(calibsummaryRow.GeneralEqmConditions)))

% logmoments: targets of the logged entries are given as logs
ntm=length(targetmomentvec); % (target entries incl. the NaN one: K.Mean, H.Mean, Gini, then N_j ages of K, N_j ages of H)
namesizes=[1,1,1,N_j,N_j];
for form=1:4
    caliboptionsL=caliboptionsC;
    if form==1
        caliboptionsL.logmoments=1;
        expected=ones(ntm,1);
        lab='scalar';
    elseif form==2
        expected=zeros(ntm,1); expected([1,3])=1; expected(3+2:3+N_j)=1; % (entry 3+1 is the NaN age-1 target)
        caliboptionsL.logmoments=expected';
        lab='one entry per target';
    elseif form==3
        caliboptionsL.logmoments=[1,0,1,1,0];
        expected=repelem([1;0;1;1;0],namesizes(:));
        lab='one entry per target name';
    else
        caliboptionsL.logmoments=struct();
        caliboptionsL.logmoments.AllStats.K.Mean=1;
        caliboptionsL.logmoments.AllStats.earnings.Gini=1;
        caliboptionsL.logmoments.AgeConditionalStats.K.Mean=1;
        expected=repelem([1;0;1;1;0],namesizes(:));
        lab='by name';
    end
    TMl=TM;
    TMl.AllStats.K.Mean=log(TM.AllStats.K.Mean)*expected(1)+TM.AllStats.K.Mean*(1-expected(1));
    TMl.AllStats.H.Mean=log(TM.AllStats.H.Mean)*expected(2)+TM.AllStats.H.Mean*(1-expected(2));
    TMl.AllStats.earnings.Gini=log(TM.AllStats.earnings.Gini)*expected(3)+TM.AllStats.earnings.Gini*(1-expected(3));
    e=expected(4:3+N_j)'; TMl.AgeConditionalStats.K.Mean=log(TM.AgeConditionalStats.K.Mean).*e+TM.AgeConditionalStats.K.Mean.*(1-e);
    e=expected(4+N_j:end)'; TMl.AgeConditionalStats.H.Mean=log(TM.AgeConditionalStats.H.Mean).*e+TM.AgeConditionalStats.H.Mean.*(1-e);
    try
        [CalibParamsL,calibsummaryL]=CalibrateOLGModel(CalibParamNames,TMl,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsC,caliboptionsL,vfoptions,simoptions);
        if length(calibsummaryL.logmoments)==ntm
            dl=max(abs(calibsummaryL.logmoments(:)-expected));
        else
            dl=1;
        end
        dm=max(abs(calibsummaryL.currentmomentvec-calibsummaryL.targetmomentvec));
        dp=max(abs([CalibParamsL.beta-Params.beta,CalibParamsL.varphi-Params.varphi]));
        failed=0;
    catch err
        fprintf('ObjectiveForms (c) logmoments %s errored: %s \n',lab,err.message)
        dl=1; dm=1; dp=1; failed=1;
    end
    fprintf('ObjectiveForms (c) logmoments %s: logmoments vector vs expected, this should be zero: %.3e; currentmomentvec vs (log) targets, this should be zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',lab,dl,dm,dp,failed)
end

%% (d) the verbose prints: one 'breakdown' line in each form
co=caliboptions; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w; co.verbose=1;
txt0=evalc('CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);');
co=caliboptions; co.vectoroutput=2; co.metric='sum_squared'; co.weights=sqrt(w); co.verbose=1;
txt2=evalc('CalibrateOLGModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);');
n0=length(strfind(txt0,'breakdown')); n2=length(strfind(txt2,'breakdown'));
fprintf('ObjectiveForms (d) verbose: number of breakdown lines minus one, scalar form, this should be zero: %.3e; vector form, this should be zero: %.3e \n',abs(n0-1),abs(n2-1))
fprintf('ObjectiveForms (d) verbose: the scalar breakdown line prints relativeGEweight (7), this should be zero: %.3e \n',1-contains(txt0,'relativeGEweight=7.000000'))

output=struct();
output.m=m; output.g=g;

end
