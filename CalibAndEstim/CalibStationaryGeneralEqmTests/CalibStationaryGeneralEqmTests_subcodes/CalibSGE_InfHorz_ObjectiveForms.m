function output=CalibSGE_InfHorz_ObjectiveForms(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 8: the forms of the BIHA calibration objective, and the option plumbing of CalibrateBIHAModel.
% The InfHorz counterpart of part 2 (see there):
% (a) the joint objective forms against their formulas (targets 5% off, weights w, relativeGEweight=7, multiGEweights=[1,2,0.5]); exact;
%     each vector form's sum of squares is its scalar form
% (a2) through CalibrateBIHAModel with lsqnonlin (sum_squared with weights; MethodOfMoments with a full W): the reported objective equals
%     the scalar formula at the solution; exact
% (b) the nested objective (which could never run: undefined Params and z_grid) equals the targets part of the joint objective at the
%     nested objective's own GE prices; exact
% (c) through CalibrateBIHAModel, runs stopped after fminsearch's initial simplex: multiGEweights as a column gives exactly what the row
%     gives; caliboptions.logmoments as a scalar, one entry per target, one per target name (the vector form referred to acsmomentnames,
%     which InfHorz does not have) and by name becomes the expected vector
% (d) verbose=1: exactly one 'breakdown' line in the scalar and in the vector form

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

%% Targets at the truth, then 5% off
[z_gridvals,pi_zJ,vfoptionsJ]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0);
simoptionsJ=simoptions;
simoptionsJ.e_gridvals=vfoptionsJ.e_gridvals;
simoptionsJ.pi_e=vfoptionsJ.pi_e;
vfoptionsJ.alreadygridvals=1;
simoptionsJ.alreadygridvals=1;
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsJ);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_zJ,simoptionsJ,Params);
simoptions_full=simoptionsJ;
simoptions_full.whichstats=ones(1,7);
AllFull=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,simoptions_full);
TM=struct();
TM.AllStats.K.Mean=gather(AllFull.K.Mean);
TM.AllStats.H.Mean=gather(AllFull.H.Mean);
TM.AllStats.earnings.Gini=gather(AllFull.earnings.Gini);
TMp=TM;
TMp.AllStats.K.Mean=1.05*TM.AllStats.K.Mean;
TMp.AllStats.H.Mean=1.05*TM.AllStats.H.Mean;
TMp.AllStats.earnings.Gini=1.05*TM.AllStats.earnings.Gini;

%% Both objectives, set up as CalibrateBIHAModel sets them up
CalibParamNames={'beta','varphi'};
calibparamsvec0=[Params.beta; Params.varphi];
caliboptionsN=struct();
caliboptionsN.verbose=0;
caliboptionsN.constrainpositive={};
caliboptionsN.constrainpositivemethod='softplus';
caliboptionsN.constrain0to1={};
caliboptionsN.constrainAtoB={};
caliboptionsN.relativeGEweight=7;
caliboptionsN.useCustomModelStats=0;
caliboptionsN.calibrateshocks=0;
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
[targetmomentvec,usingallstats,usingautocorr,usingcrosssec,usingcustomstats,allstatmomentnames,allstatcummomentsizes,AllStats_whichstats,FnsToEvaluate_AllStats,autocorrmomentnames,autocorrcummomentsizes,AutoCorrStats_whichstats,FnsToEvaluate_AutoCorrStats,crosssecmomentnames,crossseccummomentsizes,CrossSecStats_whichstats,FnsToEvaluate_CrossSecStats,cmsmomentnames,cmscummomentsizes,selectors,autocorrtimehorizons]=SetupTargetMoments_InfHorz(TMp,FnsToEvaluate,0,simoptions);
caliboptions.whichcombos=1; % (CalibrateBIHAModel's default; it hands the objective the selectors and horizons from SetupTargetMoments_InfHorz)
caliboptions.selectors=selectors;
caliboptions.autocorrtimehorizons=autocorrtimehorizons;
caliboptionsN.whichcombos=1;
caliboptionsN.selectors=selectors;
caliboptionsN.autocorrtimehorizons=autocorrtimehorizons;
caliboptions.logmoments=zeros(length(targetmomentvec),1);
caliboptionsN.logmoments=caliboptions.logmoments;
t=targetmomentvec;
nmom=length(t);
w=1+(1:nmom)'/nmom; % arbitrary positive weights
Wfull=eye(nmom)+0.3*ones(nmom,nmom)/nmom; % a full symmetric positive definite weighting matrix
statargs={usingallstats,usingautocorr,usingcrosssec,usingcustomstats,targetmomentvec,allstatmomentnames,autocorrmomentnames,crosssecmomentnames,cmsmomentnames,allstatcummomentsizes,autocorrcummomentsizes,crossseccummomentsizes,cmscummomentsizes,AllStats_whichstats,AutoCorrStats_whichstats,CrossSecStats_whichstats,FnsToEvaluate_AllStats,FnsToEvaluate_AutoCorrStats,FnsToEvaluate_CrossSecStats};
jointargs=[{[CalibParamNames,GEPriceParamNames],n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqnsCell,GEeqnNames,GeneralEqmEqnParamNames},statargs,{5,[0;1;2;3;4;5],zeros(5,1),zeros(1,1)}];
nestedargs=[{CalibParamNames,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_zJ,ReturnFn,Params,DiscountFactorParamNames,GEPriceParamNames,[],FnsToEvaluate,GeneralEqmEqns},statargs,{2,[0;1;2],zeros(2,1),zeros(1,1)}];
GEinner.heteroagentoptions=heteroagentoptions;
GEinner.heteroagentoptions.outputGEstruct=2;
GEinner.z_grid=z_grid;
GEinner.pi_z=pi_z;
GEinner.vfoptions=vfoptions;
GEinner.simoptions=simoptions;

%% (a) the joint objective forms
co=caliboptions; co.vectoroutput=1; co.metric='sum_squared'; co.weights=w;
v1=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
m=v1(1:nmom); g=v1(nmom+1:end);
GEpart=7*sum(mw(:).*g.^2);
co=caliboptions; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w;
O=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms InfHorz (a) sum_squared scalar vs formula (targets part /2 calibrated params, + 7*weighted GE), this should be zero: %.3e (value %.6e) \n',abs(O-(sum(w.*(m-t).^2)/2+GEpart)),O)
co=caliboptions; co.vectoroutput=2; co.metric='sum_squared'; co.weights=sqrt(w);
O=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
f=[sqrt(w).*(m-t)/sqrt(2); sqrt(7)*sqrt(mw(:)).*g; 0]; % (the square roots of the scalar form's weights and divisor)
if length(O)==length(f)
    d=max(abs(O(:)-f));
else
    d=1;
end
fprintf('ObjectiveForms InfHorz (a) sum_squared residual vector vs formula, this should be zero: %.3e \n',d)
fprintf('ObjectiveForms InfHorz (a) sum_squared residual vector: its sum of squares vs the scalar objective, this should be zero: %.3e \n',abs(sum(O(1:end-1).^2)-(sum(w.*(m-t).^2)/2+GEpart)))
co=caliboptions; co.vectoroutput=0; co.metric='sum_logratiosquared'; co.weights=w;
O=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms InfHorz (a) sum_logratiosquared scalar vs formula, this should be zero: %.3e \n',abs(O-(sum(w.*log(m./t).^2)/2+GEpart)))
co=caliboptions; co.vectoroutput=2; co.metric='sum_logratiosquared'; co.weights=sqrt(w);
O=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
f=[sqrt(w).*log(m./t)/sqrt(2); sqrt(7)*sqrt(mw(:)).*g; 0];
if length(O)==length(f)
    d=max(abs(O(:)-f));
else
    d=1;
end
fprintf('ObjectiveForms InfHorz (a) sum_logratiosquared residual vector vs formula, this should be zero: %.3e \n',d)
fprintf('ObjectiveForms InfHorz (a) sum_logratiosquared residual vector: its sum of squares vs the scalar objective, this should be zero: %.3e \n',abs(sum(O(1:end-1).^2)-(sum(w.*log(m./t).^2)/2+GEpart)))
co=caliboptions; co.vectoroutput=0; co.metric='MethodOfMoments'; co.weights=Wfull;
O=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms InfHorz (a) MethodOfMoments scalar (full weighting matrix) vs formula, this should be zero: %.3e \n',abs(O-((m-t)'*Wfull*(m-t)/2+GEpart)))
co=caliboptions; co.vectoroutput=2; co.metric='MethodOfMoments'; co.weights=chol(Wfull,'upper'); % (what CalibrateBIHAModel now hands the objective for MethodOfMoments under lsqnonlin)
O=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
fprintf('ObjectiveForms InfHorz (a) MethodOfMoments residual vector (upper Cholesky factor of W): its sum of squares vs the scalar objective, this should be zero: %.3e \n',abs(sum(O(1:end-1).^2)-((m-t)'*Wfull*(m-t)/2+GEpart)))

%% (a2) through CalibrateBIHAModel with lsqnonlin: the objective value it reports (lsqnonlin's sum of squares) against the scalar formula
% at the solution, from calibsummary's moments and GE conditions (default relativeGEweight=10, multiGEweights ones)
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
        [~,sR]=CalibrateBIHAModel(CalibParamNames,TMp,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,struct(),caliboptionsR,vfoptions,simoptions);
        dR=sR.currentmomentvec-sR.targetmomentvec;
        if mm==1
            formula=sum(w.*dR.^2)/2+10*sum(sR.GeneralEqmConditions.^2);
        else
            formula=dR'*Wfull*dR/2+10*sum(sR.GeneralEqmConditions.^2);
        end
        d=abs(sR.objvalue-formula);
        failed=0;
    catch err
        fprintf('ObjectiveForms InfHorz (a2) %s errored: %s \n',lab,err.message)
        d=1; formula=NaN; failed=1;
    end
    fprintf('ObjectiveForms InfHorz (a2) lsqnonlin, %s: reported objective vs the scalar formula at the solution, this should be zero: %.3e; did not run, this should be zero: %.3e (formula %.6e) \n',lab,d,failed,formula)
end

%% (b) nested vs the targets part of joint at the nested objective's GE prices
co=caliboptionsN; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w;
try
    [On,GEn]=CalibrateBIHAModel_Nested_objectivefn(calibparamsvec0,nestedargs{:},co,GEinner,vfoptionsJ,simoptionsJ);
    On=gather(On);
    calibparamsvecJn=[calibparamsvec0; GEn.prices]; % (no constraints, so the GE prices enter untransformed)
    co=caliboptions; co.vectoroutput=1; co.metric='sum_squared'; co.weights=w;
    vn=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJn,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
    co=caliboptions; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w;
    Ojn=gather(CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJn,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ));
    Ojn_targets=Ojn-7*sum(mw(:).*vn(nmom+1:end).^2);
    d=abs(On-Ojn_targets); dp=max(abs(GEn.prices-[Params.r;Params.Tr;Params.tau_c]));
    failed=0;
catch err
    fprintf('ObjectiveForms InfHorz (b) errored: %s \n',err.message)
    On=NaN; d=1; dp=1; failed=1;
end
fprintf('ObjectiveForms InfHorz (b) nested scalar objective vs the targets part of joint at the same GE prices, this should be zero: %.3e; did not run, this should be zero: %.3e (nested %.6e; nested GE prices vs equilibrium, should be near zero: %.3e) \n',d,failed,On,dp)

%% (c) through CalibrateBIHAModel: multiGEweights as a column, and the four forms of caliboptions.logmoments
caliboptionsC=struct();
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
heteroagentoptionsC=struct();
heteroagentoptionsC.multiGEweights=mw;
[CalibParamsRow,calibsummaryRow]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsC,caliboptionsC,vfoptions,simoptions);
heteroagentoptionsC.multiGEweights=mw'; % a column
try
    [CalibParamsCol,calibsummaryCol]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsC,caliboptionsC,vfoptions,simoptions);
    d=max(abs([CalibParamsCol.beta-CalibParamsRow.beta, CalibParamsCol.varphi-CalibParamsRow.varphi, CalibParamsCol.r-CalibParamsRow.r, CalibParamsCol.Tr-CalibParamsRow.Tr, CalibParamsCol.tau_c-CalibParamsRow.tau_c, calibsummaryCol.objvalue-calibsummaryRow.objvalue]));
    failed=0;
catch err
    fprintf('ObjectiveForms InfHorz (c) multiGEweights as a column errored: %s \n',err.message)
    d=1; failed=1;
end
fprintf('ObjectiveForms InfHorz (c) multiGEweights column vs row through CalibrateBIHAModel, this should be zero: %.3e; did not run, this should be zero: %.3e \n',d,failed)
fprintf('ObjectiveForms InfHorz (c) the run returned its start (the truth), this should be zero: %.3e; currentmomentvec vs targets, this should be zero: %.3e; GeneralEqmConditions at the truth, should be near zero: %.3e \n',max(abs([CalibParamsRow.beta-Params.beta,CalibParamsRow.varphi-Params.varphi])),max(abs(calibsummaryRow.currentmomentvec-calibsummaryRow.targetmomentvec)),max(abs(calibsummaryRow.GeneralEqmConditions)))

% logmoments: a Lorenz curve target (unlogged) is added so that the number of target names differs from the number of entries
nL=numel(AllFull.earnings.LorenzCurve);
for form=1:4
    caliboptionsL=caliboptionsC;
    TML=TM;
    TML.AllStats.K.Mean=log(TM.AllStats.K.Mean);
    TML.AllStats.earnings.Gini=log(TM.AllStats.earnings.Gini);
    if form==1
        % scalar: every target is logged, so no Lorenz curve here (it starts at zero) and H.Mean is given as a log too
        TML.AllStats.H.Mean=log(TM.AllStats.H.Mean);
        caliboptionsL.logmoments=1;
        expected=ones(3,1);
        lab='scalar';
    else
        TML.AllStats.earnings.LorenzCurve=gather(AllFull.earnings.LorenzCurve);
    end
    if form==2
        expected=[1;0;1;zeros(nL,1)];
        caliboptionsL.logmoments=expected';
        lab='one entry per target';
    elseif form==3
        caliboptionsL.logmoments=[1,0,1,0];
        expected=[1;0;1;zeros(nL,1)];
        lab='one entry per target name';
    elseif form==4
        caliboptionsL.logmoments=struct();
        caliboptionsL.logmoments.AllStats.K.Mean=1;
        caliboptionsL.logmoments.AllStats.earnings.Gini=1;
        expected=[1;0;1;zeros(nL,1)];
        lab='by name';
    end
    try
        [CalibParamsL,calibsummaryL]=CalibrateBIHAModel(CalibParamNames,TML,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,struct(),caliboptionsL,vfoptions,simoptions);
        if length(calibsummaryL.logmoments)==length(expected)
            dl=max(abs(calibsummaryL.logmoments(:)-expected));
        else
            dl=1;
        end
        dm=max(abs(calibsummaryL.currentmomentvec-calibsummaryL.targetmomentvec));
        dp=max(abs([CalibParamsL.beta-Params.beta,CalibParamsL.varphi-Params.varphi]));
        failed=0;
    catch err
        fprintf('ObjectiveForms InfHorz (c) logmoments %s errored: %s \n',lab,err.message)
        dl=1; dm=1; dp=1; failed=1;
    end
    fprintf('ObjectiveForms InfHorz (c) logmoments %s: logmoments vector vs expected, this should be zero: %.3e; currentmomentvec vs (log) targets, this should be zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',lab,dl,dm,dp,failed)
end

%% (d) the verbose prints: one 'breakdown' line in each form
co=caliboptions; co.vectoroutput=0; co.metric='sum_squared'; co.weights=w; co.verbose=1;
txt0=evalc('CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);');
co=caliboptions; co.vectoroutput=2; co.metric='sum_squared'; co.weights=sqrt(w); co.verbose=1;
txt2=evalc('CalibrateBIHAModel_Joint_objectivefn(calibparamsvecJ,jointargs{:},co,heteroagentoptionsJ,vfoptionsJ,simoptionsJ);');
n0=length(strfind(txt0,'breakdown')); n2=length(strfind(txt2,'breakdown'));
fprintf('ObjectiveForms InfHorz (d) verbose: number of breakdown lines minus one, scalar form, this should be zero: %.3e; vector form, this should be zero: %.3e \n',abs(n0-1),abs(n2-1))
fprintf('ObjectiveForms InfHorz (d) verbose: the scalar breakdown line prints relativeGEweight (7), this should be zero: %.3e \n',1-contains(txt0,'relativeGEweight=7.000000'))

output=struct();
output.m=m; output.g=g;

end
