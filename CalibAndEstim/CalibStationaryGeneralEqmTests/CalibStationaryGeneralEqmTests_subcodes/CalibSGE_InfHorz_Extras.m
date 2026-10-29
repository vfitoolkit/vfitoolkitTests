function output=CalibSGE_InfHorz_Extras(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 11: the options of CalibrateBIHAModel that change how the model is set up, each in joint and in
% nested mode. The InfHorz counterpart of part 5 (see there), without jequaloneDist.
% (a) vfoptions.ExogShockFn with rho_z calibrated, target the autocorrelation of earnings (AutoCorrTransProbs): at the truth (exact,
%     joint; near, nested) and from rho_z=0.8 (close), one-price model
% (b) heteroagentoptions.CustomModelStats alone (the joint objective did not keep V, so it crashed): GE conditions equal the baseline
%     (its CapitalMarket adds a flag that is 1 when CustomModelStats is handed internal vfoptions, which HeteroAgentStationaryEqm_InfHorz did)
% (c) caliboptions.CustomModelStats under CustomModelStats_usergrids=1 and =0: the hand mean of assets equals its target, the z_grid
%     handed in is the expected one, and the vfoptions handed in are the user's (they were the internal ones, alreadygridvals=1)
% (d) heteroagentoptions.intermediateEqns: GE conditions equal the baseline

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
FnsToEvaluate.earnings=@(d,aprime,a,z,r,alpha,delta,A) ((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d;
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;

%% The truth and its general eqm
p_eqm=HeteroAgentStationaryEqm_InfHorz(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params);
AllStats=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
TM=struct();
TM.AllStats.K.Mean=gather(AllStats.K.Mean);
TM.AllStats.H.Mean=gather(AllStats.H.Mean);
CalibParamNames={'beta','varphi'};

caliboptionsC=struct(); % runs that stop after fminsearch's initial simplex, so return their start
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;

% The baseline: joint and nested, started at the truth
caliboptionsJ=caliboptionsC; caliboptionsJ.jointoptimization=1;
caliboptionsN=caliboptionsC; caliboptionsN.jointoptimization=0;
[~,summJ0]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsJ,vfoptions,simoptions);
try
    [CPN0,summN0]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsN,vfoptions,simoptions);
    dgN=max(abs(summN0.GeneralEqmConditions)); failed=0;
catch err
    fprintf('Extras InfHorz baseline nested errored: %s \n',err.message)
    CPN0.r=NaN; CPN0.Tr=NaN; CPN0.tau_c=NaN; dgN=1; failed=1;
end
fprintf('Extras InfHorz baseline: joint at the truth, moments vs targets, this should be zero: %.3e; GE conditions, should be near zero: %.3e; nested did not run, this should be zero: %.3e; nested GE conditions, should be near zero: %.3e \n',max(abs(summJ0.currentmomentvec-summJ0.targetmomentvec)),max(abs(summJ0.GeneralEqmConditions)),failed,dgN)

%% (a) ExogShockFn with rho_z calibrated, one-price model
Params1=Params;
Params1.tau=0; Params1.G=0; Params1.Tr=0; Params1.tau_c=0;
GEPriceParamNames1={'r'};
GeneralEqmEqns1.CapitalMarket=GeneralEqmEqns.CapitalMarket;
vfoptionsZ=vfoptions;
vfoptionsZ.ExogShockFn=@(rho_z,sigma_z) CalibSGE_ExogShockFn(rho_z,sigma_z);
simoptionsZ=simoptions;
simoptionsZ.ExogShockFn=vfoptionsZ.ExogShockFn; % set in both, as a user would
p_eqm1=HeteroAgentStationaryEqm_InfHorz(n_d,n_a,n_z,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1,DiscountFactorParamNames,[],[],[],GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1.r=p_eqm1.r;
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params1);
AllStats1=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params1,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
AC1=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params1,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions);
TMZ=struct();
TMZ.AllStats.K.Mean=gather(AllStats1.K.Mean);
TMZ.AutoCorrTransProbs.earnings.AutoCorrelation=gather(AC1.earnings.AutoCorrelation);
for jo=1:-1:0
    caliboptionsZ=caliboptionsC; caliboptionsZ.jointoptimization=jo;
    try
        [CPz,summZ]=CalibrateBIHAModel({'rho_z'},TMZ,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],GEPriceParamNames1,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ,simoptionsZ);
        dm=max(abs(summZ.currentmomentvec-summZ.targetmomentvec)); dp=abs(CPz.rho_z-0.9);
        failed=0;
    catch err
        fprintf('Extras InfHorz (a) ExogShockFn jointoptimization=%i at the truth errored: %s \n',jo,err.message)
        dm=1; dp=1; failed=1;
    end
    if jo==1
        fprintf('Extras InfHorz (a) ExogShockFn, joint, at the truth: moments (incl. AutoCorrTransProbs) vs targets, this should be zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',dm,dp,failed)
    else
        fprintf('Extras InfHorz (a) ExogShockFn, nested, at the truth: moments vs targets, should be near zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',dm,dp,failed)
    end
    caliboptionsZ=struct(); caliboptionsZ.verbose=0;
    caliboptionsZ.fminalgo=8; caliboptionsZ.jointoptimization=jo;
    if jo==1
        caliboptionsZ.finitedifferencestepsize=[1e-2;1e-3]; % one step per unknown, [rho_z; r] (the GE price comes after the calibrated parameter): r needs 1e-3 (at 1e-2 its probes cross beta*(1+r)=1, see part 10), while rho_z needs 1e-2 (it works through the earnings autocorrelation, which over a 1e-3 step in rho_z moves less than the model noise: 1e-3 for both left rho_z 4.3e-2 off, 2026-10-21)
    end
    caliboptionsZ.toleranceparams=1e-6; caliboptionsZ.toleranceobjective=1e-10;
    Params1Start=Params1; Params1Start.rho_z=0.8;
    try
        [CPz,summZ]=CalibrateBIHAModel({'rho_z'},TMZ,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1Start,DiscountFactorParamNames,[],GEPriceParamNames1,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ,simoptionsZ);
        dp=abs(CPz.rho_z-0.9); dr=abs(CPz.r-Params1.r); dg=max(abs(summZ.GeneralEqmConditions));
        failed=0;
    catch err
        fprintf('Extras InfHorz (a) ExogShockFn jointoptimization=%i from rho_z=0.8 errored: %s \n',jo,err.message)
        dp=1; dr=1; dg=1; failed=1;
    end
    if jo==1
        % Joint struggles with this example, and its GE condition is held to the loose bar only. TMZ targets K.Mean as well as the
        % earnings autocorrelation, and in InfHorz K is very steep in r as beta*(1+r) approaches 1 (at the point it stops, +0.001 in r
        % moves K by +0.89, i.e. 10%, against +0.57 for +0.01 in rho_z, and dK/dr roughly doubles for every +0.004 in r). Joint lets r
        % move freely, so the objective is a narrow curved valley along K=target: lsqnonlin hits K exactly, then its steps along the
        % valley push K off target by more (curvature of K in r) than they gain on the autocorrelation, the trust region shrinks and it
        % stops at rho_z=0.852, r 4.3e-3 off, GE condition 4.3e-3, autocorrelation 0.818 vs 0.864 (verbose run, 111 evaluations,
        % 2026-10-22; deterministic, the [1e-2;1e-3] step was not the problem, it never got near the wall). Nested solves r inside each
        % evaluation, so its outer problem is one-dimensional in rho_z and it recovers rho_z to 8e-5. FHorz (OLG) does not have this,
        % as K stays finite in r.
        fprintf('Extras InfHorz (a) ExogShockFn, jointoptimization=%i, from rho_z=0.8: did not run, this should be zero: %.3e; rho_z vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero (loose: joint stalls in the narrow valley along the K target, see the comment above): %.3e \n',jo,failed,dp,dr,dg)
    else
        fprintf('Extras InfHorz (a) ExogShockFn, jointoptimization=%i, from rho_z=0.8: did not run, this should be zero: %.3e; rho_z vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero: %.3e \n',jo,failed,dp,dr,dg)
    end
end

%% (b) heteroagentoptions.CustomModelStats alone
GeneralEqmEqnsB=GeneralEqmEqns;
GeneralEqmEqnsB.CapitalMarket=@(r,Kcustom,N,alpha,delta,A,optionsflag) r-(alpha*A*(Kcustom^(alpha-1))*(N^(1-alpha))-delta)+optionsflag; % (optionsflag is 1 if CustomModelStats was handed the internal vfoptions)
heteroagentoptionsB=heteroagentoptions;
heteroagentoptionsB.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions) CalibSGE_InfHorz_GECustomStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions);
for jo=1:-1:0
    try
        if jo==1
            [~,summB]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsB,heteroagentoptionsB,caliboptionsJ,vfoptions,simoptions);
            d=max(abs(summB.GeneralEqmConditions-summJ0.GeneralEqmConditions));
        else
            [CPB,~]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsB,heteroagentoptionsB,caliboptionsN,vfoptions,simoptions);
            d=max(abs([CPB.r-CPN0.r,CPB.Tr-CPN0.Tr,CPB.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras InfHorz (b) heteroagentoptions.CustomModelStats jointoptimization=%i errored: %s \n',jo,err.message)
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras InfHorz (b) heteroagentoptions.CustomModelStats alone, joint: did not run, this should be zero: %.3e; GE conditions vs baseline, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras InfHorz (b) heteroagentoptions.CustomModelStats alone, nested: did not run, this should be zero: %.3e; GE prices vs baseline, should be near zero: %.3e \n',failed,d)
    end
end

%% (c) caliboptions.CustomModelStats, both usergrids settings
[z_gridvals,~,~]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions,3,0);
TMC=TM;
TMC.CustomModelStats.meanK=TM.AllStats.K.Mean;
TMC.CustomModelStats.zgridcheck=0;
TMC.CustomModelStats.optionscheck=0;
for ug=1:-1:0
    for jo=1:-1:0
        caliboptionsCC=caliboptionsC;
        caliboptionsCC.jointoptimization=jo;
        caliboptionsCC.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibSGE_InfHorz_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        caliboptionsCC.CustomModelStats_usergrids=ug;
        if ug==1
            caliboptionsCC.CalibSGE_expected_zgrid=z_grid;
        else
            caliboptionsCC.CalibSGE_expected_zgrid=gather(z_gridvals);
        end
        try
            [~,summC]=CalibrateBIHAModel(CalibParamNames,TMC,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsCC,vfoptions,simoptions);
            dK=abs(summC.currentmomentvec(end-2)-summC.targetmomentvec(end-2));
            dz=summC.currentmomentvec(end-1);
            dopt=summC.currentmomentvec(end);
            failed=0;
        catch err
            fprintf('Extras InfHorz (c) CustomModelStats usergrids=%i jointoptimization=%i errored: %s \n',ug,jo,err.message)
            dK=1; dz=1; dopt=1; failed=1;
        end
        if jo==1
            fprintf('Extras InfHorz (c) caliboptions.CustomModelStats, usergrids=%i, joint: did not run, this should be zero: %.3e; hand mean of assets vs target, this should be zero: %.3e; z_grid handed in vs expected, this should be zero: %.3e; vfoptions handed in are the internal ones, this should be zero: %.3e \n',ug,failed,dK,dz,dopt)
        else
            fprintf('Extras InfHorz (c) caliboptions.CustomModelStats, usergrids=%i, nested: did not run, this should be zero: %.3e; hand mean of assets vs target, should be near zero: %.3e; z_grid handed in vs expected, this should be zero: %.3e; vfoptions handed in are the internal ones, this should be zero: %.3e \n',ug,failed,dK,dz,dopt)
        end
    end
end

%% (d) intermediateEqns
GeneralEqmEqnsD=GeneralEqmEqns;
GeneralEqmEqnsD.GovBudget=@(tau,w,N,Tr) tau*w*N-Tr;
heteroagentoptionsD=heteroagentoptions;
heteroagentoptionsD.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
for jo=1:-1:0
    try
        if jo==1
            [~,summD]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsD,heteroagentoptionsD,caliboptionsJ,vfoptions,simoptions);
            d=max(abs(summD.GeneralEqmConditions-summJ0.GeneralEqmConditions));
        else
            [CPD,~]=CalibrateBIHAModel(CalibParamNames,TM,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsD,heteroagentoptionsD,caliboptionsN,vfoptions,simoptions);
            d=max(abs([CPD.r-CPN0.r,CPD.Tr-CPN0.Tr,CPD.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras InfHorz (d) intermediateEqns jointoptimization=%i errored: %s \n',jo,err.message)
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras InfHorz (d) intermediateEqns, joint: did not run, this should be zero: %.3e; GE conditions vs baseline, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras InfHorz (d) intermediateEqns, nested: did not run, this should be zero: %.3e; GE prices vs baseline, should be near zero: %.3e \n',failed,d)
    end
end

output=struct();
output.summJ0=summJ0;

end
