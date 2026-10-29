function output=CalibSGE_FHorz_Extras(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 5: the options of CalibrateOLGModel that change how the model is set up, each in joint and in
% nested mode.
%
% Most runs here stop after fminsearch's initial simplex (tolerances far above any change), so a run started at the truth returns the
% truth, and calibsummary holds the model at the truth: an identity against the baseline (the same model set up without the option).
% (a) vfoptions.ExogShockFn with rho_z calibrated (the shock grids are rebuilt at every evaluation), target the autocorrelation of
%     earnings by age (an AutoCorrTransProbs target, new to CalibrateOLGModel): started at the truth the moments equal the targets
%     (exact); from rho_z=0.8 lsqnonlin recovers rho_z (close), joint and nested (one-price model, for speed)
% (b) heteroagentoptions.CustomModelStats ALONE (the joint objective did not keep V, so it crashed), feeding the CapitalMarket condition
%     a hand-computed K: the GE conditions equal the baseline ones (exact, joint; near, nested)
% (c) caliboptions.CustomModelStats under CustomModelStats_usergrids=1 and =0: a hand-computed mean of assets equals its target and
%     the z_grid handed in is the user's (=1) or the internal joint grid (=0); joint and nested
% (d) heteroagentoptions.intermediateEqns (w computed once, used by GovBudget): the GE conditions equal the baseline ones (exact, joint;
%     near, nested)
% (e) jequaloneDist as a function: moments and GE conditions equal the matrix version (exact, joint; near, nested)

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
FnsToEvaluate.earnings=@(d,aprime,a,z,r,kappa_j,alpha,delta,A) ((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d;
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;

%% The truth and its general eqm
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
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
[~,summJ0]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsJ,vfoptions,simoptions);
[CPN0,summN0]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsN,vfoptions,simoptions);
fprintf('Extras baseline: joint at the truth, moments vs targets, this should be zero: %.3e; GE conditions, should be near zero: %.3e; nested GE conditions, should be near zero: %.3e \n',max(abs(summJ0.currentmomentvec-summJ0.targetmomentvec)),max(abs(summJ0.GeneralEqmConditions)),max(abs(summN0.GeneralEqmConditions)))

%% (a) ExogShockFn with rho_z calibrated, one-price model
Params1=Params;
Params1.tau=0; Params1.G=0; Params1.Tr=0; Params1.tau_c=0;
GEPriceParamNames1={'r'};
GeneralEqmEqns1.CapitalMarket=GeneralEqmEqns.CapitalMarket;
vfoptionsZ=vfoptions;
vfoptionsZ.ExogShockFn=@(rho_z,sigma_z) CalibSGE_ExogShockFn(rho_z,sigma_z);
simoptionsZ=simoptions;
simoptionsZ.ExogShockFn=vfoptionsZ.ExogShockFn; % set in both, as a user would
p_eqm1=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1,DiscountFactorParamNames,[],[],[],GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1.r=p_eqm1.r;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params1,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
AC1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params1,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
TMZ=struct();
TMZ.AllStats.K.Mean=gather(AllStats1.K.Mean);
t=gather(AC1.earnings.AutoCorrelation); t(1)=NaN; TMZ.AutoCorrTransProbs.earnings.AutoCorrelation=t;
for jo=1:-1:0
    caliboptionsZ=caliboptionsC; caliboptionsZ.jointoptimization=jo;
    try
        [CPz,summZ]=CalibrateOLGModel({'rho_z'},TMZ,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames1,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ,simoptionsZ);
        dm=max(abs(summZ.currentmomentvec-summZ.targetmomentvec)); dp=abs(CPz.rho_z-0.9);
        failed=0;
    catch err
        fprintf('Extras (a) ExogShockFn jointoptimization=%i at the truth errored: %s \n',jo,err.message)
        dm=1; dp=1; failed=1;
    end
    if jo==1
        fprintf('Extras (a) ExogShockFn, joint, at the truth: moments (incl. AutoCorrTransProbs) vs targets, this should be zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',dm,dp,failed)
    else
        fprintf('Extras (a) ExogShockFn, nested, at the truth: moments vs targets, should be near zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',dm,dp,failed)
    end
    caliboptionsZ=struct(); caliboptionsZ.verbose=0; caliboptionsZ.fminalgo=8; caliboptionsZ.jointoptimization=jo;
    caliboptionsZ.toleranceparams=1e-6; caliboptionsZ.toleranceobjective=1e-10;
    Params1Start=Params1; Params1Start.rho_z=0.8;
    try
        [CPz,summZ]=CalibrateOLGModel({'rho_z'},TMZ,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1Start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames1,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ,simoptionsZ);
        dp=abs(CPz.rho_z-0.9); dr=abs(CPz.r-Params1.r); dg=max(abs(summZ.GeneralEqmConditions));
        failed=0;
    catch err
        fprintf('Extras (a) ExogShockFn jointoptimization=%i from rho_z=0.8 errored: %s \n',jo,err.message)
        dp=1; dr=1; dg=1; failed=1;
    end
    fprintf('Extras (a) ExogShockFn, jointoptimization=%i, from rho_z=0.8: did not run, this should be zero: %.3e; rho_z vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero: %.3e \n',jo,failed,dp,dr,dg)
end

%% (b) heteroagentoptions.CustomModelStats alone
GeneralEqmEqnsB=GeneralEqmEqns;
GeneralEqmEqnsB.CapitalMarket=@(r,Kcustom,N,alpha,delta,A,optionsflag) r-(alpha*A*(Kcustom^(alpha-1))*(N^(1-alpha))-delta)+optionsflag; % (optionsflag is 1 if CustomModelStats was handed the internal vfoptions)
heteroagentoptionsB=heteroagentoptions;
heteroagentoptionsB.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions) CalibSGE_FHorz_GECustomStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions);
for jo=1:-1:0
    try
        if jo==1
            [~,summB]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsB,heteroagentoptionsB,caliboptionsJ,vfoptions,simoptions);
            d=max(abs(summB.GeneralEqmConditions-summJ0.GeneralEqmConditions));
        else
            [CPB,~]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsB,heteroagentoptionsB,caliboptionsN,vfoptions,simoptions);
            d=max(abs([CPB.r-CPN0.r,CPB.Tr-CPN0.Tr,CPB.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras (b) heteroagentoptions.CustomModelStats jointoptimization=%i errored: %s \n',jo,err.message)
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras (b) heteroagentoptions.CustomModelStats alone, joint: did not run, this should be zero: %.3e; GE conditions vs baseline, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras (b) heteroagentoptions.CustomModelStats alone, nested: did not run, this should be zero: %.3e; GE prices vs baseline, should be near zero: %.3e \n',failed,d)
    end
end

%% (c) caliboptions.CustomModelStats, both usergrids settings
[z_gridvals_J,~,~]=ExogShockSetup_FHorz(n_z,z_grid,pi_z,N_j,Params,vfoptions,3,0);
TMC=TM;
TMC.CustomModelStats.meanK=TM.AllStats.K.Mean;
TMC.CustomModelStats.zgridcheck=0;
for ug=1:-1:0
    for jo=1:-1:0
        caliboptionsCC=caliboptionsC;
        caliboptionsCC.jointoptimization=jo;
        caliboptionsCC.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibSGE_FHorz_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        caliboptionsCC.CustomModelStats_usergrids=ug;
        if ug==1
            caliboptionsCC.CalibSGE_expected_zgrid=z_grid;
        else
            caliboptionsCC.CalibSGE_expected_zgrid=gather(z_gridvals_J);
        end
        try
            [~,summC]=CalibrateOLGModel(CalibParamNames,TMC,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsCC,vfoptions,simoptions);
            dK=abs(summC.currentmomentvec(end-1)-summC.targetmomentvec(end-1));
            dz=summC.currentmomentvec(end);
            failed=0;
        catch err
            fprintf('Extras (c) CustomModelStats usergrids=%i jointoptimization=%i errored: %s \n',ug,jo,err.message)
            dK=1; dz=1; failed=1;
        end
        if jo==1
            fprintf('Extras (c) caliboptions.CustomModelStats, usergrids=%i, joint: did not run, this should be zero: %.3e; hand mean of assets vs target, this should be zero: %.3e; z_grid handed in vs expected, this should be zero: %.3e \n',ug,failed,dK,dz)
        else
            fprintf('Extras (c) caliboptions.CustomModelStats, usergrids=%i, nested: did not run, this should be zero: %.3e; hand mean of assets vs target, should be near zero: %.3e; z_grid handed in vs expected, this should be zero: %.3e \n',ug,failed,dK,dz)
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
            [~,summD]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsD,heteroagentoptionsD,caliboptionsJ,vfoptions,simoptions);
            d=max(abs(summD.GeneralEqmConditions-summJ0.GeneralEqmConditions));
        else
            [CPD,~]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqnsD,heteroagentoptionsD,caliboptionsN,vfoptions,simoptions);
            d=max(abs([CPD.r-CPN0.r,CPD.Tr-CPN0.Tr,CPD.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras (d) intermediateEqns jointoptimization=%i errored: %s \n',jo,err.message)
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras (d) intermediateEqns, joint: did not run, this should be zero: %.3e; GE conditions vs baseline, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras (d) intermediateEqns, nested: did not run, this should be zero: %.3e; GE prices vs baseline, should be near zero: %.3e \n',failed,d)
    end
end

%% (e) jequaloneDist as a function
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z) CalibSGE_jequaloneDistFn(a_grid,z_grid,n_a,n_z);
simoptionsE=simoptions;
simoptionsE.a_grid=a_grid; % (jequaloneDist as a function needs a_grid and z_grid in simoptions)
simoptionsE.z_grid=z_grid;
for jo=1:-1:0
    try
        if jo==1
            [~,summE]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsJ,vfoptions,simoptionsE);
            d=max(abs([summE.currentmomentvec-summJ0.currentmomentvec; summE.GeneralEqmConditions-summJ0.GeneralEqmConditions]));
        else
            [CPE,~]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDistFn,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsN,vfoptions,simoptionsE);
            d=max(abs([CPE.r-CPN0.r,CPE.Tr-CPN0.Tr,CPE.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras (e) jequaloneDist as a function jointoptimization=%i errored: %s \n',jo,err.message)
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras (e) jequaloneDist as a function, joint: did not run, this should be zero: %.3e; moments and GE conditions vs the matrix version, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras (e) jequaloneDist as a function, nested: did not run, this should be zero: %.3e; GE prices vs the matrix version, should be near zero: %.3e \n',failed,d)
    end
end

output=struct();
output.summJ0=summJ0;

end
