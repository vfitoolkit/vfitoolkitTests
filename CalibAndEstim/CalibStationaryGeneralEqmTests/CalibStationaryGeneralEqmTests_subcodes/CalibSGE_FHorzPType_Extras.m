function output=CalibSGE_FHorzPType_Extras(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 16: the options of CalibrateOLGModel_PType that change how the model is set up, each in joint and
% in nested mode. Part 5 with two permanent types.
%
% Most runs here stop after fminsearch's initial simplex (tolerances far above any change), so a run started at the truth returns the
% truth, and calibsummary holds the model at the truth: an identity against the baseline (the same model set up without the option).
% (a) vfoptions.ExogShockFn with rho_z calibrated (the shock grids are rebuilt at every evaluation), targets mean capital and the
%     autocorrelation of earnings by age: started at the truth the moments equal the targets (exact, joint; near, nested); from
%     rho_z=0.8 joint lsqnonlin recovers rho_z (close) (one-price model, for speed)
% (a2) the same ExogShockFn given per type (a struct over Names_i): the shock-dependence detection read it with
%     getAnonymousFnInputNames, which errors on a struct [before 2026-10-29]; started at the truth it runs and the moments equal the
%     targets (exact, joint)
% (b) heteroagentoptions.CustomModelStats ALONE, feeding the CapitalMarket condition a hand-computed K: the GE conditions equal the
%     baseline ones (exact, joint; near, nested)
% (c) caliboptions.CustomModelStats under CustomModelStats_usergrids=1 and =0: a hand-computed mean of assets equals its target, the
%     z_grid handed in is the user's (=1) or the internal joint grid (=0), and the vfoptions handed in are the user's; joint and nested
% (d) heteroagentoptions.intermediateEqns (w computed once, used by GovBudget): the GE conditions equal the baseline ones (exact, joint;
%     near, nested)
% (jequaloneDist as a function is not covered: StationaryDist_Case1_FHorz_PType, which the calibration objectives use, sums jequaloneDist
% as an array, so it does not take a function; HeteroAgentStationaryEqm_Case1_FHorz_PType does, through its own per-type loop.)

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
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
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
[~,summJ0]=CalibrateOLGModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsJ,vfoptions,simoptions);
failed=0;
try
    [CPN0,summN0]=CalibrateOLGModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsN,vfoptions,simoptions);
    dgN=max(abs(summN0.GeneralEqmConditions));
catch err
    fprintf('Extras PType baseline nested errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    CPN0.r=Params.r; CPN0.Tr=Params.Tr; CPN0.tau_c=Params.tau_c; dgN=1; failed=1;
end
fprintf('Extras PType baseline: joint at the truth, moments vs targets, this should be zero: %.3e; GE conditions, should be near zero: %.3e; nested did not run, this should be zero: %.3e; nested GE conditions, should be near zero: %.3e \n',max(abs(summJ0.currentmomentvec-summJ0.targetmomentvec)),max(abs(summJ0.GeneralEqmConditions)),failed,dgN)

%% (a) ExogShockFn with rho_z calibrated, one-price model
Params1=Params;
Params1.tau=0; Params1.G=0; Params1.Tr=0; Params1.tau_c=0;
GEPriceParamNames1={'r'};
GeneralEqmEqns1.CapitalMarket=GeneralEqmEqns.CapitalMarket;
vfoptionsZ=vfoptions;
vfoptionsZ.ExogShockFn=@(rho_z,sigma_z) CalibSGE_ExogShockFn(rho_z,sigma_z);
simoptionsZ=simoptions;
simoptionsZ.ExogShockFn=vfoptionsZ.ExogShockFn; % set in both, as a user would
p_eqm1=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1.r=p_eqm1.r;
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params1,simoptions);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params1,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AC1=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params1,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions);
TMZ=struct();
TMZ.AllStats.K.Mean=gather(AllStats1.K.Mean);
t=gather(AC1.earnings.AutoCorrelation); t(1)=NaN; TMZ.AutoCorrTransProbs.earnings.AutoCorrelation=t;
for jo=1:-1:0
    caliboptionsZ=caliboptionsC; caliboptionsZ.jointoptimization=jo;
    try
        [CPz,summZ]=CalibrateOLGModel_PType({'rho_z'},TMZ,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames1,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ,simoptionsZ);
        dm=max(abs(summZ.currentmomentvec-summZ.targetmomentvec)); dp=abs(CPz.rho_z-0.9);
        failed=0;
    catch err
        fprintf('Extras PType (a) ExogShockFn jointoptimization=%i at the truth errored: %s \n',jo,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        dm=1; dp=1; failed=1;
    end
    if jo==1
        fprintf('Extras PType (a) ExogShockFn, joint, at the truth: moments (incl. AutoCorrTransProbs) vs targets, this should be zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',dm,dp,failed)
    else
        fprintf('Extras PType (a) ExogShockFn, nested, at the truth: moments vs targets, should be near zero: %.3e; returned its start, this should be zero: %.3e; did not run, this should be zero: %.3e \n',dm,dp,failed)
    end
end
caliboptionsZ=struct(); caliboptionsZ.verbose=0; caliboptionsZ.fminalgo=8; caliboptionsZ.jointoptimization=1;
caliboptionsZ.toleranceparams=1e-6; caliboptionsZ.toleranceobjective=1e-10;
Params1Start=Params1; Params1Start.rho_z=0.8;
try
    [CPz,summZ]=CalibrateOLGModel_PType({'rho_z'},TMZ,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1Start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames1,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ,simoptionsZ);
    dp=abs(CPz.rho_z-0.9); dr=abs(CPz.r-Params1.r); dg=max(abs(summZ.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('Extras PType (a) ExogShockFn joint from rho_z=0.8 errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dp=1; dr=1; dg=1; failed=1;
end
fprintf('Extras PType (a) ExogShockFn, joint lsqnonlin, from rho_z=0.8: did not run, this should be zero: %.3e; rho_z vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero: %.3e \n',failed,dp,dr,dg)

%% (a2) the ExogShockFn given per type (a struct over Names_i)
vfoptionsZ2=vfoptions;
vfoptionsZ2.ExogShockFn.(Names_i{1})=@(rho_z,sigma_z) CalibSGE_ExogShockFn(rho_z,sigma_z);
vfoptionsZ2.ExogShockFn.(Names_i{2})=@(rho_z,sigma_z) CalibSGE_ExogShockFn(rho_z,sigma_z);
simoptionsZ2=simoptions;
simoptionsZ2.ExogShockFn=vfoptionsZ2.ExogShockFn;
caliboptionsZ=caliboptionsC; caliboptionsZ.jointoptimization=1;
try
    [CPz,summZ]=CalibrateOLGModel_PType({'rho_z'},TMZ,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames1,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsZ,vfoptionsZ2,simoptionsZ2);
    dm=max(abs(summZ.currentmomentvec-summZ.targetmomentvec)); dp=abs(CPz.rho_z-0.9);
    failed=0;
catch err
    fprintf('Extras PType (a2) ExogShockFn per type errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dm=1; dp=1; failed=1;
end
fprintf('Extras PType (a2) ExogShockFn per type (a struct over the types), joint, at the truth: did not run, this should be zero: %.3e; moments vs targets, this should be zero: %.3e; returned its start, this should be zero: %.3e \n',failed,dm,dp)

%% (b) heteroagentoptions.CustomModelStats alone
GeneralEqmEqnsB=GeneralEqmEqns;
GeneralEqmEqnsB.CapitalMarket=@(r,Kcustom,N,alpha,delta,A,optionsflag) r-(alpha*A*(Kcustom^(alpha-1))*(N^(1-alpha))-delta)+optionsflag; % (optionsflag is 1 if CustomModelStats was handed the internal vfoptions)
heteroagentoptionsB=heteroagentoptions;
heteroagentoptionsB.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions) CalibSGE_FHorzPType_GECustomStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions);
for jo=1:-1:0
    try
        if jo==1
            [~,summB]=CalibrateOLGModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsB,heteroagentoptionsB,caliboptionsJ,vfoptions,simoptions);
            d=max(abs(summB.GeneralEqmConditions-summJ0.GeneralEqmConditions));
        else
            [CPB,~]=CalibrateOLGModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsB,heteroagentoptionsB,caliboptionsN,vfoptions,simoptions);
            d=max(abs([CPB.r-CPN0.r,CPB.Tr-CPN0.Tr,CPB.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras PType (b) heteroagentoptions.CustomModelStats jointoptimization=%i errored: %s \n',jo,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras PType (b) heteroagentoptions.CustomModelStats alone, joint: did not run, this should be zero: %.3e; GE conditions vs baseline, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras PType (b) heteroagentoptions.CustomModelStats alone, nested: did not run, this should be zero: %.3e; GE prices vs baseline, should be near zero: %.3e \n',failed,d)
    end
end

%% (c) caliboptions.CustomModelStats, both usergrids settings
[z_gridvals_J,~,~]=ExogShockSetup_FHorz_PType(n_z,z_grid,pi_z,N_j,Names_i,Params,vfoptions,3,0);
TMC=TM;
TMC.CustomModelStats.meanK=TM.AllStats.K.Mean;
TMC.CustomModelStats.zgridcheck=0;
TMC.CustomModelStats.optionscheck=0;
for ug=1:-1:0
    for jo=1:-1:0
        caliboptionsCC=caliboptionsC;
        caliboptionsCC.jointoptimization=jo;
        caliboptionsCC.CustomModelStats=@(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions) CalibSGE_FHorzPType_CustomModelStats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,caliboptions,vfoptions,simoptions);
        caliboptionsCC.CustomModelStats_usergrids=ug;
        if ug==1
            caliboptionsCC.CalibSGE_expected_zgrid=z_grid;
        else
            caliboptionsCC.CalibSGE_expected_zgrid=gather(z_gridvals_J.(Names_i{1}));
        end
        try
            [~,summC]=CalibrateOLGModel_PType(CalibParamNames,TMC,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsCC,vfoptions,simoptions);
            dK=abs(summC.currentmomentvec(end-2)-summC.targetmomentvec(end-2));
            dz=summC.currentmomentvec(end-1);
            dopt=summC.currentmomentvec(end);
            failed=0;
        catch err
            fprintf('Extras PType (c) CustomModelStats usergrids=%i jointoptimization=%i errored: %s \n',ug,jo,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
            dK=1; dz=1; dopt=1; failed=1;
        end
        if jo==1
            fprintf('Extras PType (c) caliboptions.CustomModelStats, usergrids=%i, joint: did not run, this should be zero: %.3e; hand mean of assets vs target, this should be zero: %.3e; z_grid handed in vs expected, this should be zero: %.3e; vfoptions handed in are the internal ones, this should be zero: %.3e \n',ug,failed,dK,dz,dopt)
        else
            fprintf('Extras PType (c) caliboptions.CustomModelStats, usergrids=%i, nested: did not run, this should be zero: %.3e; hand mean of assets vs target, should be near zero: %.3e; z_grid handed in vs expected, this should be zero: %.3e; vfoptions handed in are the internal ones, this should be zero: %.3e \n',ug,failed,dK,dz,dopt)
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
            [~,summD]=CalibrateOLGModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsD,heteroagentoptionsD,caliboptionsJ,vfoptions,simoptions);
            d=max(abs(summD.GeneralEqmConditions-summJ0.GeneralEqmConditions));
        else
            [CPD,~]=CalibrateOLGModel_PType(CalibParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsD,heteroagentoptionsD,caliboptionsN,vfoptions,simoptions);
            d=max(abs([CPD.r-CPN0.r,CPD.Tr-CPN0.Tr,CPD.tau_c-CPN0.tau_c]));
        end
        failed=0;
    catch err
        fprintf('Extras PType (d) intermediateEqns jointoptimization=%i errored: %s \n',jo,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        d=1; failed=1;
    end
    if jo==1
        fprintf('Extras PType (d) intermediateEqns, joint: did not run, this should be zero: %.3e; GE conditions vs baseline, this should be zero: %.3e \n',failed,d)
    else
        fprintf('Extras PType (d) intermediateEqns, nested: did not run, this should be zero: %.3e; GE prices vs baseline, should be near zero: %.3e \n',failed,d)
    end
end

output=struct();
output.summJ0=summJ0;

end
