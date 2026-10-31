function output=CalibSGE_InfHorzPType_RoundTrips(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 20: round trips through CalibrateBIHAModel_PType. Part 10 with two permanent types, on the one-price
% model only (tau=G=Tr=tau_c=0, only r and CapitalMarket): two types double the cost of every solve.
%
% Targets are generated at the truth (beta=0.96, varphi=[1;1.2] per type, type masses [0.6;0.4]) and its general eqm. Every run starts
% away from the truth and from the equilibrium, and must recover both (close to zero); calibsummary's GE condition must be small (near
% zero), and a fresh HeteroAgentStationaryEqm_InfHorz_PType solve at the calibrated parameters must reproduce CalibParams.r (near).
%   (A) joint lsqnonlin, beta and varphi given as a column over the types, against mean capital and each type's mean hours
%   (B) the same, nested
%   (C) the type masses through ParametrizeParamsFn: the share of the first type from 0.5 (truth 0.6), joint lsqnonlin, against mean
%       capital and the grouped mean hours
% The joint runs set caliboptions.finitedifferencestepsize=1e-3: the InfHorz equilibrium sits close below beta*(1+r)=1, and at the
% default 1e-2 the finite-difference probes of beta and r cross that wall (see part 10).

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z; % aggregate effective labor
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
GeneralEqmEqns1.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GEPriceParamNames1={'r'};

heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;

%% The truth and its general eqm (one price)
Params1=Params;
Params1.tau=0; Params1.G=0; Params1.Tr=0; Params1.tau_c=0;
Params1.varphi=[1;1.2]; % per type, as a column over the types
p_eqm1=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1.r=p_eqm1.r;
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,Params1,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params1,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
TM1=struct();
TM1.AllStats.K.Mean=gather(AllStats.K.Mean);
TM1.AllStats.H.typeA.Mean=gather(AllStats.H.typeA.Mean);
TM1.AllStats.H.typeB.Mean=gather(AllStats.H.typeB.Mean);
truth=[Params1.beta;Params1.varphi];
Params1Start=Params1;
Params1Start.beta=0.95; Params1Start.varphi=[1.2;1]; Params1Start.r=1.1*Params1.r;

caliboptionsBase=struct();
caliboptionsBase.verbose=0;
caliboptionsBase.toleranceparams=1e-6;
caliboptionsBase.toleranceobjective=1e-10;
caliboptionsBase.fminalgo=8;

labels={'(A) joint lsqnonlin, varphi per type as a column','(B) nested lsqnonlin, varphi per type as a column'};
CP=cell(1,2);
for rr=1:2
    caliboptions=caliboptionsBase;
    caliboptions.jointoptimization=2-rr; % (A) joint, (B) nested
    if rr==1
        caliboptions.finitedifferencestepsize=1e-3; % (joint: lsqnonlin probes r, see the header)
    end
    try
        tt=tic;
        [CalibParams,calibsummary]=CalibrateBIHAModel_PType({'beta','varphi'},TM1,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1Start,DiscountFactorParamNames,GEPriceParamNames1,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptions,vfoptions,simoptions);
        time=toc(tt);
        shapebad=~(isnumeric(CalibParams.varphi) && isequal(size(CalibParams.varphi),[2,1]));
        if shapebad==0
            dpar=max(abs([CalibParams.beta;CalibParams.varphi]-truth));
        else
            dpar=1;
        end
        dpr=abs(CalibParams.r-Params1.r);
        dg=max(abs(calibsummary.GeneralEqmConditions));
        ParamsCheck=Params1Start;
        ParamsCheck.beta=CalibParams.beta; ParamsCheck.varphi=CalibParams.varphi; ParamsCheck.r=CalibParams.r;
        p_check=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,ParamsCheck,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
        dfresh=abs(p_check.r-CalibParams.r);
        CP{rr}=[CalibParams.beta;CalibParams.varphi(:);CalibParams.r];
        failed=0;
    catch err
        fprintf('RoundTrips InfHorz PType %s errored: %s \n',labels{rr},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        shapebad=1; dpar=1; dpr=1; dg=1; dfresh=1; failed=1; time=NaN; CP{rr}=nan(4,1);
    end
    fprintf('RoundTrips InfHorz PType %s: did not run, this should be zero: %.3e; varphi not returned as a column over the types, this should be zero: %.3e; (beta,varphi) vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e \n',labels{rr},failed,shapebad,dpar,dpr)
    fprintf('RoundTrips InfHorz PType %s: calibsummary GE condition, should be near zero: %.3e; fresh GE solve at the calibrated parameters vs CalibParams.r, should be near zero: %.3e (%.1f seconds) \n',labels{rr},dg,dfresh,time)
end
fprintf('RoundTrips InfHorz PType: nested lsqnonlin vs joint lsqnonlin, (beta,varphi,r), this should be close to zero: %.3e \n',max(abs(CP{1}-CP{2})))

%% (C) the type masses through ParametrizeParamsFn (the share of the first type)
TMC=struct();
TMC.AllStats.K.Mean=TM1.AllStats.K.Mean;
TMC.AllStats.H.Mean=gather(AllStats.H.Mean); % grouped: the types differ in varphi, so it pins the share
ParamsC=Params1;
ParamsC.shareA=0.5; % start (truth 0.6)
ParamsC.ptypemass=[0.5;0.5];
ParamsC.r=1.1*Params1.r;
ParametrizeWeightsFn=@(Parameters) CalibSGE_ParametrizeWeightsFn(Parameters);
caliboptions=caliboptionsBase;
caliboptions.finitedifferencestepsize=1e-3; % (joint, see the header)
try
    tt=tic;
    [CalibParamsC,calibsummaryC]=CalibrateBIHAModel_PType({'shareA'},TMC,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsC,DiscountFactorParamNames,GEPriceParamNames1,PTypeDistParamNames,ParametrizeWeightsFn,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptions,vfoptions,simoptions);
    time=toc(tt);
    dshare=abs(CalibParamsC.shareA-Params1.ptypemass(1));
    dpr=abs(CalibParamsC.r-Params1.r);
    dg=max(abs(calibsummaryC.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('RoundTrips InfHorz PType (C) type masses through ParametrizeParamsFn errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dshare=1; dpr=1; dg=1; failed=1; time=NaN;
end
fprintf('RoundTrips InfHorz PType (C) type masses through ParametrizeParamsFn, share from 0.5 (truth 0.6), joint lsqnonlin: did not run, this should be zero: %.3e; share vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero: %.3e (%.1f seconds) \n',failed,dshare,dpr,dg,time)

output=struct();
output.CP=CP;

end
