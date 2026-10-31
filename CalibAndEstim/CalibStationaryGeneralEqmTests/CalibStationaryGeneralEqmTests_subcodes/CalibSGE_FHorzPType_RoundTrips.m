function output=CalibSGE_FHorzPType_RoundTrips(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 15: round trips through CalibrateOLGModel_PType. Part 4 with two permanent types, on the one-price
% model only (tau=G=Tr=tau_c=0, only r and CapitalMarket): two types double the cost of every solve.
%
% Targets are generated at the truth (beta=0.96, varphi=[1;1.2] per type, type masses [0.6;0.4]) and its general eqm. Every run starts
% away from the truth and from the equilibrium, and must recover both (close to zero); calibsummary's GE condition must be small (near
% zero), and a fresh HeteroAgentStationaryEqm_Case1_FHorz_PType solve at the calibrated parameters must reproduce CalibParams.r (near).
%   (A) joint lsqnonlin, beta and varphi given as a column over the types (the driver calibrates it as one parameter per type and
%       returns it as the column it was given), against mean capital and each type's mean hours
%   (B) the same, nested
%   (C) the type masses through ParametrizeParamsFn: the share of the first type, calibrated from 0.5 (truth 0.6), joint lsqnonlin,
%       against mean capital and the grouped mean hours (the types differ in varphi, so the grouped mean pins the share)
%   (D) omitcalibparam on part of the row kappa_j with a row mask: the run (stopped after fminsearch's initial simplex) returns its start
%       and keeps the omitted ages [before 2026-10-29 a row omit mask was never transposed, and the check of it broadcast row against column]
% (A)-(C) run on 51 hours points: hours are a discrete choice, and on the bank's 11 the stepped objective stalls lsqnonlin (see part 4).

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
GeneralEqmEqns1.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GEPriceParamNames1={'r'};

heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;

n_d3=51;
d_grid3=linspace(0,1,n_d3)';

%% The truth and its general eqm (one price, 51 hours points)
Params1=Params;
Params1.tau=0; Params1.G=0; Params1.Tr=0; Params1.tau_c=0;
Params1.varphi=[1;1.2]; % per type, as a column over the types
p_eqm1=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d3,n_a,n_z,N_j,Names_i,0,pi_z,d_grid3,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1.r=p_eqm1.r;
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d3,n_a,n_z,N_j,Names_i,d_grid3,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d3,n_a,n_z,N_j,Names_i,pi_z,Params1,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params1,n_d3,n_a,n_z,N_j,Names_i,d_grid3,a_grid,z_grid,simoptions);
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
    try
        tt=tic;
        [CalibParams,calibsummary]=CalibrateOLGModel_PType({'beta','varphi'},TM1,n_d3,n_a,n_z,N_j,Names_i,d_grid3,a_grid,z_grid,pi_z,ReturnFn,Params1Start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames1,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptions,vfoptions,simoptions);
        time=toc(tt);
        shapebad=~(isnumeric(CalibParams.varphi) && isequal(size(CalibParams.varphi),[2,1]));
        if shapebad==0
            dpar=max(abs([CalibParams.beta;CalibParams.varphi]-truth));
        else
            dpar=1;
        end
        dpr=abs(CalibParams.r-Params1.r);
        dg=max(abs(calibsummary.GeneralEqmConditions));
        % a fresh general eqm solve at the calibrated parameters
        ParamsCheck=Params1Start;
        ParamsCheck.beta=CalibParams.beta; ParamsCheck.varphi=CalibParams.varphi; ParamsCheck.r=CalibParams.r;
        p_check=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d3,n_a,n_z,N_j,Names_i,0,pi_z,d_grid3,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,ParamsCheck,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
        dfresh=abs(p_check.r-CalibParams.r);
        CP{rr}=[CalibParams.beta;CalibParams.varphi(:);CalibParams.r];
        failed=0;
    catch err
        fprintf('RoundTrips PType %s errored: %s \n',labels{rr},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        shapebad=1; dpar=1; dpr=1; dg=1; dfresh=1; failed=1; time=NaN; CP{rr}=nan(4,1);
    end
    fprintf('RoundTrips PType %s: did not run, this should be zero: %.3e; varphi not returned as a column over the types, this should be zero: %.3e; (beta,varphi) vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e \n',labels{rr},failed,shapebad,dpar,dpr)
    fprintf('RoundTrips PType %s: calibsummary GE condition, should be near zero: %.3e; fresh GE solve at the calibrated parameters vs CalibParams.r, should be near zero: %.3e (%.1f seconds) \n',labels{rr},dg,dfresh,time)
end
fprintf('RoundTrips PType: nested lsqnonlin vs joint lsqnonlin, (beta,varphi,r), this should be close to zero: %.3e \n',max(abs(CP{1}-CP{2})))

%% (C) the type masses through ParametrizeParamsFn (the share of the first type)
TMC=struct();
TMC.AllStats.K.Mean=TM1.AllStats.K.Mean;
TMC.AllStats.H.Mean=gather(AllStats.H.Mean); % grouped: 0.6*H_A+0.4*H_B, the types differ in varphi, so it pins the share
ParamsC=Params1;
ParamsC.shareA=0.5; % start (truth 0.6)
ParamsC.ptypemass=[0.5;0.5];
ParamsC.r=1.1*Params1.r;
ParametrizeWeightsFn=@(Parameters) CalibSGE_ParametrizeWeightsFn(Parameters);
caliboptions=caliboptionsBase;
try
    tt=tic;
    [CalibParamsC,calibsummaryC]=CalibrateOLGModel_PType({'shareA'},TMC,n_d3,n_a,n_z,N_j,Names_i,d_grid3,a_grid,z_grid,pi_z,ReturnFn,ParamsC,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames1,PTypeDistParamNames,ParametrizeWeightsFn,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptions,vfoptions,simoptions);
    time=toc(tt);
    dshare=abs(CalibParamsC.shareA-Params1.ptypemass(1));
    dpr=abs(CalibParamsC.r-Params1.r);
    dg=max(abs(calibsummaryC.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('RoundTrips PType (C) type masses through ParametrizeParamsFn errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dshare=1; dpr=1; dg=1; failed=1; time=NaN;
end
fprintf('RoundTrips PType (C) type masses through ParametrizeParamsFn, share from 0.5 (truth 0.6), joint lsqnonlin: did not run, this should be zero: %.3e; share vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero: %.3e (%.1f seconds) \n',failed,dshare,dpr,dg,time)

%% (D) omitcalibparam on part of the row kappa_j, with a row omit mask (stopped after fminsearch's initial simplex, on the bank's grids)
Params1D=Params1;
Params1D.varphi=1;
p_eqmD=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,0,pi_z,d_grid,a_grid,z_grid,jequaloneDist,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1D,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1D.r=p_eqmD.r;
TMD=struct();
TMD.AllStats.K.Mean=1; % (the values do not matter: the run returns its start)
TMD.AllStats.H.Mean=1;
caliboptionsD=struct();
caliboptionsD.verbose=0;
caliboptionsD.fminalgo=1;
caliboptionsD.toleranceparams=1e3;
caliboptionsD.toleranceobjective=1e3;
omitmask=Params1D.kappa_j; omitmask(2:4)=NaN; % a row, like kappa_j: calibrate ages 2 to 4
caliboptionsD.omitcalibparam.kappa_j=omitmask;
try
    [CalibParamsD,~]=CalibrateOLGModel_PType({'beta','kappa_j'},TMD,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1D,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,GEPriceParamNames1,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptionsD,vfoptions,simoptions);
    dk=max(abs(CalibParamsD.kappa_j(:)-Params1D.kappa_j(:)));
    failed=0;
catch err
    fprintf('RoundTrips PType (D) omitcalibparam errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dk=1; failed=1;
end
fprintf('RoundTrips PType (D) omitcalibparam (row mask): did not run, this should be zero: %.3e; kappa_j returned its start (omitted ages kept), this should be zero: %.3e \n',failed,dk)

output=struct();
output.CP=CP;

end
