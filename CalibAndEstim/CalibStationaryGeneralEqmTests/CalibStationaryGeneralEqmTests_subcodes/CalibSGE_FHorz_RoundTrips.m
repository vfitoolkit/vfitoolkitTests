function output=CalibSGE_FHorz_RoundTrips(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 4: round trips through CalibrateOLGModel.
%
% Targets are generated at the truth (beta=0.96, varphi=1) and its general eqm. Every run starts away from the truth and from the
% equilibrium, and must recover both (close to zero); a fresh HeteroAgentStationaryEqm solve at the calibrated parameters must
% reproduce the GE prices in CalibParams (near zero), and calibsummary's GE conditions must be small (near zero).
%   three prices, joint, on 51 hours points (see below): (A) lsqnonlin, (C) lsqnonlin with constrain0to1 on beta, constrainpositive
%       on varphi and (heteroagentoptions) on Tr, constrainAtoB on tau_c; (D) lsqnonlin with three ages of the row kappa_j calibrated
%       under a row omit mask (the omit mask was never transposed) against the age profile of earnings
%   one price (tau=G=Tr=tau_c=0, only r and CapitalMarket), on the bank's 11 hours points: (E) nested lsqnonlin, (F) nested
%       fminsearch, (G) joint lsqnonlin, (B) joint fminsearch, and nested vs joint on the same targets (close)
% The lsqnonlin runs use the default caliboptions.finitedifferencestepsize=1e-2. Step sweep of (A) (a diagnostic block, run 2026-10-21
% and then removed; (beta,varphi) error / GE conditions): 1e-2 1.1e-5 / 4e-5; 1e-3 7.2e-2 / 6.1e-3; 1e-4 0.14 / 1.2e-2; sqrt(eps) 0.2 /
% 1.8e-2. Below 1e-2 the step drowns in model noise here. InfHorz (part 10) needs 1e-3 instead (its equilibrium sits close below
% beta*(1+r)=1, which FHorz does not have), so no single default suits both, and the toolkit default stays 1e-2.

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

% The three-price runs use 51 hours points instead of the bank's 11. Hours are a discrete choice (the grid interpolation layer is on
% assets only), so on 11 points mean hours, N and the objective move in steps, and with five unknowns both lsqnonlin and fminsearch
% stalled ("step smaller than the step tolerance", beta/varphi 5e-2 to 2e-1 off, GE conditions 4e-3 to 2e-2). The diagnostic run of
% 2026-10-13 settled it: on 51 points lsqnonlin recovers the truth to about 1e-5 with relativeGEweight 10 and 50 alike, so it is the
% grid and not the lsqnonlin weighting of the GE conditions. fminsearch with five unknowns stayed 9e-3 off even on 51 points (GE
% conditions 1.3e-3), so (B) runs on the one-price model (three unknowns) instead.
n_d3=51;
d_grid3=linspace(0,1,n_d3)';

%% The truth and its general eqm (three prices, 51 hours points)
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d3,n_a,n_z,N_j,0,pi_z,d_grid3,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end
truth=[Params.beta;Params.varphi];
peq=[Params.r;Params.Tr;Params.tau_c];
[~,Policy]=ValueFnIter_Case1_FHorz(n_d3,n_a,n_z,N_j,d_grid3,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d3,n_a,n_z,N_j,pi_z,Params,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d3,n_a,n_z,N_j,d_grid3,a_grid,z_grid,simoptions);
LCP=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d3,n_a,n_z,N_j,d_grid3,a_grid,z_grid,simoptions);
TM=struct();
TM.AllStats.K.Mean=gather(AllStats.K.Mean);
TM.AllStats.H.Mean=gather(AllStats.H.Mean);
t=gather(LCP.K.Mean); t(1)=NaN; TM.AgeConditionalStats.K.Mean=t;

ParamsStart=Params; % away from the truth and from the equilibrium
ParamsStart.beta=0.95; ParamsStart.varphi=1.2;
ParamsStart.r=1.1*Params.r; ParamsStart.Tr=0.9*Params.Tr; ParamsStart.tau_c=1.1*Params.tau_c;

CalibParamNames={'beta','varphi'};
caliboptionsBase=struct();
caliboptionsBase.verbose=0;
caliboptionsBase.toleranceparams=1e-6;
caliboptionsBase.toleranceobjective=1e-10;

labels={'(A) joint lsqnonlin','(C) joint lsqnonlin with constraints'};
for rr=1:2
    caliboptions=caliboptionsBase;
    heteroagentoptionsR=struct();
    if rr==1
        caliboptions.fminalgo=8;
    else
        caliboptions.fminalgo=8;
        caliboptions.constrain0to1={'beta'};
        caliboptions.constrainpositive={'varphi'};
        heteroagentoptionsR.constrainpositive={'Tr'};
        heteroagentoptionsR.constrainAtoB={'tau_c'};
        heteroagentoptionsR.constrainAtoBlimits.tau_c=[0,0.5];
    end
    try
        tt=tic;
        [CalibParams,calibsummary]=CalibrateOLGModel(CalibParamNames,TM,n_d3,n_a,n_z,N_j,d_grid3,a_grid,z_grid,pi_z,ReturnFn,ParamsStart,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsR,caliboptions,vfoptions,simoptions);
        time=toc(tt);
        dpar=max(abs([CalibParams.beta;CalibParams.varphi]-truth));
        dpr=max(abs([CalibParams.r;CalibParams.Tr;CalibParams.tau_c]-peq));
        dg=max(abs(calibsummary.GeneralEqmConditions));
        % a fresh general eqm solve at the calibrated parameters
        ParamsCheck=ParamsStart;
        ParamsCheck.beta=CalibParams.beta; ParamsCheck.varphi=CalibParams.varphi;
        ParamsCheck.r=CalibParams.r; ParamsCheck.Tr=CalibParams.Tr; ParamsCheck.tau_c=CalibParams.tau_c;
        p_check=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d3,n_a,n_z,N_j,0,pi_z,d_grid3,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,ParamsCheck,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
        dfresh=max(abs([p_check.r;p_check.Tr;p_check.tau_c]-[CalibParams.r;CalibParams.Tr;CalibParams.tau_c]));
        failed=0;
    catch err
        fprintf('RoundTrips %s errored: %s \n',labels{rr},err.message)
        dpar=1; dpr=1; dg=1; dfresh=1; failed=1; time=NaN;
    end
    fprintf('RoundTrips %s: did not run, this should be zero: %.3e; (beta,varphi) vs truth, this should be close to zero: %.3e; GE prices vs equilibrium, this should be close to zero: %.3e \n',labels{rr},failed,dpar,dpr)
    fprintf('RoundTrips %s: calibsummary GE conditions, should be near zero: %.3e; fresh GE solve at the calibrated parameters vs CalibParams GE prices, should be near zero: %.3e (%.1f seconds) \n',labels{rr},dg,dfresh,time)
end

%% (D) omitcalibparam on part of the row kappa_j, with a row omit mask
TMD=TM;
TMD.AgeConditionalStats.earnings.Mean=gather(LCP.earnings.Mean);
caliboptions=caliboptionsBase;
caliboptions.fminalgo=8;
omitmask=Params.kappa_j; omitmask(2:4)=NaN; % a row, like kappa_j: calibrate ages 2 to 4
caliboptions.omitcalibparam.kappa_j=omitmask;
ParamsStartD=ParamsStart;
ParamsStartD.kappa_j(2:4)=1.1*Params.kappa_j(2:4);
try
    [CalibParamsD,calibsummaryD]=CalibrateOLGModel({'beta','varphi','kappa_j'},TMD,n_d3,n_a,n_z,N_j,d_grid3,a_grid,z_grid,pi_z,ReturnFn,ParamsStartD,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,struct(),caliboptions,vfoptions,simoptions);
    dk=max(abs(CalibParamsD.kappa_j(:)-Params.kappa_j(:)));
    dfixed=max(abs(CalibParamsD.kappa_j([1,5:N_j])-Params.kappa_j([1,5:N_j])'));
    dpar=max(abs([CalibParamsD.beta;CalibParamsD.varphi]-truth));
    dg=max(abs(calibsummaryD.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('RoundTrips (D) omitcalibparam errored: %s \n',err.message)
    dk=1; dfixed=1; dpar=1; dg=1; failed=1;
end
fprintf('RoundTrips (D) omitcalibparam (row mask): did not run, this should be zero: %.3e; omitted ages of kappa_j unchanged, this should be zero: %.3e \n',failed,dfixed)
fprintf('RoundTrips (D) omitcalibparam: kappa_j vs truth, this should be close to zero: %.3e; (beta,varphi) vs truth, this should be close to zero: %.3e; GE conditions, should be near zero: %.3e \n',dk,dpar,dg)

%% The one-price model
Params1=Params;
Params1.tau=0; Params1.G=0; Params1.Tr=0; Params1.tau_c=0;
GEPriceParamNames1={'r'};
GeneralEqmEqns1.CapitalMarket=GeneralEqmEqns.CapitalMarket;
p_eqm1=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns1,Params1,DiscountFactorParamNames,[],[],[],GEPriceParamNames1,heteroagentoptions,simoptions,vfoptions);
Params1.r=p_eqm1.r;
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params1,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
TM1=struct();
TM1.AllStats.K.Mean=gather(AllStats.K.Mean);
TM1.AllStats.H.Mean=gather(AllStats.H.Mean);
Params1Start=Params1;
Params1Start.beta=0.95; Params1Start.varphi=1.2; Params1Start.r=1.1*Params1.r;

labels={'(E) nested lsqnonlin (one price)','(F) nested fminsearch (one price)','(G) joint lsqnonlin (one price)','(B) joint fminsearch (one price)'};
CP=cell(1,4);
for rr=1:4
    caliboptions=caliboptionsBase;
    if rr==1
        caliboptions.fminalgo=8; caliboptions.jointoptimization=0;
    elseif rr==2
        caliboptions.fminalgo=1; caliboptions.jointoptimization=0;
    elseif rr==3
        caliboptions.fminalgo=8; caliboptions.jointoptimization=1;
    else
        caliboptions.fminalgo=1; caliboptions.jointoptimization=1;
    end
    try
        tt=tic;
        [CalibParams,calibsummary]=CalibrateOLGModel(CalibParamNames,TM1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1Start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames1,FnsToEvaluate,GeneralEqmEqns1,heteroagentoptions,caliboptions,vfoptions,simoptions);
        time=toc(tt);
        dpar=max(abs([CalibParams.beta;CalibParams.varphi]-truth));
        dpr=abs(CalibParams.r-Params1.r);
        dg=max(abs(calibsummary.GeneralEqmConditions));
        CP{rr}=[CalibParams.beta;CalibParams.varphi;CalibParams.r];
        failed=0;
    catch err
        fprintf('RoundTrips %s errored: %s \n',labels{rr},err.message)
        dpar=1; dpr=1; dg=1; failed=1; time=NaN; CP{rr}=nan(3,1);
    end
    fprintf('RoundTrips %s: did not run, this should be zero: %.3e; (beta,varphi) vs truth, this should be close to zero: %.3e; r vs equilibrium, this should be close to zero: %.3e; GE condition, should be near zero: %.3e (%.1f seconds) \n',labels{rr},failed,dpar,dpr,dg,time)
end
fprintf('RoundTrips (one price): nested lsqnonlin vs joint lsqnonlin, (beta,varphi,r), this should be close to zero: %.3e \n',max(abs(CP{1}-CP{3})))

output=struct();
output.CP=CP;

end
