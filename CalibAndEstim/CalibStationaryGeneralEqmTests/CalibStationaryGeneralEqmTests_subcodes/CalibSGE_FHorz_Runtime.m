function output=CalibSGE_FHorz_Runtime(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,jequaloneDist,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 6: runtime of CalibrateOLGModel, informational. Each run stops after fminsearch's initial simplex
% (two calibrated parameters, so three evaluations in nested mode, six in joint, plus the one for calibsummary). Timed: joint with
% caliboptions.whichcombos=1 and =0 for a sparse and a dense target set, and nested. The only check is that whichcombos=1 and =0 give
% the same moments (exact).

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
p_eqm=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
for pp=1:length(GEPriceParamNames)
    Params.(GEPriceParamNames{pp})=p_eqm.(GEPriceParamNames{pp});
end

% Targets (values do not matter for timing): sparse = one statistic at five ages; dense = many statistics of every function
TMs=struct();
t=nan(1,N_j); t(5:9)=1; TMs.AgeConditionalStats.K.Mean=t;
TMd=struct();
TMd.AllStats.K.Mean=1; TMd.AllStats.K.Gini=1; TMd.AllStats.H.Median=1; TMd.AllStats.earnings.Variance=1; TMd.AllStats.C.QuantileCutoffs=ones(21,1); % (nquantiles+1: ventiles, min and max included)
TMd.AgeConditionalStats.K.Mean=ones(1,N_j); TMd.AgeConditionalStats.K.Median=ones(1,N_j); TMd.AgeConditionalStats.earnings.Gini=ones(1,N_j);
TMd.AgeConditionalStats.H.StdDeviation=ones(1,N_j); TMd.AgeConditionalStats.C.Mean=ones(1,N_j);

caliboptionsC=struct();
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
CalibParamNames={'beta','varphi'};
labs={'sparse','dense'};
for tt_c=1:2
    if tt_c==1
        TM=TMs;
    else
        TM=TMd;
    end
    caliboptions=caliboptionsC; caliboptions.whichcombos=1;
    tt=tic;
    [~,summ1]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptions,vfoptions,simoptions);
    time1=toc(tt);
    caliboptions=caliboptionsC; caliboptions.whichcombos=0;
    tt=tic;
    [~,summ0]=CalibrateOLGModel(CalibParamNames,TM,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptions,vfoptions,simoptions);
    time0=toc(tt);
    fprintf('Runtime %s targets: whichcombos=1 vs =0 moments, this should be zero: %.3e \n',labs{tt_c},max(abs(summ1.currentmomentvec-summ0.currentmomentvec)))
    fprintf('Runtime %s targets (informational): joint, 7 evaluations, whichcombos=1 %.2f s, whichcombos=0 %.2f s, ratio %.2f \n',labs{tt_c},time1,time0,time1/time0)
end
caliboptions=caliboptionsC; caliboptions.jointoptimization=0;
tt=tic;
CalibrateOLGModel(CalibParamNames,TMs,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],GEPriceParamNames,FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptions,vfoptions,simoptions);
timeN=toc(tt);
fprintf('Runtime (informational): nested, 4 evaluations (each a general eqm solve from the equilibrium), %.2f s \n',timeN)

output=struct();

end
