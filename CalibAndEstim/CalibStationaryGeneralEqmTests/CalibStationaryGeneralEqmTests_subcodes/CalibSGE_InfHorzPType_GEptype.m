function output=CalibSGE_InfHorzPType_GEptype(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,vfoptions,simoptions)
% CalibStationaryGeneralEqmTests part 23: general eqm conditions by permanent type (heteroagentoptions.GEptype) in CalibrateBIHAModel_PType.
% Part 22 for the InfHorz model (see there): Tr per type, GovBudget named in GEptype, so the GE conditions are [CapitalMarket;
% GovBudget_typeA; GovBudget_typeB; ConsTax].
% (a) identical types: joint with GEptype gives the economy-wide GE conditions with GovBudget twice, and the same moments (exact); Tr comes
%     back as a struct over the types
% (b) types that differ: the joint GE conditions equal HeteroAgentStationaryEqm_InfHorz_PType with GEptype and maxiter=0 at the same prices
%     (exact); the two per-type GovBudget differ (should NOT be zero)
% (c) the per-type intermediateEqn form gives the same GE conditions as (b) (exact)
% (d) nested at the GEptype general eqm: Tr returned as a struct, equal to the equilibrium (near); GE conditions near zero (loose, as in
%     part 9: the inner fminsearch over the InfHorz GE prices stalls on the discrete hours grid)
% (e) joint lsqnonlin round trip, beta from 0.95 and the four GE prices off the equilibrium, finitedifferencestepsize 1e-3 (see part 10):
%     beta and the prices recovered (close), GE conditions near zero
% (e2) Tr given as a column over the types: it comes back as that column (exact)
% (f) multiGEweights [1,2,0.5] expanded to [1,2,2,0.5]: through lsqnonlin, the reported objective equals the scalar formula (exact)

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);
FnsToEvaluate.K=@(d,aprime,a,z) a; % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption
FnsToEvaluate.H=@(d,aprime,a,z) d; % hours
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

heteroagentoptions=struct();
heteroagentoptions.fminalgo=1;
heteroagentoptions.toleranceGEprices=1e-7;
heteroagentoptions.toleranceGEcondns=1e-7;
heteroagentoptions.verbose=0;
heteroagentoptionsP=heteroagentoptions;
heteroagentoptionsP.GEptype={'GovBudget'};

caliboptionsC=struct(); % runs that stop after fminsearch's initial simplex, so return their start
caliboptionsC.verbose=0;
caliboptionsC.fminalgo=1;
caliboptionsC.toleranceparams=1e3;
caliboptionsC.toleranceobjective=1e3;
caliboptionsCJ=caliboptionsC; caliboptionsCJ.jointoptimization=1;
caliboptionsCN=caliboptionsC; caliboptionsCN.jointoptimization=0;

%% (a) identical types
% (Every run here stops after fminsearch's initial simplex, and fminsearch returns the BEST vertex of that simplex, which is the start
% only if the start minimises the objective: so the targets are the model at the start, and the start is a general eqm, of the economy-wide
% conditions in (a) and of the GEptype ones in (b)-(e2) [2026-11-01 the targets were placeholders and (b)-(c) were not at a GEptype
% equilibrium, so the runs moved to a simplex vertex and compared different points]. Each run checks that it returned its start.)
ParamsID=Params;
ParamsID.sigma=2;
p_eqm=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,ParamsID,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
ParamsID.r=p_eqm.r; ParamsID.Tr=p_eqm.Tr; ParamsID.tau_c=p_eqm.tau_c;
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsID,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,ParamsID,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,ParamsID,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
TMa=struct();
TMa.AllStats.K.Mean=gather(AllStats.K.Mean);
TMa.AllStats.H.Mean=gather(AllStats.H.Mean);
[~,summ0]=CalibrateBIHAModel_PType({'beta'},TMa,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsID,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptions,caliboptionsCJ,vfoptions,simoptions);
ParamsIDp=ParamsID;
ParamsIDp=rmfield(ParamsIDp,'Tr');
ParamsIDp.Tr.(Names_i{1})=ParamsID.Tr;
ParamsIDp.Tr.(Names_i{2})=ParamsID.Tr;
try
    [CPa,summa]=CalibrateBIHAModel_PType({'beta'},TMa,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsIDp,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsP,caliboptionsCJ,vfoptions,simoptions);
    g0=summ0.GeneralEqmConditions;
    if length(summa.GeneralEqmConditions)==4
        dg=max(abs(summa.GeneralEqmConditions-[g0(1);g0(2);g0(2);g0(3)]));
    else
        dg=1;
    end
    dm=max(abs(summa.currentmomentvec-summ0.currentmomentvec));
    if isstruct(CPa.Tr) && isfield(CPa.Tr,Names_i{1}) && isfield(CPa.Tr,Names_i{2})
        dtr=max(abs([CPa.beta-ParamsID.beta; CPa.Tr.(Names_i{1})-ParamsID.Tr; CPa.Tr.(Names_i{2})-ParamsID.Tr])); formbad=0;
    else
        dtr=1; formbad=1;
    end
    failed=0;
catch err
    fprintf('GEptype InfHorz (a) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dg=1; dm=1; dtr=1; formbad=1; failed=1;
end
fprintf('GEptype InfHorz (a) identical types, joint: did not run, this should be zero: %.3e; GE conditions vs the economy-wide ones (GovBudget twice), this should be zero: %.3e; moments vs economy-wide, this should be zero: %.3e \n',failed,dg,dm)
fprintf('GEptype InfHorz (a) identical types, joint: Tr not returned as a struct over the types, this should be zero: %.3e; returned its start (beta, Tr), this should be zero: %.3e \n',formbad,dtr)

%% The GEptype general eqm of the model with types that differ, and the targets there (used by (b)-(e2))
p_eqm=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions);
ParamsP=Params; % the economy-wide general eqm, with Tr handed in per type: the start of the GEptype solve
ParamsP.r=p_eqm.r; ParamsP.tau_c=p_eqm.tau_c;
ParamsP=rmfield(ParamsP,'Tr');
ParamsP.Tr.(Names_i{1})=p_eqm.Tr;
ParamsP.Tr.(Names_i{2})=p_eqm.Tr;
p_eqmP=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,ParamsP,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptionsP,simoptions,vfoptions);
ParamsE=ParamsP;
ParamsE.r=p_eqmP.r; ParamsE.tau_c=p_eqmP.tau_c; ParamsE.Tr=p_eqmP.Tr; % (p_eqmP.Tr is a struct over the types)
peqP=[p_eqmP.r; p_eqmP.Tr.(Names_i{1}); p_eqmP.Tr.(Names_i{2}); p_eqmP.tau_c];
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsE,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,ParamsE,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,ParamsE,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
TME=struct();
TME.AllStats.K.Mean=gather(AllStats.K.Mean);
TME.AllStats.H.Mean=gather(AllStats.H.Mean);
fprintf('GEptype InfHorz: at the GEptype general eqm the two types get different transfers (the types differ), this should NOT be zero: %.3e \n',abs(p_eqmP.Tr.(Names_i{1})-p_eqmP.Tr.(Names_i{2})))

%% (b) the joint GE conditions vs HeteroAgentStationaryEqm_InfHorz_PType with GEptype, maxiter=0, at the same prices
heteroagentoptions0=heteroagentoptionsP;
heteroagentoptions0.maxiter=0;
[~,GEc0]=HeteroAgentStationaryEqm_InfHorz_PType(n_d,n_a,n_z,Names_i,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,ParamsE,DiscountFactorParamNames,PTypeDistParamNames,GEPriceParamNames,heteroagentoptions0,simoptions,vfoptions);
g0=gather([GEc0.CapitalMarket(:); GEc0.GovBudget(:); GEc0.ConsTax(:)]);
try
    [CPb,summb]=CalibrateBIHAModel_PType({'beta'},TME,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsE,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsP,caliboptionsCJ,vfoptions,simoptions);
    if length(summb.GeneralEqmConditions)==4
        dg=max(abs(summb.GeneralEqmConditions-g0));
    else
        dg=1;
    end
    dstart=max(abs([CPb.beta; CPb.r; CPb.Tr.(Names_i{1}); CPb.Tr.(Names_i{2}); CPb.tau_c]-[ParamsE.beta; peqP]));
    failed=0;
catch err
    fprintf('GEptype InfHorz (b) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dg=1; dstart=1; failed=1;
end
fprintf('GEptype InfHorz (b) types that differ, joint: did not run, this should be zero: %.3e; GE conditions vs HeteroAgentStationaryEqm_InfHorz_PType with GEptype at the same prices, this should be zero: %.3e; returned its start, this should be zero: %.3e \n',failed,dg,dstart)

%% (c) a per-type intermediateEqn: the tax revenue by type, GovBudget comparing it to Tr
heteroagentoptionsI=heteroagentoptionsP;
heteroagentoptionsI.intermediateEqns.taxrev=@(tau,r,N,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N;
heteroagentoptionsI.intermediateEqnsptype={'taxrev'};
GeneralEqmEqnsI=GeneralEqmEqns;
GeneralEqmEqnsI.GovBudget=@(taxrev,Tr) taxrev-Tr;
try
    [CPc,summc]=CalibrateBIHAModel_PType({'beta'},TME,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsE,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqnsI,heteroagentoptionsI,caliboptionsCJ,vfoptions,simoptions);
    if length(summc.GeneralEqmConditions)==4
        dg=max(abs(summc.GeneralEqmConditions-g0));
    else
        dg=1;
    end
    dstart=max(abs([CPc.beta; CPc.r; CPc.Tr.(Names_i{1}); CPc.Tr.(Names_i{2}); CPc.tau_c]-[ParamsE.beta; peqP]));
    failed=0;
catch err
    fprintf('GEptype InfHorz (c) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dg=1; dstart=1; failed=1;
end
fprintf('GEptype InfHorz (c) per-type intermediateEqn, joint: did not run, this should be zero: %.3e; GE conditions vs (b), this should be zero: %.3e; returned its start, this should be zero: %.3e \n',failed,dg,dstart)

%% (d) nested at the GEptype general eqm
try
    [CPd,summd]=CalibrateBIHAModel_PType({'beta'},TME,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsE,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsP,caliboptionsCN,vfoptions,simoptions);
    if isstruct(CPd.Tr) && isfield(CPd.Tr,Names_i{1}) && isfield(CPd.Tr,Names_i{2})
        formbad=0; dp=max(abs([CPd.r; CPd.Tr.(Names_i{1}); CPd.Tr.(Names_i{2}); CPd.tau_c]-peqP));
    else
        formbad=1; dp=1;
    end
    if length(summd.GeneralEqmConditions)==4
        dg=max(abs(summd.GeneralEqmConditions));
    else
        dg=1;
    end
    dstart=abs(CPd.beta-ParamsE.beta);
    failed=0;
catch err
    fprintf('GEptype InfHorz (d) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    formbad=1; dp=1; dg=1; dstart=1; failed=1;
end
fprintf('GEptype InfHorz (d) nested at the GEptype general eqm: did not run, this should be zero: %.3e; Tr not returned as a struct over the types, this should be zero: %.3e; beta returned its start, this should be zero: %.3e; GE prices vs the equilibrium, should be near zero (loose: InfHorz fminsearch over several GE prices stalls on the discrete hours grid): %.3e; the four GE conditions, should be near zero (loose: same reason): %.3e \n',failed,formbad,dstart,dp,dg)

%% (e) joint lsqnonlin round trip
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsE,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,ParamsE,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,ParamsE,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
TMR=struct();
TMR.AllStats.K.Mean=gather(AllStats.K.Mean);
ParamsRStart=ParamsE;
ParamsRStart.beta=0.95; ParamsRStart.r=1.1*ParamsE.r; ParamsRStart.tau_c=1.1*ParamsE.tau_c;
ParamsRStart.Tr.(Names_i{1})=0.9*ParamsE.Tr.(Names_i{1}); ParamsRStart.Tr.(Names_i{2})=0.9*ParamsE.Tr.(Names_i{2});
caliboptionsR=struct();
caliboptionsR.verbose=0;
caliboptionsR.fminalgo=8;
caliboptionsR.toleranceparams=1e-6;
caliboptionsR.toleranceobjective=1e-10;
caliboptionsR.finitedifferencestepsize=1e-3; % (joint: lsqnonlin probes r, see part 10)
try
    tt=tic;
    [CPe,summe]=CalibrateBIHAModel_PType({'beta'},TMR,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsRStart,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsP,caliboptionsR,vfoptions,simoptions);
    time=toc(tt);
    dbeta=abs(CPe.beta-ParamsE.beta);
    dp=max(abs([CPe.r; CPe.Tr.(Names_i{1}); CPe.Tr.(Names_i{2}); CPe.tau_c]-peqP));
    dg=max(abs(summe.GeneralEqmConditions));
    failed=0;
catch err
    fprintf('GEptype InfHorz (e) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    dbeta=1; dp=1; dg=1; failed=1; time=NaN;
end
fprintf('GEptype InfHorz (e) joint lsqnonlin round trip, Tr per type: did not run, this should be zero: %.3e; beta vs truth, this should be close to zero: %.3e; GE prices (Tr per type) vs the equilibrium, this should be close to zero: %.3e; GE conditions, should be near zero: %.3e (%.1f seconds) \n',failed,dbeta,dp,dg,time)

%% (e2) Tr given as a column over the types
ParamsV=ParamsE;
ParamsV.Tr=[ParamsE.Tr.(Names_i{1}); ParamsE.Tr.(Names_i{2})];
try
    [CPv,~]=CalibrateBIHAModel_PType({'beta'},TME,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsV,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsP,caliboptionsCJ,vfoptions,simoptions);
    if isnumeric(CPv.Tr) && isequal(size(CPv.Tr),[2,1])
        formbad=0; dtr=max(abs(CPv.Tr-ParamsV.Tr));
    else
        formbad=1; dtr=1;
    end
    failed=0;
catch err
    fprintf('GEptype InfHorz (e2) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    formbad=1; dtr=1; failed=1;
end
fprintf('GEptype InfHorz (e2) Tr given as a column over the types, joint: did not run, this should be zero: %.3e; Tr not returned as that column, this should be zero: %.3e; returned its start, this should be zero: %.3e \n',failed,formbad,dtr)

%% (f) multiGEweights expanded to the GE conditions, through lsqnonlin
heteroagentoptionsW=heteroagentoptionsP;
mw=[1,2,0.5];
heteroagentoptionsW.multiGEweights=mw;
TMW=struct();
TMW.AllStats.K.Mean=1.05*TMR.AllStats.K.Mean;
caliboptionsW=struct();
caliboptionsW.verbose=0;
caliboptionsW.fminalgo=8;
caliboptionsW.toleranceparams=1e-2; % (the identity holds wherever lsqnonlin stops)
caliboptionsW.toleranceobjective=1e-4;
caliboptionsW.finitedifferencestepsize=1e-3;
try
    [~,summf]=CalibrateBIHAModel_PType({'beta'},TMW,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsE,DiscountFactorParamNames,GEPriceParamNames,PTypeDistParamNames,[],FnsToEvaluate,GeneralEqmEqns,heteroagentoptionsW,caliboptionsW,vfoptions,simoptions);
    dR=summf.currentmomentvec-summf.targetmomentvec;
    if length(summf.GeneralEqmConditions)==4
        formula=sum(dR.^2)/1+10*sum([1;2;2;0.5].*summf.GeneralEqmConditions.^2);
        d=abs(summf.objvalue-formula)/abs(formula);
    else
        formula=NaN; d=1;
    end
    failed=0;
catch err
    fprintf('GEptype InfHorz (f) errored: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    d=1; formula=NaN; failed=1;
end
fprintf('GEptype InfHorz (f) multiGEweights [1,2,0.5] expanded to [1,2,2,0.5], lsqnonlin: did not run, this should be zero: %.3e; reported objective vs the scalar formula at the solution, relative, this should be zero: %.3e (formula %.6e) \n',failed,d,formula)

output=struct();
output.p_eqmP=p_eqmP;

end
