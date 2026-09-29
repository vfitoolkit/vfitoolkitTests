function output=CoreStationaryGE_InfHorz_extraoptions(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz stationary GE: the general eqm OPTIONS that the rest of this bank never touches -
% heteroagentoptions.intermediateEqns, and shock grids that are themselves determined in general
% eqm (vfoptions.ExogShockFn taking a GE price, which sets heteroagentoptions.gridsinGE=1).
%
% HOW THESE ARE CHECKED, AND WHY IT IS CHEAP. Every one of these options is supposed to leave the
% equilibrium alone (an intermediate eqn that computes a quantity the GE eqn already computed
% inline; a grid rebuild that rebuilds the same grid), so the test is an invariance test. But it
% does NOT need an equilibrium to be found: with heteroagentoptions.maxiter=0 the toolkit just
% evaluates the general eqm conditions at the prices currently in Params and returns them, which
% is one model solve instead of the 150 a solve takes. Two runs of the same model at the same
% prices must then agree BIT-IDENTICALLY, so these are identity checks ("should be zero") rather
% than the looser "near zero" of the converged-equilibrium comparisons elsewhere in the bank.
% Cheaper AND sharper: a solve would hide a small corruption inside its own tolerance.
%
% The one thing maxiter=0 cannot see is staleness. At the initial prices, a grid frozen at the
% initial guess and a grid rebuilt at the current prices are the same grid, so a solver that
% never rebuilt would pass every check above. That needs prices to have MOVED, so the last check
% does one real solve with grids that genuinely depend on r, and then re-evaluates with the grids
% frozen at the answer: see the frozen-grid check below.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC, so it is NOT a GE price.

n_p=0;

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);

FnsToEvaluate.K=@(d,aprime,a,z) a;
FnsToEvaluate.N=@(d,aprime,a,z) d*z;
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c);

GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

% Parameters of the z process, needed because ExogShockFn is given its inputs by name out of
% Parameters. n_z_param is n_z under another name for exactly that reason (and not called n_z so
% that nothing else can mistake it for the grid size it is).
Params.rho_z=0.9;
Params.sigma_epsz=0.03;
Params.n_z_param=n_z;
Params.rgearing=0; % =0 leaves CoreStationaryGE_zgrid returning the baseline grid

% Every evaluation here uses fminalgo=1 and maxiter=0 (just evaluate the GE conditions at the
% prices in Params), except the one real solve at the end.
heteroagentoptions0=heteroagentoptions;
heteroagentoptions0.fminalgo=1;
heteroagentoptions0.maxiter=0;

%% Reference: the baseline model, GE conditions evaluated at the initial prices
[~,GEcondnsR]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% I1: the wage as an intermediate eqn
% GovBudget recomputes w from r inline; hand that expression to intermediateEqns instead and let
% the GE eqn take w as a parameter. Same arithmetic, so the same GE conditions.
GeneralEqmEqnsI1=GeneralEqmEqns;
GeneralEqmEqnsI1.GovBudget=@(tau,w,N,Tr) tau*w*N-Tr;
heteroagentoptionsI1=heteroagentoptions0;
heteroagentoptionsI1.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
[~,GEcondnsI1]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI1, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsI1, simoptions, vfoptions);

%% I2: two intermediate eqns, the second one consuming the first
% This is the check that intermediates are evaluated IN ORDER and each one is in Parameters
% before the next is evaluated: taxbase can only be computed if w is already there. Struct field
% order is insertion order, so w is declared first on purpose.
GeneralEqmEqnsI2=GeneralEqmEqns;
GeneralEqmEqnsI2.GovBudget=@(taxbase,Tr) taxbase-Tr;
heteroagentoptionsI2=heteroagentoptions0;
heteroagentoptionsI2.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
heteroagentoptionsI2.intermediateEqns.taxbase=@(tau,w,N) tau*w*N;
[~,GEcondnsI2]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI2, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsI2, simoptions, vfoptions);

%% I3: an intermediate eqn consuming an aggregate variable
% C is an AggVar, so this only works if the intermediates are evaluated AFTER the AggVars have
% been put into Parameters. Pins the ordering on that side too.
GeneralEqmEqnsI3=GeneralEqmEqns;
GeneralEqmEqnsI3.ConsTax=@(taxrev_c,G) taxrev_c-G;
heteroagentoptionsI3=heteroagentoptions0;
heteroagentoptionsI3.intermediateEqns.taxrev_c=@(tau_c,C) tau_c*C;
[~,GEcondnsI3]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI3, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsI3, simoptions, vfoptions);

%% Z1: shock grids determined in general eqm, but rebuilt to exactly the baseline grid
% r is an input of ExogShockFn, so the toolkit rebuilds the grids at every price; rgearing=0 means
% what it rebuilds is the baseline grid. Note the z_grid and pi_z inputs below are now placeholders.
vfoptionsZ1=vfoptions;
vfoptionsZ1.ExogShockFn=@(rho_z,sigma_epsz,n_z_param,r,rgearing) CoreStationaryGE_zgrid(rho_z,sigma_epsz,n_z_param,r,rgearing);
simoptionsZ1=simoptions;
simoptionsZ1.ExogShockFn=vfoptionsZ1.ExogShockFn; % set in both, as a user would
[~,GEcondnsZ1]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z2a: grids that genuinely move with r - confirm the variant actually bites
% Without this, a gearing that silently did nothing would let every check below pass vacuously.
ParamsZ=Params;
ParamsZ.rgearing=0.5; % income risk rises with r: sigma_epsz+0.5*r
[~,GEcondnsZ2a]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z2b: the staleness check - solve for real, then re-evaluate with the grids frozen at the answer
% The solve's own final evaluation of the GE conditions is made with grids rebuilt at p_eqmZ. If
% instead the solver had frozen the grids at the initial guess, those conditions would be the ones
% belonging to a different z process, and would not match this frozen-grid evaluation.
heteroagentoptionsZ=heteroagentoptions; % a real solve, so NOT maxiter=0
heteroagentoptionsZ.fminalgo=1;
[p_eqmZ,GEcondnsZ]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsZ, simoptionsZ1, vfoptionsZ1);

ParamsZfrozen=ParamsZ;
ParamsZfrozen.r=p_eqmZ.r; ParamsZfrozen.Tr=p_eqmZ.Tr; ParamsZfrozen.tau_c=p_eqmZ.tau_c;
[z_gridZ,pi_zZ]=CoreStationaryGE_zgrid(ParamsZ.rho_z,ParamsZ.sigma_epsz,n_z,p_eqmZ.r,ParamsZ.rgearing);
[~,GEcondnsZfrozen]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_zZ, d_grid, a_grid, z_gridZ, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZfrozen, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% Compare
% Each difference is the largest absolute gap across the three general eqm conditions.
dI1=max(abs([GEcondnsR.CapitalMarket-GEcondnsI1.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI1.GovBudget,GEcondnsR.ConsTax-GEcondnsI1.ConsTax]));
dI2=max(abs([GEcondnsR.CapitalMarket-GEcondnsI2.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI2.GovBudget,GEcondnsR.ConsTax-GEcondnsI2.ConsTax]));
dI3=max(abs([GEcondnsR.CapitalMarket-GEcondnsI3.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI3.GovBudget,GEcondnsR.ConsTax-GEcondnsI3.ConsTax]));
dZ1=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ1.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ1.GovBudget,GEcondnsR.ConsTax-GEcondnsZ1.ConsTax]));
dZ2a=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ2a.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ2a.GovBudget,GEcondnsR.ConsTax-GEcondnsZ2a.ConsTax]));
dZ2b=max(abs([GEcondnsZ.CapitalMarket-GEcondnsZfrozen.CapitalMarket,GEcondnsZ.GovBudget-GEcondnsZfrozen.GovBudget,GEcondnsZ.ConsTax-GEcondnsZfrozen.ConsTax]));

fprintf('\n=== InfHorz: intermediateEqns and shock grids in GE ===\n')
fprintf('GE conditions at the initial prices, baseline: CapitalMarket=%.6e GovBudget=%.6e ConsTax=%.6e \n',GEcondnsR.CapitalMarket,GEcondnsR.GovBudget,GEcondnsR.ConsTax)
fprintf('intermediateEqns, wage as an intermediate, this should be zero: %.3e \n',dI1)
fprintf('intermediateEqns, one intermediate consuming another, this should be zero: %.3e \n',dI2)
fprintf('intermediateEqns, an intermediate consuming an AggVar, this should be zero: %.3e \n',dI3)
fprintf('pi_z in GE, rebuilt to the same grid, this should be zero: %.3e \n',dZ1)
fprintf('pi_z in GE, grids geared to r, this should NOT be zero: %.3e \n',dZ2a)
fprintf('  (that solve landed at r=%.6f Tr=%.6f tau_c=%.6f) \n',p_eqmZ.r,p_eqmZ.Tr,p_eqmZ.tau_c)
fprintf('pi_z in GE, solver grids vs grids frozen at its own answer, this should be zero: %.3e \n',dZ2b)

output.GEcondnsR=GEcondnsR;
output.GEcondnsI1=GEcondnsI1;
output.GEcondnsI2=GEcondnsI2;
output.GEcondnsI3=GEcondnsI3;
output.GEcondnsZ1=GEcondnsZ1;
output.GEcondnsZ2a=GEcondnsZ2a;
output.GEcondnsZ=GEcondnsZ;
output.GEcondnsZfrozen=GEcondnsZfrozen;
output.p_eqmZ=p_eqmZ;

end
