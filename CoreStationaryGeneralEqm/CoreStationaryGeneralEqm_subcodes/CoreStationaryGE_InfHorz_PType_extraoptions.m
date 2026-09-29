function output=CoreStationaryGE_InfHorz_PType_extraoptions(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz stationary GE with permanent types (PType): the general eqm OPTIONS that the rest of this
% bank never touches - heteroagentoptions.intermediateEqns, and shock grids that are themselves
% determined in general eqm (vfoptions.ExogShockFn taking a GE price, which sets
% heteroagentoptions.gridsinGE=1, here a VECTOR with one entry per permanent type).
%
% HOW THESE ARE CHECKED, AND WHY IT IS CHEAP. See the header of
% CoreStationaryGE_InfHorz_extraoptions: each option is supposed to leave the model alone, and
% heteroagentoptions.maxiter=0 evaluates the general eqm conditions at the prices currently in
% Params without solving, so two runs of the same model must agree BIT-IDENTICALLY.
%
% The PType-specific check here is the MIXED case: the grids of one permanent type depend on a GE
% price and the other type's do not, so heteroagentoptions.gridsinGE is [0 1] and the per-type
% branching is what gets tested.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC, so it is NOT a GE price.

n_p=0;

% Permanent types: two types differing in sigma
N_i=2;
Names_i=N_i;
Params.sigma=[2.2,1.8];       % differs by permanent type
Params.ptypemass=[0.5,0.5];   % mass of each permanent type
PTypeDistParamNames={'ptypemass'};

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
[~,GEcondnsR]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% I1: the wage as an intermediate eqn
GeneralEqmEqnsI1=GeneralEqmEqns;
GeneralEqmEqnsI1.GovBudget=@(tau,w,N,Tr) tau*w*N-Tr;
heteroagentoptionsI1=heteroagentoptions0;
heteroagentoptionsI1.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
[~,GEcondnsI1]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI1, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsI1, simoptions, vfoptions);

%% I2: two intermediate eqns, the second one consuming the first
% Struct field order is insertion order, so w is declared first on purpose: taxbase can only be
% computed if w is already in Parameters.
GeneralEqmEqnsI2=GeneralEqmEqns;
GeneralEqmEqnsI2.GovBudget=@(taxbase,Tr) taxbase-Tr;
heteroagentoptionsI2=heteroagentoptions0;
heteroagentoptionsI2.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
heteroagentoptionsI2.intermediateEqns.taxbase=@(tau,w,N) tau*w*N;
[~,GEcondnsI2]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI2, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsI2, simoptions, vfoptions);

%% I3: an intermediate eqn consuming an aggregate variable
% C is an AggVar summed over the permanent types, so this also pins that the intermediates are
% evaluated after the ptype-weighted AggVars are in Parameters.
GeneralEqmEqnsI3=GeneralEqmEqns;
GeneralEqmEqnsI3.ConsTax=@(taxrev_c,G) taxrev_c-G;
heteroagentoptionsI3=heteroagentoptions0;
heteroagentoptionsI3.intermediateEqns.taxrev_c=@(tau_c,C) tau_c*C;
[~,GEcondnsI3]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI3, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsI3, simoptions, vfoptions);

%% Z1: shock grids determined in general eqm for BOTH types, rebuilt to exactly the baseline grid
vfoptionsZ1=vfoptions;
vfoptionsZ1.ExogShockFn=@(rho_z,sigma_epsz,n_z_param,r,rgearing) CoreStationaryGE_zgrid(rho_z,sigma_epsz,n_z_param,r,rgearing);
simoptionsZ1=simoptions;
simoptionsZ1.ExogShockFn=vfoptionsZ1.ExogShockFn; % set in both, as a user would
[~,GEcondnsZ1]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z4: the mixed case - only the second permanent type's grids are determined in general eqm
% heteroagentoptions.gridsinGE is then [0 1]: type 1 keeps the grids built once at setup, type 2 has
% them rebuilt every iteration. The per-type option is set the usual way, as vfoptions.ptype002
% (the names are generated as ptype001, ptype002 when Names_i is just the number of types).
vfoptionsZ4=vfoptions;
vfoptionsZ4.ptype002.ExogShockFn=vfoptionsZ1.ExogShockFn;
simoptionsZ4=simoptions;
simoptionsZ4.ptype002.ExogShockFn=vfoptionsZ1.ExogShockFn;
[~,GEcondnsZ4]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ4, vfoptionsZ4);

%% Z2a: grids that genuinely move with r - confirm the variant actually bites
ParamsZ=Params;
ParamsZ.rgearing=0.5; % income risk rises with r: sigma_epsz+0.5*r
[~,GEcondnsZ2a]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z2b: the staleness check - solve for real, then re-evaluate with the grids frozen at the answer
heteroagentoptionsZ=heteroagentoptions; % a real solve, so NOT maxiter=0
heteroagentoptionsZ.fminalgo=1;
[p_eqmZ,GEcondnsZ]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsZ, simoptionsZ1, vfoptionsZ1);

ParamsZfrozen=ParamsZ;
ParamsZfrozen.r=p_eqmZ.r; ParamsZfrozen.Tr=p_eqmZ.Tr; ParamsZfrozen.tau_c=p_eqmZ.tau_c;
[z_gridZ,pi_zZ]=CoreStationaryGE_zgrid(ParamsZ.rho_z,ParamsZ.sigma_epsz,n_z,p_eqmZ.r,ParamsZ.rgearing);
[~,GEcondnsZfrozen]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_zZ, d_grid, a_grid, z_gridZ, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZfrozen, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% Compare
dI1=max(abs([GEcondnsR.CapitalMarket-GEcondnsI1.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI1.GovBudget,GEcondnsR.ConsTax-GEcondnsI1.ConsTax]));
dI2=max(abs([GEcondnsR.CapitalMarket-GEcondnsI2.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI2.GovBudget,GEcondnsR.ConsTax-GEcondnsI2.ConsTax]));
dI3=max(abs([GEcondnsR.CapitalMarket-GEcondnsI3.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI3.GovBudget,GEcondnsR.ConsTax-GEcondnsI3.ConsTax]));
dZ1=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ1.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ1.GovBudget,GEcondnsR.ConsTax-GEcondnsZ1.ConsTax]));
dZ4=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ4.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ4.GovBudget,GEcondnsR.ConsTax-GEcondnsZ4.ConsTax]));
dZ2a=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ2a.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ2a.GovBudget,GEcondnsR.ConsTax-GEcondnsZ2a.ConsTax]));
dZ2b=max(abs([GEcondnsZ.CapitalMarket-GEcondnsZfrozen.CapitalMarket,GEcondnsZ.GovBudget-GEcondnsZfrozen.GovBudget,GEcondnsZ.ConsTax-GEcondnsZfrozen.ConsTax]));

fprintf('\n=== InfHorz PType: intermediateEqns and shock grids in GE ===\n')
fprintf('GE conditions at the initial prices, baseline: CapitalMarket=%.6e GovBudget=%.6e ConsTax=%.6e \n',GEcondnsR.CapitalMarket,GEcondnsR.GovBudget,GEcondnsR.ConsTax)
fprintf('intermediateEqns, wage as an intermediate, this should be zero: %.3e \n',dI1)
fprintf('intermediateEqns, one intermediate consuming another, this should be zero: %.3e \n',dI2)
fprintf('intermediateEqns, an intermediate consuming an AggVar, this should be zero: %.3e \n',dI3)
fprintf('pi_z in GE for both ptypes, rebuilt to the same grid, this should be zero: %.3e \n',dZ1)
fprintf('pi_z in GE for the second ptype only (gridsinGE=[0 1]), this should be zero: %.3e \n',dZ4)
fprintf('pi_z in GE, grids geared to r, this should NOT be zero: %.3e \n',dZ2a)
fprintf('  (that solve landed at r=%.6f Tr=%.6f tau_c=%.6f) \n',p_eqmZ.r,p_eqmZ.Tr,p_eqmZ.tau_c)
fprintf('pi_z in GE, solver grids vs grids frozen at its own answer, this should be zero: %.3e \n',dZ2b)

output.GEcondnsR=GEcondnsR;
output.GEcondnsI1=GEcondnsI1;
output.GEcondnsI2=GEcondnsI2;
output.GEcondnsI3=GEcondnsI3;
output.GEcondnsZ1=GEcondnsZ1;
output.GEcondnsZ4=GEcondnsZ4;
output.GEcondnsZ2a=GEcondnsZ2a;
output.GEcondnsZ=GEcondnsZ;
output.GEcondnsZfrozen=GEcondnsZfrozen;
output.p_eqmZ=p_eqmZ;

end
