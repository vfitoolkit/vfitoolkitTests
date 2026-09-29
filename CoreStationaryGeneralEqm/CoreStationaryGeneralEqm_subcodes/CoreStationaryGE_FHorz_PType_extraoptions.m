function output=CoreStationaryGE_FHorz_PType_extraoptions(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% FHorz stationary GE with permanent types (PType): the general eqm OPTIONS that the rest of this
% bank never touches - heteroagentoptions.intermediateEqns, the four forms jequaloneDist can take
% when there are permanent types, and shock grids that are themselves determined in general eqm
% (vfoptions.ExogShockFn taking a GE price, so heteroagentoptions.gridsinGE is a vector over types).
%
% HOW THESE ARE CHECKED, AND WHY IT IS CHEAP. See the header of
% CoreStationaryGE_InfHorz_extraoptions: each option is supposed to leave the model alone, and
% heteroagentoptions.maxiter=0 evaluates the general eqm conditions at the prices currently in
% Params without solving, so two runs of the same model must agree BIT-IDENTICALLY. Two checks do
% need a real solve, because maxiter=0 cannot see staleness: at the initial prices a grid frozen at
% the initial guess IS the grid rebuilt at the current prices.
%
% THE FOUR jequaloneDist FORMS. A matrix (same newborns for every type), a function, a struct with
% one function per type, and a numeric array with ptype as its last dimension. The last is the
% interesting one: its slice masses ARE the permanent type masses, overriding PTypeDistParamNames,
% and J4 below is what holds the toolkit to that - it is the check that would have caught the
% masses being computed and then discarded.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC, so it is NOT a GE price.

n_p=0;
vfoptions.divideandconquer=1; % finite-horizon tests use divide-and-conquer

% Permanent types: two types differing in sigma
N_i=2;
Names_i=N_i;
Params.sigma=[2.2,1.8];       % differs by permanent type
Params.ptypemass=[0.5,0.5];   % mass of each permanent type
PTypeDistParamNames={'ptypemass'};

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);

FnsToEvaluate.K=@(d,aprime,a,z) a;
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c);

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
% prices in Params), except the two real solves at the end.
heteroagentoptions0=heteroagentoptions;
heteroagentoptions0.fminalgo=1;
heteroagentoptions0.maxiter=0;

%% Reference: the baseline model, GE conditions evaluated at the initial prices
[~,GEcondnsR]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% I1: the wage as an intermediate eqn
GeneralEqmEqnsI1=GeneralEqmEqns;
GeneralEqmEqnsI1.GovBudget=@(tau,w,N,Tr) tau*w*N-Tr;
heteroagentoptionsI1=heteroagentoptions0;
heteroagentoptionsI1.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
[~,GEcondnsI1]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI1, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsI1, simoptions, vfoptions);

%% I2: two intermediate eqns, the second one consuming the first
% Struct field order is insertion order, so w is declared first on purpose: taxbase can only be
% computed if w is already in Parameters.
GeneralEqmEqnsI2=GeneralEqmEqns;
GeneralEqmEqnsI2.GovBudget=@(taxbase,Tr) taxbase-Tr;
heteroagentoptionsI2=heteroagentoptions0;
heteroagentoptionsI2.intermediateEqns.w=@(r,alpha,delta,A) (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
heteroagentoptionsI2.intermediateEqns.taxbase=@(tau,w,N) tau*w*N;
[~,GEcondnsI2]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI2, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsI2, simoptions, vfoptions);

%% I3: an intermediate eqn consuming an aggregate variable
% C is an AggVar summed over the permanent types, so this also pins that the intermediates are
% evaluated after the ptype-weighted AggVars are in Parameters.
GeneralEqmEqnsI3=GeneralEqmEqns;
GeneralEqmEqnsI3.ConsTax=@(taxrev_c,G) taxrev_c-G;
heteroagentoptionsI3=heteroagentoptions0;
heteroagentoptionsI3.intermediateEqns.taxrev_c=@(tau_c,C) tau_c*C;
[~,GEcondnsI3]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI3, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsI3, simoptions, vfoptions);

%% J1: jequaloneDist as a function, against the same newborn distribution as a matrix
% A jequaloneDist function is handed (a_grid,z_grid,n_a,n_z), so it must be given a_grid and z_grid
% in simoptions. CoreStationaryGE_j1dist puts a different newborn distribution on the grid from the
% bank's baseline (weights proportional to z, because a function is never given pi_z), so the
% reference here is that SAME distribution built as a matrix: the pair differs only in the route.
simoptionsJ=simoptions;
simoptionsJ.a_grid=a_grid;
simoptionsJ.z_grid=z_grid;
jequaloneDistNum=zeros([n_a,n_z],'gpuArray');
jequaloneDistNum(1,:)=shiftdim(z_grid./sum(z_grid),-1); % must match CoreStationaryGE_j1dist exactly
jequaloneDistFn=@(a_grid,z_grid,n_a,n_z) CoreStationaryGE_j1dist(a_grid,z_grid,n_a,n_z);
[~,GEcondnsJ1num]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistNum, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsJ, vfoptions);
[~,GEcondnsJ1fn]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistFn, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsJ, vfoptions);

%% J2: a struct with one function per permanent type
% Until 2026-09-20 this form did not reach a solve at all: the per-type value was left as a handle
% and StationaryDist_FHorz_Case1 then had no user z_grid to evaluate it against. Both types get the
% same function here, so it must agree with J1's single-handle run.
jequaloneDistStructFn.ptype001=jequaloneDistFn;
jequaloneDistStructFn.ptype002=jequaloneDistFn;
[~,GEcondnsJ2]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistStructFn, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsJ, vfoptions);

%% J3: ptype as a DIMENSION of jequaloneDist, with slice masses that agree with PTypeDistParamNames
% The [n_a,n_z,N_i] form is converted to the struct form internally, its slices normalized to mass
% one and their masses taken as the ptype masses. Here the slices carry 0.5 each, which is what
% ptypemass says anyway, so this must agree with the same thing passed as a struct of matrices.
jequaloneDistStructNum.ptype001=jequaloneDistNum;
jequaloneDistStructNum.ptype002=jequaloneDistNum;
jequaloneDistDim=zeros([n_a,n_z,N_i],'gpuArray');
jequaloneDistDim(:,:,1)=0.5*jequaloneDistNum;
jequaloneDistDim(:,:,2)=0.5*jequaloneDistNum;
[~,GEcondnsJ3struct]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistStructNum, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);
[~,GEcondnsJ3dim]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistDim, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% J4: ptype as a dimension, with slice masses that CONTRADICT PTypeDistParamNames
% The slices carry 0.3 and 0.7 while ptypemass still says 0.5/0.5, so the implicit masses must win
% (the toolkit warns that it is doing exactly that - the warning in the diary here is correct, not a
% failure). The reference is the same model set up the explicit way: struct of matrices, and
% ptypemass=[0.3,0.7]. If the implicit masses were computed and then dropped, these two disagree.
ParamsJ4=Params;
ParamsJ4.ptypemass=[0.3,0.7];
jequaloneDistDim4=zeros([n_a,n_z,N_i],'gpuArray');
jequaloneDistDim4(:,:,1)=0.3*jequaloneDistNum;
jequaloneDistDim4(:,:,2)=0.7*jequaloneDistNum;
[~,GEcondnsJ4struct]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistStructNum, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsJ4, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);
[~,GEcondnsJ4dim]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistDim4, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% Z1: shock grids determined in general eqm for BOTH types, rebuilt to exactly the baseline grid
vfoptionsZ1=vfoptions;
vfoptionsZ1.ExogShockFn=@(rho_z,sigma_epsz,n_z_param,r,rgearing) CoreStationaryGE_zgrid(rho_z,sigma_epsz,n_z_param,r,rgearing);
simoptionsZ1=simoptions;
simoptionsZ1.ExogShockFn=vfoptionsZ1.ExogShockFn; % set in both, as a user would
[~,GEcondnsZ1]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z3: the same, but with the grids age-dependent
% sigma_epsz as a 1-by-N_j vector makes ExogShockFn be evaluated once per age, which is a different
% branch of ExogShockSetup_FHorz (and the branch that allocates the age-deep copy of the user's own
% grid, per permanent type). Every age gets the same value, so the model is still the baseline one.
ParamsZ3=Params;
ParamsZ3.sigma_epsz=0.03*ones(1,N_j);
[~,GEcondnsZ3]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ3, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z4: the mixed case - only the second permanent type's grids are determined in general eqm
% heteroagentoptions.gridsinGE is then [0 1]: type 1 keeps the grids built once at setup, type 2 has
% them rebuilt every iteration. The per-type option is set the usual way, as vfoptions.ptype002
% (the names are generated as ptype001, ptype002 when Names_i is just the number of types).
vfoptionsZ4=vfoptions;
vfoptionsZ4.ptype002.ExogShockFn=vfoptionsZ1.ExogShockFn;
simoptionsZ4=simoptions;
simoptionsZ4.ptype002.ExogShockFn=vfoptionsZ1.ExogShockFn;
[~,GEcondnsZ4]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ4, vfoptionsZ4);

%% Z2a: grids that genuinely move with r - confirm the variant actually bites
ParamsZ=Params;
ParamsZ.rgearing=0.5; % income risk rises with r: sigma_epsz+0.5*r
[~,GEcondnsZ2a]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsZ1, vfoptionsZ1);

%% Z2b: the staleness check - solve for real, then re-evaluate with the grids frozen at the answer
heteroagentoptionsZ=heteroagentoptions; % a real solve, so NOT maxiter=0
heteroagentoptionsZ.fminalgo=1;
[p_eqmZ,GEcondnsZ]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsZ, simoptionsZ1, vfoptionsZ1);

ParamsZfrozen=ParamsZ;
ParamsZfrozen.r=p_eqmZ.r; ParamsZfrozen.Tr=p_eqmZ.Tr; ParamsZfrozen.tau_c=p_eqmZ.tau_c;
[z_gridZ,pi_zZ]=CoreStationaryGE_zgrid(ParamsZ.rho_z,ParamsZ.sigma_epsz,n_z,p_eqmZ.r,ParamsZ.rgearing);
[~,GEcondnsZfrozen]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_zZ, d_grid, a_grid, z_gridZ, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZfrozen, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptions, vfoptions);

%% C1: jequaloneDist as a function AND the grids moving with r, together
% This is the only check in the bank that reaches the line that hands each permanent type its own
% REBUILT user z_grid to evaluate jequaloneDist against: the newborn weights are read off z_grid,
% and z_grid is being rebuilt at every price. Solve for real, then re-evaluate with the grids frozen
% at the answer - if jequaloneDist had been given a grid built at the initial guess (or the internal
% joint-grid form), the newborn distribution during the solve was the wrong one and these disagree.
simoptionsC=simoptionsZ1;
simoptionsC.a_grid=a_grid;
simoptionsC.z_grid=z_grid;
[p_eqmC,GEcondnsC]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDistStructFn, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsZ, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsZ, simoptionsC, vfoptionsZ1);

ParamsCfrozen=ParamsZ;
ParamsCfrozen.r=p_eqmC.r; ParamsCfrozen.Tr=p_eqmC.Tr; ParamsCfrozen.tau_c=p_eqmC.tau_c;
[z_gridC,pi_zC]=CoreStationaryGE_zgrid(ParamsZ.rho_z,ParamsZ.sigma_epsz,n_z,p_eqmC.r,ParamsZ.rgearing);
simoptionsCfrozen=simoptions;
simoptionsCfrozen.a_grid=a_grid;
simoptionsCfrozen.z_grid=z_gridC; % the frozen grid, so jequaloneDist reads its weights off that
[~,GEcondnsCfrozen]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_zC, d_grid, a_grid, z_gridC, jequaloneDistStructFn, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsCfrozen, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions0, simoptionsCfrozen, vfoptions);

%% Compare
dI1=max(abs([GEcondnsR.CapitalMarket-GEcondnsI1.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI1.GovBudget,GEcondnsR.ConsTax-GEcondnsI1.ConsTax]));
dI2=max(abs([GEcondnsR.CapitalMarket-GEcondnsI2.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI2.GovBudget,GEcondnsR.ConsTax-GEcondnsI2.ConsTax]));
dI3=max(abs([GEcondnsR.CapitalMarket-GEcondnsI3.CapitalMarket,GEcondnsR.GovBudget-GEcondnsI3.GovBudget,GEcondnsR.ConsTax-GEcondnsI3.ConsTax]));
dJ1=max(abs([GEcondnsJ1num.CapitalMarket-GEcondnsJ1fn.CapitalMarket,GEcondnsJ1num.GovBudget-GEcondnsJ1fn.GovBudget,GEcondnsJ1num.ConsTax-GEcondnsJ1fn.ConsTax]));
dJ2=max(abs([GEcondnsJ1fn.CapitalMarket-GEcondnsJ2.CapitalMarket,GEcondnsJ1fn.GovBudget-GEcondnsJ2.GovBudget,GEcondnsJ1fn.ConsTax-GEcondnsJ2.ConsTax]));
dJ3=max(abs([GEcondnsJ3struct.CapitalMarket-GEcondnsJ3dim.CapitalMarket,GEcondnsJ3struct.GovBudget-GEcondnsJ3dim.GovBudget,GEcondnsJ3struct.ConsTax-GEcondnsJ3dim.ConsTax]));
dJ4=max(abs([GEcondnsJ4struct.CapitalMarket-GEcondnsJ4dim.CapitalMarket,GEcondnsJ4struct.GovBudget-GEcondnsJ4dim.GovBudget,GEcondnsJ4struct.ConsTax-GEcondnsJ4dim.ConsTax]));
dZ1=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ1.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ1.GovBudget,GEcondnsR.ConsTax-GEcondnsZ1.ConsTax]));
dZ3=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ3.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ3.GovBudget,GEcondnsR.ConsTax-GEcondnsZ3.ConsTax]));
dZ4=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ4.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ4.GovBudget,GEcondnsR.ConsTax-GEcondnsZ4.ConsTax]));
dZ2a=max(abs([GEcondnsR.CapitalMarket-GEcondnsZ2a.CapitalMarket,GEcondnsR.GovBudget-GEcondnsZ2a.GovBudget,GEcondnsR.ConsTax-GEcondnsZ2a.ConsTax]));
dZ2b=max(abs([GEcondnsZ.CapitalMarket-GEcondnsZfrozen.CapitalMarket,GEcondnsZ.GovBudget-GEcondnsZfrozen.GovBudget,GEcondnsZ.ConsTax-GEcondnsZfrozen.ConsTax]));
dC1=max(abs([GEcondnsC.CapitalMarket-GEcondnsCfrozen.CapitalMarket,GEcondnsC.GovBudget-GEcondnsCfrozen.GovBudget,GEcondnsC.ConsTax-GEcondnsCfrozen.ConsTax]));

fprintf('\n=== FHorz PType: intermediateEqns, the four jequaloneDist forms, shock grids in GE ===\n')
fprintf('GE conditions at the initial prices, baseline: CapitalMarket=%.6e GovBudget=%.6e ConsTax=%.6e \n',GEcondnsR.CapitalMarket,GEcondnsR.GovBudget,GEcondnsR.ConsTax)
fprintf('intermediateEqns, wage as an intermediate, this should be zero: %.3e \n',dI1)
fprintf('intermediateEqns, one intermediate consuming another, this should be zero: %.3e \n',dI2)
fprintf('intermediateEqns, an intermediate consuming an AggVar, this should be zero: %.3e \n',dI3)
fprintf('jequaloneDist as a function vs the same distribution as a matrix, this should be zero: %.3e \n',dJ1)
fprintf('jequaloneDist as a struct of one function per ptype vs one function for all, this should be zero: %.3e \n',dJ2)
fprintf('jequaloneDist with ptype as a dimension vs as a struct, masses agreeing, this should be zero: %.3e \n',dJ3)
fprintf('jequaloneDist with ptype as a dimension, its slice masses overriding PTypeDistParamNames, this should be zero: %.3e \n',dJ4)
fprintf('pi_z in GE for both ptypes, rebuilt to the same grid, this should be zero: %.3e \n',dZ1)
fprintf('pi_z in GE, age-dependent, rebuilt to the same grid, this should be zero: %.3e \n',dZ3)
fprintf('pi_z in GE for the second ptype only (gridsinGE=[0 1]), this should be zero: %.3e \n',dZ4)
fprintf('pi_z in GE, grids geared to r, this should NOT be zero: %.3e \n',dZ2a)
fprintf('  (that solve landed at r=%.6f Tr=%.6f tau_c=%.6f) \n',p_eqmZ.r,p_eqmZ.Tr,p_eqmZ.tau_c)
fprintf('pi_z in GE, solver grids vs grids frozen at its own answer, this should be zero: %.3e \n',dZ2b)
fprintf('  (the jequaloneDist-as-a-function solve landed at r=%.6f Tr=%.6f tau_c=%.6f) \n',p_eqmC.r,p_eqmC.Tr,p_eqmC.tau_c)
fprintf('jequaloneDist as a function with the grids geared to r, frozen-grid check, this should be zero: %.3e \n',dC1)

output.GEcondnsR=GEcondnsR;
output.GEcondnsI1=GEcondnsI1;
output.GEcondnsI2=GEcondnsI2;
output.GEcondnsI3=GEcondnsI3;
output.GEcondnsJ1num=GEcondnsJ1num;
output.GEcondnsJ1fn=GEcondnsJ1fn;
output.GEcondnsJ2=GEcondnsJ2;
output.GEcondnsJ3struct=GEcondnsJ3struct;
output.GEcondnsJ3dim=GEcondnsJ3dim;
output.GEcondnsJ4struct=GEcondnsJ4struct;
output.GEcondnsJ4dim=GEcondnsJ4dim;
output.GEcondnsZ1=GEcondnsZ1;
output.GEcondnsZ3=GEcondnsZ3;
output.GEcondnsZ4=GEcondnsZ4;
output.GEcondnsZ2a=GEcondnsZ2a;
output.GEcondnsZ=GEcondnsZ;
output.GEcondnsZfrozen=GEcondnsZfrozen;
output.GEcondnsC=GEcondnsC;
output.GEcondnsCfrozen=GEcondnsCfrozen;
output.p_eqmZ=p_eqmZ;
output.p_eqmC=p_eqmC;

end
