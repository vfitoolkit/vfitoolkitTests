function output=CoreStationaryGE_FHorz_PType_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% FHorz stationary GE with permanent types (PType): same 3-price/3-eqn model as
% CoreStationaryGE_FHorz_fminalgo, but with N_i=2 permanent types that differ in
% the value of sigma (2.2 and 1.8). Solve with fminalgo=1,5,8,4 and confirm agreement.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% These are economy-wide (not ptype-dependent). w is hardcoded from r.

n_p=0;

% Permanent types: two types differing in sigma
N_i=2;
Names_i=N_i;
Params.sigma=[2.2,1.8];       % differs by permanent type
Params.ptypemass=[0.5,0.5];   % mass of each permanent type
PTypeDistParamNames={'ptypemass'};

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);

% Aggregates needed by the GE eqns (economy-wide, aggregated across ptypes)
FnsToEvaluate.K=@(d,aprime,a,z) a;                   % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption

% GE eqns: written as LHS-RHS so that =0 at equilibrium (w hardcoded from r)
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

%% fminalgo=1 (fminsearch)
heteroagentoptions1=heteroagentoptions;
heteroagentoptions1.fminalgo=1;
[p_eqm1,GEcondns1]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions1, simoptions, vfoptions);

% At the fminalgo=1 equilibrium: compute V and StationaryDist, then plot the (economy-wide) cumulative asset distribution
Params1=Params;
Params1.r=p_eqm1.r; Params1.Tr=p_eqm1.Tr; Params1.tau_c=p_eqm1.tau_c;
[V,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params1,simoptions);
% economy-wide marginal over assets = sum over ptypes of ptweight*(asset marginal of that type)
fns=fieldnames(StationaryDist); fns=fns(~strcmp(fns,'ptweights'));
assetdist=zeros(n_a,1);
for ii=1:numel(fns)
    assetdist=assetdist+StationaryDist.ptweights(ii)*sum(gather(StationaryDist.(fns{ii})),[2,3]);
end
figure;
plot(a_grid,cumsum(assetdist))
title('FHorz PType: cumulative distribution over asset grid'); xlabel('assets'); ylabel('cumulative probability')

%% fminalgo=5 (shooting; one row per GE eqn: {GEeqnName, PriceName, add, factor})
heteroagentoptions5=heteroagentoptions;
heteroagentoptions5.fminalgo=5;
heteroagentoptions5.fminalgo5.howtoupdate={...
    'CapitalMarket','r',0,0.01;    % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.1;        % Tr_new = Tr + factor*(tau*w*N-Tr)
    'ConsTax','tau_c',0,0.1};      % tau_c_new = tau_c - factor*(tau_c*C-G)
heteroagentoptions5.maxiter=1000;
[p_eqm5,GEcondns5]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions5, simoptions, vfoptions);

%% fminalgo=8 (lsqnonlin)
heteroagentoptions8=heteroagentoptions;
heteroagentoptions8.fminalgo=8;
[p_eqm8,GEcondns8]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions8, simoptions, vfoptions);

%% fminalgo=4 (CMA-ES; slow but very globally robust)
heteroagentoptions4=heteroagentoptions;
heteroagentoptions4.fminalgo=4;
[p_eqm4,GEcondns4]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions4, simoptions, vfoptions);

%% Compare
fprintf('\n=== FHorz PType: fminalgo agreement (r,Tr,tau_c) ===\n')
fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm1.r,p_eqm1.Tr,p_eqm1.tau_c)
fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm5.r,p_eqm5.Tr,p_eqm5.tau_c)
fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm8.r,p_eqm8.Tr,p_eqm8.tau_c)
fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm4.r,p_eqm4.Tr,p_eqm4.tau_c)
d15=max(abs([p_eqm1.r-p_eqm5.r,p_eqm1.Tr-p_eqm5.Tr,p_eqm1.tau_c-p_eqm5.tau_c]));
d18=max(abs([p_eqm1.r-p_eqm8.r,p_eqm1.Tr-p_eqm8.Tr,p_eqm1.tau_c-p_eqm8.tau_c]));
d14=max(abs([p_eqm1.r-p_eqm4.r,p_eqm1.Tr-p_eqm4.Tr,p_eqm1.tau_c-p_eqm4.tau_c]));
fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
fprintf('fminalgo 1 vs 4 (CMA-ES, lower accuracy), this should be small: %.8f \n',d14)

output.p_eqm1=p_eqm1;
output.p_eqm5=p_eqm5;
output.p_eqm8=p_eqm8;
output.p_eqm4=p_eqm4;

output.GEcondns1=GEcondns1;
output.GEcondns5=GEcondns5;
output.GEcondns8=GEcondns8;
output.GEcondns4=GEcondns4;

end
