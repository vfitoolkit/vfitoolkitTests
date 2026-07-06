function output=CoreStationaryGE_InfHorz_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz stationary GE: solve the same 3-price/3-eqn model with
% fminalgo=1 (fminsearch), 5 (shooting), 8 (lsqnonlin), 4 (CMA-ES),
% 9 (Anderson acceleration of the shooting map) and confirm they agree
% (CMA-ES is stochastic/lower-accuracy, so a looser match).
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC (computed inside the ReturnFn, and
% inline wherever else it is needed), so it is NOT a GE price.

n_p=0;
heteroagentoptions.countGEsolves=1; % record number of model solves (GE-condn evaluations) for each fminalgo

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);

% Aggregates needed by the GE eqns
FnsToEvaluate.K=@(d,aprime,a,z) a;      % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z) d*z;    % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption

% GE eqns: written as LHS-RHS so that =0 at equilibrium (w hardcoded from r)
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

%% fminalgo=1 (fminsearch)
heteroagentoptions1=heteroagentoptions;
heteroagentoptions1.fminalgo=1;
StationaryGeneralEqm_subcode_GEsolvecounter('reset'); tt=tic;
[p_eqm1,GEcondns1]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions1, simoptions, vfoptions);
time1=toc(tt); nsolves1=StationaryGeneralEqm_subcode_GEsolvecounter('get');

% At the fminalgo=1 equilibrium: compute V and StationaryDist, then plot the cumulative asset distribution
Params1=Params;
Params1.r=p_eqm1.r; Params1.Tr=p_eqm1.Tr; Params1.tau_c=p_eqm1.tau_c;
[V,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,simoptions,Params1,[]);
assetdist=sum(StationaryDist,2); % marginal distribution over assets (sum out z)
figure;
plot(a_grid,cumsum(assetdist))
title('InfHorz: cumulative distribution over asset grid'); xlabel('assets'); ylabel('cumulative probability')

%% fminalgo=5 (shooting; need update rules, one row per GE eqn: {GEeqnName, PriceName, add, factor})
heteroagentoptions5=heteroagentoptions;
heteroagentoptions5.fminalgo=5;
heteroagentoptions5.fminalgo5.howtoupdate={...
    'CapitalMarket','r',0,0.005;   % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.05;        % Tr_new = Tr + factor*(tau*w*N-Tr)
    'ConsTax','tau_c',0,0.05};      % tau_c_new = tau_c - factor*(tau_c*C-G)
heteroagentoptions5.maxiter=1e4;
StationaryGeneralEqm_subcode_GEsolvecounter('reset'); tt=tic;
[p_eqm5,GEcondns5]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions5, simoptions, vfoptions);
time5=toc(tt); nsolves5=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=8 (lsqnonlin)
heteroagentoptions8=heteroagentoptions;
heteroagentoptions8.fminalgo=8;
heteroagentoptions8.toleranceGEcondns=heteroagentoptions.toleranceGEcondns/10; % lsqnonlin() is much looser on accuracy
heteroagentoptions8.toleranceGEprices=heteroagentoptions.toleranceGEprices/10;
StationaryGeneralEqm_subcode_GEsolvecounter('reset'); tt=tic;
[p_eqm8,GEcondns8]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions8, simoptions, vfoptions);
time8=toc(tt); nsolves8=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=4 (CMA-ES; slow but very globally robust)
heteroagentoptions4=heteroagentoptions;
heteroagentoptions4.fminalgo=4;
StationaryGeneralEqm_subcode_GEsolvecounter('reset'); tt=tic;
[p_eqm4,GEcondns4]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions4, simoptions, vfoptions);
time4=toc(tt); nsolves4=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=9 (Anderson acceleration of the shooting fixed-point map)
% Accelerates the same map as fminalgo=5, using the same howtoupdate rules.
heteroagentoptions9=heteroagentoptions;
heteroagentoptions9.fminalgo=9;
heteroagentoptions9.fminalgo9.howtoupdate={...
    'CapitalMarket','r',0,0.005;   % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.05;        % Tr_new = Tr + factor*(tau*w*N-Tr)
    'ConsTax','tau_c',0,0.05};      % tau_c_new = tau_c - factor*(tau_c*C-G)
heteroagentoptions9.anderson.maxiter=1e4;
StationaryGeneralEqm_subcode_GEsolvecounter('reset'); tt=tic;
[p_eqm9,GEcondns9]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions9, simoptions, vfoptions);
time9=toc(tt); nsolves9=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% Compare
fprintf('\n=== InfHorz: fminalgo agreement (r,Tr,tau_c) ===\n')
fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm1.r,p_eqm1.Tr,p_eqm1.tau_c)
fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm5.r,p_eqm5.Tr,p_eqm5.tau_c)
fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm8.r,p_eqm8.Tr,p_eqm8.tau_c)
fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm4.r,p_eqm4.Tr,p_eqm4.tau_c)
fprintf('fminalgo=9: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm9.r,p_eqm9.Tr,p_eqm9.tau_c)
d15=max(abs([p_eqm1.r-p_eqm5.r,p_eqm1.Tr-p_eqm5.Tr,p_eqm1.tau_c-p_eqm5.tau_c]));
d18=max(abs([p_eqm1.r-p_eqm8.r,p_eqm1.Tr-p_eqm8.Tr,p_eqm1.tau_c-p_eqm8.tau_c]));
d14=max(abs([p_eqm1.r-p_eqm4.r,p_eqm1.Tr-p_eqm4.Tr,p_eqm1.tau_c-p_eqm4.tau_c]));
d19=max(abs([p_eqm1.r-p_eqm9.r,p_eqm1.Tr-p_eqm9.Tr,p_eqm1.tau_c-p_eqm9.tau_c]));
fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
fprintf('fminalgo 1 vs 4, this should be near zero: %.8f \n',d14)
fprintf('fminalgo 1 vs 9, this should be near zero: %.8f \n',d19)

fprintf('\n=== InfHorz: fminalgo runtime & model solves ===\n')
fprintf('fminalgo=1 (fminsearch): time=%7.2fs, model solves=%d \n',time1,nsolves1)
fprintf('fminalgo=5 (shooting):   time=%7.2fs, model solves=%d \n',time5,nsolves5)
fprintf('fminalgo=8 (lsqnonlin):  time=%7.2fs, model solves=%d \n',time8,nsolves8)
fprintf('fminalgo=4 (CMA-ES):     time=%7.2fs, model solves=%d \n',time4,nsolves4)
fprintf('fminalgo=9 (Anderson):   time=%7.2fs, model solves=%d \n',time9,nsolves9)

output.p_eqm1=p_eqm1;
output.p_eqm5=p_eqm5;
output.p_eqm8=p_eqm8;
output.p_eqm4=p_eqm4;
output.p_eqm9=p_eqm9;

output.time1=time1;   output.nsolves1=nsolves1;
output.time5=time5;   output.nsolves5=nsolves5;
output.time8=time8;   output.nsolves8=nsolves8;
output.time4=time4;   output.nsolves4=nsolves4;
output.time9=time9;   output.nsolves9=nsolves9;

output.GEcondns1=GEcondns1;
output.GEcondns5=GEcondns5;
output.GEcondns8=GEcondns8;
output.GEcondns4=GEcondns4;
output.GEcondns9=GEcondns9;

end
