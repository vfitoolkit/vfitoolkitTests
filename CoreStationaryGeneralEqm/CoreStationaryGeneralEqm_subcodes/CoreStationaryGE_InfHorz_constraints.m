function output=CoreStationaryGE_InfHorz_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz stationary GE: solve unconstrained, then re-solve applying each of
% the three kinds of parameter constraint (each to a different price):
%   constrain0to1     on r
%   constrainpositive on Tr
%   constrainAtoB     on tau_c (limits [0,1])
% and then all three at once. All equilibria are interior, so the constrained
% solves must reproduce the unconstrained answer.
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

% All solves here use fminalgo=1 (so any difference is due to the constraint transform only)
heteroagentoptions.fminalgo=1;

%% Unconstrained reference
[p_eqm0,GEcondns0]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

%% constrain0to1 on r
heteroagentoptionsA=heteroagentoptions;
heteroagentoptionsA.constrain0to1={'r'};
[p_eqmA,GEcondnsA]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsA, simoptions, vfoptions);

%% constrainpositive on Tr
heteroagentoptionsB=heteroagentoptions;
heteroagentoptionsB.constrainpositive={'Tr'};
[p_eqmB,GEcondnsB]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsB, simoptions, vfoptions);

%% constrainAtoB on tau_c (limits [0,1])
heteroagentoptionsC=heteroagentoptions;
heteroagentoptionsC.constrainAtoB={'tau_c'};
heteroagentoptionsC.constrainAtoBlimits.tau_c=[0,1];
[p_eqmC,GEcondnsC]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsC, simoptions, vfoptions);

%% all three constraints at once
heteroagentoptionsD=heteroagentoptions;
heteroagentoptionsD.constrain0to1={'r'};
heteroagentoptionsD.constrainpositive={'Tr'};
heteroagentoptionsD.constrainAtoB={'tau_c'};
heteroagentoptionsD.constrainAtoBlimits.tau_c=[0,1];
[p_eqmD,GEcondnsD]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsD, simoptions, vfoptions);

%% Compare each constrained solve to the unconstrained reference
fprintf('\n=== InfHorz: parameter-constraint invariance ===\n')
fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm0.r,p_eqm0.Tr,p_eqm0.tau_c)
fprintf('constrain0to1 on r:   r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmA.r,p_eqmA.Tr,p_eqmA.tau_c)
fprintf('constrainpositive Tr: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB.r,p_eqmB.Tr,p_eqmB.tau_c)
fprintf('constrainAtoB tau_c:  r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmC.r,p_eqmC.Tr,p_eqmC.tau_c)
fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmD.r,p_eqmD.Tr,p_eqmD.tau_c)
dA=max(abs([p_eqm0.r-p_eqmA.r,p_eqm0.Tr-p_eqmA.Tr,p_eqm0.tau_c-p_eqmA.tau_c]));
dB=max(abs([p_eqm0.r-p_eqmB.r,p_eqm0.Tr-p_eqmB.Tr,p_eqm0.tau_c-p_eqmB.tau_c]));
dC=max(abs([p_eqm0.r-p_eqmC.r,p_eqm0.Tr-p_eqmC.Tr,p_eqm0.tau_c-p_eqmC.tau_c]));
dD=max(abs([p_eqm0.r-p_eqmD.r,p_eqm0.Tr-p_eqmD.Tr,p_eqm0.tau_c-p_eqmD.tau_c]));
fprintf('constrain0to1 on r, this should be near zero: %.8f \n',dA)
fprintf('constrainpositive on Tr, this should be near zero: %.8f \n',dB)
fprintf('constrainAtoB on tau_c, this should be near zero: %.8f \n',dC)
fprintf('constrain all three, this should be near zero: %.8f \n',dD)

output.p_eqm0=p_eqm0;
output.p_eqmA=p_eqmA;
output.p_eqmB=p_eqmB;
output.p_eqmC=p_eqmC;
output.p_eqmD=p_eqmD;

output.GEcondns0=GEcondns0;
output.GEcondnsA=GEcondnsA;
output.GEcondnsB=GEcondnsB;
output.GEcondnsC=GEcondnsC;
output.GEcondnsD=GEcondnsD;

end
