function output=CoreStationaryGE_FHorz_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% FHorz stationary GE: solve unconstrained, then re-solve applying each of
% the three kinds of parameter constraint (each to a different price):
%   constrain0to1     on r
%   constrainpositive on Tr (three times: constrainpositivemethod left unset, set
%                     explicitly to 'softplus', and set to 'log'. Since 2026-09-12
%                     unset IS softplus, so the first two are the same transform and
%                     must agree exactly; 'log' is the other map onto (0,infty), so it
%                     must reach the same equilibrium, though less accurately here)
%   constrainAtoB     on tau_c (limits [0,1])
% and then all three at once. All equilibria are interior, so the constrained
% solves must reproduce the unconstrained answer.
%
% SINCE 2026-10-09 the single constrain0to1-on-r and constrainAtoB-on-tau_c solves are not run: 'all
% three at once' puts both transforms through a full solve, so they are covered there (if that check ever
% fails, putting the singles back says which transform broke). This cut the part from 7 converged solves
% to 4, as part of getting the bank under a day.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% w is hardcoded from r via the firm FOC, so it is NOT a GE price.

n_p=0;
vfoptions.divideandconquer=1; % finite-horizon tests use divide-and-conquer

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);

FnsToEvaluate.K=@(d,aprime,a,z) a;
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d;
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c);

GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;

% All solves here use fminalgo=1 (so any difference is due to the constraint transform only)
heteroagentoptions.fminalgo=1;

%% Unconstrained reference
[p_eqm0,GEcondns0]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

%% constrainpositive on Tr (constrainpositivemethod left unset, so this is the default)
heteroagentoptionsB=heteroagentoptions;
heteroagentoptionsB.constrainpositive={'Tr'};
[p_eqmB,GEcondnsB]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsB, simoptions, vfoptions);

%% constrainpositive on Tr: default vs softplus requested explicitly, both stopped by maxiter=20
% softplus has been the default since 2026-09-12, so case B (option unset) is already softplus, and the
% two must agree EXACTLY, not merely to tolerance: what this checks is that the default in this solver's
% own option handling has not silently reverted to log. An exact identity holds along the whole path,
% not only at the answer, so neither solve needs to converge: both are stopped by maxiter=20 (for fminsearch that is 20 model
% solves, MaxIter=MaxFunEvals=maxiter) instead of two full solves. (Until 2026-10-09 the softplus run went to convergence
% and was compared only to the unconstrained answer, never directly to case B.)
heteroagentoptionsB1=heteroagentoptions;
heteroagentoptionsB1.constrainpositive={'Tr'};
heteroagentoptionsB1.maxiter=20;
[p_eqmB1,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsB1, simoptions, vfoptions);
heteroagentoptionsB2=heteroagentoptionsB1;
heteroagentoptionsB2.constrainpositivemethod='softplus';
[p_eqmB2,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsB2, simoptions, vfoptions);

%% constrainpositive on Tr, with log requested explicitly
% uparam=log(cparam), cparam=exp(uparam). The other map onto (0,infty), and the
% default before 2026-09-12. It must find the same equilibrium as case B, but it
% is known to get there less accurately on this model: when the default was
% switched, this very check came out at 5.3e-2 under log against 2.6e-5 under
% softplus. So the threshold on it is deliberately loose, and what the bank
% watches is the ratio of the two errors. This is also the only case that puts
% a non-default constrainpositivemethod through a full GE solve, which is why
% the print block checks that it DIFFERS from the softplus case: an exact zero
% there means the option never reached the solver.
heteroagentoptionsB3=heteroagentoptions;
heteroagentoptionsB3.constrainpositive={'Tr'};
heteroagentoptionsB3.constrainpositivemethod='log';
[p_eqmB3,GEcondnsB3]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsB3, simoptions, vfoptions);

%% all three constraints at once
heteroagentoptionsD=heteroagentoptions;
heteroagentoptionsD.constrain0to1={'r'};
heteroagentoptionsD.constrainpositive={'Tr'};
heteroagentoptionsD.constrainAtoB={'tau_c'};
heteroagentoptionsD.constrainAtoBlimits.tau_c=[0,1];
[p_eqmD,GEcondnsD]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptionsD, simoptions, vfoptions);

%% Compare each constrained solve to the unconstrained reference
fprintf('\n=== FHorz: parameter-constraint invariance ===\n')
fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm0.r,p_eqm0.Tr,p_eqm0.tau_c)
fprintf('constrainpos default: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB.r,p_eqmB.Tr,p_eqmB.tau_c)
fprintf('constrainpos default,  stopped by maxiter=20: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB1.r,p_eqmB1.Tr,p_eqmB1.tau_c)
fprintf('constrainpos softplus, stopped by maxiter=20: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB2.r,p_eqmB2.Tr,p_eqmB2.tau_c)
fprintf('constrainpos log:     r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB3.r,p_eqmB3.Tr,p_eqmB3.tau_c)
fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmD.r,p_eqmD.Tr,p_eqmD.tau_c)
dB=max(abs([p_eqm0.r-p_eqmB.r,p_eqm0.Tr-p_eqmB.Tr,p_eqm0.tau_c-p_eqmB.tau_c]));
dB1B2=max(abs([p_eqmB1.r-p_eqmB2.r,p_eqmB1.Tr-p_eqmB2.Tr,p_eqmB1.tau_c-p_eqmB2.tau_c]));
dB3=max(abs([p_eqm0.r-p_eqmB3.r,p_eqm0.Tr-p_eqmB3.Tr,p_eqm0.tau_c-p_eqmB3.tau_c]));
dBB3=max(abs([p_eqmB.r-p_eqmB3.r,p_eqmB.Tr-p_eqmB3.Tr,p_eqmB.tau_c-p_eqmB3.tau_c]));
dD=max(abs([p_eqm0.r-p_eqmD.r,p_eqm0.Tr-p_eqmD.Tr,p_eqm0.tau_c-p_eqmD.tau_c]));
fprintf('constrainpositive on Tr, default method, this should be near zero: %.8f \n',dB)
fprintf('constrainpositive on Tr, default vs softplus requested explicitly (both stopped by maxiter=20), this should be zero: %.3e \n',dB1B2)
fprintf('constrainpositive on Tr with log, this should be near zero (loose: log is the worse transform on this model): %.8f \n',dB3)
fprintf('  ratio of the log error to the softplus error, log is expected to be the larger: %.1f \n',dB3/dB)
fprintf('log vs softplus (the default), this should NOT be zero: %.3e \n',dBB3)
fprintf('constrain all three, this should be near zero: %.8f \n',dD)

output.p_eqm0=p_eqm0;
output.p_eqmB=p_eqmB;
output.p_eqmB1=p_eqmB1;
output.p_eqmB2=p_eqmB2;
output.p_eqmB3=p_eqmB3;
output.p_eqmD=p_eqmD;

output.GEcondns0=GEcondns0;
output.GEcondnsB=GEcondnsB;
output.GEcondnsB3=GEcondnsB3;
output.GEcondnsD=GEcondnsD;

end
