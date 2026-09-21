function output=CoreStationaryGE_InfHorz_PType_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz stationary GE with permanent types (PType): same as
% CoreStationaryGE_InfHorz_constraints, but with N_i=2 permanent types that
% differ in the value of sigma (2.2 and 1.8). Solve unconstrained, then re-solve
% applying each of the three kinds of parameter constraint (each to a different
% price), and then all three at once:
%   constrain0to1     on r
%   constrainpositive on Tr (three times: constrainpositivemethod left unset, set
%                     explicitly to 'softplus', and set to 'log'. Since 2026-09-12
%                     unset IS softplus, so the first two are the same transform and
%                     must agree exactly; 'log' is the other map onto (0,infty), so it
%                     must reach the same equilibrium, though less accurately here)
%   constrainAtoB     on tau_c (limits [0,1])
% All equilibria are interior, so the constrained solves must reproduce the
% unconstrained answer.

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

% All solves here use fminalgo=1 (so any difference is due to the constraint transform only)
heteroagentoptions.fminalgo=1;

%% Unconstrained reference
[p_eqm0,GEcondns0]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions, simoptions, vfoptions);

%% constrain0to1 on r
heteroagentoptionsA=heteroagentoptions;
heteroagentoptionsA.constrain0to1={'r'};
[p_eqmA,GEcondnsA]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsA, simoptions, vfoptions);

%% constrainpositive on Tr (constrainpositivemethod left unset, so this is the default)
heteroagentoptionsB=heteroagentoptions;
heteroagentoptionsB.constrainpositive={'Tr'};
[p_eqmB,GEcondnsB]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsB, simoptions, vfoptions);

%% constrainpositive on Tr, with softplus requested explicitly
% uparam=log(exp(cparam)-1), cparam=log(1+exp(uparam)). Sends the real line to
% (0,infty) just as log does, so it must find the SAME equilibrium as case B.
% softplus has been the default since 2026-09-12, so case B (option unset) is
% already softplus: this run must reproduce case B EXACTLY, not merely to
% tolerance, and what it checks is that the default in this solver's own option
% handling has not silently reverted to log.
heteroagentoptionsB2=heteroagentoptions;
heteroagentoptionsB2.constrainpositive={'Tr'};
heteroagentoptionsB2.constrainpositivemethod='softplus';
[p_eqmB2,GEcondnsB2]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsB2, simoptions, vfoptions);

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
[p_eqmB3,GEcondnsB3]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsB3, simoptions, vfoptions);

%% constrainAtoB on tau_c (limits [0,1])
heteroagentoptionsC=heteroagentoptions;
heteroagentoptionsC.constrainAtoB={'tau_c'};
heteroagentoptionsC.constrainAtoBlimits.tau_c=[0,1];
[p_eqmC,GEcondnsC]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsC, simoptions, vfoptions);

%% all three constraints at once
heteroagentoptionsD=heteroagentoptions;
heteroagentoptionsD.constrain0to1={'r'};
heteroagentoptionsD.constrainpositive={'Tr'};
heteroagentoptionsD.constrainAtoB={'tau_c'};
heteroagentoptionsD.constrainAtoBlimits.tau_c=[0,1];
[p_eqmD,GEcondnsD]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptionsD, simoptions, vfoptions);

%% Compare each constrained solve to the unconstrained reference
fprintf('\n=== InfHorz PType: parameter-constraint invariance ===\n')
fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm0.r,p_eqm0.Tr,p_eqm0.tau_c)
fprintf('constrain0to1 on r:   r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmA.r,p_eqmA.Tr,p_eqmA.tau_c)
fprintf('constrainpos default: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB.r,p_eqmB.Tr,p_eqmB.tau_c)
fprintf('constrainpos softplus:r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB2.r,p_eqmB2.Tr,p_eqmB2.tau_c)
fprintf('constrainpos log:     r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmB3.r,p_eqmB3.Tr,p_eqmB3.tau_c)
fprintf('constrainAtoB tau_c:  r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmC.r,p_eqmC.Tr,p_eqmC.tau_c)
fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqmD.r,p_eqmD.Tr,p_eqmD.tau_c)
dA=max(abs([p_eqm0.r-p_eqmA.r,p_eqm0.Tr-p_eqmA.Tr,p_eqm0.tau_c-p_eqmA.tau_c]));
dB=max(abs([p_eqm0.r-p_eqmB.r,p_eqm0.Tr-p_eqmB.Tr,p_eqm0.tau_c-p_eqmB.tau_c]));
dB2=max(abs([p_eqm0.r-p_eqmB2.r,p_eqm0.Tr-p_eqmB2.Tr,p_eqm0.tau_c-p_eqmB2.tau_c]));
dB3=max(abs([p_eqm0.r-p_eqmB3.r,p_eqm0.Tr-p_eqmB3.Tr,p_eqm0.tau_c-p_eqmB3.tau_c]));
dB2B3=max(abs([p_eqmB2.r-p_eqmB3.r,p_eqmB2.Tr-p_eqmB3.Tr,p_eqmB2.tau_c-p_eqmB3.tau_c]));
dC=max(abs([p_eqm0.r-p_eqmC.r,p_eqm0.Tr-p_eqmC.Tr,p_eqm0.tau_c-p_eqmC.tau_c]));
dD=max(abs([p_eqm0.r-p_eqmD.r,p_eqm0.Tr-p_eqmD.Tr,p_eqm0.tau_c-p_eqmD.tau_c]));
fprintf('constrain0to1 on r, this should be near zero: %.8f \n',dA)
fprintf('constrainpositive on Tr, default method, this should be near zero: %.8f \n',dB)
fprintf('constrainpositive on Tr with softplus, this should be near zero: %.8f \n',dB2)
fprintf('constrainpositive on Tr with log, this should be near zero (loose: log is the worse transform on this model): %.8f \n',dB3)
fprintf('  ratio of the log error to the softplus error, log is expected to be the larger: %.1f \n',dB3/dB2)
fprintf('log vs softplus, this should NOT be zero: %.3e \n',dB2B3)
fprintf('constrainAtoB on tau_c, this should be near zero: %.8f \n',dC)
fprintf('constrain all three, this should be near zero: %.8f \n',dD)

output.p_eqm0=p_eqm0;
output.p_eqmA=p_eqmA;
output.p_eqmB=p_eqmB;
output.p_eqmB2=p_eqmB2;
output.p_eqmB3=p_eqmB3;
output.p_eqmC=p_eqmC;
output.p_eqmD=p_eqmD;

output.GEcondns0=GEcondns0;
output.GEcondnsA=GEcondnsA;
output.GEcondnsB=GEcondnsB;
output.GEcondnsB2=GEcondnsB2;
output.GEcondnsB3=GEcondnsB3;
output.GEcondnsC=GEcondnsC;
output.GEcondnsD=GEcondnsD;

end
