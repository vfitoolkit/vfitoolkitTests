function output=CoreStationaryGE_FHorz_lsqnonlinweights(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% FHorz stationary GE: heteroagentoptions.multiGEweights must reach lsqnonlin (fminalgo=8), as sqrt(weights).
%
% Until 2026-10-07 they did not. The command took sqrt(multiGEweights) before building the objective,
% but the subfn's vector output ignored the weights entirely, so lsqnonlin always solved the
% unweighted problem and the sqrt was dead code. fminalgo=8 now uses outputGEform=3, which returns
% sqrt(multiGEweights).*GEcondns, so lsqnonlin's sum of squares is sum(multiGEweights.*GEcondns.^2):
% the square of the multiGEcriterion=1 scalar that fminsearch minimizes.
%
% WHY THE MODEL HAS A FOURTH, CONFLICTING GENERAL EQM CONDITION. The first version of this part (run
% 2026-10-08) used the bank's 3-price/3-condition model and compared weights [1,1,1] with [25,1,1]. The two
% answers differed by 4e-16 - not zero, so the weights did reach lsqnonlin, but nothing a test can see.
% That is not a fluke: with as many conditions as prices the Jacobian is square, and the Gauss-Newton
% step solves J*d=-F, which row-scaling by the weights does not change. Weights only matter when the
% conditions cannot all hold at once. So here a fourth condition, rTarget: r-rtarget=0, is added with
% rtarget=0.07 against the model's own r of about 0.080, and the solvers have to trade CapitalMarket
% off against it. Then:
%  - weighted vs unweighted lsqnonlin must NOT agree (the weights reached it), and the weighted answer
%    must sit closer to rtarget (they pulled the right way);
%  - weighted lsqnonlin and weighted fminsearch minimize the same sum(w.*GEcondns.^2), so they must find
%    the same prices. That is the check that the weights reach lsqnonlin as sqrt(w): had they arrived as
%    w, lsqnonlin would be minimizing sum(w.^2.*GEcondns.^2) and land somewhere else.
% FHorz rather than InfHorz because lsqnonlin is known to stop short on InfHorz (doPart(2)).
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).

n_p=0;
vfoptions.divideandconquer=1; % finite-horizon tests use divide-and-conquer
heteroagentoptions.countGEsolves=1; % record number of model solves (GE-condn evaluations)

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_FHorz(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A,sigma,eta,varphi);

% Aggregates needed by the GE eqns
FnsToEvaluate.K=@(d,aprime,a,z) a;                   % aggregate capital (assets)
FnsToEvaluate.N=@(d,aprime,a,z,kappa_j) kappa_j*z*d; % aggregate effective labor
FnsToEvaluate.C=@(d,aprime,a,z,r,tau,Tr,tau_c,kappa_j,alpha,delta,A) ((1+r)*a+(1-tau)*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*kappa_j*z*d+Tr-aprime)/(1+tau_c); % aggregate consumption

% GE eqns: written as LHS-RHS so that =0 at equilibrium (w hardcoded from r), plus the conflicting fourth
Params.rtarget=0.07;
GeneralEqmEqns.CapitalMarket=@(r,K,N,alpha,delta,A) r-(alpha*A*(K^(alpha-1))*(N^(1-alpha))-delta);
GeneralEqmEqns.GovBudget=@(tau,r,N,Tr,alpha,delta,A) tau*((1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)))*N-Tr;
GeneralEqmEqns.ConsTax=@(tau_c,C,G) tau_c*C-G;
GeneralEqmEqns.rTarget=@(r,rtarget) r-rtarget;

%% fminalgo=8 (lsqnonlin), unweighted
heteroagentoptions8=heteroagentoptions;
heteroagentoptions8.fminalgo=8;
heteroagentoptions8.multiGEweights=[1,1,1,1];
tt=tic;
[p_eqm8,GEcondns8]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions8, simoptions, vfoptions);
time8=toc(tt); nsolves8=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=8 (lsqnonlin), weight 100 on rTarget (so its residual is scaled by 10)
heteroagentoptions8w=heteroagentoptions;
heteroagentoptions8w.fminalgo=8;
heteroagentoptions8w.multiGEweights=[1,1,1,100];
tt=tic;
[p_eqm8w,GEcondns8w]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions8w, simoptions, vfoptions);
time8w=toc(tt); nsolves8w=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=1 (fminsearch), the same weights, minimizing sqrt(sum(w.*GEcondns.^2)) (multiGEcriterion=1)
heteroagentoptions1w=heteroagentoptions;
heteroagentoptions1w.fminalgo=1;
heteroagentoptions1w.multiGEweights=[1,1,1,100];
tt=tic;
[p_eqm1w,GEcondns1w]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], GEPriceParamNames,heteroagentoptions1w, simoptions, vfoptions);
time1w=toc(tt); nsolves1w=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% Compare
fprintf('\n=== FHorz: multiGEweights reach lsqnonlin (fminalgo=8), 4 conditions in 3 prices ===\n')
fprintf('lsqnonlin  [1,1,1,1]:   r=%.8f Tr=%.8f tau_c=%.8f | CapitalMarket=%.3e rTarget=%.3e \n',p_eqm8.r,p_eqm8.Tr,p_eqm8.tau_c,GEcondns8.CapitalMarket,GEcondns8.rTarget)
fprintf('lsqnonlin  [1,1,1,100]: r=%.8f Tr=%.8f tau_c=%.8f | CapitalMarket=%.3e rTarget=%.3e \n',p_eqm8w.r,p_eqm8w.Tr,p_eqm8w.tau_c,GEcondns8w.CapitalMarket,GEcondns8w.rTarget)
fprintf('fminsearch [1,1,1,100]: r=%.8f Tr=%.8f tau_c=%.8f | CapitalMarket=%.3e rTarget=%.3e \n',p_eqm1w.r,p_eqm1w.Tr,p_eqm1w.tau_c,GEcondns1w.CapitalMarket,GEcondns1w.rTarget)
d88w=max(abs([p_eqm8.r-p_eqm8w.r,p_eqm8.Tr-p_eqm8w.Tr,p_eqm8.tau_c-p_eqm8w.tau_c]));
d8w1w=max(abs([p_eqm8w.r-p_eqm1w.r,p_eqm8w.Tr-p_eqm1w.Tr,p_eqm8w.tau_c-p_eqm1w.tau_c]));
fprintf('weighted vs unweighted lsqnonlin, this should NOT be zero: %.3e \n',d88w)
fprintf('weighted lsqnonlin closer to rtarget than unweighted, this should be one: %d \n',abs(p_eqm8w.r-Params.rtarget)<abs(p_eqm8.r-Params.rtarget))
fprintf('weighted lsqnonlin vs weighted fminsearch (same objective), this should be near zero: %.8f \n',d8w1w)

fprintf('\n=== FHorz: lsqnonlin weights runtime & model solves ===\n')
fprintf('lsqnonlin unweighted: time=%7.2fs, model solves=%d \n',time8,nsolves8)
fprintf('lsqnonlin weighted:   time=%7.2fs, model solves=%d \n',time8w,nsolves8w)
fprintf('fminsearch weighted:  time=%7.2fs, model solves=%d \n',time1w,nsolves1w)

output.p_eqm8=p_eqm8;
output.p_eqm8w=p_eqm8w;
output.p_eqm1w=p_eqm1w;
output.GEcondns8=GEcondns8;
output.GEcondns8w=GEcondns8w;
output.GEcondns1w=GEcondns1w;
output.time8=time8;   output.nsolves8=nsolves8;
output.time8w=time8w; output.nsolves8w=nsolves8w;
output.time1w=time1w; output.nsolves1w=nsolves1w;

end
