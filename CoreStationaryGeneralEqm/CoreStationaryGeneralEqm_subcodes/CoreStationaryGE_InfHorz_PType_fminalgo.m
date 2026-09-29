function output=CoreStationaryGE_InfHorz_PType_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptions,simoptions,vfoptions)
% InfHorz stationary GE with permanent types (PType): same 3-price/3-eqn model
% as CoreStationaryGE_InfHorz_fminalgo, but with N_i=2 permanent types that
% differ in the value of sigma (2.2 and 1.8). Solve with fminalgo=1,5,8,4,9 and
% confirm they agree.
%
% GE prices: r (capital market), Tr (gov budget), tau_c (consumption-tax budget).
% These are economy-wide (not ptype-dependent). w is hardcoded from r.

n_p=0;
heteroagentoptions.countGEsolves=1; % record number of model solves (GE-condn evaluations) for each fminalgo

% Permanent types: two types differing in sigma
N_i=2;
Names_i=N_i;
Params.sigma=[2.2,1.8];       % differs by permanent type
Params.ptypemass=[0.5,0.5];   % mass of each permanent type
PTypeDistParamNames={'ptypemass'};

ReturnFn=@(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi) ...
    ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi);

% Aggregates needed by the GE eqns (economy-wide, aggregated across ptypes)
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
tt=tic;
[p_eqm1,GEcondns1]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions1, simoptions, vfoptions);
time1=toc(tt); nsolves1=StationaryGeneralEqm_subcode_GEsolvecounter('get');

% At the fminalgo=1 equilibrium: compute V and StationaryDist, then plot the (economy-wide) cumulative asset distribution
Params1=Params;
Params1.r=p_eqm1.r; Params1.Tr=p_eqm1.Tr; Params1.tau_c=p_eqm1.tau_c;
[V,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,Params1,simoptions);
% economy-wide marginal over assets = sum over ptypes of ptweight*(asset marginal of that type)
fns=fieldnames(StationaryDist); fns=fns(~strcmp(fns,'ptweights'));
assetdist=zeros(n_a,1);
for ii=1:numel(fns)
    assetdist=assetdist+StationaryDist.ptweights(ii)*sum(gather(StationaryDist.(fns{ii})),2);
end
figure;
plot(a_grid,cumsum(assetdist))
title('InfHorz PType: cumulative distribution over asset grid'); xlabel('assets'); ylabel('cumulative probability')

%% fminalgo=5 (shooting; need update rules, one row per GE eqn: {GEeqnName, PriceName, add, factor})
heteroagentoptions5=heteroagentoptions;
heteroagentoptions5.fminalgo=5;
% The r factor is 0.001. Measured, not guessed: at the previous value the shooting map's multiplier
% |1-factor*dCapitalMarket/dr| was at or past its stability boundary, so the solve reached about
% 1e-05 and then sat in a period-2 limit cycle (the condition alternating sign every iteration with
% constant amplitude ~2.4e-05, above toleranceGEcondns) and ran to maxiter instead of exiting. The
% 2026-09-23 diary of doPart(17) shows it directly; that fixed dCapitalMarket/dr at about 1000 for
% this model, which is why 0.001 puts the multiplier near zero. dc/dTr is -1 and dc/dtau_c is C, so
% those two rows are nowhere near the boundary and keep their factors.
heteroagentoptions5.fminalgo5.howtoupdate={...
    'CapitalMarket','r',0,0.001;   % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.02;        % Tr_new = Tr + factor*(tau*w*N-Tr)
    'ConsTax','tau_c',0,0.02};      % tau_c_new = tau_c - factor*(tau_c*C-G)
heteroagentoptions5.maxiter=3000; % backstop: with the factor above this converges in tens of iterations
tt=tic;
[p_eqm5,GEcondns5]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions5, simoptions, vfoptions);
time5=toc(tt); nsolves5=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=8 (lsqnonlin)
heteroagentoptions8=heteroagentoptions;
heteroagentoptions8.fminalgo=8;
tt=tic;
[p_eqm8,GEcondns8]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions8, simoptions, vfoptions);
time8=toc(tt); nsolves8=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=4 (CMA-ES; slow but very globally robust)
heteroagentoptions4=heteroagentoptions;
heteroagentoptions4.fminalgo=4;
% insigma and MaxFunEvals are set here rather than left to the defaults. The default insigma is 30% of
% |p0|, which at p0=[0.04,0.2,0.1] is a 1-sigma spread on r of +-0.012: CMA-ES then routinely samples r
% near zero and negative, and most of its evaluations go on shrinking that ellipsoid rather than
% locating the equilibrium (4,595 model solves, 1.4 hours, in the 2026-09-20 run of this part - and it
% is what the search looks like when watched, prices jumping back and forth). 5% keeps it a stochastic
% search without the excursions. MaxFunEvals then bounds the part: the default is Inf, with only
% StopFitness and stagnation detection to stop it. CMA-ES returns its best-so-far, which is all the
% 1-vs-4 comparison below needs - and that comparison is deliberately the loosest in the file.
heteroagentoptions4.insigma=0.05*abs([Params.r;Params.Tr;Params.tau_c]); % in GEPriceParamNames order
heteroagentoptions4.inopts.MaxFunEvals=1500;
tt=tic;
[p_eqm4,GEcondns4]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions4, simoptions, vfoptions);
time4=toc(tt); nsolves4=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=9 (Anderson acceleration of the shooting fixed-point map)
% Accelerates the same map as fminalgo=5, using the same howtoupdate rules.
heteroagentoptions9=heteroagentoptions;
heteroagentoptions9.fminalgo=9;
heteroagentoptions9.fminalgo9.howtoupdate={... % same rules as fminalgo=5 above, including the r factor
    'CapitalMarket','r',0,0.001;   % r_new = r - factor*(r-MPK)
    'GovBudget','Tr',1,0.02;        % Tr_new = Tr + factor*(tau*w*N-Tr)
    'ConsTax','tau_c',0,0.02};      % tau_c_new = tau_c - factor*(tau_c*C-G)
heteroagentoptions9.anderson.maxiter=1e4;
% A shorter memory than the default of 5. With the default settings Anderson extrapolated r clean out
% of the feasible region (r below -delta, so the MPK inversion took a fractional power of a negative
% number and the model threw); fewer past iterates in the mix shrinks the extrapolation, and the
% safeguard now rejects such a step rather than the model ending the run.
% regularization is left at its default ON PURPOSE. It was 1e-6 for the 2026-09-26 run, which was a
% mistake: the term was absolute back then, and 1e-6 against normal equations of size ~1e-19 meant the
% Anderson step was pure regularization, i.e. exactly the plain shooting step it should accelerate
% (fminalgo=9 took the same iterations as fminalgo=5 and paid the safeguard's second solve, so it was
% strictly slower). The term is relative from 2026-09-27, so the default is the right thing to run.
heteroagentoptions9.anderson.memory=3;
tt=tic;
[p_eqm9,GEcondns9]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions9, simoptions, vfoptions);
time9=toc(tt); nsolves9=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% fminalgo=9 with Anderson Type-I (Zhang, O'Donoghue & Boyd 2020; globally convergent)
% The Type-I step x-H*g is unbounded, and on this model it overshoots: on
% 2026-09-17 it drove r below -delta at iteration 108, where the wage the
% ReturnFn builds from the firm FOC, (1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1)),
% has a negative base and a fractional exponent. That is a complex result, which
% gpuArray/arrayfun cannot return, so the whole test bank died there. Note that
% AndersonAcceleration already halves a Type-I step back toward the previous
% iterate when the GE conditions come back non-finite, but that guard tests
% isfinite() and this failure THROWS instead, so it never fires.
% constrainpositive on r is what keeps it out of that region: Anderson then
% iterates on uparam over the whole real line while the solver only ever sees
% r=softplus(uparam)>0, so r+delta can never go negative however long the step.
heteroagentoptions9I=heteroagentoptions9;   % same fminalgo9.howtoupdate as the Type-II run
heteroagentoptions9I.anderson.type='I';
heteroagentoptions9I.constrainpositive={'r'};
tt=tic;
[p_eqm9I,GEcondns9I]=HeteroAgentStationaryEqm_InfHorz_PType(n_d, n_a, n_z, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, PTypeDistParamNames, GEPriceParamNames,heteroagentoptions9I, simoptions, vfoptions);
time9I=toc(tt); nsolves9I=StationaryGeneralEqm_subcode_GEsolvecounter('get');

%% Compare
% Note that the 9I run is the one solve here that carries a parameter constraint
% (constrainpositive on r, see above). The equilibrium r is interior to (0,infty),
% so the constraint is not binding and 9I must still agree with the rest; it is
% only the path the optimizer takes to get there that differs.
fprintf('\n=== InfHorz PType: fminalgo agreement (r,Tr,tau_c) ===\n')
fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm1.r,p_eqm1.Tr,p_eqm1.tau_c)
fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm5.r,p_eqm5.Tr,p_eqm5.tau_c)
fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm8.r,p_eqm8.Tr,p_eqm8.tau_c)
fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm4.r,p_eqm4.Tr,p_eqm4.tau_c)
fprintf('fminalgo=9: r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm9.r,p_eqm9.Tr,p_eqm9.tau_c)
fprintf('fminalgo=9I (Type-I): r=%.6f Tr=%.6f tau_c=%.6f \n',p_eqm9I.r,p_eqm9I.Tr,p_eqm9I.tau_c)
d15=max(abs([p_eqm1.r-p_eqm5.r,p_eqm1.Tr-p_eqm5.Tr,p_eqm1.tau_c-p_eqm5.tau_c]));
d18=max(abs([p_eqm1.r-p_eqm8.r,p_eqm1.Tr-p_eqm8.Tr,p_eqm1.tau_c-p_eqm8.tau_c]));
d14=max(abs([p_eqm1.r-p_eqm4.r,p_eqm1.Tr-p_eqm4.Tr,p_eqm1.tau_c-p_eqm4.tau_c]));
d19=max(abs([p_eqm1.r-p_eqm9.r,p_eqm1.Tr-p_eqm9.Tr,p_eqm1.tau_c-p_eqm9.tau_c]));
d19I=max(abs([p_eqm1.r-p_eqm9I.r,p_eqm1.Tr-p_eqm9I.Tr,p_eqm1.tau_c-p_eqm9I.tau_c]));
fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
fprintf('fminalgo 1 vs 4 (CMA-ES, lower accuracy), this should be small: %.8f \n',d14)
fprintf('fminalgo 1 vs 9, this should be near zero: %.8f \n',d19)
fprintf('fminalgo 1 vs 9I (Type-I), this should be near zero: %.8f \n',d19I)

fprintf('\n=== InfHorz PType: fminalgo runtime & model solves ===\n')
fprintf('fminalgo=1 (fminsearch): time=%7.2fs, model solves=%d \n',time1,nsolves1)
fprintf('fminalgo=5 (shooting):   time=%7.2fs, model solves=%d \n',time5,nsolves5)
fprintf('fminalgo=8 (lsqnonlin):  time=%7.2fs, model solves=%d \n',time8,nsolves8)
fprintf('fminalgo=4 (CMA-ES):     time=%7.2fs, model solves=%d \n',time4,nsolves4)
fprintf('fminalgo=9 (Anderson):   time=%7.2fs, model solves=%d \n',time9,nsolves9)
fprintf('fminalgo=9 Type-I (Anders.I): time=%7.2fs, model solves=%d \n',time9I,nsolves9I)

output.p_eqm1=p_eqm1;
output.p_eqm5=p_eqm5;
output.p_eqm8=p_eqm8;
output.p_eqm4=p_eqm4;
output.p_eqm9=p_eqm9;
output.p_eqm9I=p_eqm9I;

output.time1=time1;   output.nsolves1=nsolves1;
output.time5=time5;   output.nsolves5=nsolves5;
output.time8=time8;   output.nsolves8=nsolves8;
output.time4=time4;   output.nsolves4=nsolves4;
output.time9=time9;   output.nsolves9=nsolves9;
output.time9I=time9I; output.nsolves9I=nsolves9I;

output.GEcondns1=GEcondns1;
output.GEcondns5=GEcondns5;
output.GEcondns8=GEcondns8;
output.GEcondns4=GEcondns4;
output.GEcondns9=GEcondns9;
output.GEcondns9I=GEcondns9I;

end
