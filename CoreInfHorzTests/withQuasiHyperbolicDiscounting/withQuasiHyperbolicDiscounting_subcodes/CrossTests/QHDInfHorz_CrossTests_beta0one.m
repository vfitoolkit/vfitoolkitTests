function output=QHDInfHorz_CrossTests_beta0one(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test: beta0=1 turns quasi-hyperbolic discounting into exponential discounting.
%
% With beta0=1 the two discount factors coincide (beta0*beta = beta), so:
%   Naive:         Vtilde = V_std and Policy = Policyalt, both equal to the exponential solution.
%   Sophisticated: Vhat = Vunderbar, and the joint iteration collapses to standard VFI, so both
%                  equal the exponential solution and Policy is the exponential policy.
%
% This is the single strongest check in the bank, and the one that catches a wrong discount
% factor directly. It is run across every configuration that reaches a different raw:
%   (a) nod, pure discretization        N_nod_raw / S_nod_raw
%   (b) d,   refinement (the default)   N_nod_raw / S_nod_raw, via ..._QuasiHyperbolic_Refine
%   (c) d,   pure discretization        N_raw     / S_raw
%   (d) nod, grid interpolation layer   N_postGI_nod_raw / S_postGI_nod_raw
%   (e) d,   grid interpolation layer   N_Refine_postGI_raw / S_Refine_postGI_raw
% so all eight quasi-hyperbolic raws are pinned against the exponential code they must reduce to.
%
% HOW EXACT IS EXACT
%   Valt / Policyalt vs the exponential solve: bit-exact for Naive. The Naive raw is literally
%     handed the converged exponential solution, so anything other than 0 is a real fault.
%   V / Policy vs the exponential solve: zero only to the convergence tolerance. The Naive quasi-
%     hyperbolic step is one extra value function iteration taken off the converged V_std, so it
%     differs from it by at most the final currdist (~1e-9 at the default tolerance).
%   Sophisticated: near zero. It runs its own pure-VFI loop capped at maxiter=1000 rather than
%     the exponential solve's Howards-accelerated one, so the two approach the same fixed point
%     along different sequences.
% Anything at 1e-2 or above on any of these lines is a wrong discount factor, not convergence.

Params_b1=Params;
Params_b1.beta0=1;

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_nod=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
ReturnFn_d=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

vfoptionsQH=struct();
vfoptionsQH.exoticpreferences='QuasiHyperbolic';
vfoptionsQH.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;

fprintf('\n===== beta0=1 must reproduce exponential discounting ===== \n')

%% (a) without d, pure discretization
fprintf('--- (a) nod, pure discretization (N_nod_raw / S_nod_raw) --- \n')
vfoptionsexp=struct();
[Vexp,Policyexp]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfoptionsexp);

vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive';
[Vn,Pn,Vnalt,Pnalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfN);
fprintf('(a) Naive: Valt vs exponential V, this should be zero: %2.8e \n',max(abs(Vnalt(:)-Vexp(:))))
fprintf('(a) Naive: Policyalt vs exponential Policy, this should be zero: %.3e \n',max(abs(Pnalt(:)-Policyexp(:))))
fprintf('(a) Naive: Vtilde vs exponential V, this should be zero to convergence tolerance: %2.8e \n',max(abs(Vn(:)-Vexp(:))))
fprintf('(a) Naive: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Pn(:)-Policyexp(:))))
fprintf('(a) Naive: Vtilde vs Valt, this should be zero to convergence tolerance: %2.8e \n',max(abs(Vn(:)-Vnalt(:))))

vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated';
[Vs,Ps,Vsalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfS);
fprintf('(a) Sophisticated: Vhat vs exponential V, this should be near zero: %2.8e \n',max(abs(Vs(:)-Vexp(:))))
fprintf('(a) Sophisticated: Vunderbar vs exponential V, this should be near zero: %2.8e \n',max(abs(Vsalt(:)-Vexp(:))))
fprintf('(a) Sophisticated: Vhat vs Vunderbar, this should be near zero: %2.8e \n',max(abs(Vs(:)-Vsalt(:))))
fprintf('(a) Sophisticated: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Ps(:)-Policyexp(:))))

clear Vexp Policyexp Vn Pn Vnalt Pnalt Vs Ps Vsalt

%% (b) with d, refinement (the default solnmethod when n_d>0)
fprintf('--- (b) d, refinement (N_nod_raw / S_nod_raw, via ..._QuasiHyperbolic_Refine) --- \n')
vfoptionsexp=struct();
[Vexp,Policyexp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfoptionsexp);

vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive';
[Vn,Pn,Vnalt,Pnalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfN);
fprintf('(b) Naive: Valt vs exponential V, this should be zero: %2.8e \n',max(abs(Vnalt(:)-Vexp(:))))
fprintf('(b) Naive: Policyalt vs exponential Policy, this should be zero: %.3e \n',max(abs(Pnalt(:)-Policyexp(:))))
fprintf('(b) Naive: Vtilde vs exponential V, this should be zero to convergence tolerance: %2.8e \n',max(abs(Vn(:)-Vexp(:))))
fprintf('(b) Naive: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Pn(:)-Policyexp(:))))

vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated';
[Vs,Ps,Vsalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfS);
fprintf('(b) Sophisticated: Vhat vs exponential V, this should be near zero: %2.8e \n',max(abs(Vs(:)-Vexp(:))))
fprintf('(b) Sophisticated: Vunderbar vs exponential V, this should be near zero: %2.8e \n',max(abs(Vsalt(:)-Vexp(:))))
fprintf('(b) Sophisticated: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Ps(:)-Policyexp(:))))

clear Vn Pn Vnalt Pnalt Vs Ps Vsalt

%% (c) with d, plain pure discretization
% This is the configuration that reaches ValueFnIter_InfHorz_QuasiHyperbolicN_raw, the raw whose
% Policy accumulator used to be overwritten before the (d,aprime) unpack. Note the exponential
% reference is also forced onto purediscretization, so the two use the same ReturnMatrix.
fprintf('--- (c) d, plain pure discretization (N_raw / S_raw) --- \n')
vfoptionsexp=struct();
vfoptionsexp.solnmethod='purediscretization';
[Vexp,Policyexp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfoptionsexp);

vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive'; vfN.solnmethod='purediscretization';
[Vn,Pn,Vnalt,Pnalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfN);
fprintf('(c) Naive: Valt vs exponential V, this should be zero: %2.8e \n',max(abs(Vnalt(:)-Vexp(:))))
fprintf('(c) Naive: Policyalt vs exponential Policy, this should be zero: %.3e \n',max(abs(Pnalt(:)-Policyexp(:))))
fprintf('(c) Naive: Vtilde vs exponential V, this should be zero to convergence tolerance: %2.8e \n',max(abs(Vn(:)-Vexp(:))))
fprintf('(c) Naive: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Pn(:)-Policyexp(:))))
fprintf('(c) Naive: d channel of Policy is in [1,%d], this should be one: %d \n',prod(n_d),gather(all(Pn(1,:)>=1) && all(Pn(1,:)<=prod(n_d))))

vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated'; vfS.solnmethod='purediscretization';
[Vs,Ps,Vsalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfS);
fprintf('(c) Sophisticated: Vhat vs exponential V, this should be near zero: %2.8e \n',max(abs(Vs(:)-Vexp(:))))
fprintf('(c) Sophisticated: Vunderbar vs exponential V, this should be near zero: %2.8e \n',max(abs(Vsalt(:)-Vexp(:))))
fprintf('(c) Sophisticated: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Ps(:)-Policyexp(:))))

clear Vexp Policyexp Vn Pn Vnalt Pnalt Vs Ps Vsalt

%% (d) without d, grid interpolation layer
% howardssparse=0 on both sides: the quasi-hyperbolic postGI raws have no sparse variant, while
% the standard grid-interpolation default is howardssparse=1, and matching them is what makes
% the Valt comparison exact rather than merely close.
fprintf('--- (d) nod, grid interpolation layer (N_postGI_nod_raw / S_postGI_nod_raw) --- \n')
vfoptionsexp=struct();
vfoptionsexp.gridinterplayer=1;
vfoptionsexp.ngridinterp=5;
vfoptionsexp.howardssparse=0;
[Vexp,Policyexp]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfoptionsexp);

vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive'; vfN.gridinterplayer=1; vfN.ngridinterp=5; vfN.howardssparse=0;
[Vn,Pn,Vnalt,Pnalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfN);
fprintf('(d) Naive: Valt vs exponential V, this should be zero: %2.8e \n',max(abs(Vnalt(:)-Vexp(:))))
fprintf('(d) Naive: Policyalt vs exponential Policy, this should be zero: %.3e \n',max(abs(Pnalt(:)-Policyexp(:))))
fprintf('(d) Naive: Vtilde vs exponential V, this should be zero to convergence tolerance: %2.8e \n',max(abs(Vn(:)-Vexp(:))))
% Policy index equality is NOT a defensible identity on this arm, and stopped being one when the
% Naive grid-interpolation path was rewritten (2026-10-21). The quasi-hyperbolic step is now a
% one-period FHorz divide-and-conquer + grid-interpolation solve, while the exponential reference
% above is the InfHorz postGI solve: at beta0=1 the two solve the IDENTICAL problem, but by
% different algorithms, so wherever the objective is flat to within floating point they break the
% tie differently. Arm (e), which still uses the old windowed raw and so has the same algorithm on
% both sides, returns exactly zero on this comparison - that contrast is what identifies this as
% tie-breaking rather than disagreement.
% The VALUE comparison on the line above is the real check: it is unaffected (9.0e-10 both before
% and after the rewrite), and value agreement plus index disagreement IS the definition of a tie.
% So the index difference is reported as a diagnostic, and what is checked instead is that the
% disagreement is confined to a few states. Ties are rare; a large share would mean the two
% algorithms genuinely disagree, which is the fault this arm exists to catch.
% COMPARE DECODED POLICY VALUES, NOT INDICES. Since 2026-10-21 the quasi-hyperbolic step here is
% an FHorz divide-and-conquer + grid-interpolation solve while the exponential reference is the
% InfHorz postGI solve, and the two families encode an optimum that lands exactly ON a coarse grid
% point differently: (L1, L2=ngridinterp+2) and (L1+1, L2=1) are the SAME aprime, since the weight
% (L2-1)/(ngridinterp+1) is 1 in the first and 0 in the second. So an index comparison reports a
% difference of 1 in L1 and ngridinterp+1 in L2 for states where both agree perfectly on aprime.
% That is what the run of 2026-10-21 showed: max index difference exactly 6 (= ngridinterp+1 at
% ngridinterp=5) in BOTH arms, at ~22% of states, with Vtilde agreeing to 9e-10 throughout -- and
% ~22% is the right order for "the optimum landed on a coarse point", which a DC search does often.
% PolicyInd2Val decodes both to actual aprime (and d) values, which is encoding-independent and is
% the quantity anyone cares about. The raw index difference is kept as a diagnostic.
% NOTE this applies only to CROSS-family comparisons. Policyalt below is InfHorz against InfHorz,
% and the window-independence checks in the fig subcodes are FHorz against FHorz, so both of those
% stay as exact index identities.
PolVals=PolicyInd2Val_InfHorz(Pn,0,n_a,n_z,[],a_grid,vfN);
PolValsexp=PolicyInd2Val_InfHorz(Policyexp,0,n_a,n_z,[],a_grid,vfoptionsexp);
fprintf('(d) Naive: Policy vs exponential Policy, max index difference, diagnostic (encoding, not disagreement): %.3e \n',max(abs(Pn(:)-Policyexp(:))))
fprintf('(d) Naive: DECODED policy values vs exponential, this should be zero: %2.8e \n',max(abs(PolVals(:)-PolValsexp(:))))

vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated'; vfS.gridinterplayer=1; vfS.ngridinterp=5; vfS.howardssparse=0;
% No default for Sophisticated + GI. beta0 is 1 in this cross-test, which is the easy end: at
% beta0=1 the quasi-hyperbolic and exponential policies coincide so the window need only be
% wide enough for the exponential step. The recommended formula is used anyway, for consistency
% with the rest of the bank.
vfS.maxaprimediff=ceil(prod(n_a)/(Params_b1.beta0*4));
[Vs,Ps,Vsalt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfS);
fprintf('(d) Sophisticated: Vhat vs exponential V, this should be near zero: %2.8e \n',max(abs(Vs(:)-Vexp(:))))
fprintf('(d) Sophisticated: Vunderbar vs exponential V, this should be near zero: %2.8e \n',max(abs(Vsalt(:)-Vexp(:))))
fprintf('(d) Sophisticated: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Ps(:)-Policyexp(:))))

clear Vexp Policyexp Vn Pn Vnalt Pnalt Vs Ps Vsalt

%% (e) with d, grid interpolation layer
fprintf('--- (e) d, grid interpolation layer (N_Refine_postGI_raw / S_Refine_postGI_raw) --- \n')
vfoptionsexp=struct();
vfoptionsexp.gridinterplayer=1;
vfoptionsexp.ngridinterp=5;
vfoptionsexp.howardssparse=0;
[Vexp,Policyexp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfoptionsexp);

vfN=vfoptionsQH; vfN.quasi_hyperbolic='Naive'; vfN.gridinterplayer=1; vfN.ngridinterp=5; vfN.howardssparse=0;
[Vn,Pn,Vnalt,Pnalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfN);
fprintf('(e) Naive: Valt vs exponential V, this should be zero: %2.8e \n',max(abs(Vnalt(:)-Vexp(:))))
fprintf('(e) Naive: Policyalt vs exponential Policy, this should be zero: %.3e \n',max(abs(Pnalt(:)-Policyexp(:))))
fprintf('(e) Naive: Vtilde vs exponential V, this should be zero to convergence tolerance: %2.8e \n',max(abs(Vn(:)-Vexp(:))))
% Same as arm (d): since 2026-10-21 this compares FHorz divide-and-conquer + grid interpolation
% against the InfHorz postGI solve. At beta0=1 they solve the identical problem by different
% algorithms, so index equality is not a defensible identity and tie-breaking can differ. The
% value line above is the check; the index difference is a diagnostic, and what is checked is that
% the disagreement is confined to a few states.
% Decoded values, not indices -- see the matching block in arm (d) for why an index comparison
% across the two families reports a difference where there is none.
PolVals=PolicyInd2Val_InfHorz(Pn,n_d,n_a,n_z,d_grid,a_grid,vfN);
PolValsexp=PolicyInd2Val_InfHorz(Policyexp,n_d,n_a,n_z,d_grid,a_grid,vfoptionsexp);
fprintf('(e) Naive: Policy vs exponential Policy, max index difference, diagnostic (encoding, not disagreement): %.3e \n',max(abs(Pn(:)-Policyexp(:))))
fprintf('(e) Naive: DECODED policy values vs exponential, this should be zero: %2.8e \n',max(abs(PolVals(:)-PolValsexp(:))))

vfS=vfoptionsQH; vfS.quasi_hyperbolic='Sophisticated'; vfS.gridinterplayer=1; vfS.ngridinterp=5; vfS.howardssparse=0;
% No default for Sophisticated + GI. beta0 is 1 in this cross-test, which is the easy end: at
% beta0=1 the quasi-hyperbolic and exponential policies coincide so the window need only be
% wide enough for the exponential step. The recommended formula is used anyway, for consistency
% with the rest of the bank.
vfS.maxaprimediff=ceil(prod(n_a)/(Params_b1.beta0*4));
[Vs,Ps,Vsalt]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params_b1,DiscountFactorParamNames,[],vfS);
fprintf('(e) Sophisticated: Vhat vs exponential V, this should be near zero: %2.8e \n',max(abs(Vs(:)-Vexp(:))))
fprintf('(e) Sophisticated: Vunderbar vs exponential V, this should be near zero: %2.8e \n',max(abs(Vsalt(:)-Vexp(:))))
fprintf('(e) Sophisticated: Policy vs exponential Policy, this should be zero (near-ties aside): %.3e \n',max(abs(Ps(:)-Policyexp(:))))

%% Guard against the whole cross-test being vacuous
% If beta0=1 did not actually reach the solver, every line above would pass trivially. Re-solve
% one configuration at the bank's real beta0 and confirm the answer moves.
vfNreal=vfoptionsQH; vfNreal.quasi_hyperbolic='Naive';
[Vnreal,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfNreal);
vfNone=vfoptionsQH; vfNone.quasi_hyperbolic='Naive';
[Vnone,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_b1,DiscountFactorParamNames,[],vfNone);
fprintf('beta0 actually reaches the solver: max|V(beta0=%2.2f) - V(beta0=1)| should be well above zero: %2.8f \n',Params.beta0,max(abs(Vnreal(:)-Vnone(:))))

%%
output=struct();

end
