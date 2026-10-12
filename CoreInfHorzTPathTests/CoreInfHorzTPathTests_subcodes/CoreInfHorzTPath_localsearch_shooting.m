function output=CoreInfHorzTPath_localsearch_shooting(T,n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Checks vfoptions.localsearch=1 on a REAL general eqm transition path, solved by shooting.
% doPart(4) and doPart(5) test the restricted value fn step against a given price path. This is the
% first part where the restriction is inside a general eqm solve and so can affect the answer.
%
% THE SCHEME. Iteration 1 is a standard sweep, whatever divideandconquer and gridinterplayer are set
% to, and its Policy becomes the reference. Later iterations restrict aprime to a window around the
% previous iteration's answer, refreshing the reference as they go. When those converge, ONE standard
% sweep re-evaluates the SAME price path: still converged means done, and moved means its Policy
% becomes the new reference and local search resumes.
% WHAT THE VERIFICATION SWEEP DOES AND DOES NOT CERTIFY. It makes the general eqm conditions hold
% under an UNRESTRICTED value function solve, so the prices and aggregates are certified to
% toleranceGEcondns. It does NOT certify that the policy equals the unrestricted policy state by
% state, because it is accepted on PricePathDist and GEcondnPathDist and never on whether the policy
% moved. Two gaps follow, and neither is closed here:
%   - The reference iteration can have a spurious fixed point. Where the objective in aprime is
%     multi-peaked -- which is what makes a policy function non-monotone -- the restricted sweep can
%     sit on a local peak while the global one lies outside the window, and the new reference is that
%     local peak, so the next restricted sweep sits there too. An unrestricted sweep would move. If
%     the states affected are few, or carry little mass, the aggregate residual moves by less than
%     toleranceGEcondns and the solve terminates with a policy that is wrong at those states.
%   - With multiple equilibria the restriction can steer the price trajectory into a different basin,
%     and both endpoints pass verification.
% So the restriction does not merely change how fast the answer is reached. Neither gap can occur
% with a single-peaked objective in aprime, which covers concave problems on convex feasible sets,
% and the checks below are run on such a model. On a model with non-convexities, fixed costs, or
% sophisticated quasi-hyperbolic discounting, local search is not a free speed-up.
%
% SHOOTING ONLY, and the refusal that enforces it is checked below. Anderson hands a PURE residual
% function to AndersonAcceleration, which is shared with the stationary general eqm code, so the
% reference cannot be carried through it without changing that shared contract; quasi-Newton perturbs
% prices to build a Jacobian, where the reference has to be the baseline path's and not the previous
% perturbation's. Were they allowed through, they would run every iteration against the stay-put
% default reference and never take a verification sweep, and return a wrong path in silence.
%
% NLOCALSEARCH IS A STARTING WINDOW, NOT THE WINDOW. The shooting loop ratchets it up by
% vfoptions.nlocalsearchup (default 2) whenever the answer sat on a window edge anywhere, and down by
% nlocalsearchdown (default 1) otherwise, resetting to the input value after every standard sweep.
% Edges at the GRID's own ends do not count: those are corner solutions, and counting them would
% ratchet the window to its cap forever on any model with a top-of-grid saver or a constrained
% household -- which is to say, on all of them.
% So the three starts below (1, 3 and 25) test that the controller reaches a workable window from
% BOTH directions: 1 has to ratchet up, 25 has to ratchet down. The trajectory each one took is in
% the Local search report lines, and is the real diagnostic -- it is the window this application
% needed, found by the solver instead of guessed.
%
% WHAT THIS PART DOES NOT MEASURE: wall clock. At n_a_GE=201 the sweep is about 2e5 elements, far
% below the kernel launch floor, so every solve here is launch bound and a timing would say nothing.
% The speed of the restricted sweep is measured in doPart(5) on the big grid. What this part is for
% is whether the scheme CONVERGES TO THE SAME PLACE, how many verification sweeps it costs, and the
% window trajectory the ratchet settles on. Check (7) adds the measurements that bear on whether the
% Anderson algorithm could be given local search at all; the edge-hit RATE is not that number, since
% the ratchet drives it towards zero by construction and so it says nothing about accuracy.

vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
n_d=0; d_grid=[]; % the no-d half of this part; the with-d half is sections (9) to (11) below

gpudev=gpuDevice;

%% The model, taken verbatim from the general eqm block of CoreInfHorzTPath_nod_z_noe_nosemiz
% firm* names because the 2A entrepreneur ReturnFn already uses alpha and delta with its own meanings
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(aprime,a,z) a;
FnsToEvaluateGE.N=@(aprime,a,z) z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct();
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a_GE, n_z, 0, pi_z, d_grid, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
fprintf('localsearch shooting: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)

ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d,n_a_GE,n_z,d_grid,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d,n_a_GE,n_z,pi_z,simoptions,ParamsGE,[]);

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing changes, this is a null reform

% A bumped start, so the solver has to travel. Both prices bumped by 1 to 5 per cent, by a different
% amount in each period, over the first two thirds of the path. A path that is already at the answer
% converges in one iteration and would never exercise the scheme at all.
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+bumpr);
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+bumpw);

transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.GEnewprice=3; % shooting
transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                           'LabourMarket','w',0,0.1};
transpathoptionsS.maxiter=1000; % shooting at a damping factor of 0.1 needs far more iterations than Anderson
transpathoptionsS.verbose=0;
% Set explicitly rather than left to the defaults inside TransitionPath_InfHorz, because the checks
% below print them as the scale to read the price comparisons against, and a bare transpathoptions
% does not carry them.
% toleranceGEcondns is a MULTIPLE OF THE MEASURED STATIONARY RESIDUAL, not a fixed number. The path
% cannot satisfy the general eqm conditions any more tightly than the stationary problem does at the
% same prices, and the default fminalgo leaves CapitalMarket about 6e-4 from zero on this calibration,
% so a fixed 1e-4 is unreachable and the shooting solve would run to maxiter and return nothing
% useful. Tying it to the floor the stationary solve actually achieved means a run on a different
% calibration does not fail for a reason that has nothing to do with local search.
% toleranceGEprices stays fixed: the shooting update is factor*residual, so the per-iteration price
% change is about a tenth of the residual and 1e-4 is comfortably reachable.
GEfloor=max(abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket));
transpathoptionsS.toleranceGEprices=10^(-4);
transpathoptionsS.toleranceGEcondns=5*GEfloor;
fprintf('localsearch shooting: the stationary solve leaves a residual of %.3e, so toleranceGEcondns for the path is set to 5 times that, %.3e \n',GEfloor,transpathoptionsS.toleranceGEcondns)

%% (A) standard, (B) local search nlocalsearch=1, (C) local search nlocalsearch=3
tic;
[PricePathA,GEcondnPathA]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptions, []);
wait(gpudev)
tA=toc;

vfoptionsB=vfoptions;
vfoptionsB.localsearch=1;
vfoptionsB.nlocalsearch=1;
tic;
[PricePathB,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsB, []);
wait(gpudev)
tB=toc;

vfoptionsC=vfoptions;
vfoptionsC.localsearch=1;
vfoptionsC.nlocalsearch=3;
tic;
[PricePathC,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsC, []);
wait(gpudev)
tC=toc;

fprintf('localsearch shooting: standard solve %6.2fs, nlocalsearch=1 %6.2fs, nlocalsearch=3 %6.2fs (launch bound at this grid size, so read the iteration counts above, not these) \n',tA,tB,tC)
fprintf('localsearch shooting: max deviation of the standard solve from the stationary eqm, r %.3e, w %.3e, against a stationary CapitalMarket residual of %.3e \n',max(abs(gather(PricePathA.r)-p_eqm.r)),max(abs(gather(PricePathA.w)-p_eqm.w)),abs(GEcondns.CapitalMarket))

%% (1) Shooting still works: the standard shooting solve agrees with Anderson on the same path
% The regression check on the loop surgery. Every new line in the shooting loop is gated on
% vfoptions.localsearch, so a standard solve should be untouched, and Anderson is an independent
% route to the same equilibrium path.
transpathoptionsAnd=transpathoptionsbaseline;
transpathoptionsAnd.GEnewprice=2;
transpathoptionsAnd.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
transpathoptionsAnd.maxiter=200;
transpathoptionsAnd.verbose=0;
[PricePathAnd,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAnd, simoptions, vfoptions, []);
fprintf('localsearch shooting: (1) standard shooting vs Anderson, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathA.r)-gather(PricePathAnd.r))))
fprintf('localsearch shooting: (1) standard shooting vs Anderson, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathA.w)-gather(PricePathAnd.w))))

%% (2)-(3) The local search solves land where the standard one lands
% Not an identity, and deliberately not held to one: the two take different iteration sequences and
% can stop at different points inside the same tolerance ball. The scale to read these against is
% transpathoptions.toleranceGEprices, printed beside them.
fprintf('localsearch shooting: (2) nlocalsearch starting at 1 vs the standard solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathB.r)-gather(PricePathA.r))))
fprintf('localsearch shooting: (2) nlocalsearch starting at 1 vs the standard solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathB.w)-gather(PricePathA.w))))
fprintf('localsearch shooting: (3) nlocalsearch starting at 3 vs the standard solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathC.r)-gather(PricePathA.r))))
fprintf('localsearch shooting: (3) nlocalsearch starting at 3 vs the standard solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathC.w)-gather(PricePathA.w))))

%% (3b) A deliberately TOO WIDE start, which the ratchet has to walk down
% The complement of the nlocalsearch=1 start above. 25 is far wider than this model needs, so the
% window should ratchet down towards whatever it actually requires, and the answer must be unaffected
% either way. Together the two starts say the controller recovers from a bad guess in both directions,
% which is what makes nlocalsearch a starting value a user does not have to get right.
vfoptionsD=vfoptions;
vfoptionsD.localsearch=1;
vfoptionsD.nlocalsearch=25;
[PricePathD,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsD, []);
fprintf('localsearch shooting: (3b) nlocalsearch starting at 25 vs the standard solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathD.r)-gather(PricePathA.r))))
fprintf('localsearch shooting: (3b) nlocalsearch starting at 25 vs the standard solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathD.w)-gather(PricePathA.w))))
fprintf('localsearch shooting: for scale, toleranceGEprices is %.3e and toleranceGEcondns is %.3e \n',transpathoptionsS.toleranceGEprices,transpathoptionsS.toleranceGEcondns)

%% (4) The answer is certified from outside
% Re-evaluate the local search answer with an EXACT sweep and check the general eqm conditions. This
% does not take the solver's word for it: factor=0 and maxiter=1 make the call a pure evaluation of
% the path handed in, with vfoptions.localsearch off, so what is tested is whether the path the
% restricted solve returned satisfies the conditions under a full unrestricted value fn solve.
transpathoptionsV=transpathoptionsbaseline;
transpathoptionsV.GEnewprice=3;
transpathoptionsV.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0; ...
                                           'LabourMarket','w',0,0};
transpathoptionsV.maxiter=1;
transpathoptionsV.verbose=0;
[~,GEcondnPathV]=TransitionPath_InfHorz(PricePathB, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsV, simoptions, vfoptions, []);
GEdistV=max(max(abs(gather(GEcondnPathV.CapitalMarket))),max(abs(gather(GEcondnPathV.LabourMarket))));
fprintf('localsearch shooting: (4) the nlocalsearch=1 answer re-evaluated with an exact sweep, its general eqm conditions are within tolerance, should be one: %i \n',GEdistV<=transpathoptionsS.toleranceGEcondns)
fprintf('localsearch shooting: (4) that distance is %.3e against toleranceGEcondns %.3e \n',GEdistV,transpathoptionsS.toleranceGEcondns)

%% (5) The refusals that keep the scheme honest
err=zeros(1,3);
% QUASI-NEWTON is the one still refused. It perturbs prices to build a Jacobian, and the reference
% for a perturbed path has to be the BASELINE path's and not the previous perturbation's, or every
% column of the Jacobian is computed against a different restriction. Anderson USED to be refused
% here too and no longer is -- see section (12) -- so this check moved rather than being deleted.
transpathoptionsQN=transpathoptionsbaseline;
transpathoptionsQN.GEnewprice=1;
transpathoptionsQN.maxiter=2;
transpathoptionsQN.verbose=0;
vfoptionsE=vfoptionsB;
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsQN, simoptions, vfoptionsE, []);
catch ME
    err(1)=contains(ME.message,'GEnewprice=3 (shooting) or 2 (Anderson)');
end
vfoptionsE=vfoptionsB; vfoptionsE.n_e=3; vfoptionsE.e_grid=linspace(0.9,1.1,3)'; vfoptionsE.pi_e=ones(3,1)/3;
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsE, []);
catch ME
    err(2)=contains(ME.message,'not implemented with an e variable');
end
vfoptionsE=vfoptionsB; vfoptionsE.experienceasset=1;
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsE, []);
catch ME
    err(3)=contains(ME.message,'not implemented with an experience asset');
end
fprintf('localsearch shooting: (5) local search with quasi-Newton is refused by name, should be one: %i \n',err(1))
fprintf('localsearch shooting: (5) local search with an e variable is refused by name, should be one: %i \n',err(2))
fprintf('localsearch shooting: (5) local search with an experience asset is refused by name, should be one: %i \n',err(3))

%% (6) The ratchet options are validated
err6=0;
vfoptionsE=vfoptionsB; vfoptionsE.nlocalsearchup=0;
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsE, []);
catch ME
    err6=contains(ME.message,'must both be at least 1');
end
fprintf('localsearch shooting: (6) a ratchet step below one is refused by name, should be one: %i \n',err6)

%% (7) THE GRID INTERPOLATION LAYER, which is the GE scheme's untested half
% doPart(4) and doPart(5) cover the restricted value fn STEP with GI. What they cannot cover is the
% two GI paths that only exist inside the general eqm scheme, and that nothing had ever run before
% this block:
%   - Step1 converting a STANDARD sweep's Policy into a coarse reference, which under GI means
%     compressing (lower coarse point, L2) to the nearest coarse point.
%   - The GI branch of the EDGE TEST that drives the ratchet. Choosing a fine point strictly inside
%     the window's end coarse points must NOT ratchet; only the two extreme coarse points count.
% Both are statically sound, and both are the same kind of code that has twice been statically sound
% and wrong in this programme -- the linidx_upper overrun, and the GI refusal that outlived the
% feature it refused. So they get run.
%
% The comparison is GI-LOCALSEARCH AGAINST GI-STANDARD, not against the no-GI standard solve of (A).
% GI changes the discretisation and so moves the answer by a little on its own; comparing across that
% would mix the discretisation difference into a check that is supposed to isolate the restriction.
% Terminal V and initial distribution are likewise recomputed with GI, so the pair solves one problem.
ngi=5; % the value the rest of this bank uses for ngridinterp
vfoptionsGI=vfoptions;
vfoptionsGI.gridinterplayer=1;
vfoptionsGI.ngridinterp=ngi;
simoptionsGI=simoptions;
simoptionsGI.gridinterplayer=vfoptionsGI.gridinterplayer;
simoptionsGI.ngridinterp=vfoptionsGI.ngridinterp;
[V_finalGE_GI,Policy_finalGE_GI]=ValueFnIter_InfHorz(n_d,n_a_GE,n_z,d_grid,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptionsGI);
AgentDist_initialGE_GI=StationaryDist_InfHorz(Policy_finalGE_GI,n_d,n_a_GE,n_z,pi_z,simoptionsGI,ParamsGE,[]);

[PricePathGIstd,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptionsGI, vfoptionsGI, []);

vfoptionsGIls=vfoptionsGI;
vfoptionsGIls.localsearch=1;
vfoptionsGIls.nlocalsearch=3;
[PricePathGIls,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptionsGI, vfoptionsGIls, []);

fprintf('localsearch shooting: (7) GI, nlocalsearch starting at 3 vs the standard GI solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathGIls.r)-gather(PricePathGIstd.r))))
fprintf('localsearch shooting: (7) GI, nlocalsearch starting at 3 vs the standard GI solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathGIls.w)-gather(PricePathGIstd.w))))

%% (8) The Anderson diagnostic, measured on the nlocalsearch=3 solve
% transpathoptions.localsearchdiagnostic=1 turns on two measurements inside the shooting loop, both
% aimed at ONE question: could Anderson be given local search, and if so with what reference.
%   - CUMULATIVE movement. How far the policy ends up from the reference the FIRST standard sweep
%     produced, against the per-iteration movement. A FROZEN reference is what Anderson could have
%     without forking AndersonAcceleration, because a pure residual closure can capture one by value.
%     It has to span the cumulative figure, so if that is far larger than the per-iteration one then
%     a frozen reference needs a window wide enough to give most of the speed-up back.
%   - RESIDUAL CONTAMINATION, per iteration: the restricted sweep's general eqm residual against the
%     residual an EXACT sweep gives at the SAME prices, divided by the genuine change in residual
%     from the previous iteration. That is precisely the signal-to-noise of what Anderson would store
%     in DeltaF. Small means a per-iteration reference is safe to hand it; of order one means the
%     differences Anderson extrapolates from are mostly the error the restriction introduced.
% The ratio is taken per iteration and NOT as one maximum over another, because the genuine change is
% largest early while contamination matters most late, so two separately-maximised numbers would be
% drawn from different iterations and read far too favourably.
%
% The measurement costs a second full iteration every restricted sweep, which is why it is off by
% default and why this runs on one solve rather than all of them.
transpathoptionsDiag=transpathoptionsS;
transpathoptionsDiag.localsearchdiagnostic=1;
[PricePathDiag,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsDiag, simoptions, vfoptionsC, []);

% The diagnostic must be PURE OBSERVATION. Its extra sweep takes only the first output, so nothing it
% computes is written back to the reference path, the Policy path or the loop's own state, and the
% solve it watches must therefore follow a bit-identical trajectory to the same solve without it.
% Held as an IDENTITY, not a tolerance: anything other than exact agreement means the measurement is
% changing what it measures, which would make every number above describe a different solve.
fprintf('localsearch shooting: (8) the diagnostic is pure observation, its solve vs the same solve without it, r, this should be zero: %2.8f \n',max(abs(gather(PricePathDiag.r)-gather(PricePathC.r))))
fprintf('localsearch shooting: (8) the diagnostic is pure observation, its solve vs the same solve without it, w, this should be zero: %2.8f \n',max(abs(gather(PricePathDiag.w)-gather(PricePathC.w))))

%% (9)-(11) WITH A d, which needs a SEED SWEEP and is the point of this half of the part
% A standard sweep returns one (dstar,aprimestar) per state, so it pins the reference only for the d
% it chose and says nothing about the other N_d-1 values. The next restricted sweep searches a window
% PER d, so it needs a centre per d. That gap is what kept local search out of a general eqm path
% with a d, and vfoptions.localsearchseed is what closes it: on iteration 1 and on every verification
% sweep the dispatcher runs an unrestricted sweep that ALSO returns the coarse aprime argmax per
% (d,a,z). Without the interpolation layer that comes from a dedicated raw, LSseed_raw, which
% unfolds (d,aprime) and maxes over aprime first; with it, the standard GI raw already computes
% exactly that argmax as its midpoint and now returns it.
%
% Seed sweeps HONOUR divideandconquer. DC is a memory option as much as a speed one here: the seed
% sweep is the unrestricted one, about 2 GB with a d at n_a=1001 against about 63 MB under DC, so a
% model that needs DC to solve at all can still use local search. The restricted sweeps still refuse
% DC, and still ignore lowmemory.
ReturnFnd=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
FnsToEvaluateGEd.K=@(d,aprime,a,z) a;
FnsToEvaluateGEd.N=@(d,aprime,a,z) d*z;

% THE WITH-d STATIONARY PROBLEM HAS TO BE SOLVED WITH GI FIRST, AND THE ANSWER HANDED TO THE NO-GI
% SOLVE AS ITS STARTING GUESS. Solved the obvious way -- no interpolation layer, default fminalgo,
% default start -- it does not reach an equilibrium at all: it STALLS at r=0.05309198, w=0.68065688
% with a CapitalMarket residual of 5.553e-02. K=E[a] and N=E[d*z] are step functions of the prices on
% coarse grids, so the general eqm objective is piecewise-constant and fminsearch stops where it
% stops. The setup file records the same mechanism and is why n_d_GE is 51 rather than 9, but 51
% points only cut the residual (from about 170x the no-d model's to 94x) without stopping the stall.
% WITH THE INTERPOLATION LAYER the chosen aprime moves continuously in the prices, so K is smooth and
% the search finds r=0.05219468, w=0.43766914 at a residual of -9.983e-04. That is a factor of 0.02,
% and the w is essentially the no-d model's 0.43777959 -- which the same firm technology should give.
% So 0.68 was never an approximate equilibrium, it was a stall, and everything in this block used to
% be built on it. That is why its tolerance had to be 2.8e-01 and why its agreements came out at
% 1e-03 where the no-d ones are 1e-06: transition paths around a point that is not a steady state
% are not wrong exactly, but they test much less.
% WHAT TO READ: whether the no-GI residual is now small. If it is, the stall was a starting-point
% problem and the no-GI half of this block is on its own discretisation's equilibrium. If it is still
% about 5e-02, the no-GI discretisation genuinely cannot do better and the tolerance derived from it
% below is correct for it rather than merely loose.
vfoptionsGEd_GI=vfoptions;
vfoptionsGEd_GI.gridinterplayer=1;
vfoptionsGEd_GI.ngridinterp=ngi;
simoptionsGEd_GI=simoptions;
simoptionsGEd_GI.gridinterplayer=1;
simoptionsGEd_GI.ngridinterp=ngi;
[p_eqmd_GI,GEcondnsd_GI]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptionsGEd_GI, vfoptionsGEd_GI);
fprintf('localsearch shooting d: stationary GE WITH GI, r=%2.8f, w=%2.8f, CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqmd_GI.r,p_eqmd_GI.w,GEcondnsd_GI.CapitalMarket,GEcondnsd_GI.LabourMarket)

% The no-GI solve, started from the GI answer instead of the default. HeteroAgentStationaryEqm_InfHorz
% takes its initial guess from Parameters.r and Parameters.w.
ParamsStartd=Params;
ParamsStartd.r=p_eqmd_GI.r;
ParamsStartd.w=p_eqmd_GI.w;
[p_eqmd,GEcondnsd]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsStartd, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions, vfoptions);
fprintf('localsearch shooting d: stationary GE no GI, started from the GI answer, r=%2.8f, w=%2.8f, CapitalMarket=%.3e LabourMarket=%.3e \n',p_eqmd.r,p_eqmd.w,GEcondnsd.CapitalMarket,GEcondnsd.LabourMarket)
fprintf('localsearch shooting d: the same no-GI solve from the DEFAULT start stalled at w=0.68065688 with CapitalMarket 5.553e-02, which is the measurement this restructure rests on \n')

ParamsGEd=Params;
ParamsGEd.r=p_eqmd.r;
ParamsGEd.w=p_eqmd.w;
[V_finalGEd,Policy_finalGEd]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFnd,ParamsGEd,DiscountFactorParamNames,[],vfoptions);
AgentDist_initialGEd=StationaryDist_InfHorz(Policy_finalGEd,n_d_GE,n_a_GE,n_z,pi_z,simoptions,ParamsGEd,[]);

PricePathBumpedd.r=p_eqmd.r*ones(1,T);
PricePathBumpedd.w=p_eqmd.w*ones(1,T);
PricePathBumpedd.r(1:Tbump)=PricePathBumpedd.r(1:Tbump).*(1+bumpr);
PricePathBumpedd.w(1:Tbump)=PricePathBumpedd.w(1:Tbump).*(1+bumpw);

% Tolerance tied to the stationary residual for this calibration, as in the no-d half: a path cannot
% satisfy the general eqm conditions more tightly than the stationary problem does at the same prices.
transpathoptionsSd=transpathoptionsS;
GEfloord=max(abs(GEcondnsd.CapitalMarket),abs(GEcondnsd.LabourMarket));
transpathoptionsSd.toleranceGEcondns=5*GEfloord;
fprintf('localsearch shooting d: the stationary solve leaves a residual of %.3e, so toleranceGEcondns is %.3e \n',GEfloord,transpathoptionsSd.toleranceGEcondns)

%% (9) with d, no interpolation layer
[PricePathAd,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd, AgentDist_initialGEd, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsSd, simoptions, vfoptions, []);
vfoptionsBd=vfoptions;
vfoptionsBd.localsearch=1;
vfoptionsBd.nlocalsearch=3;
[PricePathBd,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd, AgentDist_initialGEd, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsSd, simoptions, vfoptionsBd, []);
fprintf('localsearch shooting d: (9) with d, nlocalsearch starting at 3 vs the standard solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathBd.r)-gather(PricePathAd.r))))
fprintf('localsearch shooting d: (9) with d, nlocalsearch starting at 3 vs the standard solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathBd.w)-gather(PricePathAd.w))))

%% (10) with d AND the interpolation layer, which is the fourth corner of the table
% Compared against the standard GI solve, not the no-GI one: GI moves the answer on its own and
% comparing across that would mix the discretisation difference into a check about the restriction.
vfoptionsGId=vfoptions;
vfoptionsGId.gridinterplayer=1;
vfoptionsGId.ngridinterp=ngi;
simoptionsGId=simoptions;
simoptionsGId.gridinterplayer=1;
simoptionsGId.ngridinterp=ngi;
[V_finalGEd_GI,Policy_finalGEd_GI]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFnd,ParamsGEd,DiscountFactorParamNames,[],vfoptionsGId);
AgentDist_initialGEd_GI=StationaryDist_InfHorz(Policy_finalGEd_GI,n_d_GE,n_a_GE,n_z,pi_z,simoptionsGId,ParamsGEd,[]);
[PricePathGIdstd,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd_GI, AgentDist_initialGEd_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsSd, simoptionsGId, vfoptionsGId, []);
vfoptionsGIdls=vfoptionsGId;
vfoptionsGIdls.localsearch=1;
vfoptionsGIdls.nlocalsearch=3;
[PricePathGIdls,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd_GI, AgentDist_initialGEd_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsSd, simoptionsGId, vfoptionsGIdls, []);
fprintf('localsearch shooting d: (10) with d and GI, nlocalsearch starting at 3 vs the standard GI solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathGIdls.r)-gather(PricePathGIdstd.r))))
fprintf('localsearch shooting d: (10) with d and GI, nlocalsearch starting at 3 vs the standard GI solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathGIdls.w)-gather(PricePathGIdstd.w))))

%% (11) The SEED RAW itself, held to the folded standard raw it replaces
% LSseed_raw searches the same set as ValueFnIter_InfHorz_TPath_SingleStep_raw -- every (d,aprime) --
% and differs only in taking the max over aprime per d and then over d, rather than one max over a
% folded (d,aprime) dimension. So V is an IDENTITY. Policy can differ only at exact ties, where the
% two orders disagree: folded order runs d fastest and so picks the lowest aprime then the lowest d
% within it, while the unfolded version picks the lowest d then that d's lowest aprime. The tie count
% is reported rather than checked, because it is a property of the model and not of the code.
% This is the check on the one piece of genuinely new arithmetic in the with-d general eqm path. The
% solves above exercise it too, but only through the general eqm residual, which would hide a small
% error.
ReturnFnParamNamesd=ReturnFnParamNamesFn(ReturnFnd,n_d_GE,n_a_GE,n_z,0,vfoptions,ParamsGEd);
% The same branch doPart(5) uses: with a single d the dispatcher wants the grid itself, and
% CreateGridvals is only for a joint grid over several d. Passing the wrong one of these is a
% kind-mismatch that no arity check would catch.
if isscalar(n_d_GE)
    d_gridvalsd=gpuArray(d_grid_GE);
else
    d_gridvalsd=CreateGridvals(n_d_GE,gpuArray(d_grid_GE),1);
end
z_gridvalsd=CreateGridvals(n_z,gpuArray(z_grid),1);
VnextSd=gpuArray(reshape(V_finalGEd,[prod(n_a_GE),prod(n_z)]));
vfoptionsSd=vfoptions;
vfoptionsSd.exoticpreferences='None';
vfoptionsSd.lowmemory=0;
vfoptionsSd.gridinterplayer=0;
vfoptionsSd.divideandconquer=0;
vfoptionsSd.localsearch=0;
vfoptionsSd.localsearchseed=0;
[Vfold,Pfold]=ValueFnIter_InfHorz_TPath_SingleStep(VnextSd,n_d_GE,n_a_GE,n_z,d_gridvalsd, gpuArray(a_grid_GE), z_gridvalsd, gpuArray(pi_z), ReturnFnd, ParamsGEd, DiscountFactorParamNames, ReturnFnParamNamesd, [], vfoptionsSd);
vfoptionsSeedd=vfoptionsSd;
vfoptionsSeedd.localsearchseed=1;
[Vseed,Pseed,refseed]=ValueFnIter_InfHorz_TPath_SingleStep(VnextSd,n_d_GE,n_a_GE,n_z,d_gridvalsd, gpuArray(a_grid_GE), z_gridvalsd, gpuArray(pi_z), ReturnFnd, ParamsGEd, DiscountFactorParamNames, ReturnFnParamNamesd, [], vfoptionsSeedd);
fprintf('localsearch shooting d: (11) the seed raw searches the same set as the folded standard raw, V, this should be zero: %.3e \n',max(abs(gather(Vseed)-gather(Vfold)),[],'all'))
fprintf('localsearch shooting d: (11) for information, Policy differs at %i of %i states, which is exact ties broken in a different order \n',sum(gather(Pseed(1,:,:))~=gather(Pfold(1,:,:)),'all'),prod(n_a_GE)*prod(n_z))
% The reference is a coarse aprime index per (d,a,z), and at the CHOSEN d it must be the aprime the
% Policy reports -- the one place the per-d lookup could be wired to the wrong index.
N_dd=prod(n_d_GE); N_ad=prod(n_a_GE); N_zd=prod(n_z);
refseed=gather(refseed);
dstarseed=gather(Pseed(1,:,:));
linseed=dstarseed+N_dd*(0:1:N_ad-1)+N_dd*N_ad*shiftdim((0:1:N_zd-1),-1);
fprintf('localsearch shooting d: (11) the reference read at the chosen d vs the Policy aprime, this should be zero: %.3e \n',max(abs(refseed(linseed)-gather(Pseed(2,:,:))),[],'all'))

%% (12) ANDERSON WITH LOCAL SEARCH, against Anderson without it and against shooting
% Anderson used to be refused with local search because AndersonAcceleration takes a PURE residual
% function and is shared with the stationary general eqm code, so the reference could not be threaded
% through it. AndersonAccelerationTPath is a fork that threads it, and certifies convergence with an
% unrestricted VerifyFn before believing it. The fork changes no line of the algorithm -- gamma, the
% step, the history update, the Tikhonov scaling and the safeguard are identical -- only the contract.
%
% WHY THE VERIFICATION CANNOT BE AN ORDINARY ITERATION, which is the one structural difference from
% shooting. DeltaX and DeltaF hold successive DIFFERENCES, so re-evaluating at the same prices would
% contribute a zero DeltaX against a non-zero DeltaF: a degenerate secant pair, exactly what the
% restart guards exist to discard. So it is a side evaluation outside the history, and a failure
% restarts the history with the reference reseeded from the exact sweep.
%
% THE COMPARISON THAT MATTERS IS AGAINST ANDERSON WITHOUT LOCAL SEARCH, not against shooting: it
% isolates what the restriction does to the accelerator. Against shooting is also reported, because
% what a user wants is the fastest configuration and Anderson is meant to be that.
%
% RUNTIME HERE MEANS LITTLE. At n_a_GE=201 every solve is launch bound, which is why the shooting
% timings above show local search LOSING. Read the iteration counts instead. The relevant arithmetic
% is that local search always pays two full sweeps -- iteration 1 and each verification -- so its
% gain is k/(2+(k-2)/s) for k iterations and a per-sweep speed-up s. Anderson exists to cut k, which
% is what that gain is proportional to, so its multiplier is SMALLER than shooting's while its
% absolute cost is lower. The iteration counts below are what settles that.
transpathoptionsAndV=transpathoptionsAnd;
transpathoptionsAndV.verbose=1; % so the Anderson iteration count reaches the diary
tic;
[PricePathAndP,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAndV, simoptions, vfoptions, []);
wait(gpudev)
tAndPlain=toc;
vfoptionsAndLS=vfoptions;
vfoptionsAndLS.localsearch=1;
vfoptionsAndLS.nlocalsearch=3;
tic;
[PricePathAndLS,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAndV, simoptions, vfoptionsAndLS, []);
wait(gpudev)
tAndLS=toc;
fprintf('localsearch Anderson: shooting %6.2fs, shooting+localsearch %6.2fs, Anderson %6.2fs, Anderson+localsearch %6.2fs (all launch bound at n_a=%i, so read the iteration counts) \n',tA,tC,tAndPlain,tAndLS,prod(n_a_GE))
fprintf('localsearch Anderson: (12) Anderson+localsearch vs Anderson without it, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathAndLS.r)-gather(PricePathAndP.r))))
fprintf('localsearch Anderson: (12) Anderson+localsearch vs Anderson without it, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathAndLS.w)-gather(PricePathAndP.w))))
fprintf('localsearch Anderson: (12) Anderson+localsearch vs the standard SHOOTING solve, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathAndLS.r)-gather(PricePathA.r))))
fprintf('localsearch Anderson: (12) Anderson+localsearch vs the standard SHOOTING solve, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathAndLS.w)-gather(PricePathA.w))))

%% (13) The same with the grid interpolation layer, so the fork is exercised on a coarse reference
% Against the standard GI Anderson solve, not the no-GI one: GI moves the answer on its own.
tic;
[PricePathAndGI,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAndV, simoptionsGI, vfoptionsGI, []);
wait(gpudev)
tAndGI=toc;
vfoptionsAndGILS=vfoptionsGI;
vfoptionsAndGILS.localsearch=1;
vfoptionsAndGILS.nlocalsearch=3;
tic;
[PricePathAndGILS,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE_GI, AgentDist_initialGE_GI, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAndV, simoptionsGI, vfoptionsAndGILS, []);
wait(gpudev)
tAndGILS=toc;
fprintf('localsearch Anderson: GI, Anderson %6.2fs, Anderson+localsearch %6.2fs \n',tAndGI,tAndGILS)
fprintf('localsearch Anderson: (13) GI, Anderson+localsearch vs Anderson without it, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathAndGILS.r)-gather(PricePathAndGI.r))))
fprintf('localsearch Anderson: (13) GI, Anderson+localsearch vs Anderson without it, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathAndGILS.w)-gather(PricePathAndGI.w))))

%% (14) anderson.type='I' with local search is refused by name
% The Type-I implementation was not forked: its step is unbounded and backtracks toward the previous
% iterate, so a backtracked trial is not cleanly separable from an accepted one and the rule that the
% reference commits only on an accepted step has nowhere to live.
err14=0;
transpathoptionsAndI=transpathoptionsAnd;
transpathoptionsAndI.anderson.type='I';
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAndI, simoptions, vfoptionsAndLS, []);
catch ME
    err14=contains(ME.message,'not implemented with vfoptions.localsearch on a transition path');
end
fprintf('localsearch Anderson: (14) anderson.type=''I'' with local search is refused by name, should be one: %i \n',err14)

%% (15) WITH A d, ANDERSON WITH LOCAL SEARCH -- the fourth and last corner
% Nothing in the fork or the seed route is algorithm-specific: the seed sweep is produced by the
% dispatcher and the reference is carried by whichever solver holds it, so with-d Anderson should
% work with no further plumbing. "Should" is exactly what a test is for, and two things could have
% gone wrong without one. The reference is [N_d,N_a,N_z] here rather than [1,N_a,N_z], so the shim
% has to carry an object of a different shape than it was written against; and the seed sweep with a
% d routes to LSseed_raw rather than to the standard raw, which is a path Anderson had never taken.
%
% THE TOLERANCE HAS TO BE TIGHT ENOUGH THAT ANDERSON ACTUALLY ITERATES, which is the trap this
% section fell into on its first run. The with-d stationary solve leaves a residual of about 5.6e-02,
% two orders worse than the no-d model's 5.9e-04, and tying the path tolerance to it (5 times, so
% 2.8e-01) made that tolerance looser than the residual of the BUMPED STARTING PATH. Anderson stops
% as soon as the residual is inside tolerance, so all four solves "converged in 1 iterations",
% returned their input unchanged, and the comparisons below came back as exactly 0.000e+00 four times
% over -- four checks passing while testing nothing at all. Shooting does not short-circuit the same
% way, which is why (9) and (10) were unaffected and agreed to 3.8e-06 and 8.4e-08.
% So the tolerance is a twentieth of that. The iteration budget then has to be large enough for both
% halves to CONVERGE: at maxiter=30 the no-GI pair converged in 15 and 16 iterations but the GI pair
% both ran out, and comparing two NON-CONVERGED paths is only a statement that two unfinished
% trajectories had not yet diverged much -- it passed the 0.1 bar on a 2.2e-02 disagreement, which is
% larger than the tolerance they were aiming for, so the bar was carrying the check rather than the
% code. Two converged answers is the only version that tests anything.
%
% WHY THE TOLERANCE HERE IS LOOSE, AND WHY THAT IS PROBABLY NOT FIXABLE BY TIGHTENING ANYTHING. The
% with-d stationary solve leaves 5.6e-02, so every tolerance in this block is wide and two converged
% solves may stop anywhere inside it: the with-d agreements below are 1e-03-ish where the no-d ones
% in (12) are 1e-06-ish, purely because the ball is two orders wider.
% It is tempting to call that an under-converged stationary solve and tighten it. It is not. K=E[a]
% and N=E[d*z] are step functions of the prices on coarse grids, so the general eqm objective is
% piecewise-constant and the search stalls -- no tolerance or iteration count reaches past that. The
% setup file records the same thing and is why n_d_GE is 51 rather than 9.
% The diagnostic beside the stationary solve above measures whether the remaining floor is the a grid
% (in which case the interpolation layer smooths it and the with-d block should be built on the GI
% equilibrium) or something else (in which case these tolerances are simply right for this model).
transpathoptionsAndd=transpathoptionsAnd;
transpathoptionsAndd.toleranceGEcondns=transpathoptionsSd.toleranceGEcondns/20;
transpathoptionsAndd.maxiter=120;
transpathoptionsAndd.verbose=1; % so the iteration counts reach the diary

[PricePathAnddP,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd, AgentDist_initialGEd, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsAndd, simoptions, vfoptions, []);
vfoptionsAnddLS=vfoptions;
vfoptionsAnddLS.localsearch=1;
vfoptionsAnddLS.nlocalsearch=3;
[PricePathAnddLS,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd, AgentDist_initialGEd, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsAndd, simoptions, vfoptionsAnddLS, []);
fprintf('localsearch Anderson d: (15) with d, Anderson+localsearch vs Anderson without it, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathAnddLS.r)-gather(PricePathAnddP.r))))
fprintf('localsearch Anderson d: (15) with d, Anderson+localsearch vs Anderson without it, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathAnddLS.w)-gather(PricePathAnddP.w))))
% THE ANTI-VACUITY GUARD, which is what the first run of this section lacked. An agreement check
% between two solves that both stopped at iteration 1 is satisfied by them both returning the price
% path handed in, and reads as a perfect 0.000e+00. So require that the solve MOVED: if the restricted
% Anderson solve is still sitting on its starting path then whatever the agreement check says, it
% says nothing about local search.
fprintf('localsearch Anderson d: (15) the restricted solve moved off its starting path, w, should NOT be zero: %.3e \n',max(abs(gather(PricePathAnddLS.w)-gather(PricePathBumpedd.w))))
% NO COMPARISON AGAINST THE SHOOTING ANSWER. On this calibration toleranceGEcondns permits price
% paths that differ by tens of per cent -- w is about 0.68 and the first run of this section reported
% a 3.5e-01 gap -- so comparing two different ALGORITHMS across that ball tests the width of the ball
% and not the code. (9) already certifies the with-d shooting answer against its own standard solve,
% and the Anderson-against-Anderson check above is what isolates local search.

% And with the interpolation layer, which is the corner where the reference is a COARSE index per d
% and the seed comes from the standard GI raw's midpoint rather than from LSseed_raw.
[PricePathAnddGI,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd_GI, AgentDist_initialGEd_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsAndd, simoptionsGId, vfoptionsGId, []);
vfoptionsAnddGILS=vfoptionsGId;
vfoptionsAnddGILS.localsearch=1;
vfoptionsAnddGILS.nlocalsearch=3;
[PricePathAnddGILS,~]=TransitionPath_InfHorz(PricePathBumpedd, ParamPathGE, T, V_finalGEd_GI, AgentDist_initialGEd_GI, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFnd, FnsToEvaluateGEd, GeneralEqmEqnsGE, ParamsGEd, DiscountFactorParamNames, transpathoptionsAndd, simoptionsGId, vfoptionsAnddGILS, []);
fprintf('localsearch Anderson d: (15) with d and GI, Anderson+localsearch vs Anderson without it, r, this should be close to zero: %.3e \n',max(abs(gather(PricePathAnddGILS.r)-gather(PricePathAnddGI.r))))
fprintf('localsearch Anderson d: (15) with d and GI, Anderson+localsearch vs Anderson without it, w, this should be close to zero: %.3e \n',max(abs(gather(PricePathAnddGILS.w)-gather(PricePathAnddGI.w))))
fprintf('localsearch Anderson d: (15) with GI, the restricted solve moved off its starting path, w, should NOT be zero: %.3e \n',max(abs(gather(PricePathAnddGILS.w)-gather(PricePathBumpedd.w))))
fprintf('localsearch Anderson d: (15) for scale, toleranceGEcondns here is %.3e, against %.3e for the with-d shooting solves \n',transpathoptionsAndd.toleranceGEcondns,transpathoptionsSd.toleranceGEcondns)

output=struct();

end
