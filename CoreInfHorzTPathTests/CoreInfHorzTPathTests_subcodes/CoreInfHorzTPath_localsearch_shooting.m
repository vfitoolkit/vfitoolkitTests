function output=CoreInfHorzTPath_localsearch_shooting(T,n_a_GE,n_z,a_grid_GE,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Checks vfoptions.localsearch=1 on a REAL general eqm transition path, solved by shooting.
% doPart(4) and doPart(5) test the restricted value fn step against a given price path. This is the
% first part where the restriction is inside a general eqm solve and so can affect the answer.
%
% THE SCHEME. Iteration 1 is a standard sweep, whatever divideandconquer and gridinterplayer are set
% to, and its Policy becomes the reference. Later iterations restrict aprime to a window around the
% previous iteration's answer, refreshing the reference as they go. When those converge, ONE standard
% sweep re-evaluates the SAME price path: still converged means done, and moved means its Policy
% becomes the new reference and local search resumes.
% So the restriction cannot change the answer, only how fast it is reached, and the certification
% rests on nothing -- no claim about the problem being single-peaked, no per-step condition. That is
% what makes this testable as a real check rather than a tolerance.
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
n_d=0; d_grid=[]; % local search with a d is doPart(5); the general eqm scheme is no-d for now

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
vfoptionsE=vfoptionsB;
try
    TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d, n_a_GE, n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAnd, simoptions, vfoptionsE, []);
catch ME
    err(1)=contains(ME.message,'requires transpathoptions.GEnewprice=3');
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
fprintf('localsearch shooting: (5) local search with an algorithm other than shooting is refused by name, should be one: %i \n',err(1))
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

output=struct();

end
