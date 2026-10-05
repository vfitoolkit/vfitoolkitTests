function output=CoreInfHorzTPath_tminus1tplus1(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Checks that _tminus1 and _tplus1 inputs do what they are meant to on an InfHorz transition path.
% The InfHorz counterpart of CoreFHorzTPath_tminus1tplus1.
%
% THE METHOD. Each lag or lead is checked against an ORACLE that does the same shift by hand: the
% shifted series is put on the ParamPath under its own name (rlag=[initialvalue, r(1:T-1)]), and the
% same function is written once with r_tminus1 and once with rlag. Both are evaluated at the same
% price path, and the general eqm conditions (the second output of TransitionPath_InfHorz) must be
% identical. Each check also reports, for scale, what the UNLAGGED version gives, so that a zero is
% not mistaken for a comparison of two things that were never different.
%
% EVALUATING WITHOUT MOVING THE PATH. The probe uses factor=0 in howtoupdate, so the update rule
% cannot move any price whatever the algorithm does, and GEcondnPath is therefore the conditions at
% exactly the path handed in. That is what makes it safe to use GEnewprice=2 here: this bank runs
% Anderson throughout (the algorithms are compared in with_CoreInfHorzTPathAlgoTests, not here), and
% with factor=0 Anderson, shooting and quasi-Newton would all return the same thing. maxiter=1 on top
% so nothing is paid for iterations that cannot change anything.
%
% Semantics being checked (TransitionPath_InfHorz_singlepathiter, EvalFnOnTransPath_AggVars_InfHorz):
%   X_tminus1 is X in period t-1, and transpathoptions.initialvalues.X in period 1. X can be a price,
%     a ParamPath parameter, or an aggregate (a FnsToEvaluate).
%   p_tplus1 is the price in period t+1 (only prices). At the end of the path it is period T.
%   They can be used in GeneralEqmEqns, FnsToEvaluate and transpathoptions.intermediateEqns.
%
% (1)-(7) GE eqns, FnsToEvaluate and intermediateEqns, each against its oracle
% (8) a real general eqm path whose condition uses lags, solved by Anderson and by shooting; and the
%     converged path re-evaluated with the oracle must give back the conditions the solver returned
% (9) EvalFnOnTransPath_AggVars_InfHorz with lagged/led FnsToEvaluate against the same aggregates
%     computed inside the transition path
% (10) SimPanelValues_TransPath_InfHorz with a lagged FnsToEvaluate: the per-period panel mean must
%     track the same aggregate (Monte Carlo, so close rather than zero)
%
% No fastOLG: that is an OLG-specific optimization and does not exist for InfHorz, so where the FHorz
% version runs (1)-(7) twice this runs them once.

vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
n_d=0; d_grid=[]; % this check does not need d; the lag machinery is independent of it

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;

% V_final and the initial distribution: the stationary model at the period-1 parameters. Any would do,
% the checks are identities, but this keeps the path somewhere sensible.
[V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_final,n_d,n_a,n_z,pi_z,simoptions,Params,[]);

%% The initial values, deliberately not equal to period 1 of anything, so that period 1 tests them
initialvalues.r=0.04;      % PricePath.r(1) is 0.03
initialvalues.w=1.2;       % ParamPath.w(1) is 1.1
initialvalues.assets=0.5;

%% The oracles: the shifted series, as ParamPath entries of their own
ParamPathO=ParamPath;
ParamPathO.rlag=[initialvalues.r, PricePath.r(1:T-1)];
ParamPathO.rlead=[PricePath.r(2:T), PricePath.r(T)]; % the last entry is never read (the path runs over t=1,...,T-1)
ParamPathO.wlag=[initialvalues.w, ParamPath.w(1:T-1)];
Params.rlag=ParamPathO.rlag(1);
Params.rlead=ParamPathO.rlead(1);
Params.wlag=ParamPathO.wlag(1);
Params.Alag=initialvalues.assets;

% The evaluate-only probe. factor=0, so no price can move and the conditions come back at the path
% that went in.
transpathoptionsE=transpathoptionsbaseline;
transpathoptionsE.maxiter=1;
transpathoptionsE.GEnewprice=2;
transpathoptionsE.GEnewprice2.howtoupdate={'g','r',0,0};
transpathoptionsE.initialvalues=initialvalues;
transpathoptionsE.verbose=0;

% The probe runs the solver at maxiter=1, so every one of the calls below emits an
% 'AndersonAcceleration: reached maxiter (1) without convergence' warning. That is expected and is not
% a failure: the probe is an evaluator, not a solve, and with factor=0 it could never converge. Said
% once here so the warnings in the diary are not mistaken for results. They cannot be suppressed from
% here without a blanket warning('off'), which would hide real ones.
fprintf('NOTE: the checks below run the solver at maxiter=1 purely to evaluate the GE conditions at a \n')
fprintf('      given path, so each emits an expected "reached maxiter (1) without convergence" warning. \n')

%% The assets path, exposed, and its lag (needed by the oracles of (4), (6) and (7))
GeneralEqmEqnsA.g=@(assets) assets;
[~,GEcondnPathA]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsA, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
Apath=gather(GEcondnPathA.g); % periods 1,...,T-1
ParamPathO.Alag=[initialvalues.assets, Apath(1:T-2), Apath(T-1)]; % the last entry is never read

%% (1) _tminus1 of a price, in a GE eqn
GeneralEqmEqns1.g=@(r_tminus1,assets) r_tminus1*assets;
[~,GEcondnPath1]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns1, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
GeneralEqmEqns1O.g=@(rlag,assets) rlag*assets;
[~,GEcondnPath1O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns1O, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
GeneralEqmEqns1U.g=@(r,assets) r*assets; % unlagged: only to show the comparison above is not empty
[~,GEcondnPath1U]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns1U, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tminus1: for scale, r_tminus1 vs r in the same GE eqn differ by %.3e (must not be zero, or (1) below is empty) \n',max(abs(GEcondnPath1.g-GEcondnPath1U.g)))
fprintf('_tminus1: (1) price r_tminus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',max(abs(GEcondnPath1.g-GEcondnPath1O.g)))

%% (2) _tplus1 of a price, in a GE eqn
GeneralEqmEqns2.g=@(r_tplus1,assets) r_tplus1*assets;
[~,GEcondnPath2]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns2, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
GeneralEqmEqns2O.g=@(rlead,assets) rlead*assets;
[~,GEcondnPath2O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns2O, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tplus1: for scale, r_tplus1 vs r_tminus1 in the same GE eqn differ by %.3e (must not be zero) \n',max(abs(GEcondnPath2.g-GEcondnPath1.g)))
fprintf('_tplus1: (2) price r_tplus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',max(abs(GEcondnPath2.g-GEcondnPath2O.g)))

%% (3) _tminus1 of a ParamPath parameter, in a GE eqn
GeneralEqmEqns3.g=@(w_tminus1,earnings) w_tminus1*earnings;
[~,GEcondnPath3]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns3, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
GeneralEqmEqns3O.g=@(wlag,earnings) wlag*earnings;
[~,GEcondnPath3O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns3O, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tminus1: (3) parameter w_tminus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',max(abs(GEcondnPath3.g-GEcondnPath3O.g)))

%% (4) _tminus1 of an aggregate, in a GE eqn
GeneralEqmEqns4.g=@(assets_tminus1,r) assets_tminus1*r;
[~,GEcondnPath4]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns4, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
GeneralEqmEqns4O.g=@(Alag,r) Alag*r;
[~,GEcondnPath4O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns4O, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tminus1: (4) aggregate assets_tminus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',max(abs(GEcondnPath4.g-GEcondnPath4O.g)))

%% (5) _tminus1 and _tplus1 of prices and a parameter, inside a FnsToEvaluate
FnsToEvaluate5=FnsToEvaluate;
FnsToEvaluate5.x=@(aprime,a,z,r_tminus1,r_tplus1,w_tminus1) a*r_tminus1+aprime*r_tplus1+w_tminus1*z;
FnsToEvaluate5O=FnsToEvaluate;
FnsToEvaluate5O.x=@(aprime,a,z,rlag,rlead,wlag) a*rlag+aprime*rlead+wlag*z;
GeneralEqmEqns5.g=@(x) x;
[~,GEcondnPath5]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate5, GeneralEqmEqns5, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
[~,GEcondnPath5O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate5O, GeneralEqmEqns5, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tminus1/_tplus1: (5) r_tminus1, r_tplus1, w_tminus1 in a FnsToEvaluate vs the shifted paths, this should be zero: %.3e \n',max(abs(GEcondnPath5.g-GEcondnPath5O.g)))

%% (6) _tminus1 of an aggregate, inside another FnsToEvaluate
FnsToEvaluate6=FnsToEvaluate;
FnsToEvaluate6.y=@(aprime,a,z,assets_tminus1) a*assets_tminus1;
FnsToEvaluate6O=FnsToEvaluate;
FnsToEvaluate6O.y=@(aprime,a,z,Alag) a*Alag;
GeneralEqmEqns6.g=@(y) y;
[~,GEcondnPath6]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate6, GeneralEqmEqns6, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
[~,GEcondnPath6O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate6O, GeneralEqmEqns6, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tminus1: (6) aggregate assets_tminus1 in a FnsToEvaluate vs the shifted path, this should be zero: %.3e \n',max(abs(GEcondnPath6.g-GEcondnPath6O.g)))

%% (7) _tminus1 and _tplus1 inside an intermediateEqn (and nowhere else)
% This is also the only place in the InfHorz banks where intermediateEqns is checked for CORRECTNESS
% rather than just run: elsewhere it is set and the path solves, which would not catch a wrong value.
transpathoptions7=transpathoptionsE;
transpathoptions7.intermediateEqns.m=@(r_tminus1,r_tplus1,assets_tminus1,assets) r_tminus1*assets+r_tplus1*assets_tminus1;
GeneralEqmEqns7.g=@(m) m;
[~,GEcondnPath7]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns7, Params, DiscountFactorParamNames, transpathoptions7, simoptions, vfoptions, []);
transpathoptions7.intermediateEqns.m=@(rlag,rlead,Alag,assets) rlag*assets+rlead*Alag;
[~,GEcondnPath7O]=TransitionPath_InfHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns7, Params, DiscountFactorParamNames, transpathoptions7, simoptions, vfoptions, []);
fprintf('_tminus1/_tplus1: (7) r_tminus1, r_tplus1, assets_tminus1 in an intermediateEqn vs the shifted paths, this should be zero: %.3e \n',max(abs(GEcondnPath7.g-GEcondnPath7O.g)))

%% (8) A general eqm path whose condition uses lags, solved for real
% r_t = 0.5*r_{t-1} + 0.5*0.05*(1+A_{t-1})^(-0.5). The lag on r makes it a contraction, so both
% solvers can find it. On a moving path the shifted-series oracle cannot be built in advance (it
% depends on the answer), so instead: Anderson and shooting must agree on the answer, and the answer,
% re-evaluated with the oracle built from it, must give back exactly the conditions Anderson returned.
GeneralEqmEqns8.rmarket=@(r,r_tminus1,assets_tminus1) r-(0.5*r_tminus1+0.5*0.05*(1+assets_tminus1)^(-0.5));

transpathoptions8A=transpathoptionsbaseline;
transpathoptions8A.GEnewprice=2;
transpathoptions8A.GEnewprice2.howtoupdate={'rmarket','r',0,0.5};
transpathoptions8A.maxiter=200;
transpathoptions8A.initialvalues=initialvalues;
transpathoptions8A.verbose=0;
[PricePath8A,GEcondnPath8A]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns8, Params, DiscountFactorParamNames, transpathoptions8A, simoptions, vfoptions, []);

transpathoptions8B=rmfield(transpathoptions8A,'GEnewprice2');
transpathoptions8B.GEnewprice=3;
transpathoptions8B.GEnewprice3.howtoupdate={'rmarket','r',0,0.5};
transpathoptions8B.maxiter=500;
[PricePath8B,GEcondnPath8B]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns8, Params, DiscountFactorParamNames, transpathoptions8B, simoptions, vfoptions, []);

fprintf('_tminus1 GE path: how far Anderson moved r from the starting guess: %.3e \n',max(abs(PricePath8A.r-PricePath.r)))
fprintf('_tminus1 GE path, Anderson, GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',max(abs(GEcondnPath8A.rmarket)),1e-4)
fprintf('_tminus1 GE path, shooting, GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',max(abs(GEcondnPath8B.rmarket)),1e-4)
fprintf('_tminus1 GE path, Anderson vs shooting, same eqm path, should be small: %.3e \n',max(abs(PricePath8A.r-PricePath8B.r)))

% Re-evaluate Anderson's answer with the oracle built from it
GeneralEqmEqnsA8.g=@(assets) assets;
[~,GEcondnPathA8]=TransitionPath_InfHorz(PricePath8A, ParamPath, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsA8, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
Apath8=gather(GEcondnPathA8.g);
ParamPath8O=ParamPath;
ParamPath8O.rlag=[initialvalues.r, gather(PricePath8A.r(1:T-1))];
ParamPath8O.Alag=[initialvalues.assets, Apath8(1:T-2), Apath8(T-1)];
GeneralEqmEqns8O.g=@(r,rlag,Alag) r-(0.5*rlag+0.5*0.05*(1+Alag)^(-0.5));
[~,GEcondnPath8O]=TransitionPath_InfHorz(PricePath8A, ParamPath8O, T, V_final, AgentDist_initial, n_d, n_a, n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns8O, Params, DiscountFactorParamNames, transpathoptionsE, simoptions, vfoptions, []);
fprintf('_tminus1 GE path: conditions Anderson returned vs its answer re-evaluated with the shifted paths, this should be zero: %.3e \n',max(abs(GEcondnPath8A.rmarket-GEcondnPath8O.g)))

%% (9) EvalFnOnTransPath_AggVars_InfHorz with lagged and led FnsToEvaluate
% The same aggregates as (5) and (6), computed after the fact from the policy and agent distribution
% paths, must equal the ones the transition path computed inside itself. Compared over t=1,...,T-1,
% the periods the transition path evaluates. This is the check that needed
% EvalFnOnTransPath_AggVars_InfHorz to take transpathoptions: without it there is no initialvalues
% for period 1 and the lag cannot be formed at all.
transpathoptions9=transpathoptionsbaseline;
transpathoptions9.initialvalues=initialvalues;
[~,PolicyPath9]=ValueFnOnTransPath_InfHorz(PricePath, ParamPathO, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions9, vfoptions);
AgentDistPath9=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPathO, PolicyPath9, n_d,n_a,n_z, pi_z, T, Params, transpathoptions9, simoptions);
FnsToEvaluate9=FnsToEvaluate;
FnsToEvaluate9.x=FnsToEvaluate5.x;
FnsToEvaluate9.y=FnsToEvaluate6.y;
AggVarsPath9=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate9, AgentDistPath9, PolicyPath9, PricePath, ParamPathO, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid, transpathoptions9, simoptions);
fprintf('EvalFnOnTransPath_AggVars: (9) assets vs inside the transition path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath9.assets.Mean(1:T-1))-Apath)))
fprintf('EvalFnOnTransPath_AggVars: (9) x with r_tminus1, r_tplus1, w_tminus1 vs inside the transition path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath9.x.Mean(1:T-1))-gather(GEcondnPath5.g))))
fprintf('EvalFnOnTransPath_AggVars: (9) y with assets_tminus1 vs inside the transition path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath9.y.Mean(1:T-1))-gather(GEcondnPath6.g))))

%% (10) SimPanelValues_TransPath_InfHorz with a lagged FnsToEvaluate
% The other command that needed transpathoptions. A Monte Carlo panel, so this is a 'close to zero'
% check and not an identity: the per-period panel mean of x estimates the same aggregate (9) computes
% exactly. A lag that was silently wrong would show up as a bias far larger than the sampling error,
% which is why this is worth having even though it cannot be tight.
% x only, NOT y: SimPanelValues_TransPath_InfHorz refuses the _tminus1 of an AGGREGATE by design, and
% rightly so, since a simulated panel gives individual agents and not the population aggregate that
% assets_tminus1 would need. Price and parameter lags are fine, because those are path objects every
% agent shares. (9) is where the aggregate lag gets its identity check.
FnsToEvaluate10=FnsToEvaluate;
FnsToEvaluate10.x=FnsToEvaluate5.x;
simoptions10=simoptions;
simoptions10.numbersims=10^4; % the same count the bank's other SimPanel checks use
simoptions10.simperiods=T;
SimPanel10=SimPanelValues_TransPath_InfHorz(PolicyPath9, PricePath, ParamPathO, T, AgentDist_initial, n_d, n_a, n_z, pi_z, d_grid,a_grid,z_grid, FnsToEvaluate10, Params, transpathoptions9, simoptions10);
panelmean=mean(SimPanel10.x,2)';
fprintf('SimPanel along TPath: (10) x with r_tminus1, r_tplus1, w_tminus1, panel mean vs AggVars, should be close to zero (Monte Carlo): %.3e \n',max(abs(gather(panelmean(1:T-1))-gather(AggVarsPath9.x.Mean(1:T-1)))))

output=struct();

end
