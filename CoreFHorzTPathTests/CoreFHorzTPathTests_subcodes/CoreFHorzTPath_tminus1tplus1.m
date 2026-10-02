function output=CoreFHorzTPath_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% Checks that _tminus1 and _tplus1 inputs do what they are meant to on an FHorz transition path.
%
% THE METHOD. Each lag or lead is checked against an oracle that does the same shift by hand: the
% shifted series is put on the ParamPath under its own name (e.g. rlag=[initialvalue, r(1:T-1)]), and
% the same function is written once with r_tminus1 and once with rlag. Both are evaluated at the same
% price path, and the general eqm conditions (the second output of TransitionPath_Case1_FHorz) must be
% identical. To evaluate without moving the path, the shooting algorithm is run with maxiter=1: it
% returns the conditions at exactly the path it was handed (Anderson would take a step first). Where
% an aggregate path is needed, a general eqm eqn that just returns it, @(assets) assets, exposes it.
%
% Semantics being checked (TransitionPath_FHorz_singlepathiter, EvalFnOnTransPath_AggVars_Case1_FHorz):
%   X_tminus1 is X in period t-1, and transpathoptions.initialvalues.X in period 1. X can be a price,
%     a ParamPath parameter, or an aggregate (a FnsToEvaluate).
%   p_tplus1 is the price in period t+1 (only prices). In period T-1 it is period T of the price path.
%   They can be used in GeneralEqmEqns, FnsToEvaluate and transpathoptions.intermediateEqns.
%
% (1)-(7) GE eqns, FnsToEvaluate and intermediateEqns, each against its oracle, with fastOLG=1 and 0
% (8) a general eqm path whose condition uses lags, solved by Anderson and by shooting; and the
%     converged path re-evaluated with the oracle must give back the conditions the solver returned
% (9) EvalFnOnTransPath_AggVars_Case1_FHorz with lagged/led FnsToEvaluate against the same aggregates
%     inside the transition path

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];

% zero assets, mid point for the shock
jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

% V_final and the initial distribution: the stationary model at the period-1 parameters. Any would do,
% the checks are identities, but this keeps the path somewhere sensible.
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_final,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

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

fastOLGopts=[1,0];
for ff=1:length(fastOLGopts)
    % Evaluate-only: one shooting iteration returns the conditions at the path it was given
    transpathoptionsE=transpathoptionsbaseline;
    transpathoptionsE.fastOLG=fastOLGopts(ff);
    transpathoptionsE.maxiter=1;
    transpathoptionsE.GEnewprice=3;
    transpathoptionsE.GEnewprice3.howtoupdate={'g','r',0,0.1};
    transpathoptionsE.initialvalues=initialvalues;
    transpathoptionsE.verbose=0;

    %% The assets path, exposed, and its lag (needed by the oracles of (4), (6) and (7))
    GeneralEqmEqnsA.g=@(assets) assets;
    [~,GEcondnPathA]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsA, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    Apath=gather(GEcondnPathA.g); % periods 1,...,T-1
    ParamPathO.Alag=[initialvalues.assets, Apath(1:T-2), Apath(T-1)]; % the last entry is never read

    %% (1) _tminus1 of a price, in a GE eqn
    GeneralEqmEqns1.g=@(r_tminus1,assets) r_tminus1*assets;
    [~,GEcondnPath1]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns1, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    GeneralEqmEqns1O.g=@(rlag,assets) rlag*assets;
    [~,GEcondnPath1O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns1O, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    GeneralEqmEqns1U.g=@(r,assets) r*assets; % unlagged: only to show the comparison above is not empty
    [~,GEcondnPath1U]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns1U, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('_tminus1 (fastOLG=%i): for scale, r_tminus1 vs r in the same GE eqn differ by %.3e (must not be zero, or (1) below is empty) \n',fastOLGopts(ff),max(abs(GEcondnPath1.g-GEcondnPath1U.g)))
    fprintf('_tminus1 (fastOLG=%i): (1) price r_tminus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath1.g-GEcondnPath1O.g)))

    %% (2) _tplus1 of a price, in a GE eqn
    GeneralEqmEqns2.g=@(r_tplus1,assets) r_tplus1*assets;
    [~,GEcondnPath2]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns2, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    GeneralEqmEqns2O.g=@(rlead,assets) rlead*assets;
    [~,GEcondnPath2O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns2O, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('_tplus1 (fastOLG=%i): (2) price r_tplus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath2.g-GEcondnPath2O.g)))

    %% (3) _tminus1 of a ParamPath parameter, in a GE eqn
    GeneralEqmEqns3.g=@(w_tminus1,earnings) w_tminus1*earnings;
    [~,GEcondnPath3]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns3, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    GeneralEqmEqns3O.g=@(wlag,earnings) wlag*earnings;
    [~,GEcondnPath3O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns3O, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('_tminus1 (fastOLG=%i): (3) parameter w_tminus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath3.g-GEcondnPath3O.g)))

    %% (4) _tminus1 of an aggregate, in a GE eqn
    GeneralEqmEqns4.g=@(assets_tminus1,r) assets_tminus1*r;
    [~,GEcondnPath4]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns4, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    GeneralEqmEqns4O.g=@(Alag,r) Alag*r;
    [~,GEcondnPath4O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns4O, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('_tminus1 (fastOLG=%i): (4) aggregate assets_tminus1 in a GE eqn vs the shifted path, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath4.g-GEcondnPath4O.g)))

    %% (5) _tminus1 and _tplus1 of prices and a parameter, inside a FnsToEvaluate
    FnsToEvaluate5=FnsToEvaluate;
    FnsToEvaluate5.x=@(aprime,a,z,r_tminus1,r_tplus1,w_tminus1) a*r_tminus1+aprime*r_tplus1+w_tminus1*z;
    FnsToEvaluate5O=FnsToEvaluate;
    FnsToEvaluate5O.x=@(aprime,a,z,rlag,rlead,wlag) a*rlag+aprime*rlead+wlag*z;
    GeneralEqmEqns5.g=@(x) x;
    [~,GEcondnPath5]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate5, GeneralEqmEqns5, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath5O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate5O, GeneralEqmEqns5, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('_tminus1/_tplus1 (fastOLG=%i): (5) r_tminus1, r_tplus1, w_tminus1 in a FnsToEvaluate vs the shifted paths, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath5.g-GEcondnPath5O.g)))

    %% (6) _tminus1 of an aggregate, inside another FnsToEvaluate
    FnsToEvaluate6=FnsToEvaluate;
    FnsToEvaluate6.y=@(aprime,a,z,assets_tminus1) a*assets_tminus1;
    FnsToEvaluate6O=FnsToEvaluate;
    FnsToEvaluate6O.y=@(aprime,a,z,Alag) a*Alag;
    GeneralEqmEqns6.g=@(y) y;
    [~,GEcondnPath6]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate6, GeneralEqmEqns6, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    [~,GEcondnPath6O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate6O, GeneralEqmEqns6, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
    fprintf('_tminus1 (fastOLG=%i): (6) aggregate assets_tminus1 in a FnsToEvaluate vs the shifted path, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath6.g-GEcondnPath6O.g)))

    %% (7) _tminus1 and _tplus1 inside an intermediateEqn (and nowhere else)
    % Before the lags were looked for in intermediateEqns, this one could not run at all.
    transpathoptions7=transpathoptionsE;
    transpathoptions7.intermediateEqns.m=@(r_tminus1,r_tplus1,assets_tminus1,assets) r_tminus1*assets+r_tplus1*assets_tminus1;
    GeneralEqmEqns7.g=@(m) m;
    [~,GEcondnPath7]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns7, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions7, simoptions, vfoptions);
    transpathoptions7.intermediateEqns.m=@(rlag,rlead,Alag,assets) rlag*assets+rlead*Alag;
    [~,GEcondnPath7O]=TransitionPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns7, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions7, simoptions, vfoptions);
    fprintf('_tminus1/_tplus1 (fastOLG=%i): (7) r_tminus1, r_tplus1, assets_tminus1 in an intermediateEqn vs the shifted paths, this should be zero: %.3e \n',fastOLGopts(ff),max(abs(GEcondnPath7.g-GEcondnPath7O.g)))
end


%% (8) A general eqm path whose condition uses lags, solved for real
% r_t = 0.5*r_{t-1} + 0.5*0.05*(1+A_{t-1})^(-0.5). The lag on r makes it a contraction, so both solvers
% can find it. On a moving path the shifted-series oracle cannot be built in advance (it depends on the
% answer), so instead: Anderson and shooting must agree on the answer, and the answer, re-evaluated with
% the oracle built from it, must give back exactly the conditions Anderson returned with it.
GeneralEqmEqns8.rmarket=@(r,r_tminus1,assets_tminus1) r-(0.5*r_tminus1+0.5*0.05*(1+assets_tminus1)^(-0.5));

transpathoptions8A=transpathoptionsbaseline;
transpathoptions8A.GEnewprice=2;
transpathoptions8A.GEnewprice2.howtoupdate={'rmarket','r',0,0.5};
transpathoptions8A.maxiter=200;
transpathoptions8A.initialvalues=initialvalues;
transpathoptions8A.verbose=0;
[PricePath8A,GEcondnPath8A]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns8, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions8A, simoptions, vfoptions);

transpathoptions8B=rmfield(transpathoptions8A,'GEnewprice2');
transpathoptions8B.GEnewprice=3;
transpathoptions8B.GEnewprice3.howtoupdate={'rmarket','r',0,0.5};
transpathoptions8B.maxiter=500;
[PricePath8B,GEcondnPath8B]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns8, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions8B, simoptions, vfoptions);

fprintf('_tminus1 GE path: how far Anderson moved r from the starting guess: %.3e \n',max(abs(PricePath8A.r-PricePath.r)))
fprintf('_tminus1 GE path, Anderson, GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',max(abs(GEcondnPath8A.rmarket)),1e-4)
fprintf('_tminus1 GE path, shooting, GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',max(abs(GEcondnPath8B.rmarket)),1e-4)
fprintf('_tminus1 GE path, Anderson vs shooting, same eqm path, should be small: %.3e \n',max(abs(PricePath8A.r-PricePath8B.r)))

% Re-evaluate Anderson's answer with the oracle built from it
transpathoptionsE=transpathoptionsbaseline;
transpathoptionsE.maxiter=1;
transpathoptionsE.GEnewprice=3;
transpathoptionsE.GEnewprice3.howtoupdate={'g','r',0,0.1};
transpathoptionsE.initialvalues=initialvalues;
transpathoptionsE.verbose=0;
GeneralEqmEqnsA.g=@(assets) assets;
[~,GEcondnPathA8]=TransitionPath_Case1_FHorz(PricePath8A, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsA, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
Apath8=gather(GEcondnPathA8.g);
ParamPath8O=ParamPath;
ParamPath8O.rlag=[initialvalues.r, gather(PricePath8A.r(1:T-1))];
ParamPath8O.Alag=[initialvalues.assets, Apath8(1:T-2), Apath8(T-1)];
GeneralEqmEqns8O.g=@(r,rlag,Alag) r-(0.5*rlag+0.5*0.05*(1+Alag)^(-0.5));
[~,GEcondnPath8O]=TransitionPath_Case1_FHorz(PricePath8A, ParamPath8O, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns8O, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsE, simoptions, vfoptions);
fprintf('_tminus1 GE path: conditions Anderson returned vs its answer re-evaluated with the shifted paths, this should be zero: %.3e \n',max(abs(GEcondnPath8A.rmarket-GEcondnPath8O.g)))


%% (9) EvalFnOnTransPath_AggVars_Case1_FHorz with lagged and led FnsToEvaluate
% The same aggregates as (5) and (6), computed after the fact from the policy and agent distribution
% paths, must equal the ones the transition path computed inside itself. Compared over t=1,...,T-1,
% the periods the transition path evaluates.
transpathoptions9=transpathoptionsbaseline;
transpathoptions9.initialvalues=initialvalues;
[~,PolicyPath9]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPathO, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions9, vfoptions);
AgentDistPath9=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPathO, PolicyPath9, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions9, simoptions);
FnsToEvaluate9=FnsToEvaluate;
FnsToEvaluate9.x=FnsToEvaluate5.x;
FnsToEvaluate9.y=FnsToEvaluate6.y;
AggVarsPath9=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate9, AgentDistPath9, PolicyPath9, PricePath, ParamPathO, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions9, simoptions);
% (fastOLG=0 was the last pass of the loop above, and gives the same paths as fastOLG=1, as (1)-(7) check)
fprintf('EvalFnOnTransPath_AggVars: (9) assets vs inside the transition path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath9.assets.Mean(1:T-1))-Apath)))
fprintf('EvalFnOnTransPath_AggVars: (9) x with r_tminus1, r_tplus1, w_tminus1 vs inside the transition path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath9.x.Mean(1:T-1))-gather(GEcondnPath5.g))))
fprintf('EvalFnOnTransPath_AggVars: (9) y with assets_tminus1 vs inside the transition path, this should be zero: %.3e \n',max(abs(gather(AggVarsPath9.y.Mean(1:T-1))-gather(GEcondnPath6.g))))

%% (10) LifeCycleProfiles_TransPath_FHorz_Case1 with lagged and led FnsToEvaluate
% The age-conditional means of x and y must be the same whether they use the lags or the shifted
% paths, for the ages alive at the start of the transition and for those born during it. Both outputs
% are NaN where an age-period combination does not exist, in the same places, so those are skipped.
FnsToEvaluate10O=FnsToEvaluate;
FnsToEvaluate10O.x=FnsToEvaluate5O.x;
FnsToEvaluate10O.y=FnsToEvaluate6O.y;
LCP10=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate9, AgentDistPath9, PolicyPath9, PricePath, ParamPathO, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, transpathoptions9, simoptions);
LCP10O=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate10O, AgentDistPath9, PolicyPath9, PricePath, ParamPathO, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, transpathoptions9, simoptions);
fprintf('LifeCycleProfiles_TransPath: (10) x with r_tminus1, r_tplus1, w_tminus1, age-conditional means vs the shifted paths, this should be zero: %.3e \n',max([max(abs(LCP10.x.initialages.Mean(:)-LCP10O.x.initialages.Mean(:)),[],'omitnan'), max(abs(LCP10.x.bornduringtranstion.Mean(:)-LCP10O.x.bornduringtranstion.Mean(:)),[],'omitnan')]))
fprintf('LifeCycleProfiles_TransPath: (10) y with assets_tminus1, age-conditional means vs the shifted path, this should be zero: %.3e \n',max([max(abs(LCP10.y.initialages.Mean(:)-LCP10O.y.initialages.Mean(:)),[],'omitnan'), max(abs(LCP10.y.bornduringtranstion.Mean(:)-LCP10O.y.bornduringtranstion.Mean(:)),[],'omitnan')]))

output=struct();

end
