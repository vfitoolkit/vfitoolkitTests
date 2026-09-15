function output=CoreInfHorzTPathAlgo_jacobiandiagnostic(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c)
% Diagnostic, not a test: is the Jacobian badly conditioned because the decision variable makes the
% residual jumpy, or is it intrinsic to the price path?
%
% With a decision variable the optimal d flips discretely as a price moves, so the general eqm
% residual is a step function of the prices and a finite difference may catch one flip or none. Two
% things would fix that if it is the story: a finer d grid, which makes each flip smaller, and a
% larger epsprice, which averages over more flips. This builds one Jacobian for each combination and
% reports rcond, at maxiter=1, so it costs one Jacobian per cell rather than a whole Newton solve.
%
% If rcond improves markedly, jumpiness is the story and a full run at the winning configuration is
% worth paying for. If it stays around 1e-6 everywhere, the ill-conditioning is intrinsic (neighbouring
% columns of the Jacobian are nearly parallel because a price change at period k and at period k+1
% have almost the same effect) and no grid refinement will help.


% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d,aprime,a,z) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*d*z;
FnsToEvaluate.nextassets=@(d,aprime,a,z) aprime; % aprime-dependent, so it tests the policy decode directly
FnsToEvaluate.zval=@(d,aprime,a,z) z; % the shock itself, so the path of its mean/std deviation tracks the z process


%% Period-0 VFI: gives the final-step V/Policy (used as both V_final for TPath and the steady state to compare against)
vfoptions1=vfoptions;
simoptions1=simoptions;

%% General equilibrium: solve the stationary GE first, then run a null-reform transition path
% The only part of this file that exercises the GE machinery itself (the shooting algorithm, the
% price update rule, the convergence test). Everything else either runs a constant path through
% the non-GE commands, or calls TransitionPath_InfHorz with a dummy GE eqn and maxiter=1, so
% nothing ever updates.
%
% A Cobb-Douglas firm supplies BOTH prices from its first order conditions:
%     r = firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha)) - firmdelta   (MPK less depreciation)
%     w = (1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha))                (MPL)
% so w is a general eqm price here, not a ParamPath entry as it is elsewhere in this file.
% K is aggregate assets and N aggregate labour supply.
%
% Both solvers are given the SAME residual-form GeneralEqmEqns. That works because the path uses
% transpathoptions.GEnewprice=3 rather than the default 1: GEnewprice=1 treats the eqns as price
% UPDATING FORMULAE (updatePricePathNew_TPath_tt sets PricePathNew_tt=p_i directly), whereas
% GEnewprice=3 treats them as residuals and updates via howtoupdate, which is the same convention
% HeteroAgentStationaryEqm_InfHorz uses. One definition of the firm, not two.
%
% firm* names are used because the 2A entrepreneur ReturnFn already uses alpha and delta with its
% own (Kitao) meanings.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluateGE.K=@(d,aprime,a,z) a;
FnsToEvaluateGE.N=@(d,aprime,a,z) d*z;

GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w

[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions1, vfoptions1);
fprintf('Stationary GE: r=%2.8f, w=%2.8f \n',p_eqm.r,p_eqm.w)
fprintf('Stationary GE conditions (these should be close to zero): CapitalMarket=%.3e, LabourMarket=%.3e \n',GEcondns.CapitalMarket,GEcondns.LabourMarket)

% The path needs the EQUILIBRIUM V_final and initial dist, else the null reform is not exact
ParamsGE=Params;
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions1);
AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d_GE,n_a_GE,n_z,pi_z,simoptions1,ParamsGE,[]);

PricePathGE.r=p_eqm.r*ones(1,T);
PricePathGE.w=p_eqm.w*ones(1,T);
ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: nothing actually changes, this is a null reform

transpathoptionsGE=transpathoptionsbaseline;
transpathoptionsGE.maxiter=25; % the no-change solve starts at the answer, so it needs few iterations
transpathoptionsGE.verbose=0;
% GEnewprice=3: treat the GeneralEqmEqns as residuals and update each price by a damped fraction
% of its own residual. howtoupdate columns are {GEcondn name, price name, add, factor}, and
% updatePricePathNew_TPath_tt does new = old + add*factor*residual - (1-add)*factor*residual.
% Both residuals are of the form (price - its firm FOC), so a POSITIVE residual means the price is
% too high and must come down: hence add=0 (subtract) for both. The rows have to be in the same
% order as the PricePath fields (r then w), because add and factor are applied to the prices
% positionally by row. setupGEnewprice3_shooting has a reorder meant to lift that requirement,
% but it matches on the wrong column of howtoupdate and so never fires.
transpathoptionsGE.GEnewprice=3;
transpathoptionsGE.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; ...
                                            'LabourMarket','w',0,0.1};

% (i) Start from the equilibrium path. The update formulae return the same prices, so the solver
% must leave the path alone.
%% One Jacobian per (n_d_GE, epsprice) cell. maxiter=1 builds the Jacobian, prints rcond, takes one
% step and stops, so each cell costs 198 path solves rather than a full solve.
n_d_GE_vals=[51,200];
epsprice_vals=[1e-5,1e-3];
for nn=1:length(n_d_GE_vals)
    n_d_GE=n_d_GE_vals(nn);
    d_grid_GE=linspace(0,1,n_d_GE)';
    % The stationary eqm has to be re-solved on each d grid, since p_eqm depends on it
    [p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d_GE, n_a_GE, n_z, 0, pi_z, d_grid_GE, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptionsGE, simoptions1, vfoptions1);
    ParamsGE=Params; ParamsGE.r=p_eqm.r; ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_d_GE,n_a_GE,n_z,d_grid_GE,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions1);
    AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_d_GE,n_a_GE,n_z,pi_z,simoptions1,ParamsGE,[]);
    PricePathGE.r=p_eqm.r*ones(1,T); PricePathGE.w=p_eqm.w*ones(1,T);
    Tbump=floor(2*T/3);
    PricePathBumped=PricePathGE;
    PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
    PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));
    for ee=1:length(epsprice_vals)
        transpathoptionsJ=transpathoptionsGE;
        transpathoptionsJ.GEnewprice=1;
        transpathoptionsJ.maxiter=1; % one Jacobian is the whole diagnostic
        transpathoptionsJ.verbose=1; % this is what prints the rcond line
        transpathoptionsJ.epsprice=epsprice_vals(ee);
        fprintf('--- Jacobian diagnostic: n_d_GE=%i, epsprice=%g --- \n',n_d_GE,epsprice_vals(ee))
        tic;
        PricePathJdiag=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsJ, simoptions1, vfoptions1, []);
        fprintf('    that cell took %2.2f seconds \n',toc)
    end
end

output=1;

end
