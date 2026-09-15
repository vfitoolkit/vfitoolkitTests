function output=CoreInfHorzTPathAlgo_d_z_noe_nosemiz_with2A(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE,figure_c)
% Compares the transition-path algorithms on one model: GEnewprice=1 (quasi-Newton with Broyden) and
% GEnewprice=2 (Anderson acceleration) against GEnewprice=3 (shooting), the last with and without
% the additionalfactor ramp. All of them solve the same
% null-reform general eqm path from the same bumped starting guess, and all three are timed.
% The model and its stationary general eqm are lifted from CoreInfHorzTPath_d_z_noe_nosemiz_with2A.

% Two-endogenous-state InfHorz TPath test, mirroring the fig-4 (d, z) battery.
% Model: simplified Kitao (2008) with Bruggemann (2021) endogenous labor supply, no taxes --
% the same model CoreInfHorzTests uses for its d 2A subcode.
%   Endogenous states: asset, occupation (0=worker, 1=entrepreneur)
%   Decision variable:  l (labor supply; entrepreneurs must supply exactly lbar)
%   Exogenous markov z: eta (labor productivity), theta (entrepreneurial ability)
% n_a(2)=2 is what triggers the DC2A/GI2A code paths; the decision variable additionally
% routes through the Refine variants of them.
%
% Two things differ from the 1A subcodes and are worth knowing before reading the checks:
%
% (1) V is -Inf at infeasible states (an entrepreneur with near-zero assets has no feasible
%     positive consumption). -Inf minus -Inf is NaN, and max() omits NaN, so the comparisons
%     below silently ignore states where both sides agree on -Inf -- which is what we want. A
%     state where only ONE side is -Inf gives +-Inf, so max() reports Inf and the check still
%     fails loudly. So no special handling is needed, but these lines are not doing quite what
%     they appear to at a glance.
%
% (2) The plain tier (no DC, no GI) works for 2A as it stands: pure discretization is generic
%     in N_a=prod(n_a), so it does not care how many endogenous states there are. The
%     GI-without-DC tier (vfoptions3) does NOT. ValueFnIter_InfHorz_TPath_SingleStep_GI1_nod_raw
%     assumes N_a==length(a_grid), which is false for two endogenous states: here N_a=202 while
%     the stacked a_grid has 103 entries, so its interp1 on line 17 errors. The dispatcher sends
%     gridinterplayer==1 with divideandconquer==0 to that 1A-only raw with no length(n_a) check.
%     A GI2A-without-DC raw is the missing piece -- the stationary solver already has one
%     (ValueFnIter_InfHorz_preGI2A_* / _postGI2A_*), so the gap is in the TPath tier only.

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();

vfoptions.verbose_advice=0; % 2-endo GI would otherwise sound the postGI advice on every solve

% Terminal parameters are the end of the path (same convention as the 1A subcodes)
Params.r=PricePath.r(T);
Params.w=ParamPath.w(T);

ReturnFn=@(l,aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi,xi,sigma2,lbar) ...
    ReturnFn_d_z_noe_nosemiz_with2A(l,aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi,xi,sigma2,lbar);

% Setup some FnsToEvaluate (functions of (l,aprime,eprime,a,e,eta,theta) -- decision variable first)
FnsToEvaluate.assets=@(l,aprime,eprime,a,e,eta,theta) a;
FnsToEvaluate.entrepreneur=@(l,aprime,eprime,a,e,eta,theta) e; % fraction who are entrepreneurs
FnsToEvaluate.nextassets=@(l,aprime,eprime,a,e,eta,theta) aprime; % aprime-dependent, so it tests the policy decode directly
FnsToEvaluate.laborsupply=@(l,aprime,eprime,a,e,eta,theta) l; % d-dependent, so it tests the decision-variable decode
FnsToEvaluate.etaval=@(l,aprime,eprime,a,e,eta,theta) eta; % first markov shock, this is the one KFTT makes period-dependent
FnsToEvaluate.thetaval=@(l,aprime,eprime,a,e,eta,theta) theta; % second markov shock, left alone, so it acts as a control


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
% Corporate sector, alongside the entrepreneurs. Kitao/Bruggemann structure: households' total
% assets are split between capital tied up in entrepreneurs' own businesses and capital rented to
% the corporate sector, and the labour entrepreneurs hire comes out of the same pool the corporate
% sector draws on. So the corporate factor inputs are RESIDUALS, and it is those -- not total
% assets and total labour -- that enter the firm's first order conditions.
% firmdelta is set equal to the entrepreneurs' depreciation rate so there is one depreciation
% rate in the model.
Params.firmalpha=0.36;
Params.firmdelta=Params.delta;
Params.firmA=1;

FnsToEvaluateGE.A=@(l,aprime,eprime,a,e,eta,theta) a; % total household assets
FnsToEvaluateGE.K_noncorp=@(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi) EntreCapital_d_with2A(l,aprime,eprime,a,e,theta,r,w,delta,upsilon1,upsilon2,leverage,phi); % capital inside entrepreneurs' own businesses
FnsToEvaluateGE.L_workers=@(l,aprime,eprime,a,e,eta,theta) (e==0)*l*eta; % labour supplied to the market (entrepreneurs' own eta goes to their own firm)
FnsToEvaluateGE.N_hired=@(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi) EntreHiredLabour_d_with2A(l,aprime,eprime,a,e,eta,theta,r,w,delta,upsilon1,upsilon2,leverage,phi); % labour entrepreneurs hire in

% Intermediate eqns are evaluated after the AggVars are written into Parameters and before the
% GeneralEqmEqns, with each result written back as a named parameter (see the useintermediateEqns
% block in TransitionPath_InfHorz_shooting / HeteroAgentStationaryEqm_InfHorz). That is exactly
% what is needed to form the corporate residuals once, and use them in both GE conditions.
intermediateEqnsGE.CorporateCapital=@(A,K_noncorp) A-K_noncorp;
intermediateEqnsGE.CorporateLabour=@(L_workers,N_hired) L_workers-N_hired;

GeneralEqmEqnsGE.CapitalMarket=@(r,CorporateCapital,CorporateLabour,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(CorporateCapital^(firmalpha-1))*(CorporateLabour^(1-firmalpha))-firmdelta);
GeneralEqmEqnsGE.LabourMarket=@(w,CorporateCapital,CorporateLabour,firmalpha,firmA) w-((1-firmalpha)*firmA*(CorporateCapital^firmalpha)*(CorporateLabour^(-firmalpha)));

heteroagentoptionsGE=struct(); % default fminalgo, and no constraints on r or w
heteroagentoptionsGE.intermediateEqns=intermediateEqnsGE;

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
transpathoptionsGE.intermediateEqns=intermediateEqnsGE; % same corporate residuals as the stationary solve

% (i) Start from the equilibrium path. The update formulae return the same prices, so the solver
% must leave the path alone.
%% A bumped starting guess, so all three algorithms have to actually travel
Tbump=floor(2*T/3);
bumpr=0.01+0.04*linspace(0,1,Tbump);
bumpw=0.05-0.04*linspace(0,1,Tbump);
PricePathBumped=PricePathGE;
PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+bumpr);
PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+bumpw);
fprintf('Bumped starting guess, max initial deviation from p_eqm: r %2.8f, w %2.8f \n',max(abs(PricePathBumped.r-p_eqm.r)),max(abs(PricePathBumped.w-p_eqm.w)))

%% (B) GEnewprice=3: the shooting algorithm, no ramp
transpathoptionsB=transpathoptionsGE;
transpathoptionsB.maxiter=250;
tic;
[PricePathB,GEcondnPathB]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsB, simoptions1, vfoptions1, []);
time_B=toc;

%% (C) GEnewprice=3 again, this time with the additionalfactor ramp
transpathoptionsC=transpathoptionsB;
% f_add=0.5 rather than the 3 the models without d use. On this model f_add=3 ramps the step up until
% the path leaves the region where the model can be solved, and the general eqm conditions come back
% NaN (which TransitionPath_InfHorz_shooting now errors on rather than returning). f_add<1 ramps the
% step the other way, damping it down over iterations, so it can never reach anywhere plain shooting
% (B) does not, and B survives its full maxiter here. The ramp code path is exercised either way.
% Note that the with/without-ramp comparison below says little on this model while B itself does not
% converge: both stop at maxiter, so they stop at different points rather than at the same answer.
transpathoptionsC.GEnewprice3.additionalfactor=[0.5,10,30];
tic;
[PricePathC,GEcondnPathC]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsC, simoptions1, vfoptions1, []);
time_C=toc;

%% (F) GEnewprice=2: Anderson acceleration of the same shooting map as (B)
% Anderson mixes the last m iterates of the shooting map to take a quasi-Newton-like step, without
% ever building a Jacobian. So an iteration costs the same as a shooting iteration and the hope is
% simply that there are fewer of them. howtoupdate is taken from (B) verbatim: it is the same map,
% all that differs is how the iterates are combined.
transpathoptionsF=rmfield(transpathoptionsB,'GEnewprice3');
transpathoptionsF.GEnewprice=2;
transpathoptionsF.GEnewprice2.howtoupdate=transpathoptionsB.GEnewprice3.howtoupdate;

% (F1) memory=0 empties the history as fast as it fills it, so every step is a plain shooting step
% and Anderson is then exactly the shooting algorithm: same map, same oldpathweight=0, no ramp, same
% convergence distance. This is what checks the GEnewprice=2 plumbing (howtoupdate parsed through the
% new shooting field, the permute and keepold in updatePricePathNew_TPath_T, the flattening and
% unflattening of the path, the terminal row staying put), independently of whether Anderson
% accelerates anything.
% Both are run for a fixed small number of iterations, deliberately short of convergence, so that
% each does the same number of evaluate-then-step cycles and the two paths must then agree to machine
% precision. Comparing the two converged solves instead does NOT give zero, and that is not a bug in
% either: TransitionPath_InfHorz_shooting tests its distance at the top of the loop but still applies
% one more update at the bottom before exiting, so it returns a path one shooting step past the one it
% evaluated (and returns that earlier path's GEcondnPath alongside it), whereas Anderson returns the
% point it actually evaluated. The gap is exactly factor*(final GE condn), which is what the full-run
% comparison showed.
transpathoptionsF1=transpathoptionsF;
transpathoptionsF1.maxiter=10;
transpathoptionsF1.anderson.memory=0;
[PricePathF0,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsF1, simoptions1, vfoptions1, []);
transpathoptionsB1=transpathoptionsB;
transpathoptionsB1.maxiter=10;
[PricePathB1,~]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsB1, simoptions1, vfoptions1, []);

% (F2) And now with a history to work with. memory is the option that really changes the algorithm.
% safeguard=1 buys the fallback to a plain shooting step whenever an Anderson step makes the
% distance worse, at the cost of a second path solve on every Anderson step.
% type='I' is the stabilized Type-I of Zhang, O'Donoghue & Boyd (2020), a different algorithm rather
% than a different setting, and this is the only place it gets exercised. On the first run it
% extrapolated to a price path where the GE conditions come back NaN and errored out, taking the rest
% of the bank with it: its own safeguard is about the decay rate of the fixed-point residual, not
% about staying where the model can be solved. AndersonAcceleration now halves such a step back
% toward the previous iterate (anderson.maxbacktrack) and reinitialises H rather than dividing by
% zero in the rank-one update, so this row is what checks that.
memopts=[3,5,10,5,5];
safeopts=[0,0,0,1,0];
typeopts={'II','II','II','II','I'};
nanderson=length(memopts);
time_F=zeros(1,nanderson); distF=zeros(1,nanderson);
devF_r=zeros(1,nanderson); devF_w=zeros(1,nanderson);
andersonname=cell(1,nanderson);
for mm=1:nanderson
    andersonname{mm}=sprintf('memory=%2i, safeguard=%i, type=%s',memopts(mm),safeopts(mm),typeopts{mm});
    transpathoptionsF.anderson.memory=memopts(mm);
    transpathoptionsF.anderson.safeguard=safeopts(mm);
    transpathoptionsF.anderson.type=typeopts{mm};
    fprintf('--- Anderson acceleration, %s --- \n',andersonname{mm})
    tic;
    [PricePathF,GEcondnPathF]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsF, simoptions1, vfoptions1, []);
    time_F(mm)=toc;
    distF(mm)=max(sqrt(GEcondnPathF.CapitalMarket.^2+GEcondnPathF.LabourMarket.^2));
    devF_r(mm)=max(abs(PricePathF.r-PricePathB.r));
    devF_w(mm)=max(abs(PricePathF.w-PricePathB.w));
end

%% (A) GEnewprice=1: Newton on the whole price path, run once for each way of regularising the step.
% The Jacobian of a price path is badly conditioned, because neighbouring columns are nearly parallel,
% so a plain J\f divides the least determined direction by a tiny singular value and the step is
% dominated by noise. The three options differ in what they do about that: minimumnorm drops the
% directions the Jacobian does not resolve, TikhonovRegularisation fades them out smoothly, and
% stepcap keeps the plain step but limits how far it may move the path. Run last, after the cheap
% shooting solves, so their answers are already in the diary if a Newton solve fails.
jacopts={'LudwigPath','LudwigStationary','FullJacobian','LudwigSSJ'}; % Ludwig (2007) GSQN, Omega kron I at nPrices path solves, against the brute-force finite-difference Jacobian at (T-1)*nPrices
regopts={'minimumnorm','TikhonovRegularisation','stepcap'};
ncombo=length(jacopts)*length(regopts);
distA=zeros(1,ncombo); time_A=zeros(1,ncombo); devA_r=zeros(1,ncombo); devA_w=zeros(1,ncombo); comboname=cell(1,ncombo);
cc=0;
for jj=1:length(jacopts)
    for rr=1:length(regopts)
        cc=cc+1;
        comboname{cc}=[jacopts{jj},' / ',regopts{rr}];
        transpathoptionsA=transpathoptionsGE;
        transpathoptionsA.maxiter=100;
        transpathoptionsA.GEnewprice=1;
        transpathoptionsA.verbose=1; % quasi-Newton prints a few lines per iteration. Do NOT turn this
        % on for the shooting solves above: with GEnewprice=3 verbose dumps the whole
        % [PricePathOld,PricePathNew] matrix on every iteration, which would be hundreds of screens.
        transpathoptionsA.GEnewprice1.Jacobianmethod=jacopts{jj};
        transpathoptionsA.GEnewprice1.BroydenRegularisation=regopts{rr};
        fprintf('--- quasi-Newton, Jacobianmethod=%s, BroydenRegularisation=%s --- \n',jacopts{jj},regopts{rr})
        tic;
        [PricePathA,GEcondnPathA]=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsA, simoptions1, vfoptions1, []);
        time_A(cc)=toc;
        distA(cc)=max(sqrt(GEcondnPathA.CapitalMarket.^2+GEcondnPathA.LabourMarket.^2));
        devA_r(cc)=max(abs(PricePathA.r-PricePathB.r));
        devA_w(cc)=max(abs(PricePathA.w-PricePathB.w));
    end
end

%% Each algorithm must reach its own stated tolerance
distB=max(sqrt(GEcondnPathB.CapitalMarket.^2+GEcondnPathB.LabourMarket.^2));
distC=max(sqrt(GEcondnPathC.CapitalMarket.^2+GEcondnPathC.LabourMarket.^2));
for cc=1:ncombo
    fprintf('quasi-Newton (%-34s), GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',comboname{cc},distA(cc),1e-4)
end
fprintf('Shooting,                              GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',distB,1e-4)
fprintf('Shooting+ramp,                         GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',distC,1e-4)
for mm=1:nanderson
    fprintf('Anderson (%-31s), GE condn distance: %2.10f, should be at or below toleranceGEcondns=%2.10f \n',andersonname{mm},distF(mm),1e-4)
end

%% And they must agree on where the equilibrium path is. Not exactly: each stops the first iteration
% its own residual drops under the tolerance, so they land at different points inside it. Expect
% differences of the order of toleranceGEcondns in r, and larger in w, which moves only because K does.
for cc=1:ncombo
    fprintf('quasi-Newton (%-34s) vs shooting, same eqm path, should be small, r: %2.10f, w: %2.10f \n',comboname{cc},devA_r(cc),devA_w(cc))
end
fprintf('Shooting with/without ramp, same eqm path, should be small, r: %2.10f, w: %2.10f \n',max(abs(PricePathC.r-PricePathB.r)),max(abs(PricePathC.w-PricePathB.w)))
fprintf('Anderson with memory=0 vs shooting, both at maxiter=10, the same algorithm, should be zero, r: %.3e, w: %.3e \n',max(abs(PricePathF0.r-PricePathB1.r)),max(abs(PricePathF0.w-PricePathB1.w)))
for mm=1:nanderson
    fprintf('Anderson (%-31s) vs shooting, same eqm path, should be small, r: %2.10f, w: %2.10f \n',andersonname{mm},devF_r(mm),devF_w(mm))
end
for cc=1:ncombo
    fprintf('Runtime, quasi-Newton (%-34s): %2.2f seconds \n',comboname{cc},time_A(cc))
end
fprintf('Runtime, shooting %2.2f seconds, shooting+ramp %2.2f seconds \n',time_B,time_C)
for mm=1:nanderson
    fprintf('Runtime, Anderson (%-31s): %2.2f seconds \n',andersonname{mm},time_F(mm))
end

%% (D) The order of the howtoupdate rows must not matter
% The rows above are given in PricePath order (r then w). setupGEnewprice3_shooting reorders them into
% that order itself, so giving them in the other order has to produce the same path. One iteration is
% enough to catch it: if the reorder does not work, the very first update pairs each price with the
% other price's residual. Both rows use the same add and factor, so row order changes nothing except
% which general eqm condition drives which price.
transpathoptionsD=transpathoptionsGE;
transpathoptionsD.maxiter=1;
PricePathD1=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsD, simoptions1, vfoptions1, []);
transpathoptionsD.GEnewprice3.howtoupdate={'LabourMarket','w',0,0.1; ...
                                          'CapitalMarket','r',0,0.1}; % the other order
PricePathD2=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_d_GE, n_a_GE, n_z, d_grid_GE,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsD, simoptions1, vfoptions1, []);
fprintf('howtoupdate rows in price order vs eqn order, one iteration, this should be zero, r: %2.10f, w: %2.10f \n',max(abs(PricePathD2.r-PricePathD1.r)),max(abs(PricePathD2.w-PricePathD1.w)))

output=1;

end
