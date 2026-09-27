function output=EEInfHorz_generaleqm(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% General equilibrium with endogenous entry and endogenous exit, no d, GPU.
%
% Exit imposes no general equilibrium condition; ENTRY does, and it is a non-standard one
% because it involves the distribution of ENTRANTS rather than the distribution of existing
% agents. That is what heteroagentoptions.specialgeneqmcondn={0,'entry'} is for: the first
% condition is standard, the second is the free-entry condition.
%
% Two general eqm prices:
%   p   the output price, pinned down by the goods market condition A/Y=p
%   Ne  the mass of new entrants, pinned down by the free entry condition EValueFn=p*ce
%
% This subcode is the reason the bank covers GE at all: HeteroAgentStationaryEqm_InfHorz_EntryExit
% and its four siblings call the value fn iteration many times, from many different price
% vectors, and a value function that carries NaN in its infeasible corners will steer the
% root-finder rather than merely being wrong in a corner. A GE that converges is therefore a
% much stronger statement about the exit arithmetic than a single solve.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,p,alpha,tau,cf,empcap) ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap);

% With entry-exit the mass of the distribution matters, and FnsToEvaluate may take 'agentmass'
% as an extra input placed after the action space and before the parameters.
FnsToEvaluate.Y=@(aprime,a,z,agentmass,alpha) z*(aprime^alpha);
FnsToEvaluate.employment=@(aprime,a,z,agentmass) aprime;

vfoptions1=vfoptionsbaseline;
simoptions1=simoptionsbaseline;
simoptions1.verbose=1; % the distribution iteration then prints its distance every 50 blocks, so a
                       % stationary dist that is not converging is visible while it happens

%% A solve at the initial guess, to have something to compare the eqm against
[V0,Policy0,ExitPolicy0]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
V0=gather(V0); ExitPolicy0=gather(ExitPolicy0);
fprintf('at the initial guess, V census: %i finite, %i -Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V0(:))),sum(V0(:)==-Inf),sum(isnan(V0(:))),numel(V0));
fprintf('NaN in V at the initial guess, this should be zero: %i \n',sum(isnan(V0(:))))
% Survival is the endogenous exit decision TIMES an exogenous death hazard. Without the
% exogenous piece the stationary mass is Ne/(1-mean survival), which diverges wherever no firm
% wants to exit -- and the free-entry root sits at a price high enough that almost none do. The
% first run of the scan hung there: StationaryDist_InfHorz_Iteration_EntryExit_raw is bounded by
% simoptions.maxit, so it does not loop forever, it just grinds the whole budget without ever
% converging. deathrate>0 bounds the mass at Ne/deathrate at every price the GE can reach.
Params.deathrate=0.1;
Params.zeta=(1-ExitPolicy0)*(1-Params.deathrate);

%% Bracket the free-entry condition, and start the solver inside the bracket
% The first run of this subcode wandered to p=-0.409 and Ne=2e14 and hit MaxFunEvals. Positivity
% is now enforced (constrainpositive was wired through the EntryExit path, which had been
% ignoring it), but that alone does not say whether a root EXISTS. The scan does, and it also
% gives the solver a start inside the bracket rather than a guess.
%
% Free entry is EValueFn-p*ce. Firm profit is p*z*n^alpha-n-cf, so the value of a firm rises
% faster than linearly in p while the entry cost rises linearly: the condition is negative at low
% p and positive at high p, and the root is where entry just breaks even.
fprintf('\nscanning the free-entry condition in p (it must change sign for the GE to have a root) \n')
pscan=[0.5 1 1.5 2 3 4 6];
entryscan=nan(size(pscan));
% The exit share is printed alongside, because it decides whether endogenousexit=1 is usable here
% at all. Under =1 survival is exactly 1-ExitPolicy, so the stationary mass is Ne/(1-mean
% survival): if exit collapses as p rises, the mass diverges (or converges far too slowly to be
% usable) exactly where the free-entry root sits, and the =1 GE has no workable equilibrium. Two
% forces pull against each other and the sign is not obvious a priori -- higher p makes firms more
% profitable, which discourages exit, but it also makes them hire more, which pushes them into the
% a>empcap*z region where staying is infeasible and exit is FORCED. This measures which wins.
exitscan=nan(size(pscan));
for ii=1:length(pscan)
    Pscan=Params; Pscan.p=pscan(ii);
    [Vscan,~,ExitPolicyscan]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Pscan,DiscountFactorParamNames,[],vfoptions1);
    entryscan(ii)=gather(sum(Params.upsilon(:).*Vscan(:)))-pscan(ii)*Params.ce;
    exitscan(ii)=mean(gather(ExitPolicyscan(:)));
    fprintf('  p=%5.2f   EValueFn-p*ce = %+.4f   share of states that exit = %.4f \n',pscan(ii),entryscan(ii),exitscan(ii))
end
fprintf('  (the share above floors at %i/%i=%.4f, which is exactly the a>empcap*z region where exit is\n',sum(sum(repmat(a_grid,1,N_z)>Params.empcap*repmat(gather(z_grid)',N_a,1))),N_a*N_z,sum(sum(repmat(a_grid,1,N_z)>Params.empcap*repmat(gather(z_grid)',N_a,1)))/(N_a*N_z))
fprintf('   FORCED: above p=2 voluntary exit has vanished and only the mechanical kind is left) \n')
fprintf('the free-entry condition has no NaN across the scan, this should be zero: %i \n',sum(isnan(entryscan)))
iib=find(entryscan(1:end-1)<0 & entryscan(2:end)>0,1);
fprintf('the free-entry condition changes sign in the scanned range, this should be one: %i \n',~isempty(iib))
if ~isempty(iib)
    Params.p=0.5*(pscan(iib)+pscan(iib+1));
    fprintf('  sign change between p=%g and p=%g, so starting the solver at p=%g \n',pscan(iib),pscan(iib+1),Params.p)
    % Now pick Ne to satisfy the goods market at that p. Aggregate Y is proportional to the mass
    % of firms, so solve once at Ne=1 and rescale: A/Y=p wants Ne = A/(p*Y_at_Ne1).
    Pne=Params; Pne.Ne=1;
    [Vne,Policyne,ExitPolicyne]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Pne,DiscountFactorParamNames,[],vfoptions1); %#ok<ASGLU>
    Pne.zeta=(1-gather(ExitPolicyne))*(1-Params.deathrate);
    Dne=StationaryDist_InfHorz(Policyne,0,n_a,n_z,pi_z,simoptions1,Pne,EntryExitParamNames);
    AVne=EvalFnOnAgentDist_AggVars_InfHorz(Dne,Policyne,FnsToEvaluate,Pne,[],0,n_a,n_z,[],a_grid,z_grid,simoptions1,EntryExitParamNames);
    Params.Ne=Params.A/(Params.p*gather(AVne.Y.Aggregate));
    fprintf('  and at Ne=1 that gives Y=%g, so starting Ne=%g \n',gather(AVne.Y.Aggregate),Params.Ne)

    % What actually governs whether endogenousexit=1 is usable here is the exit rate BY MASS, not
    % the share of STATES that can exit. Stationary mass balance is Ne = mass*(exit rate by mass),
    % so the distribution iteration needs roughly 1/(exit rate) periods; StationaryDist_InfHorz
    % defaults to maxit=10^6 outer blocks of 100 steps, so a tiny exit rate means it grinds
    % essentially forever rather than stopping. Measured with BINARY survival, i.e. what the GE
    % subfn actually uses (it recomputes CondlProbOfSurvival=1-ExitPolicy every evaluation and
    % discards any exogenous hazard set here).
    simoptionsdiag=simoptions1; simoptionsdiag.maxit=20000; % maxit counts ITERATIONS, so 20,000 periods
    Pbin=Pne; Pbin.zeta=1-gather(ExitPolicyne);
    Dbin=StationaryDist_InfHorz(Policyne,0,n_a,n_z,pi_z,simoptionsdiag,Pbin,EntryExitParamNames);
    exitbymass=gather(sum(Dbin.pdf(:).*reshape(gather(ExitPolicyne),[],1)));
    fprintf('  DIAGNOSTIC at p=%g with BINARY survival (no death rate): exit rate by mass = %.3e \n',Params.p,exitbymass)
    % Note the non-contraction guard inside StationaryDist_InfHorz_Iteration_EntryExit_raw will
    % normally stop this early (after ~1000 periods), so the mass below is whatever it had reached
    % by then. With Ne=1 and zero outflow the mass just equals the number of periods elapsed.
    fprintf('    implied periods to converge ~ %.3g, mass reached before stopping = %g \n',1/max(exitbymass,eps),gather(Dbin.mass))
    fprintf('    (for comparison, with the death rate the same solve converged and gave Y=%g) \n',gather(AVne.Y.Aggregate))
end

%% General equilibrium
GEPriceParamNames={'p','Ne'};

heteroagentoptions=struct();
heteroagentoptions.specialgeneqmcondn={0,'entry'};
heteroagentoptions.verbose=2; % 2 prints the prices BEFORE each evaluation, so a slow one is visible while it runs
% Both prices are economically positive, and the first run of this subcode wandered to p=-0.409.
% The EntryExit path ignored constrainpositive until it was wired through; this is the check that
% it now bites, alongside the 'general eqm prices are positive' line below.
heteroagentoptions.constrainpositive={'p','Ne'};

GeneralEqmEqns.RealOutput=@(Y,p,A) A/Y-p;            % goods market: price pinned by demand
GeneralEqmEqns.Entry=@(EValueFn,ce,p) EValueFn-p*ce; % free entry: expected value of entry is zero

n_p=0;
fprintf('\nsolving for the stationary general equilibrium \n')
[p_eqm,p_eqm_index,GeneralEqmCondition]=HeteroAgentStationaryEqm_InfHorz(0,n_a,n_z,n_p,pi_z,[],a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],GEPriceParamNames,heteroagentoptions,simoptions1,vfoptions1,EntryExitParamNames); %#ok<ASGLU>

fprintf('eqm p  = %g \n',p_eqm.p)
fprintf('eqm Ne = %g \n',p_eqm.Ne)
fprintf('general eqm prices are finite, this should be zero: %i \n',sum(~isfinite([p_eqm.p,p_eqm.Ne])))
fprintf('general eqm prices are positive, this should be zero: %i \n',sum([p_eqm.p,p_eqm.Ne]<=0))
% GeneralEqmCondition comes back as a vector in some versions and a struct in others, so flatten
% it before reporting rather than assuming one shape and aborting the bank on the other.
if isstruct(GeneralEqmCondition)
    GEC=cell2mat(struct2cell(GeneralEqmCondition));
else
    GEC=GeneralEqmCondition;
end
GEC=gather(GEC(:));
fprintf('general eqm conditions, these should be close to zero: %s \n',mat2str(GEC',4))
fprintf('max abs general eqm condition, this should be close to zero: %.3e \n',max(abs(GEC)))

%% Re-solve at the equilibrium prices and check everything is consistent
Params.p=p_eqm.p;
Params.Ne=p_eqm.Ne;
[Ve,Policye,ExitPolicye]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
Ve=gather(Ve); ExitPolicye=gather(ExitPolicye);
fprintf('\nat the eqm, V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(Ve(:))),sum(Ve(:)==-Inf),sum(Ve(:)==Inf),sum(isnan(Ve(:))),numel(Ve));
fprintf('NaN in V at the eqm, this should be zero: %i \n',sum(isnan(Ve(:))))
fprintf('ExitPolicy is binary at the eqm, this should be zero: %i \n',sum(~(ExitPolicye(:)==0 | ExitPolicye(:)==1)))

Params.zeta=(1-ExitPolicye)*(1-Params.deathrate);
StationaryDiste=StationaryDist_InfHorz(Policye,0,n_a,n_z,pi_z,simoptions1,Params,EntryExitParamNames);
fprintf('StationaryDist.pdf sums to one, this should be zero: %.3e \n',abs(sum(gather(StationaryDiste.pdf(:)))-1))
fprintf('StationaryDist.pdf has no NaN, this should be zero: %i \n',sum(isnan(gather(StationaryDiste.pdf(:)))))
fprintf('mass of existing agents at the eqm: %g \n',gather(StationaryDiste.mass))

AggVarse=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDiste,Policye,FnsToEvaluate,Params,[],0,n_a,n_z,[],a_grid,z_grid,simoptions1,EntryExitParamNames);
fprintf('eqm Y = %g, eqm employment = %g \n',AggVarse.Y.Aggregate,AggVarse.employment.Aggregate)

% The goods market condition, checked directly rather than through the root-finder's own report
fprintf('goods market A/Y-p at the eqm, this should be close to zero: %.3e \n',abs(Params.A/AggVarse.Y.Aggregate-Params.p))

%% Plot the exit decision at the equilibrium
fig=figure(figure_c); %#ok<NASGU>
surf(ExitPolicye)
title('Exit decision at the general eqm (1 indicates exit)'); xlabel('z'); ylabel('a')

%%
output=struct();
output.p_eqm=p_eqm;

end
