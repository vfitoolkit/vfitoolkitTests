function output=EEInfHorz_nod(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% endogenousexit=1, no d, GPU. Reaches ValueFnIter_InfHorz_EndogExit_nod_Par2_raw.
%
% This is the fullest subcode in the bank: the census and the two live-region counts are set up
% here and the later subcodes reuse the pattern more briefly.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,p,alpha,tau,cf,empcap) ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap);

FnsToEvaluate.employment=@(aprime,a,z) aprime;
FnsToEvaluate.output=@(aprime,a,z,alpha) z*(aprime^alpha);

%% Baseline
vfoptions1=vfoptionsbaseline;
simoptions1=simoptionsbaseline;
[V1,Policy1,ExitPolicy1]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
V1=gather(V1); Policy1=gather(Policy1); ExitPolicy1=gather(ExitPolicy1);

%% Are the two infeasible regions actually present?
% Every check below is vacuous if they are not: the exit arithmetic only misbehaves when a zero
% weight lands on a -Inf, and these are the only two ways that happens.
avals=repmat(a_grid,1,N_z);
zvals=repmat(gather(z_grid)',N_a,1);
n_stayinfeasible=sum(sum(avals>Params.empcap*zvals));
n_exitforbidden=sum(sum(avals<Params.minexit));
fprintf('states where staying is infeasible (a>empcap*z), this should be non-zero: %i of %i \n',n_stayinfeasible,N_a*N_z)
fprintf('states where exit is forbidden (a<minexit), this should be non-zero: %i of %i \n',n_exitforbidden,N_a*N_z)

%% Census of V
% max(abs(A-B)) IGNORES NaN, so every 'this should be zero' below is silently vacuous wherever
% a NaN sits. The census is the only instrument that sees the NaN class, so it comes first.
fprintf('V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V1(:))),sum(V1(:)==-Inf),sum(V1(:)==Inf),sum(isnan(V1(:))),numel(V1));
fprintf('NaN in V, this should be zero: %i \n',sum(isnan(V1(:))))
fprintf('ExitPolicy is binary, this should be zero: %i \n',sum(~(ExitPolicy1(:)==0 | ExitPolicy1(:)==1)))
fprintf('NaN in ExitPolicy, this should be zero: %i \n',sum(isnan(ExitPolicy1(:))))

%% Howards on/off
% Howards improvement is just an accelerator, so turning it off must give the same answer. This
% is the ONLY check in the bank that reaches the Howards branch, and that branch carries two of
% the three exit-weighted expressions in the raw (the Ftemp recombination and the
% beta*(1-ExitPolicy).*EVKrontemp continuation). Without this check most of the file is
% unexercised.
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1noH,Policy1noH,ExitPolicy1noH]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('howards=0 (pure VFI), this should be zero: V %.3e \n',max(abs(V1(:)-gather(V1noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: Policy %.3e \n',max(abs(Policy1(:)-gather(Policy1noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: ExitPolicy %.3e \n',max(abs(ExitPolicy1(:)-gather(ExitPolicy1noH(:)))))
fprintf('howards=0 V census: %i finite, %i -Inf, %i NaN (of %i) \n', ...
    sum(isfinite(gather(V1noH(:)))),sum(gather(V1noH(:))==-Inf),sum(isnan(gather(V1noH(:)))),numel(V1noH));

%% V0-independence
% The iteration is a contraction, so the converged V cannot depend on the starting guess.
vfoptions1_V0=vfoptions1;
vfoptions1_V0.V0=-10*ones(N_a,N_z);
[V1_V0,~,ExitPolicy1_V0]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_V0);
fprintf('different V0 gives the same V, this should be zero: %.3e \n',max(abs(V1(:)-gather(V1_V0(:)))))
fprintf('different V0 gives the same ExitPolicy, this should be zero: %.3e \n',max(abs(ExitPolicy1(:)-gather(ExitPolicy1_V0(:)))))

%% keeppolicyonexit
% The default (0) zeroes Policy wherever the firm exits; keeppolicyonexit=1 keeps the policy
% that would have been chosen. The two must agree everywhere the firm does NOT exit, and the
% default must be zero everywhere it does.
vfoptions1_keep=vfoptions1;
vfoptions1_keep.keeppolicyonexit=1;
[~,Policy1keep,ExitPolicy1keep]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_keep);
Policy1keep=gather(Policy1keep); ExitPolicy1keep=gather(ExitPolicy1keep);
stay=(ExitPolicy1keep==0);
P1=reshape(Policy1,[N_a,N_z]); P1k=reshape(Policy1keep,[N_a,N_z]);
fprintf('keeppolicyonexit agrees where the firm stays, this should be zero: %.3e \n',max(abs(P1(stay)-P1k(stay))))
fprintf('default Policy is zero where the firm exits, this should be zero: %.3e \n',max(abs(P1(~stay))))
fprintf('keeppolicyonexit does not change ExitPolicy, this should be zero: %.3e \n',max(abs(ExitPolicy1(:)-ExitPolicy1keep(:))))

%% Stationary distribution with entry and exit
Params.zeta=1-ExitPolicy1;
StationaryDist1=StationaryDist_InfHorz(Policy1,0,n_a,n_z,pi_z,simoptions1,Params,EntryExitParamNames);
fprintf('StationaryDist.pdf sums to one, this should be zero: %.3e \n',abs(sum(gather(StationaryDist1.pdf(:)))-1))
fprintf('StationaryDist.pdf has no negative mass, this should be zero: %i \n',sum(gather(StationaryDist1.pdf(:))<0))
fprintf('StationaryDist.pdf has no NaN, this should be zero: %i \n',sum(isnan(gather(StationaryDist1.pdf(:)))))
fprintf('StationaryDist.mass is finite and positive, this should be one: %i \n',isfinite(gather(StationaryDist1.mass)) && gather(StationaryDist1.mass)>0)
fprintf('mass of existing agents: %g \n',gather(StationaryDist1.mass))

%% Statistics
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1,FnsToEvaluate,Params,[],0,n_a,n_z,[],a_grid,z_grid,simoptions1,EntryExitParamNames);
fprintf('AggVars employment %g, output %g \n',AggVars1.employment.Aggregate,AggVars1.output.Aggregate)
fprintf('AggVars has no NaN, this should be zero: %i \n',sum(isnan([AggVars1.employment.Aggregate,AggVars1.output.Aggregate])))

% NOT TESTED HERE: AllStats/AutoCorrTransProbs/CrossSectionCovarCorr on an entry-exit model.
% EvalFnOnAgentDist_AllStats_InfHorz has no EntryExitParamNames argument and no _Mass subcode, so
% there is nothing to call: it cannot mass-weight, and it cannot cope with the Policy==0 exit
% sentinel either. Only AggVars, ProbDensityFn and ValuesOnGrid have entry-exit versions.
fprintf('NOT TESTED: AllStats/AutoCorr/CrossSectionCovarCorr have no entry-exit support \n')

% What AggVars actually returns on the entry-exit path, read off the code rather than assumed:
% EvalFnOnAgentDist_AggVars_InfHorz_Mass zeroes the pdf at the exiting agents, forms
% sum(Values.*pdf), and then its last line multiplies the lot by StationaryDist.mass. So the
% field is named .Aggregate (not .Mean, as on the non-entry-exit path) and equals
%     mass * survivorshare * (mean over survivors).
% Undo both factors before comparing with anything drawn from a panel of survivors.
survivorshare=gather(sum(StationaryDist1.pdf(:).*Params.zeta(:)));
fprintf('survivor share of the stationary distribution: %g \n',survivorshare)
percapita=@(agg) agg/(gather(StationaryDist1.mass)*survivorshare);

%% SimPanelValues
% entryinpanel=0 follows the entry-exit convention: the panel then follows survivors only.
% Monte Carlo, so these are roughly-equal checks and never an exact zero.
simoptionsPanel=simoptions1;
simoptionsPanel.entryinpanel=0;
simoptionsPanel.numbersims=10^4;
simoptionsPanel.simperiods=20;
SimPanel1=SimPanelValues_InfHorz(StationaryDist1,Policy1,FnsToEvaluate,[],Params,0,n_a,n_z,[],a_grid,z_grid,pi_z,simoptionsPanel,EntryExitParamNames);
fprintf('SimPanelValues ran; period-1 mean employment vs AggVars per-survivor mean (Monte Carlo) \n')
[percapita(AggVars1.employment.Aggregate), mean(SimPanel1.employment(1,:),'omitnan')]

%% Plot the exit decision
fig=figure(figure_c); %#ok<NASGU>
surf(ExitPolicy1)
title('Exit decision, no d (1 indicates exit)'); xlabel('z'); ylabel('a')

%%
output=struct();

end
