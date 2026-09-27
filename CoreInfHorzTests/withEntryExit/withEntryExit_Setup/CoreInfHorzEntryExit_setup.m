% Setup for the Core InfHorz entry-and-exit tests.
%
% Unlike the quasi-hyperbolic sub-bank, this one cannot reuse the parent bank's setup: entry and
% exit is a FIRM problem (a is previous employment, aprime is current employment, exit means
% shutting down and paying firing costs), so it needs its own grids, parameters and return
% functions. The shape follows Hopenhayn-Rogerson (1993).
%
% THE MODEL IS DELIBERATELY BUILT TO HAVE TWO INFEASIBLE REGIONS
% The exit decision is a BINARY weight:
%     ExitPolicy=((ReturnToExitMatrix-Vtemp)>0);
%     VKron=ExitPolicy.*ReturnToExitMatrix+(1-ExitPolicy).*Vtemp;
% so every entry is 1*x+0*y, and the only question is whether the zero-weighted side can be
% -Inf. It can, on both branches, and only in states that actually exist in the model:
%
%   (i)  staying is infeasible: a>empcap*z makes the return -Inf for EVERY aprime, so Vtemp=-Inf.
%        That is exactly why ReturnToExit-(-Inf)>0 fires and ExitPolicy=1 there. The losing
%        term is then 0*(-Inf).
%   (ii) exit is forbidden: a<minexit makes ReturnToExit=-Inf, which is exactly why the test
%        returns false and ExitPolicy=0 there. The losing term is again 0*(-Inf).
%
% If the model had neither region, every check in this bank would pass identically before and
% after any fix to that arithmetic, and the bank would prove nothing. Each subcode therefore
% prints the size of both regions, so the diary shows the checks are live.

n_z=9;   % firm productivity
n_a=101; % previous employment
n_d=5;   % utilisation (only used by the 'with d' subcodes)

a_grid=linspace(0,20,n_a)';
d_grid=linspace(0.6,1.4,n_d)';

% setup z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.1,n_z);
z_grid=exp(z_grid);

%% Parameters
Params.beta=0.96;
DiscountFactorParamNames={'beta'};

Params.alpha=0.64;  % span of control
Params.p=1;         % output price (a general eqm price in the GE subcode)
Params.tau=0.2;     % firing cost per worker
Params.cf=0.5;      % fixed cost of operating
Params.psi=0.5;     % cost of moving utilisation away from 1

% The two infeasible regions (see the note above). Both must be non-empty.
Params.empcap=15;   % staying is infeasible when a>empcap*z
Params.minexit=1;   % exit is forbidden when a<minexit

%% Entry and exit
% Distribution of new agents: entrants arrive small and unproductive. Must be a pdf.
Params.upsilon=zeros(n_a,n_z);
Params.upsilon(1:5,1:floor(0.65*n_z))=1; % small a, lower part of the z range
Params.upsilon=Params.upsilon/sum(Params.upsilon(:));

Params.Ne=0.5;      % mass of new agents (a general eqm price in the GE subcode)
Params.ce=1;        % fixed cost of entry (a general eqm price in the GE subcode)
Params.A=1;         % demand shifter, for the GE output condition
Params.zeta=[];     % conditional prob of survival; filled in as 1-ExitPolicy once V is solved

EntryExitParamNames.DistOfNewAgents={'upsilon'};
EntryExitParamNames.MassOfNewAgents={'Ne'};
EntryExitParamNames.CondlProbOfSurvival={'zeta'};

%% endogenousexit=2 extras
% exitprobabilities names the [endogenous, exogenous] legs; the toolkit prepends the residual
% no-exit leg as 1-sum(). So exitprob_endog=1 with exitprob_exog=0 gives [0,1,0], i.e. pure
% endogenous exit, which is what the exit2-vs-exit1 cross-test needs. Note that this makes two
% of the three weights EXACTLY zero, which is the point: a scalar zero times a -Inf value is
% NaN just as surely as an array one.
Params.exitprob_endog=1;
Params.exitprob_exog=0;
Params.continuationcost=0;

%% vfoptions/simoptions baselines
% Only the GPU raws are covered. The parallel=0 and parallel=1 raws exist (6 of the 11) but are
% deliberately out of scope for this bank.
vfoptionsbaseline=struct();
vfoptionsbaseline.parallel=2;
vfoptionsbaseline.endogenousexit=1;
vfoptionsbaseline.ReturnToExitFn=@(a,z,tau,minexit) ReturnToExitFn_EE(a,z,tau,minexit);

simoptionsbaseline=struct();
simoptionsbaseline.parallel=2;
simoptionsbaseline.agententryandexit=1;
simoptionsbaseline.endogenousexit=1;
