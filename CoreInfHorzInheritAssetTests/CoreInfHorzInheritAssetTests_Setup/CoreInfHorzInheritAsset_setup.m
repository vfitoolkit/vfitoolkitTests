% Setup for the Core InfHorz InheritanceAsset tests
%
% Uses the same d, a, z in every model that needs them, so results are comparable.
%
% THE MODEL SHAPE IS FORCED BY THE TOOLKIT
% ValueFnIter_InfHorz_InheritAsset has exactly one live branch; every other combination hits an
% 'Have not yet implemented' error. That branch is: d1 present, no a1, z present, no e. So there
% is no tier cube here, and n_d must have at least two elements (d1 and d2) and n_a exactly one.
%
% The inheritance asset is a Case-2 style object: next period's a2 is a function of the decision
% d2 and of the shock TRANSITION (z,zprime), and does not depend on this period's a2. So the
% return function takes (d1,d2,a,z) with NO aprime, and aprimeFn takes (d2,z,zprime,params).
%
% TWO DELIBERATE DESIGN CHOICES, both so that checks are not vacuous:
%
% 1. n_d2 (11) is NOT equal to n_a (101). The refine step stores d1star as [N_d2,N_a,N_z] and
%    then recovers d1 by linear indexing, so any confusion between the two sizes is invisible
%    whenever they happen to be equal. Keeping them different is what makes the check bite.
%
% 2. The aprimeFn is built to push a2prime off BOTH ends of a2_grid. Off the top, the toolkit
%    sets the lower-point probability to exactly 0; off the bottom, to exactly 1. Those exact
%    0/1 weights are the only thing that reaches the zero-weight guard in the raw
%    (a zero weight against a -Inf node used to give 0*(-Inf)=NaN). If a2prime never left the
%    grid, the guard would never fire and the NaN census below would prove nothing.
%    The setup prints the count of exact 0s and exact 1s so the diary shows the checks are live.

n_d1=5;    % labour supply
n_d2=11;   % resources placed into the inheritance asset (deliberately ~= n_a, see above)
n_d=[n_d1,n_d2];
n_a=101;   % the inheritance asset (a2). There is no a1: the toolkit does not support one here.
n_z=5;

d1_grid=linspace(0,1,n_d1)';        % labour supply, in [0,1]
d2_grid=linspace(0,5,n_d2)';        % same span as a_grid, so the decision can reach the whole asset grid
d_grid=[d1_grid; d2_grid];

a_grid=5*linspace(0,1,n_a)'.^3;     % the inheritance asset, in [0,5]

% setup z
[z_grid,pi_z]=discretizeAR1_FarmerToda(0,0.9,0.03,n_z);
z_grid=exp(z_grid);

%% Parameters
Params.beta=0.95; % discount factor
DiscountFactorParamNames={'beta'};

% Preferences
Params.sigma=2;    % CES utility param for consumption
Params.eta=1.5;    % curvature of leisure
Params.varphi=0.8; % relative weight of leisure in utility

% Prices
Params.w=1;
Params.r=0.05;

% Inheritance asset
% a2prime = d2*(zprime/z)^inheritrisk
% With inheritrisk=1 and z spanning roughly exp(+-0.2), the ratio zprime/z reaches about
% [0.67,1.5], so d2 near the top of d2_grid lands above a_grid(end) and d2=0 lands at
% a_grid(1). Both ends are reached, which is the point (see note 2 above).
Params.inheritrisk=1;

%% aprimeFn
% Inputs are (d2,z,zprime) followed by parameters. a2 is NOT an input: that is what makes this
% an inheritance asset rather than an experience asset.
aprimeFn=@(d2,z,zprime,inheritrisk) d2*(zprime/z)^inheritrisk;

%% vfoptions/simoptions baselines
% inheritanceasset=1 plus the aprimeFn is the whole interface on the value fn side. The
% distribution side additionally needs a_grid, d_grid and z_grid handed to it on simoptions.
vfoptionsbaseline=struct();
vfoptionsbaseline.inheritanceasset=1;
vfoptionsbaseline.aprimeFn=aprimeFn;

simoptionsbaseline=struct();
simoptionsbaseline.inheritanceasset=1;
simoptionsbaseline.aprimeFn=aprimeFn;
simoptionsbaseline.a_grid=a_grid;
simoptionsbaseline.d_grid=d_grid;
simoptionsbaseline.z_grid=z_grid;
