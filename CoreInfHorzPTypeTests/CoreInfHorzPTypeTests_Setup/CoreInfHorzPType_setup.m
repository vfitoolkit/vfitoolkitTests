% PType setup for CoreInfHorzPTypeTests, on top of CoreInfHorz_setup (the d, a, z grids and the parameters of CoreInfHorzTests).
%
% The models of this bank: two types that differ by a per-type wage w_pt (given as a vector with one entry per type, which
% PType reads as a type-dependent parameter), with ptweights 0.6/0.4. Each part sets its own return function and
% FnsToEvaluate. The stationary distributions are solved to a tight tolerance (simoptionsPT.tolerance), so that the
% identities that rely on the distribution being stationary hold to rounding error rather than to the default tolerance.

PTypeDistParamNames={'ptypemass'};
Params.ptypemass=[0.6;0.4];
Params.w_pt=[1;1.4];
simoptionsPT=struct();
simoptionsPT.tolerance=10^(-12);
