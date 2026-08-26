function WG=EZRiskyWarmGlowFn_cons(aprime,wg1,wg2)
% De Nardi luxury-good warm-glow of bequests, CONSUMPTION-UNITS Epstein-Zin case
% (vfoptions.EZutils=0): the warm-glow is a consumption-equivalent, strictly positive
% for all aprime>=0 (curvature comes from the EZ preferences, like the cons-units ReturnFn).
% In the riskyasset banks the warm-glow argument is a2prime (the risky-asset holding
% realized via the (d,u) lottery); the toolkit evaluates it on the a2 grid.
% wg1: strength of the bequest motive (De Nardi's theta); wg2: luxury-good shifter (kappa).

WG=wg1*(1+aprime/wg2);

end
