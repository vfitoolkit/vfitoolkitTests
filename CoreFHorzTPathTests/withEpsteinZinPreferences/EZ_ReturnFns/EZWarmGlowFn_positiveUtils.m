function WG=EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3)
% De Nardi luxury-good warm-glow of bequests, UTILITY-UNITS POSITIVE case
% (vfoptions.EZutils=1, EZpositiveutility=1): the LifeCycleModel12 form adapted with the
% same shift trick as the positiveUtils ReturnFns — ((2+aprime/wg2)^(1-wg3)-1)/(1-wg3) is
% strictly positive and increasing in aprime for wg3>1, and nonzero at aprime=0
% (=(2^(1-wg3)-1)/(1-wg3)>0), avoiding the WG==0 mask conventions.
% wg1: strength of the bequest motive (theta); wg2: luxury shifter (kappa); wg3: curvature.

WG=wg1*(((2+aprime/wg2)^(1-wg3))-1)/(1-wg3);

end
