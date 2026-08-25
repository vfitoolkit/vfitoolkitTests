function WG=EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3)
% De Nardi luxury-good warm-glow of bequests, UTILITY-UNITS NEGATIVE case
% (vfoptions.EZutils=1, EZpositiveutility=0): the exact form of LifeCycleModel12 of
% IntroToLifeCycleModels — strictly negative for all aprime>=0 with wg3>1 (the (1+aprime/wg2)
% shifter means it is nonzero at aprime=0, avoiding the WG==0 mask conventions).
% wg1: strength of the bequest motive (theta); wg2: luxury shifter (kappa); wg3: curvature
% (set equal to ezsigma so the warm-glow 'outputs the same thing as the return fn').

WG=wg1*((1+aprime/wg2)^(1-wg3))/(1-wg3);

end
