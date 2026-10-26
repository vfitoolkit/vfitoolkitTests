function F=EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma)
% Utility-units, NEGATIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=0).
% No -1 term: F=(x^(1-ezsigma))/(1-ezsigma) is strictly negative for ezsigma>1
% [the baseline (x^(1-ezsigma)-1)/(1-ezsigma) is mixed-sign, and the Epstein-Zin raws
% treat a return of exactly zero as infeasible].
% Also doubles as the vNM reference return fn for the gamma=1/phi collapse test
% (run with ezsigma=ezgamma).
% InfHorz form of EZReturnFn_negativeUtils_nod_z_noe_nosemiz (no ages, no retirement).

F=-Inf;

c=(1+r)*a+w*z-aprime;

if c>0
    x=c;
    F=(x^(1-ezsigma))/(1-ezsigma);
end


end
