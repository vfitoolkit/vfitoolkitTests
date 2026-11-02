function F=EZReturnFn_negativeUtils_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension)
% Utility-units, NEGATIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=0).
% No -1 term: F=(x^(1-ezsigma))/(1-ezsigma) is strictly negative for ezsigma>1
% [the baseline (x^(1-ezsigma)-1)/(1-ezsigma) is mixed-sign].
% Also doubles as the vNM reference return fn for the gamma=1/phi collapse test
% (run with ezsigma=ezgamma).

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*e-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    x=c;
    F=(x^(1-ezsigma))/(1-ezsigma);
end


end
