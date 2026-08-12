function F=EZReturnFn_negativeUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)
% Utility-units, NEGATIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=0).
% No -1 term: F=(x^(1-ezsigma))/(1-ezsigma) is strictly negative for ezsigma>1
% [the baseline (x^(1-ezsigma)-1)/(1-ezsigma) is mixed-sign].
% Also doubles as the vNM reference return fn for the gamma=1/phi collapse test
% (run with ezsigma=ezgamma).

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z*e*d-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0 && d<1
    x=(c^varphi)*((1-d)^(1-varphi));
    F=(x^(1-ezsigma))/(1-ezsigma);
end


end
