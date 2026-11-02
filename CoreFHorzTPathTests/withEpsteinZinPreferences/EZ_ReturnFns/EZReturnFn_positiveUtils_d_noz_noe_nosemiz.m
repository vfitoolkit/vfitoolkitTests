function F=EZReturnFn_positiveUtils_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)
% Utility-units, POSITIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=1).
% The (1+x) shift keeps F strictly positive for all x>0 [a plain
% (x^(1-ezsigma)-1)/(1-ezsigma) would change sign at x=1].

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*d-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0 && d<1
    x=(c^varphi)*((1-d)^(1-varphi));
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end


end
