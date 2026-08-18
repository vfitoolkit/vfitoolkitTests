function F=EZRiskyReturnFn_positiveUtils_d1_noz_e_nosemiz(h,savings,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)
% RiskyAsset, with d1 (h=labour), no z, e, no semiz. UTILITY-UNITS Epstein-Zin, POSITIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=1): utility of the composite consumption-leisure
% good x=(c^varphi)*((1-h)^(1-varphi)); the (1+x) shift keeps utility strictly positive
% for all x>0 despite ezsigma>1.

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*e + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    x=(c^varphi)*((1-h)^(1-varphi));
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end

end
