function F=EZRiskyReturnFn_positiveUtils_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,ezsigma,agej,Jr,pension)
% RiskyAsset, no d1, no z, e, no semiz. UTILITY-UNITS Epstein-Zin, POSITIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=1): the (1+c) shift keeps utility strictly positive
% for all c>0 despite ezsigma>1.

F=-Inf;

if agej<Jr
    c=w*kappa_j*e + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=((1+c)^(1-ezsigma)-1)/(1-ezsigma);
end

end
