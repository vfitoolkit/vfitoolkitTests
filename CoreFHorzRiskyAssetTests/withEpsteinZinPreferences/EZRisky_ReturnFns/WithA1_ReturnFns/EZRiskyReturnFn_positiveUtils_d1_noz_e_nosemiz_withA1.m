function F=EZRiskyReturnFn_positiveUtils_d1_noz_e_nosemiz_withA1(h,savings,a1prime,a1,a2,e,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension) %#ok<INUSD>
% RiskyAsset with a1 (safe, earns r_a1) + a2 (risky, return realised via aprimeFn); d1, no z, e, no semiz.
% UTILITY-UNITS Epstein-Zin, POSITIVE-valued utility fn (vfoptions.EZutils=1,
% EZpositiveutility=1): the (1+x) shift keeps utility strictly positive for all x>0
% despite ezsigma>1.
% (r is unused here but kept in the signature for consistency with the vNM withA1 ReturnFns)

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*e + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0 && h<1
    x=(c^varphi)*((1-h)^(1-varphi));
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end

end
