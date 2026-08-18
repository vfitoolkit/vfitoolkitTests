function F=EZRiskyReturnFn_positiveUtils_nod1_z_e_nosemiz_withA1(savings,a1prime,a1,a2,z,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension) %#ok<INUSD>
% RiskyAsset with a1 (safe, earns r_a1) + a2 (risky, return realised via aprimeFn); no d1, z, e, no semiz.
% UTILITY-UNITS Epstein-Zin, POSITIVE-valued utility fn (vfoptions.EZutils=1,
% EZpositiveutility=1): the (1+x) shift keeps utility strictly positive for all x>0
% despite ezsigma>1.
% (r is unused here but kept in the signature for consistency with the vNM withA1 ReturnFns)

F=-Inf;

if agej<Jr
    c=w*kappa_j*z*e + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0
    F=((1+c)^(1-ezsigma)-1)/(1-ezsigma);
end

end
