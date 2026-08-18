function F=EZRiskyReturnFn_cons_d1_noz_e_nosemiz_withA1(h,savings,a1prime,a1,a2,e,r,w,kappa_j,varphi,r_a1,agej,Jr,pension) %#ok<INUSD>
% RiskyAsset with a1 (safe, earns r_a1) + a2 (risky, return realised via aprimeFn); d1, no z, e, no semiz.
% CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0): the return fn is the composite
% consumption good x=(c^varphi)*((1-h)^(1-varphi)) itself (curvature comes from the EZ
% preferences).
% (r is unused here but kept in the signature for consistency with the vNM withA1 ReturnFns)

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*e + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0 && h<1
    x=(c^varphi)*((1-h)^(1-varphi));
    F=x;
end

end
