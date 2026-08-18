function F=EZRiskyReturnFn_cons_nod1_noz_noe_nosemiz_withA1(savings,a1prime,a1,a2,r,w,kappa_j,r_a1,agej,Jr,pension) %#ok<INUSD>
% RiskyAsset with a1 (safe, earns r_a1) + a2 (risky, return realised via aprimeFn); no d1, no z, no e, no semiz.
% CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0): the return fn is the consumption
% good itself (curvature comes from the EZ preferences).
% (r is unused here but kept in the signature for consistency with the vNM withA1 ReturnFns)

F=-Inf;

if agej<Jr
    c=w*kappa_j + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0
    F=c;
end

end
