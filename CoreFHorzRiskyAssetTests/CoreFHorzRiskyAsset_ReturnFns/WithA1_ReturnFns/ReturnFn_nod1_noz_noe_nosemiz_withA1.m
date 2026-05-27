function F=ReturnFn_nod1_noz_noe_nosemiz_withA1(savings,a1prime,a1,a2,r,w,kappa_j,sigma,r_a1,agej,Jr,pension)
% RiskyAsset with a1 (standard safe asset) + a2 (risky asset).
% Variant: nod1_noz_noe_nosemiz. a1 earns safe return r_a1; a2 enters directly (risky return realised via aprimeFn).

F=-Inf;

if agej<Jr
    c=w*kappa_j + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end

end
