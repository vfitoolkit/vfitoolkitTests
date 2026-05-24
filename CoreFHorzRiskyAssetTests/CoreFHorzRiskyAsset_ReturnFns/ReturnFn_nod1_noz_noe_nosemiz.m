function F=ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension)
% RiskyAsset, no d1, no z, no e, no semiz.
% Under refine_d=[0,1,1], only d3 (savings) is in ReturnFn.
% Note: asset returns are realised via aprimeFn, so 'a' enters the budget directly (no (1+r)*a here).

F=-Inf;

if agej<Jr
    c=w*kappa_j + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end

end
