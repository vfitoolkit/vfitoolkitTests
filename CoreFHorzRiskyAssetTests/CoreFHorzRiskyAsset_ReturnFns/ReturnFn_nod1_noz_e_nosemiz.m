function F=ReturnFn_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension)
% RiskyAsset, no d1, no z, e, no semiz.

F=-Inf;

if agej<Jr
    c=w*kappa_j*e + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end

end
