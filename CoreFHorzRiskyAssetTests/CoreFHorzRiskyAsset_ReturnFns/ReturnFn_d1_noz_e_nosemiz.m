function F=ReturnFn_d1_noz_e_nosemiz(h,savings,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)
% RiskyAsset, with d1 (h=labour), no z, e, no semiz.

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*e + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    F=(c^(1-sigma)-1)/(1-sigma) + varphi*((1-h)^(1-eta)-1)/(1-eta);
end

end
