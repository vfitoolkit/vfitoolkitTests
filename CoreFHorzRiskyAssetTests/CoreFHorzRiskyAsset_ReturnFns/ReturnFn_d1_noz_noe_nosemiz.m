function F=ReturnFn_d1_noz_noe_nosemiz(h,savings,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)
% RiskyAsset, with d1 (h=labour), no z, no e, no semiz.
% Under refine_d=[1,1,1]: ReturnFn takes (d1, d3, ...).

F=-Inf;

if agej<Jr
    c=w*kappa_j*h + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    F=(c^(1-sigma)-1)/(1-sigma) + varphi*((1-h)^(1-eta)-1)/(1-eta);
end

end
