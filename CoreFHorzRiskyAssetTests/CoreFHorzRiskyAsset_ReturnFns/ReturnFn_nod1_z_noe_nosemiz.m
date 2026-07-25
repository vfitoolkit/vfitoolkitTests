function F=ReturnFn_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension)
% RiskyAsset, no d1, z, no e, no semiz.

F=-Inf;

if agej<Jr
    c=w*kappa_j*z + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end

end
