function F=ReturnFn_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension)
% RiskyAsset, no d1, z, e, no semiz.

F=-Inf;

if agej<Jr
    c=w*kappa_j*z*e + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end

end
