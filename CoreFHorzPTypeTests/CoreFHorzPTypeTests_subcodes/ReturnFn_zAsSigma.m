function F=ReturnFn_zAsSigma(aprime,a,z,r,w,kappa_j,agej,Jr,pension)
% z is interpreted as the curvature (sigma). Used by CoreFHorzPType_NivsZidentity.

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=(c^(1-z)-1)/(1-z);
end

end
