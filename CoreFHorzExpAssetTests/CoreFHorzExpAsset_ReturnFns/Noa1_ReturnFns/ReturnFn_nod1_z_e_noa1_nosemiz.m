function F=ReturnFn_nod1_z_e_noa1_nosemiz(d2,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d2*a*z*e;
else
    c=pension;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
