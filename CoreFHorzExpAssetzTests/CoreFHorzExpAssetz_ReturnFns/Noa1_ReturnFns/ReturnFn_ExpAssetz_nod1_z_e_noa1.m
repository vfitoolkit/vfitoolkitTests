function F=ReturnFn_ExpAssetz_nod1_z_e_noa1(d2,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension)
% noa1 version of ReturnFn_ExpAssetz_nod1_z_e: drop a1prime and a1 (a is the experience asset a2)

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
