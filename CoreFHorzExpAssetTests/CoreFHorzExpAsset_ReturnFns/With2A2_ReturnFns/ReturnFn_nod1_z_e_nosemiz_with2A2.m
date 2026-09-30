function F=ReturnFn_nod1_z_e_nosemiz_with2A2(d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2_1*z*e-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end