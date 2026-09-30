function F=ReturnFn_d1_z_noe_nosemiz_with2A2(d1,d2,a1prime,a1,a2_1,a2_2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,pensionrate)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2_1*z-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end