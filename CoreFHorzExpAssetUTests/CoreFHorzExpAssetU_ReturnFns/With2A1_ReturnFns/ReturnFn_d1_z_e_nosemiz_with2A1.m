function F=ReturnFn_d1_z_e_nosemiz_with2A1(d1,d2,a1prime,a1_2prime,a1,a1_2,a2,z,e,r,r2,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+(1+r2)*a1_2+w*kappa_j*d1*d2*a2*z*e-a1prime-a1_2prime;
else
    c=(1+r)*a1+(1+r2)*a1_2+pension-a1prime-a1_2prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end