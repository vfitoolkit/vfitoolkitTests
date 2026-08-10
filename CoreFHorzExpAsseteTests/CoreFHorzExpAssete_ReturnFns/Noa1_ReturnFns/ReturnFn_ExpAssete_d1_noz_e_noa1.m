function F=ReturnFn_ExpAssete_d1_noz_e_noa1(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*d2*a*e;
else
    c=pension;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end
