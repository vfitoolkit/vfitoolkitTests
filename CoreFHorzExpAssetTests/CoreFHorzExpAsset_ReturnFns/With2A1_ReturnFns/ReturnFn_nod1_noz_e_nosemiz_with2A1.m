function F=ReturnFn_nod1_noz_e_nosemiz_with2A1(d2,a1prime,a1_2prime,a1,a1_2,a2,e,r,r2,w,kappa_j,sigma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+(1+r2)*a1_2+w*kappa_j*d2*a2*e-a1prime-a1_2prime;
else
    c=(1+r)*a1+(1+r2)*a1_2+pension-a1prime-a1_2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end