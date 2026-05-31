function F=ReturnFn_nod1_noz_e_nosemiz(d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2*e-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end