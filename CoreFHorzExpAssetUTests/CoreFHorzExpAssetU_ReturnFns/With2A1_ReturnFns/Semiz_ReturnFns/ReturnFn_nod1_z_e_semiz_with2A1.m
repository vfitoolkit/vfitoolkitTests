function F=ReturnFn_nod1_z_e_semiz_with2A1(d2,d3,a1prime,a1_2prime,a1,a1_2,a2,semiz,z,e,r,r2,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+(1+r2)*a1_2+w*kappa_j*d2*a2*semiz*z*e+uempbenefit*(1-semiz)-a1prime-a1_2prime;
else
    c=(1+r)*a1+(1+r2)*a1_2+pension-a1prime-a1_2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end