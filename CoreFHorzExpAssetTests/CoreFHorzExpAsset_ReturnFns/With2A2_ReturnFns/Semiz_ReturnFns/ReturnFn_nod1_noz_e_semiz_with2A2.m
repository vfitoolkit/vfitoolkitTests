function F=ReturnFn_nod1_noz_e_semiz_with2A2(d2,d3,a1prime,a1,a2_1,a2_2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2_1*semiz*e+uempbenefit*(1-semiz)-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end