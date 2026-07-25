function F=ReturnFn_ExpAssetsemiz_nod1_noz_e(d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2*semiz*e+uempbenefit*(1-semiz)-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
