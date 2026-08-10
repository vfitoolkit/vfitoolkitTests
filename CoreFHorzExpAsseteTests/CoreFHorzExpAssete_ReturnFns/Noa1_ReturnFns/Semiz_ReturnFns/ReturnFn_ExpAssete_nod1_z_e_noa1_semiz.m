function F=ReturnFn_ExpAssete_nod1_z_e_noa1_semiz(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d2*a*semiz*z*e+uempbenefit*(1-semiz);
else
    c=pension;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
