function F=ReturnFn_ExpAssetz_nod1_z_noe_semiz_noa1(d2,d3,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% noa1 version of ReturnFn_ExpAssetz_nod1_z_noe_semiz: drop a1prime and a1 (a is the experience asset a2)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d2*a*semiz*z+uempbenefit*(1-semiz);
else
    c=pension;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
