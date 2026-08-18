function F=ReturnFn_ExpAssetsemiz_nod1_z_e_noa1(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% noa1: the experienceassetsemiz a2 is the only endogenous state, so there is no a1/a1prime.
% r is kept in the signature (unused) so the noa1 and withA1 ReturnFns line up.

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
