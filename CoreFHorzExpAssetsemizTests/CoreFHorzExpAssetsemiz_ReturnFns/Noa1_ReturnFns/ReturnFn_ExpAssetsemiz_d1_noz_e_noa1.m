function F=ReturnFn_ExpAssetsemiz_d1_noz_e_noa1(d1,d2,d3,a,semiz,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)
% noa1: the experienceassetsemiz a2 is the only endogenous state, so there is no a1/a1prime.
% r is kept in the signature (unused) so the noa1 and withA1 ReturnFns line up.

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*d2*a*semiz*e+uempbenefit*(1-semiz);
else
    c=pension;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d3;
end


end
