function F=ReturnFn_ExpAssetsemiz_CrossTest_2A2sym_d1(d1,d2,d3,a1prime,a1,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)
% Cross-test ReturnFn: SYMMETRIC in the two experience assets (they enter only through their
% average). With identical grids and identical laws of motion, V, Policy and the agent
% distribution must then be invariant to swapping the two a2 dimensions.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*(a2_1+a2_2)/2*semiz*z*e+uempbenefit*(1-semiz)-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end

end
