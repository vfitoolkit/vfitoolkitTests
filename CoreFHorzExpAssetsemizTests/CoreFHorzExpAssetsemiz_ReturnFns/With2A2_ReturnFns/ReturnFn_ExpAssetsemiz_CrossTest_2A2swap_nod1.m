function F=ReturnFn_ExpAssetsemiz_CrossTest_2A2swap_nod1(d2,d3,a1prime,a1,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% Cross-test ReturnFn: identical to ReturnFn_ExpAssetsemiz_nod1_z_e except that earnings are driven by the
% SECOND experience asset a2_2 (a2_1 does not appear). Used by CrossTests7, where a2_1
% is the inert dimension and a2_2 carries the real law of motion.
% checkcode flags a2_1 as an unused input (INUSD). That is the POINT of this function: a2_1
% must not enter the return, so the inert-first-asset reduction in CrossTests7 is exact. Leave it.
F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2_2*semiz*z*e+uempbenefit*(1-semiz)-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
