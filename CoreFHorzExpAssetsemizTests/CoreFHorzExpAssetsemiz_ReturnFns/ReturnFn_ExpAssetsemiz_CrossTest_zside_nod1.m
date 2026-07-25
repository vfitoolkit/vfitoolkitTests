function F=ReturnFn_ExpAssetsemiz_CrossTest_zside_nod1(d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit)
% Comparison side for CrossTest1: plain experienceassetz, with z playing the
% role of semiz. Same payoff as ReturnFn_ExpAssetsemiz_nod1_noz_noe when d3=0.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2*z+uempbenefit*(1-z)-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
