function F=ReturnFn_ExpAssetsemiz_d1_noz_e_with2A2(d1,d2,d3,a1prime,a1,a2_1,a2_2,semiz,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost,pensionrate)
% with2A2: TWO experience assets, both driven by semiz. a2_1 carries the law of motion of
% the one-experience-asset tier, so with a2_2 inert this reduces to that tier exactly.
% a2_2 raises the retirement pension (pensionrate), which is what makes V genuinely depend
% on it rather than it being a payoff-irrelevant passenger.


F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2_1*semiz*e+uempbenefit*(1-semiz)-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d3;
end


end
