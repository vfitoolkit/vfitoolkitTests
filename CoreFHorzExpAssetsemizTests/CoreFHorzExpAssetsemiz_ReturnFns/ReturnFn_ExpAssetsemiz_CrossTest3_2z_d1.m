function F=ReturnFn_ExpAssetsemiz_CrossTest3_2z_d1(d1,d2,a1prime,a1,a2,z1,z2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit)
% Comparison side for CrossTest3 (with d1): experienceassetz with combined z=[semiz,z].
% z1 plays semiz (fast index of bothz), z2 the ordinary z.
% Same payoff as ReturnFn_ExpAssetsemiz_d1_z_noe when d3=0.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2*z1*z2+uempbenefit*(1-z1)-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end
