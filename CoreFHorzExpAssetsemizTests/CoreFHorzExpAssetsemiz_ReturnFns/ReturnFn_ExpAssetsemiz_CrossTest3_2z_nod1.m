function F=ReturnFn_ExpAssetsemiz_CrossTest3_2z_nod1(d2,a1prime,a1,a2,z1,z2,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit)
% Comparison side for CrossTest3: experienceassetz with a combined z=[semiz,z].
% z1 plays semiz (fast index of bothz), z2 the ordinary z.
% Same payoff as ReturnFn_ExpAssetsemiz_nod1_z_noe when d3=0.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2*z1*z2+uempbenefit*(1-z1)-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
