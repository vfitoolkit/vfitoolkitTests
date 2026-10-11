function F=ReturnFn_ExpAssetz_d1_z_noe_semiz_with2A2(d1, d2, d3, a1prime, a1, a2_1, a2_2, semiz, z, r, w, kappa_j, sigma, varphi, eta, agej, Jr, pension, uempbenefit, searcheffortcost, pensionrate)
% with2A2 version of ReturnFn_ExpAssetz_d1_z_noe_semiz: the experience asset splits into a2_1 (human capital, same law
% as the 1A tier) and a2_2 (cumulated lifetime earnings), which pays a pension at pensionrate.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2_1*semiz*z+uempbenefit*(1-semiz)-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d3;
end


end
