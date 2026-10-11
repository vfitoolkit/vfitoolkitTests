function F=ReturnFn_ExpAssetz_nod1_z_e_semiz_noa1_with2A2(d2, d3, a2_1, a2_2, semiz, z, e, r, w, kappa_j, sigma, agej, Jr, pension, uempbenefit, searcheffortcost, pensionrate)
% with2A2 version of ReturnFn_ExpAssetz_nod1_z_e_semiz_noa1: the experience asset splits into a2_1 (human capital, same law
% as the 1A tier) and a2_2 (cumulated lifetime earnings), which pays a pension at pensionrate.

% noa1 version of ReturnFn_ExpAssetz_nod1_z_e_semiz: drop a1prime and a1 (a2_1 and a2_2 are the two experience assets)

F=-Inf;

if agej<Jr
    c=w*kappa_j*d2*a2_1*semiz*z*e+uempbenefit*(1-semiz);
else
    c=pension+pensionrate*a2_2;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
