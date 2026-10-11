function F=ReturnFn_ExpAssetz_nod1_z_noe_with2A2(d2, a1prime, a1, a2_1, a2_2, z, r, w, kappa_j, sigma, agej, Jr, pension, pensionrate)
% with2A2 version of ReturnFn_ExpAssetz_nod1_z_noe: the experience asset splits into a2_1 (human capital, same law
% as the 1A tier) and a2_2 (cumulated lifetime earnings), which pays a pension at pensionrate.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2_1*z-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
