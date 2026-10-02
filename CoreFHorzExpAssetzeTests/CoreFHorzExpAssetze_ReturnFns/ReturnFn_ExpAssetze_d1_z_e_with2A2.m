function F=ReturnFn_ExpAssetze_d1_z_e_with2A2(d1,d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,pensionrate)
% with2A2: derived from the one-experience-asset version above by taking a2 -> a2_1 and adding
% the second experience asset a2_2, which raises the retirement pension via pensionrate.
% Without that term V would not depend on a2_2 at all and every check in this tier would be
% vacuous. pensionrate is appended at the END of the parameter list.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2_1*z*e-a1prime;
else
    c=(1+r)*a1+pension+pensionrate*a2_2-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end