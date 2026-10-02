function F=ReturnFn_ExpAssetze_d1_z_e_semiz_noa1_with2A2(d1,d2,d3,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost,pensionrate)
% with2A2: derived from the one-experience-asset version above by taking a2 -> a2_1 and adding
% the second experience asset a2_2, which raises the retirement pension via pensionrate.
% Without that term V would not depend on a2_2 at all and every check in this tier would be
% vacuous. pensionrate is appended at the END of the parameter list.

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*d2*a2_1*semiz*z*e+uempbenefit*(1-semiz);
else
    c=pension+pensionrate*a2_2;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d3;
end


end
