function F=ReturnFn_ExpAssetze_d1_z_e_2A2sym(d1,d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)
% Cross-test ReturnFn: SYMMETRIC in the two experience assets (they enter only through their
% average). Used by CrossTests9: with symmetric preferences, identical grids and identical laws
% of motion, V and Policy must be invariant to swapping the two a2 dimensions.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*(a2_1+a2_2)/2*z*e-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end
