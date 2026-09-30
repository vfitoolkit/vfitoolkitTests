function F=ReturnFn_nod1_noz_noe_nosemiz_2A2sym(d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,agej,Jr,pension)
% Cross-test ReturnFn: SYMMETRIC in the two experience assets (they enter only through their
% average). Used by CrossTests7: with symmetric preferences, identical grids and identical laws
% of motion, V and Policy must be invariant to swapping the two a2 dimensions.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*(a2_1+a2_2)/2-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
