function F=ReturnFn_d1_noz_noe_nosemiz_2A2swap(d1,d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension)
% Cross-test ReturnFn: identical to ReturnFn_d1_noz_noe_nosemiz except that earnings are
% driven by the SECOND experience asset a2_2 (a2_1 does not appear). Used by CrossTests6, where
% a2_1 is the inert dimension and a2_2 carries the real law of motion.

F=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2_2-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta);
end


end
