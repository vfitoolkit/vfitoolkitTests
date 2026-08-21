function F=ReturnFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,r,r2,w,sigma,eta,varphi)
% Two endogenous states: two assets, paying r and r2, plus an endogenous labour supply d.
% Deliberately plain, the point of this model is to exercise the two-endogenous-state code paths.

F=-Inf;

c=(1+r)*a1+(1+r2)*a2+w*d*z-a1prime-a2prime;

if c>0 && d<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
end


end
