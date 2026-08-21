function F=ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,r2,w,sigma)
% Two endogenous states: two assets, paying r and r2. Deliberately plain, the point of this model is
% to exercise the two-endogenous-state code paths, not to be interesting economics.

F=-Inf;

c=(1+r)*a1+(1+r2)*a2+w*z-a1prime-a2prime;

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
