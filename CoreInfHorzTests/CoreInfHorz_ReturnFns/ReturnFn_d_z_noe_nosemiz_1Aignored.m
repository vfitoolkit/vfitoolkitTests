function F=ReturnFn_d_z_noe_nosemiz_1Aignored(d,a1prime,a2prime,a1,a2,z,r,w,sigma,eta,varphi)
% Same as ReturnFn_d_z_noe_nosemiz, but with two endogenous states where the FIRST one (a1)
% is completely ignored: neither a1 nor a1prime appear anywhere in F. The second endogenous
% state (a2) is the actual asset.
% So the model is economically identical to the one-endogenous-state model.

F=-Inf;

c=(1+r)*a2+w*d*z-a2prime;

if c>0 && d<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
end


end
