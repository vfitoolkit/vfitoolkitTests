function F=ReturnFn_d_z_noe_nosemiz_2Aignored(d,aprime,a2prime,a,a2,z,r,w,sigma,eta,varphi)
% Same as ReturnFn_d_z_noe_nosemiz, but with a second endogenous state (a2) that is
% completely ignored: neither a2 nor a2prime appear anywhere in F.
% So the model is economically identical to the one-endogenous-state model.

F=-Inf;

c=(1+r)*a+w*d*z-aprime;

if c>0 && d<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
end


end
