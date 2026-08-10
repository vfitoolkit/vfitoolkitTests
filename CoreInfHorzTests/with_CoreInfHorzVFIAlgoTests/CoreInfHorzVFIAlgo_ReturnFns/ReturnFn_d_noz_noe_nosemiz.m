function F=ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,sigma,eta,varphi)

F=-Inf;

c=(1+r)*a+w*d-aprime;

if c>0 && d<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d)^(1-eta)-1)/(1-eta);
end


end
