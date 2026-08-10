function F=ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma)

F=-Inf;

c=(1+r)*a+w-aprime;

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
