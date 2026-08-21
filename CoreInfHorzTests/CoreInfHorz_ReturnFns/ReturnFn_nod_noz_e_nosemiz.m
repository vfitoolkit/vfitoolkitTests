function F=ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma)

F=-Inf;

c=(1+r)*a+w*e-aprime;

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
