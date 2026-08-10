function F=ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma)

F=-Inf;

c=(1+r)*a+w*z*e-aprime;

if c>0
    F=(c^(1-sigma)-1)/(1-sigma);
end


end
