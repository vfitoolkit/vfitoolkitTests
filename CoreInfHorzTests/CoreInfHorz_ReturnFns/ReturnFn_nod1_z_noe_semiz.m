function F=ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost)

F=-Inf;

c=(1+r)*a+w*z*semiz+uempbenefit*(1-semiz)-aprime;

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d2;
end


end
