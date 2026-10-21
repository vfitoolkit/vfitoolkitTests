function F=ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost)

F=-Inf;

c=(1+r)*a+w*d1*semiz*e+uempbenefit*(1-semiz)-aprime;

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d2;
end


end
