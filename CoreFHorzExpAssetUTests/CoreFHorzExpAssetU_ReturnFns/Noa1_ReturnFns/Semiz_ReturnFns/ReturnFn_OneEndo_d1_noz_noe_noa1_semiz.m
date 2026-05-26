function F=ReturnFn_OneEndo_d1_noz_noe_noa1_semiz(d1,d3,aprime,a,semiz,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost)
% Standard 1-endo Case1 with d1+semiz (no d2). Used as side A of CrossTests3 noa1 semiz d1.

F=-Inf;

if agej<Jr
    c=w*kappa_j*d1*aprime*a*semiz+uempbenefit*(1-semiz);
else
    c=pension;
end

if c>0 && d1<1
    F=(c^(1-sigma)-1)/(1-sigma)+varphi*((1-d1)^(1-eta)-1)/(1-eta)-searcheffortcost*d3;
end


end
