function F=ReturnFn_OneEndo_nod1_noz_noe_noa1_semiz(d3,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% Standard 1-endo Case1 with semiz (no d2 since aprime is the choice). Used as side A of
% CrossTests3 noa1 semiz: aprime plays the role of d2 in the ExpAsset side (where aprimeFn=@(d2,a2) d2).

F=-Inf;

if agej<Jr
    c=w*kappa_j*aprime*a*semiz+uempbenefit*(1-semiz);
else
    c=pension;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
