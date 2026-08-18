function F=ReturnFn_ExpAssetsemiz_nod1_noz_noe_with2A1(d2,d3,a1prime,a1_2prime,a1,a1_2,a2,semiz,r,r2,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% Two standard endogenous assets (triggers DC2A/GI2A):
%  a1   = liquid asset, return r        (this is the one that is divide-conquered)
%  a1_2 = binary capped asset, return r2>r, holdings in {0,1} (this is the one that is folded)
%  a2   = experienceassetsemiz (unchanged)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+(1+r2)*a1_2+w*kappa_j*d2*a2*semiz+uempbenefit*(1-semiz)-a1prime-a1_2prime;
else
    c=(1+r)*a1+(1+r2)*a1_2+pension-a1prime-a1_2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
