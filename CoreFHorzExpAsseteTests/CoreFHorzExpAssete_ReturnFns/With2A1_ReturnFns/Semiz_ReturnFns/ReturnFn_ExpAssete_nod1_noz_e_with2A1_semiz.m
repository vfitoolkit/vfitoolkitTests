function F=ReturnFn_ExpAssete_nod1_noz_e_with2A1_semiz(d2,d3,a1prime,a1_2prime,a1,a1_2,a2,semiz,e,r,r2,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% Two standard endogenous assets (triggers DC2A/GI2A), plus semiz:
%  a1   = liquid asset, return r        (this is the one that is divide-conquered)
%  a1_2 = binary capped asset, return r2>r, holdings in {0,1} (this is the one that is folded)
%  a2   = experienceassete (unchanged)
%  semiz = employment state (semi-exogenous, driven by d3 search effort)

F=-Inf;

if agej<Jr
    c=(1+r)*a1+(1+r2)*a1_2+w*kappa_j*d2*a2*semiz*e+uempbenefit*(1-semiz)-a1prime-a1_2prime;
else
    c=(1+r)*a1+(1+r2)*a1_2+pension-a1prime-a1_2prime;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma)-searcheffortcost*d3;
end


end
