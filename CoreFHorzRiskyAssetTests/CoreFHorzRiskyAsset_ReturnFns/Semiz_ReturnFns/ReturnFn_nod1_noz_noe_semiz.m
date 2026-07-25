function F=ReturnFn_nod1_noz_noe_semiz(savings,dsemiz,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, no d1, no z, no e.
% Under refine_d=[0,1,1,1]: ReturnFn takes (d3, d4, a, semiz, ...).
%   d3 = savings, d4 = dsemiz.

F=-Inf;

if agej<Jr
    c=w*kappa_j*semiz + uempbenefit*(1-semiz) + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma) - searcheffortcost*dsemiz;
end

end
