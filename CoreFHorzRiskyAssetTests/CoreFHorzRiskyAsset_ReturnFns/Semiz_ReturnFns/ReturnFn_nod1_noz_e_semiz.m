function F=ReturnFn_nod1_noz_e_semiz(savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, no d1, no z, e.

F=-Inf;

if agej<Jr
    c=w*kappa_j*e*semiz + uempbenefit*(1-semiz) + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma) - searcheffortcost*dsemiz;
end

end
