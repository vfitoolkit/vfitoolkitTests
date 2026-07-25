function F=ReturnFn_nod1_z_noe_semiz(savings,dsemiz,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, no d1, z, no e.

F=-Inf;

if agej<Jr
    c=w*kappa_j*z*semiz + uempbenefit*(1-semiz) + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-sigma)-1)/(1-sigma) - searcheffortcost*dsemiz;
end

end
