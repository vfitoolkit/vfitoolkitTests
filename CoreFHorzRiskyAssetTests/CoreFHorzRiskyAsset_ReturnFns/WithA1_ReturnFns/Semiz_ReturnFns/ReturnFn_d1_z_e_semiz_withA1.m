function F=ReturnFn_d1_z_e_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset with a1 (standard safe asset) + a2 (risky asset).
% Variant: d1_z_e_semiz. a1 earns safe return r_a1; a2 enters directly (risky return realised via aprimeFn).

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*z*e*semiz + uempbenefit*(1-semiz) + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0 && h<1
    F=(c^(1-sigma)-1)/(1-sigma) + varphi*((1-h)^(1-eta)-1)/(1-eta) - searcheffortcost*dsemiz;
end

end
