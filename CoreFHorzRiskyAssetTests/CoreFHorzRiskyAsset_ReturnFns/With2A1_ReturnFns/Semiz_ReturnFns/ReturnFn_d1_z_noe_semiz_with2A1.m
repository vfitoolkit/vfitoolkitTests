function F=ReturnFn_d1_z_noe_semiz_with2A1(h,savings,dsemiz,a1prime,a1_2prime,a1,a1_2,a2,semiz,z,r,w,kappa_j,sigma,eta,varphi,r_a1,r_a1_2,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset with TWO standard endogenous assets + a2 (risky asset):
%  a1_1 = liquid safe asset, return r_a1      (this is the one that is divide-conquered)
%  a1_2 = binary capped asset, return r_a1_2>r_a1, holdings in {0,1} (this is the one that is folded)
% Variant: d1_z_noe_semiz. a1_1 earns r_a1 and a1_2 earns r_a1_2; a2 enters directly (risky return realised via aprimeFn).

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*z*semiz + uempbenefit*(1-semiz) + (1+r_a1)*a1 + (1+r_a1_2)*a1_2 + a2 - a1prime - a1_2prime - savings;
else
    c=pension + (1+r_a1)*a1 + (1+r_a1_2)*a1_2 + a2 - a1prime - a1_2prime - savings;
end

if c>0 && h<1
    F=(c^(1-sigma)-1)/(1-sigma) + varphi*((1-h)^(1-eta)-1)/(1-eta) - searcheffortcost*dsemiz;
end

end
