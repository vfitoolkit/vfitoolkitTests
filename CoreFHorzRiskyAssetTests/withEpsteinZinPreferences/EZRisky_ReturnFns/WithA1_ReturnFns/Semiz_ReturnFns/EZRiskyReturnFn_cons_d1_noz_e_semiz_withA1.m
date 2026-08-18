function F=EZRiskyReturnFn_cons_d1_noz_e_semiz_withA1(h,savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,varphi,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) %#ok<INUSD>
% RiskyAsset with a1 (standard safe asset, earns r_a1) + a2 (risky asset; risky return
% realised via aprimeFn), and semi-exogenous shock semiz. Variant: d1_noz_e_semiz.
% CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0): F is the composite consumption
% -leisure good (strictly positive whenever feasible); curvature comes from the EZ preferences.
% Note: requires searcheffortcost<1 so the (1-searcheffortcost*dsemiz) factor stays positive.
% (r is unused: a1 earns the safe return r_a1, a2's risky return is realised via aprimeFn)

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*e*semiz + uempbenefit*(1-semiz) + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0 && h<1
    F=((c^varphi)*((1-h)^(1-varphi)))*(1-searcheffortcost*dsemiz);
end

end
