function F=EZRiskyReturnFn_positiveUtils_nod1_noz_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) %#ok<INUSD>
% RiskyAsset with a1 (standard safe asset, earns r_a1) + a2 (risky asset; risky return
% realised via aprimeFn), and semi-exogenous shock semiz. Variant: nod1_noz_noe_semiz.
% Utility-units, POSITIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=1). The (1+x) shift keeps F strictly
% positive for all x>0 despite ezsigma>1.
% Note: requires searcheffortcost<1 so the (1-searcheffortcost*dsemiz) factor stays positive.
% (r is unused: a1 earns the safe return r_a1, a2's risky return is realised via aprimeFn)

F=-Inf;

if agej<Jr
    c=w*kappa_j*semiz + uempbenefit*(1-semiz) + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0
    x=c*(1-searcheffortcost*dsemiz);
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end

end
