function F=EZRiskyReturnFn_negativeUtils_nod1_z_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) %#ok<INUSD>
% RiskyAsset with a1 (standard safe asset, earns r_a1) + a2 (risky asset; risky return
% realised via aprimeFn), and semi-exogenous shock semiz. Variant: nod1_z_noe_semiz.
% Utility-units, NEGATIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=0): CRRA without the -1 shift, so
% strictly negative for ezsigma>1. Doubles as the vNM reference return fn for the
% gamma=1/phi collapse test (run with ezsigma=ezgamma).
% Note: requires searcheffortcost<1 so the (1-searcheffortcost*dsemiz) factor stays positive.
% (r is unused: a1 earns the safe return r_a1, a2's risky return is realised via aprimeFn)

F=-Inf;

if agej<Jr
    c=w*kappa_j*z*semiz + uempbenefit*(1-semiz) + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0
    x=c*(1-searcheffortcost*dsemiz);
    F=(x^(1-ezsigma))/(1-ezsigma);
end

end
