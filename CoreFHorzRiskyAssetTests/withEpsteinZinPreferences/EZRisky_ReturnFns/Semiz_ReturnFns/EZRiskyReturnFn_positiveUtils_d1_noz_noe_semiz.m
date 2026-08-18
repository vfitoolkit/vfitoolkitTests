function F=EZRiskyReturnFn_positiveUtils_d1_noz_noe_semiz(h,savings,dsemiz,a,semiz,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, with d1 (h=labour), no z, no e. UTILITY-UNITS Epstein-Zin, POSITIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=1): the (1+x) shift keeps utility strictly positive
% for all x>0 despite ezsigma>1.
% Search effort enters multiplicatively [the (1-searcheffortcost*dsemiz) factor] rather than
% additively as in the vNM bank, so the composite keeps a single sign (as EZ requires);
% requires searcheffortcost<1 so the factor stays positive.
% Under refine_d=[1,1,1,1]: ReturnFn takes (d1, d3, d4, a, semiz, ...).

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*semiz + uempbenefit*(1-semiz) + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    x=((c^varphi)*((1-h)^(1-varphi)))*(1-searcheffortcost*dsemiz);
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end

end
