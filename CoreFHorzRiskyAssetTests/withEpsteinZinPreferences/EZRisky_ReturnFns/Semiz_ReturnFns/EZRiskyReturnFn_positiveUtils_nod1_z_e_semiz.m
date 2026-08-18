function F=EZRiskyReturnFn_positiveUtils_nod1_z_e_semiz(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, no d1, z, e. UTILITY-UNITS Epstein-Zin, POSITIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=1): the (1+x) shift keeps utility strictly positive
% for all x>0 despite ezsigma>1.
% Search effort enters multiplicatively [the (1-searcheffortcost*dsemiz) factor] rather than
% additively as in the vNM bank, so the composite keeps a single sign (as EZ requires);
% requires searcheffortcost<1 so the factor stays positive.

F=-Inf;

if agej<Jr
    c=w*kappa_j*z*e*semiz + uempbenefit*(1-semiz) + a - savings;
else
    c=pension + a - savings;
end

if c>0
    x=c*(1-searcheffortcost*dsemiz);
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end

end
