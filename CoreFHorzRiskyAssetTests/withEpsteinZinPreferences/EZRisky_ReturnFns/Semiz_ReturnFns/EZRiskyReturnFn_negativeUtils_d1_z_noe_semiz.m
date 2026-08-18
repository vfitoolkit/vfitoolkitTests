function F=EZRiskyReturnFn_negativeUtils_d1_z_noe_semiz(h,savings,dsemiz,a,semiz,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, with d1 (h=labour), z, no e. UTILITY-UNITS Epstein-Zin, NEGATIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=0): CRRA without the -1 shift, so strictly negative
% for ezsigma>1. Doubles as the vNM reference utility fn for the gamma=1/phi collapse test
% (at ezsigma=ezgamma).
% Search effort enters multiplicatively [the (1-searcheffortcost*dsemiz) factor] rather than
% additively as in the vNM bank, so the composite keeps a single sign (as EZ requires);
% requires searcheffortcost<1 so the factor stays positive.

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*z*semiz + uempbenefit*(1-semiz) + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    x=((c^varphi)*((1-h)^(1-varphi)))*(1-searcheffortcost*dsemiz);
    F=(x^(1-ezsigma))/(1-ezsigma);
end

end
