function F=EZRiskyReturnFn_negativeUtils_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,ezsigma,agej,Jr,pension)
% RiskyAsset, no d1, no z, no e, no semiz. UTILITY-UNITS Epstein-Zin, NEGATIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=0): CRRA without the -1 shift, so strictly negative
% for ezsigma>1. Doubles as the vNM reference utility fn for the gamma=1/phi collapse test
% (at ezsigma=ezgamma).

F=-Inf;

if agej<Jr
    c=w*kappa_j + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=(c^(1-ezsigma))/(1-ezsigma);
end

end
