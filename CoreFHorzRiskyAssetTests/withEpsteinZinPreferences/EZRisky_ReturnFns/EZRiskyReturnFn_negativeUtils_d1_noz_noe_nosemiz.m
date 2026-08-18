function F=EZRiskyReturnFn_negativeUtils_d1_noz_noe_nosemiz(h,savings,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)
% RiskyAsset, with d1 (h=labour), no z, no e, no semiz. UTILITY-UNITS Epstein-Zin, NEGATIVE-valued utility fn
% (vfoptions.EZutils=1, EZpositiveutility=0): CRRA of the composite consumption-leisure
% good x=(c^varphi)*((1-h)^(1-varphi)), without the -1 shift, so strictly negative for
% ezsigma>1. Doubles as the vNM reference utility fn for the gamma=1/phi collapse test
% (at ezsigma=ezgamma).

F=-Inf;

if agej<Jr
    c=w*kappa_j*h + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    x=(c^varphi)*((1-h)^(1-varphi));
    F=(x^(1-ezsigma))/(1-ezsigma);
end

end
