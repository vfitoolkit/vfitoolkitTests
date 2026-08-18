function F=EZRiskyReturnFn_negativeUtils_d1_z_noe_nosemiz_withA1(h,savings,a1prime,a1,a2,z,r,w,kappa_j,ezsigma,varphi,r_a1,agej,Jr,pension) %#ok<INUSD>
% RiskyAsset with a1 (safe, earns r_a1) + a2 (risky, return realised via aprimeFn); d1, z, no e, no semiz.
% UTILITY-UNITS Epstein-Zin, NEGATIVE-valued utility fn (vfoptions.EZutils=1,
% EZpositiveutility=0): CRRA without the -1 shift, so strictly negative for ezsigma>1.
% Doubles as the vNM reference utility fn for the gamma=1/phi collapse test (at
% ezsigma=ezgamma).
% (r is unused here but kept in the signature for consistency with the vNM withA1 ReturnFns)

F=-Inf;

if agej<Jr
    c=w*kappa_j*h*z + (1+r_a1)*a1 + a2 - a1prime - savings;
else
    c=pension + (1+r_a1)*a1 + a2 - a1prime - savings;
end

if c>0 && h<1
    x=(c^varphi)*((1-h)^(1-varphi));
    F=(x^(1-ezsigma))/(1-ezsigma);
end

end
