function F=EZRiskyReturnFn_cons_d1_noz_noe_nosemiz(h,savings,a,r,w,kappa_j,varphi,agej,Jr,pension)
% RiskyAsset, with d1 (h=labour), no z, no e, no semiz. CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0):
% the return fn is the composite consumption-leisure good x=(c^varphi)*((1-h)^(1-varphi))
% itself (curvature comes from the EZ preferences). Keeps the riskyasset budget: asset
% returns are realised via aprimeFn, so 'a' enters the budget directly (no (1+r)*a term).

F=-Inf;

if agej<Jr
    c=w*kappa_j*h + a - savings;
else
    c=pension + a - savings;
end

if c>0 && h<1
    F=(c^varphi)*((1-h)^(1-varphi));
end

end
