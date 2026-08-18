function F=EZRiskyReturnFn_cons_d1_noz_noe_semiz(h,savings,dsemiz,a,semiz,r,w,kappa_j,varphi,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, with d1 (h=labour), no z, no e. CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0):
% the return fn is the composite consumption good (curvature comes from the EZ preferences).
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
    F=((c^varphi)*((1-h)^(1-varphi)))*(1-searcheffortcost*dsemiz);
end

end
