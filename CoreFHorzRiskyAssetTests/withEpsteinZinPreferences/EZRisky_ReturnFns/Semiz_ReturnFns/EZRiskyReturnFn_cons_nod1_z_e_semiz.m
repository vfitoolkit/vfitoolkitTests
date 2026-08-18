function F=EZRiskyReturnFn_cons_nod1_z_e_semiz(savings,dsemiz,a,semiz,z,e,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)
% RiskyAsset + semiz, no d1, z, e. CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0):
% the return fn is the composite consumption good (curvature comes from the EZ preferences).
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
    F=c*(1-searcheffortcost*dsemiz);
end

end
