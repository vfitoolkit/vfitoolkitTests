function F=EZRiskyReturnFn_cons_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,agej,Jr,pension)
% RiskyAsset, no d1, z, e, no semiz. CONSUMPTION-UNITS Epstein-Zin (vfoptions.EZutils=0):
% the return fn is consumption itself (curvature comes from the EZ preferences).

F=-Inf;

if agej<Jr
    c=w*kappa_j*z*e + a - savings;
else
    c=pension + a - savings;
end

if c>0
    F=c;
end

end
