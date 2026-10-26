function F=EZReturnFn_negativeUtils_nod_z_floor(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,kappa,floorval)
% EZReturnFn_negativeUtils_nod_z_noe_nosemiz shifted up by kappa and SATIATED at zero: utility is
% u(c)+kappa while that is negative, and floorval (0, or a tiny negative number for the continuity
% comparison) once it would reach zero. Used by the feasibility cross-test: a return of EXACTLY zero
% at a feasible choice, here at the TOP of the utility range (negative-valued utility fns have zero
% as their supremum, so this is the only place a zero can sit).

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    u=(c^(1-ezsigma))/(1-ezsigma)+kappa;
    if u<0
        F=u;
    else
        F=floorval;
    end
end


end
