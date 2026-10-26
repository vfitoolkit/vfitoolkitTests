function F=EZReturnFn_negativeUtils_nod_z_floor(aprime,a,z,r,w,ezsigma,kappa,floorval)
% EZReturnFn_negativeUtils_nod_z shifted up by kappa and SATIATED at zero: utility is u(c)+kappa
% while that is negative, and floorval (0, or a tiny negative number for the continuity comparison)
% once it would reach zero. Used by the feasibility cross-test: a return of EXACTLY zero at a
% feasible choice, here at the TOP of the utility range (negative-valued utility fns have zero as
% their supremum, so this is the only place a zero can sit).

F=-Inf;

c=(1+r)*a+w*z-aprime;

if c>0
    u=(c^(1-ezsigma))/(1-ezsigma)+kappa;
    if u<0
        F=u;
    else
        F=floorval;
    end
end


end
