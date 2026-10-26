function F=EZReturnFn_cons_nod_z_floor(aprime,a,z,r,w,kappa_j,agej,Jr,pension,xbar,floorval)
% EZReturnFn_cons_nod_z_noe_nosemiz with a subsistence floor: consumption below xbar is FEASIBLE but
% yields the composite good floorval (0, or a tiny positive number for the continuity comparison).
% Used by the feasibility cross-test: a return of EXACTLY zero at a feasible choice.

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    if c>=xbar
        F=c;
    else
        F=floorval;
    end
end


end
