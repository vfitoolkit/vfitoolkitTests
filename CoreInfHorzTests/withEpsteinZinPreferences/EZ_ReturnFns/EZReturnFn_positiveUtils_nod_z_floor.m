function F=EZReturnFn_positiveUtils_nod_z_floor(aprime,a,z,r,w,ezsigma,xbar,floorval)
% EZReturnFn_positiveUtils_nod_z with a subsistence floor: consumption below xbar is FEASIBLE but
% yields utility floorval (0, or a tiny positive number for the continuity comparison).
% Used by the feasibility cross-test: a return of EXACTLY zero at a feasible choice.

F=-Inf;

c=(1+r)*a+w*z-aprime;

if c>0
    if c>=xbar
        F=((1+c)^(1-ezsigma)-1)/(1-ezsigma);
    else
        F=floorval;
    end
end


end
