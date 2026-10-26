function F=EZReturnFn_cons_nod_z(aprime,a,z,r,w)
% Consumption-units composite for traditional Epstein-Zin (vfoptions.EZutils=0):
% F is the composite consumption good (strictly positive whenever feasible).
% InfHorz form of EZReturnFn_cons_nod_z_noe_nosemiz (no ages, no retirement).

F=-Inf;

c=(1+r)*a+w*z-aprime;

if c>0
    F=c;
end


end
