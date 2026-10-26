function F=EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi)
% Consumption-units composite for traditional Epstein-Zin (vfoptions.EZutils=0):
% F is the composite consumption-leisure good (strictly positive whenever feasible).
% InfHorz form of EZReturnFn_cons_d_z_noe_nosemiz (no ages, no retirement).

F=-Inf;

c=(1+r)*a+w*d*z-aprime;

if c>0 && d<1
    F=(c^varphi)*((1-d)^(1-varphi));
end


end
