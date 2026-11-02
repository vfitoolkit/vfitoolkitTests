function F=EZReturnFn_cons_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension)
% Consumption-units composite for traditional Epstein-Zin (vfoptions.EZutils=0):
% F is the composite consumption good (strictly positive whenever feasible).

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z*e*d-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0 && d<1
    F=(c^varphi)*((1-d)^(1-varphi));
end


end
