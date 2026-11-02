function F=EZReturnFn_cons_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,agej,Jr,pension)
% Consumption-units composite for traditional Epstein-Zin (vfoptions.EZutils=0):
% F is the composite consumption good (strictly positive whenever feasible).

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*e-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=c;
end


end
