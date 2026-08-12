function F=EZReturnFn_cons_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,kappa_j,agej,Jr,pension,varphi,uempbenefit,searcheffortcost)
% Consumption-units composite for traditional Epstein-Zin (vfoptions.EZutils=0):
% F is the composite consumption good (strictly positive whenever feasible).
% Note: requires searcheffortcost<1 so the (1-searcheffortcost*d2) factor stays positive.

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*e*d1*semiz+uempbenefit*(1-semiz)-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0 && d1<1
    F=((c^varphi)*((1-d1)^(1-varphi)))*(1-searcheffortcost*d2);
end


end
