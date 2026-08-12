function F=EZReturnFn_cons_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,agej,Jr,pension,uempbenefit,searcheffortcost)
% Consumption-units composite for traditional Epstein-Zin (vfoptions.EZutils=0):
% F is the composite consumption good (strictly positive whenever feasible).
% Note: requires searcheffortcost<1 so the (1-searcheffortcost*d2) factor stays positive.

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*semiz+uempbenefit*(1-semiz)-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    F=c*(1-searcheffortcost*d2);
end


end
