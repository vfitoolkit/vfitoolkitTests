function F=EZReturnFn_positiveUtils_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,ezsigma,agej,Jr,pension,uempbenefit,searcheffortcost)
% Utility-units, POSITIVE-valued utility fn for Epstein-Zin
% (vfoptions.EZutils=1, vfoptions.EZpositiveutility=1).
% The (1+x) shift keeps F strictly positive for all x>0 [a plain
% (x^(1-ezsigma)-1)/(1-ezsigma) would change sign at x=1].
% Note: requires searcheffortcost<1 so the (1-searcheffortcost*d2) factor stays positive.

F=-Inf;

if agej<Jr
    c=(1+r)*a+w*kappa_j*z*e*semiz+uempbenefit*(1-semiz)-aprime;
else
    c=(1+r)*a+pension-aprime;
end

if c>0
    x=c*(1-searcheffortcost*d2);
    F=((1+x)^(1-ezsigma)-1)/(1-ezsigma);
end


end
