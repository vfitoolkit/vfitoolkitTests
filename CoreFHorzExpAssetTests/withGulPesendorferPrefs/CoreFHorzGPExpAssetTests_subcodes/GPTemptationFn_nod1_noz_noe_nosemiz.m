function T=GPTemptationFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP.
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d2*a2-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
