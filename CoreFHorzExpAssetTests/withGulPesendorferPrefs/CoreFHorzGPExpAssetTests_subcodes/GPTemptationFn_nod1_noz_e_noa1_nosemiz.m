function T=GPTemptationFn_nod1_noz_e_noa1_nosemiz(d2,a,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP.
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.

T=-Inf;

if agej<Jr
    c=w*kappa_j*d2*a*e;
else
    c=pension;
end

if c>0
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
