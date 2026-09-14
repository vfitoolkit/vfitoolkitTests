function T=GPTemptationFn_d1_noz_e_nosemiz(d1,d2,a1prime,a1,a2,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension)
% Temptation utility: consumption is tempting, v = lambdaGP*u_c(c) + shiftGP (no leisure
% term, so the temptation params deliberately differ from the ReturnFn params: no eta/varphi).
% Feasibility (-Inf) is identical to the matching ReturnFn; the budget is the same.

T=-Inf;

if agej<Jr
    c=(1+r)*a1+w*kappa_j*d1*d2*a2*e-a1prime;
else
    c=(1+r)*a1+pension-a1prime;
end

if c>0 && d1<1
    T=lambdaGP*(c^(1-sigma)-1)/(1-sigma)+shiftGP;
end

end
