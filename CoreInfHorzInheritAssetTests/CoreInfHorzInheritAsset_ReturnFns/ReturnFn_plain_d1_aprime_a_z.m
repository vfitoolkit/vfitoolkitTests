function F=ReturnFn_plain_d1_aprime_a_z(d1,aprime,a,z,r,w,sigma,eta,varphi)
% The plain (Case-1) twin of ReturnFn_d1d2_a_z, for the cross-test.
%
% In the cross-test the inheritance asset is de-risked (aprimeFn returns d2 itself, ignoring the
% shock transition) and d2_grid is set equal to a_grid. Then choosing d2 IS choosing aprime, and
% the inheritance-asset model is literally an ordinary one-asset model. This return function is
% the same as the inheritance one with d2 renamed to aprime, so the two must agree exactly.

F=-Inf;

c=w*z*d1+(1+r)*a-aprime;

if c>0
    F=(c^(1-sigma))/(1-sigma) - varphi*(d1^(1+eta))/(1+eta);
end

end
