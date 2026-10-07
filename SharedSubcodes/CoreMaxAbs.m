function m=CoreMaxAbs(x)
% The largest absolute entry of x, with NaN entries counted. max() ignores NaN, so a check written as max(abs(a-b)) passes
% silently when a or b has a NaN where the other does not, or where both do (the 2026-10-07 EstimFHorzPTypeTests run passed its
% Jacobian checks against NaN moments this way). Here every NaN entry adds one, so a NaN in a difference fails a "this should
% be zero" check as a finite number, and the integer part says how many entries were NaN. An empty x gives zero.
% Used only where the two sides are NaN-free by design (moment vectors with the NaN targets dropped, Jacobians, parameters).
x=abs(x(:));
m=max([0;x(~isnan(x))])+sum(isnan(x));
end
