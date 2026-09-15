function F=CoreInfHorzTPathAlgo_jointdgrid_ReturnFn(d1,d2,aprime,a,z,r,w,sigma)
% Two decision variables, with d1+d2>1 declared infeasible by returning -Inf. Those combinations can
% therefore never be chosen, whatever the prices, which is what lets the joint-grid cross test drop
% them and still expect an identical answer.
F=-Inf;
if d1+d2<=1
    c=w*z*(d1+d2)+(1+r)*a-aprime;
    if c>0
        F=(c^(1-sigma))/(1-sigma);
    end
end
end
