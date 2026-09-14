function F=ReturnFn_d1d2_a_z(d1,d2,a,z,r,w,sigma,eta,varphi)
% Inheritance asset return function. Case-2 style: there is NO aprime input, because next
% period's asset is determined by d2 and the shock transition, not chosen directly.
%   d1 is labour supply
%   d2 is the resources placed into the inheritance asset
%   a  is the inheritance asset carried in

F=-Inf;

c=w*z*d1+(1+r)*a-d2;

if c>0
    F=(c^(1-sigma))/(1-sigma) - varphi*(d1^(1+eta))/(1+eta);
end

end
