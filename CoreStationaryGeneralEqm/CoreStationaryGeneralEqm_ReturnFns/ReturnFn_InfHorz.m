function F=ReturnFn_InfHorz(d,aprime,a,z,r,tau,Tr,tau_c,alpha,delta,A,sigma,eta,varphi)
% InfHorz baseline GE model: endogenous labor d, asset a, productivity z.
% w is hardcoded from r via the firm's FOC (so w is NOT a separate GE price).
w=(1-alpha)*A*((r+delta)/(alpha*A))^(alpha/(alpha-1));
% Budget (tau_c is a consumption tax): (1+tau_c)*c+aprime=(1+r)*a+(1-tau)*w*z*d+Tr
c=((1+r)*a+(1-tau)*w*z*d+Tr-aprime)/(1+tau_c);

F=-Inf;
if c>0
    F=(c^(1-sigma))/(1-sigma)-varphi*(d^(1+eta))/(1+eta);
end

end
