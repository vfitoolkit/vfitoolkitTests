function F=ReturnFn_EE_d(d,aprime,a,z,p,alpha,tau,cf,empcap,psi)
% The 'with d' twin of ReturnFn_EE_nod. d is utilisation, costly to move away from 1.
%
% At d=1 the utilisation cost is zero and this is EXACTLY ReturnFn_EE_nod, which is what lets
% the d-vs-nod cross-test compare the two code paths as a machine zero.

F=-Inf;

if a<=empcap*z
    F=p*z*((d*aprime)^alpha) - aprime - tau*max(0,a-aprime) - cf - psi*((d-1)^2);
end

end
