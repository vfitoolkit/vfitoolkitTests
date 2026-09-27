function F=ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap)
% Firm profit. a is previous employment, aprime is employment chosen for this period.
%
% The a>empcap*z gate makes EVERY aprime infeasible at those states, so the max over aprime is
% -Inf there. That is deliberate: it is what forces the exit decision to weigh a finite exit
% value against a -Inf continuation, which is the case the exit arithmetic has to get right.

F=-Inf;

if a<=empcap*z
    F=p*z*(aprime^alpha) - aprime - tau*max(0,a-aprime) - cf;
end

end
