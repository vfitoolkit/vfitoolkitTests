function F=ReturnToExitFn_EE(a,z,tau,minexit)
% Value of shutting down: pay the firing cost on every remaining worker.
% Inputs are (a,z) followed by parameters; this is the endogenousexit=1 signature.
%
% The a<minexit gate makes exit itself infeasible at those states. That is deliberate: it is
% the mirror of the empcap gate in the return function, and it is what puts a -Inf on the OTHER
% side of the exit decision so that both branches of the binary weight are exercised.

F=-Inf;

if a>=minexit
    F=-tau*a;
end

end
