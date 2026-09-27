function F=ReturnToExitFn_EE_aprime(aprime,a,z,tau,minexit)
% The endogenousexit=2 signature for the return to exit is (d,aprime,a,z), NOT the (a,z) of
% endogenousexit=1: that path builds it with CreateReturnFnMatrix_Disc rather than
% CreateReturnToExitFnMatrix_Case1_Disc_Par2, and then maximises over (d,aprime) to get
% PolicyWhenExit. This subcode has no d, and CreateReturnFnMatrix_Disc drops the d input
% entirely when n_d=0, so here the signature is (aprime,a,z).
%
% This version ignores aprime and returns exactly ReturnToExitFn_EE, so the maximisation over it
% is a tie and the value matches endogenousexit=1. That is what makes the exit2-vs-exit1
% cross-test an exact comparison rather than an approximate one.

F=-Inf;

if a>=minexit
    F=-tau*a+0*aprime;
end

end
