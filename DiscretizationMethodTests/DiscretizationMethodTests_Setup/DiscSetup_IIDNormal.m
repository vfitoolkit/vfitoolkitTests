% Setup for P1: iid normal, e~N(mew,sigma^2)

calibIID=struct();
calibIID.mew=0.5;      % deliberately non-zero, see below
calibIID.sigma=0.3;

% mew=0.5 rather than 0 on purpose. At mew=0 the mean check passes for free, and the symmetry
% check (grid symmetric about the mean, probabilities centrosymmetric) degenerates to symmetry
% about zero, where a command that centred its grid on the wrong thing would be indistinguishable
% from a correct one. That is how the Rouwenhorst grid-centring bug survived as long as it did.

% ONE znum sweep, serving both the accuracy comparison and the runtime table: the timed calls are
% the calls whose output the accuracy analysis reads. Discretization is deterministic, so the
% repeats all return identical (grid,pi) and the last is simply kept.
znums=[5,9,15,31,51,101];

% Timing protocol (identical in every method subcode in this bank; it is inlined rather than
% shared, per the toolkit's no-helper-functions rule):
%   one warm-up call, discarded, then nreps timed calls, report the MEDIAN
%   nreps drops to 1 for any config whose warm-up exceeded timethreshold seconds
calibIID.nreps=5;
calibIID.timethreshold=0.5;
