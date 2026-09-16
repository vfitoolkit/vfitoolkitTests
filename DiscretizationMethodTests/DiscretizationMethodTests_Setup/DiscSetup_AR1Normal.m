% Setup for P2: AR(1) with normal innovations
%
%    z' = mew + rho*z + e,   e ~ N(0,sigma^2)
%
% Convention: mew is the INTERCEPT, so E(z)=mew/(1-rho) and Var(z)=sigma^2/(1-rho^2). All four
% stationary AR(1) commands follow this after the grid-centring fix; discretizeAR1wGM_FarmerToda
% is the one command still on the other reading (mew = the mean), and P2's cross-tests are what
% keep asking about it.

calibAR1=struct();

%% Three calibrations, all run
% moderate: where Farmer-Toda should beat Tauchen
calibAR1.moderate.mew=0;    calibAR1.moderate.rho=0.6;  calibAR1.moderate.sigma=0.2;
% persistent: where Rouwenhorst should win (Farmer & Toda 2017, last para of pg 678, which
% discretizeAR1_FarmerToda itself prints a comment about for rho>=0.99)
calibAR1.persistent.mew=0;  calibAR1.persistent.rho=0.99; calibAR1.persistent.sigma=0.1;
% drift: the mew-convention probe. At mew=0 a command that centred its grid on mew rather than on
% mew/(1-rho) is indistinguishable from a correct one, which is how that bug survived. Here
% E(z)=0.1/0.1=1.0, so the two conventions are an order of magnitude apart.
calibAR1.drift.mew=0.1;     calibAR1.drift.rho=0.9;     calibAR1.drift.sigma=0.2;

calibAR1.names={'moderate','persistent','drift'};

%% ONE znum sweep, serving both accuracy and runtime
znums=[5,9,15,31,51,101];

%% Timing protocol (inlined in every method subcode, per the no-helper-functions rule)
calibAR1.nreps=5;
calibAR1.timethreshold=0.5;

%% Tolerances that the methods actually promise, measured in P1 rather than assumed
% Rouwenhorst and Tauchen-Hussey reach their moments algebraically or by quadrature exactness.
calibAR1.exacttol=10^(-12);
% Farmer-Toda reaches its moments through an entropy solve, whose accuracy degrades with the size
% of the problem: measured in P1, 3e-17 at znum=5 rising to 1.3e-09 at znum=101.
calibAR1.entropytol=10^(-7);

%% Floden (2008) Table 1, as an external oracle
% Floden, Martin (2008), "A note on the accuracy of Markov-chain approximations to highly
% persistent AR(1) processes", Economics Letters 99(3), 516-520. doi:10.1016/j.econlet.2007.09.040
%
% His three calibrations, and the Tauchen column of his Table 1 (n = 5, 9, 15).
% Floden spaces the Tauchen nodes over +-1.2*sigma_z*log(n), NOT the toolkit's +-Tauchen_q*sigma_z,
% so the toolkit has to be called with Tauchen_q=1.2*log(znum) to reproduce him. His z^n/sigma_z
% row is the built-in check that this is right: it prints 1.9313, 2.6367, 3.2497, which is exactly
% 1.2*log(n).
%
% TRANSCRIPTION STATUS, and this matters for how the checks below are written:
%   sigma_e and sigma_z rows  - CONFIRMED. Transcribed from the paper AND independently recomputed
%                               from Floden's description of his method; the two agree to all four
%                               printed decimals in all nine cells.
%   rho row                   - NOT CONFIRMED. The independent recomputation is systematically
%                               lower than the transcription, by 0.0004 to 0.0248, shrinking as n
%                               grows. Either the transcription is wrong or Floden computes the
%                               implied autocorrelation by something other than
%                               cov(z,z')/var(z) under the stationary distribution. So rho is
%                               REPORTED here, not asserted, until someone checks the paper.
calibAR1.floden.name={'Aiyagari','HubbardSkinnerZeldes','StoreslettenTelmerYaron'};
calibAR1.floden.rho=[0.60,0.95,0.98];
calibAR1.floden.sigmasq_e=[0.013,0.030,0.020];
calibAR1.floden.nlist=[5,9,15];
% Tauchen column, rows sigma_e and sigma_z, indexed (process, n)
% Note on the (2,1) entry, Hubbard-Skinner-Zeldes at n=5: the PDF render reads 0.1843, but BOTH
% independent computations - the toolkit, and a from-scratch reimplementation of Floden's method
% written before the toolkit was run - give 0.18441. Two computations against one hand-transcription
% from a four-decimal table, so 0.1844 is recorded here. Worth a human check against the paper.
calibAR1.floden.tauchen_sigma_e=[0.1167,0.1165,0.1155; 0.1844,0.1982,0.1883; 0.0838,0.1466,0.1634];
calibAR1.floden.tauchen_sigma_z=[0.1430,0.1451,0.1443; 0.6037,0.6205,0.5995; 0.7938,0.8448,0.8306];
% rho, transcribed but NOT confirmed - reported only
calibAR1.floden.tauchen_rho_unconfirmed=[0.5844,0.6000,0.5998; 0.9770,0.9503,0.9499; 0.9952,0.9861,0.9810];
% the z^n/sigma_z row, which is the transcription self-check
calibAR1.floden.znoversigmaz=[1.9313,2.6367,3.2497];
calibAR1.floden.tol=10^(-4); % the table is printed to four decimals
% Floden's z^n/sigma_z row is printed to four decimals too, so comparing it against the exact
% 1.2*log(n) can only ever agree to half a unit in the last place. 5e-5, not a T1 tolerance.
calibAR1.floden.roundtol=5*10^(-5);
