% Setup for P0: validating MarkovChainMoments and MarkovChainMoments_FHorz
%
% P0 builds its own chains rather than reusing later blocks', so that it can run first.
% Three kinds are needed:
%   (i)  a chain whose moments are known in CLOSED FORM, so the check is exact-ish rather than
%        an approximation of anything. Two of these: a hand-built two-state chain, and a
%        Rouwenhorst chain (whose mean, variance and autocorrelation are exact by construction).
%  (ii)  a hand-built AGE-DEPENDENT chain, for the same reason, for the _FHorz instrument.
% (iii)  representative chains from real discretization commands, for the panel validation.

calibP0=struct();

%% (i) A two-state Markov chain, everything in closed form
% pi=[1-p, p; q, 1-q] on z=[z1,z2]
calibP0.two.p=0.2;
calibP0.two.q=0.3;
calibP0.two.z1=-0.4;
calibP0.two.z2=0.9;
% Deliberately asymmetric (p~=q, and the grid not centred on zero) so that a wrong stationary
% distribution cannot be hidden by symmetry.

%% (i) A Rouwenhorst chain: mean, variance and autocorrelation are matched exactly for any znum
calibP0.rouw.mew=0.1;   % note: mew is the INTERCEPT, so E(z)=mew/(1-rho)
calibP0.rouw.rho=0.85;
calibP0.rouw.sigma=0.2;
calibP0.rouw.znum=9;
% mew~=0 deliberately: at mew=0 a command that centred its grid on mew rather than on
% mew/(1-rho) would be indistinguishable from a correct one.

%% (ii) A hand-built age-dependent chain, for MarkovChainMoments_FHorz
% Same 3-point grid at every age, a different (known) transition matrix at each age, and an
% initial distribution that is not uniform. Everything is recomputed by an explicit loop in the
% subcode, so the check does not go through any toolkit code twice.
calibP0.agedep.znum=3;
calibP0.agedep.N_j=6;
calibP0.agedep.z_grid=[-1;0.25;1.5];
calibP0.agedep.jequaloneDistz=[0.5;0.3;0.2];

%% (iii) Representative chains for the panel validation
% A Farmer-Toda AR(1), a Rouwenhorst AR(1) at high persistence, and a deliberately asymmetric one
calibP0.panel.mew=0;
calibP0.panel.rho=0.7;
calibP0.panel.sigma=0.25;
calibP0.panel.znum=7;
calibP0.panel.rho_persistent=0.95;
% Life-cycle chain for the FHorz half
calibP0.panelJ.N_j=20;
calibP0.panelJ.znum=7;
calibP0.panelJ.mew=zeros(1,calibP0.panelJ.N_j);
calibP0.panelJ.rho=0.9*ones(1,calibP0.panelJ.N_j);
calibP0.panelJ.sigma=0.2+0.1*sin(pi*(1:calibP0.panelJ.N_j)/calibP0.panelJ.N_j); % hump in the innovation std dev

%% Panel size (D8)
% se(sigmahat)=sigma/sqrt(2N), so N=10^5 resolves ~0.9% of sigma_z at the 4*se band. Every check
% prints its achieved se, so the diary records what the block was actually able to detect.
calibP0.numbersims=10^5;
calibP0.simperiods=50;   % InfHorz panel length per agent
calibP0.kse=4;           % width of the band, in standard errors

%% The trivial n_a=1 model used by the panel validation (and, with more asset points, by the
%% downstream subcodes). With a single asset point there is nothing to choose, so a simulated
%% panel of z is a pure Monte Carlo draw from the discretized chain.
calibP0.Params.beta=0.9;
