% Setup for P6: life-cycle (age-dependent) AR(1) with normally distributed innovations
%
%    z(j) = mew(j) + rho(j)*z(j-1) + e(j),   e(j) ~ N(0,sigma(j)^2)
%    with z(0)=0 by default, so z(1)=mew(1)+e(1)
%
% Two calibrations. The AGE-VARYING one has all three parameters genuinely varying with age, so a
% command that quietly froze any of them would fail. The FROZEN one has all three constant, and is
% what the frozen-life-cycle identity in the cross-tests is built on. The frozen one is DRIFTLESS,
% because discretizeLifeCycleAR1_FellaGallipoliPan takes no mew argument at all and the three-way
% identity cannot otherwise be stated.

calibLC=struct();
calibLC.J=41;                  % a working life, ages 1 to 41
J=calibLC.J;
agevec=(1:J); % used only to build the age profiles below

%% The age-varying calibration
calibLC.vary.rho=0.98-0.08*(agevec-1)/(J-1);            % persistence declining 0.98 to 0.90
calibLC.vary.sigma=0.15+0.10*sin(pi*(agevec-1)/(J-1));  % hump-shaped innovation std dev
calibLC.vary.mew=0.05-0.10*(agevec-1)/(J-1);            % a drift that changes sign mid-life
calibLC.vary.name='age-varying';

%% The frozen calibration (constant parameters, no drift)
calibLC.frozen.rho=0.90*ones(1,J);
calibLC.frozen.sigma=0.20*ones(1,J);
calibLC.frozen.mew=zeros(1,J);
calibLC.frozen.name='frozen';
% The stationary distribution of the frozen process, which is what the initialj0 options have to be
% set to for the identity to hold at EVERY age rather than only in the limit.
calibLC.frozen.statsigmaz=0.20/sqrt(1-0.90^2);
calibLC.frozen.statmewz=0;

%% Analytic truth: the age recursion
% mewz(j)   = mew(j) + rho(j)*mewz(j-1)
% sigmaz(j)^2 = rho(j)^2*sigmaz(j-1)^2 + sigma(j)^2
% autocorr between j-1 and j = rho(j)*sigmaz(j-1)/sigmaz(j)
%
% This is recomputed inside every subcode that uses it, NOT read from otheroutputs.sigma_z. The
% commands return otheroutputs.sigma_z, and comparing the two is a real check; using the command's
% own answer as the truth would make it a tautology.
cnamesLC={'vary','frozen'};
for c_c=1:2
    cn=cnamesLC{c_c};
    mew=calibLC.(cn).mew; rho=calibLC.(cn).rho; sigma=calibLC.(cn).sigma;
    mewz=zeros(1,J); sigmaz=zeros(1,J);
    mewz(1)=mew(1); sigmaz(1)=sigma(1);   % from z(0)=0, the default initial condition
    for j_c=2:J
        mewz(j_c)=mew(j_c)+rho(j_c)*mewz(j_c-1);
        sigmaz(j_c)=sqrt(rho(j_c)^2*sigmaz(j_c-1)^2+sigma(j_c)^2);
    end
    aclc=nan(1,J);
    for j_c=2:J
        aclc(j_c)=rho(j_c)*sigmaz(j_c-1)/sigmaz(j_c);
    end
    calibLC.(cn).mewz=mewz;
    calibLC.(cn).sigmaz=sigmaz;
    calibLC.(cn).autocorr=aclc;
    fprintf('\nP6 calibration %s: sigmaz runs %2.4f at age 1 to %2.4f at age %i, peak %2.4f \n',cn,sigmaz(1),sigmaz(J),J,max(sigmaz))
    fprintf('P6 calibration %s: mewz runs %+2.4f at age 1 to %+2.4f at age %i \n',cn,mewz(1),mewz(J),J)
end

% The frozen calibration's own recursion must approach the stationary standard deviation, which is
% what makes the transient form of the identity work at the last age. Checked rather than assumed.
fprintf('P6 frozen calibration: sigmaz(J)=%2.6f against the stationary %2.6f, difference %2.3e \n',calibLC.frozen.sigmaz(J),calibLC.frozen.statsigmaz,abs(calibLC.frozen.sigmaz(J)-calibLC.frozen.statsigmaz))
if abs(calibLC.frozen.sigmaz(J)-calibLC.frozen.statsigmaz)>10^(-4)
    warning('P6: J is not large enough for the frozen calibration to have converged; the transient form of the identity will be loose')
end

% The age-varying calibration must actually vary, or the whole block is vacuous
if max(calibLC.vary.sigmaz)/min(calibLC.vary.sigmaz)<1.5
    warning('P6 age-varying calibration barely varies with age; the age-dependent methods will look no better than a stationary one')
end

%% Grid sizes
% znums straddles 17 deliberately. discretizeAR1_Rouwenhorst uses a grid half-width of exactly
% sqrt(znum-1)*sigmaz, with no cap, and that exact spread is what makes Rouwenhorst match the
% variance to machine precision. discretizeLifeCycleAR1_FellaGallipoliPan defaults nSigmas to
% min(sqrt(znum-1),4), so the cap binds for znum>17. Whether that costs it the exact variance match
% is measured in DiscP6_LCAR1_FGP; sqrt(znum-1)=4 exactly at znum=17, so 17 is the last grid size
% where the two agree by construction.
znums=[5,9,15,17,31,51];
znumsEven=[6,10,16,32];   % B15: the FGP methods must work for even znum too
calibLC.nreps=3;
calibLC.timethreshold=0.5;
calibLC.entropytol=10^(-7);
calibLC.Jbig=200;         % for the transient form of the frozen identity
