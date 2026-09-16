% Setup for P8: life-cycle (age-dependent) VAR(1) with normally distributed innovations
%
%    Z(j) = Mew(:,j) + Rho(:,:,j)*Z(j-1) + e(j),   e(j) ~ N(0,SigmaSq(:,:,j))
%    with Z(0)=0 by default, so Z(1) ~ N(Mew(:,1),SigmaSq(:,:,1))
%
% P8 has ONE command, discretizeLifeCycleVAR1_Tauchen, which is the last user-facing command in
% DiscretizationMethods with no coverage at all. That shapes the block: with no sibling method there
% is no accuracy race to run, so the weight falls on the analytic age profiles and on the cross-test
% identities, where discretizeVAR1_Tauchen and discretizeLifeCycleAR1_FellaGallipoliPanTauchen act
% as independent implementations of special cases this command must reproduce.
%
% THREE CALIBRATIONS, all M=2.
%   vary      every one of Mew, Rho and SigmaSq genuinely age-dependent, and Rho has non-zero
%             off-diagonals at every age, so a command that quietly froze a parameter or dropped the
%             cross terms would fail. This is the one the accuracy sweep runs on.
%   frozen    constant parameters, no drift. The frozen-life-cycle identity is built on this: with
%             constant parameters the life-cycle command must collapse onto discretizeVAR1_Tauchen.
%   diagonal  Rho and SigmaSq both diagonal, so the two variables are independent and the joint
%             chain must equal the kron of two univariate life-cycle chains. This is the check that
%             the joint transition matrix is assembled in the right order; P5 found that the two
%             stationary VAR commands disagree about that ordering, so it is not a formality.

calibLCV=struct();
calibLCV.J=25;                % shorter than P6's 41: the joint matrix is prod(znum) square
J=calibLCV.J;
calibLCV.M=2;
M=calibLCV.M;
tvec=(0:J-1)/(J-1);           % 0 at age 1, 1 at age J

%% The age-varying calibration
Mew_J=zeros(M,J); Rho_J=zeros(M,M,J); SigmaSq_J=zeros(M,M,J);
for j_c=1:J
    t=tvec(j_c);
    Mew_J(:,j_c)=[0.04-0.08*t; -0.02+0.05*t];              % drifts that change sign mid-life
    Rho_J(:,:,j_c)=[0.95-0.10*t, 0.05*t; 0.02, 0.88+0.04*t]; % off-diagonals non-zero, and one of them grows from zero
    s1=0.15+0.05*sin(pi*t); s2=0.12+0.04*t;                 % one hump-shaped, one monotone
    SigmaSq_J(:,:,j_c)=[s1^2, 0.3*s1*s2; 0.3*s1*s2, s2^2];  % correlated innovations
end
calibLCV.vary.Mew_J=Mew_J; calibLCV.vary.Rho_J=Rho_J; calibLCV.vary.SigmaSq_J=SigmaSq_J;
calibLCV.vary.name='age-varying';

%% The frozen calibration
Rho_f=[0.90, 0.05; 0.00, 0.85];
SigmaSq_f=[0.04, 0.01; 0.01, 0.03];
calibLCV.frozen.Mew_J=zeros(M,J);
calibLCV.frozen.Rho_J=repmat(Rho_f,[1,1,J]);
calibLCV.frozen.SigmaSq_J=repmat(SigmaSq_f,[1,1,J]);
calibLCV.frozen.Rho=Rho_f; calibLCV.frozen.SigmaSq=SigmaSq_f;
calibLCV.frozen.name='frozen';
% The stationary covariance solves Sigma = Rho*Sigma*Rho' + SigmaSq, i.e.
% vec(Sigma) = (I-kron(Rho,Rho))\vec(SigmaSq). This is what initialj1SigmaSqz has to be set to for
% the identity to hold at EVERY age rather than only asymptotically.
calibLCV.frozen.statSigmaSqz=reshape((eye(M^2)-kron(Rho_f,Rho_f))\SigmaSq_f(:),[M,M]);
calibLCV.frozen.statmewz=zeros(M,1); % Mew is zero, so the stationary mean is too

%% The diagonal calibration
rho_d=[0.92; 0.80]; sig_d=[0.18; 0.14];
calibLCV.diag.Mew_J=zeros(M,J);
calibLCV.diag.Rho_J=repmat(diag(rho_d),[1,1,J]);
calibLCV.diag.SigmaSq_J=repmat(diag(sig_d.^2),[1,1,J]);
calibLCV.diag.rho=rho_d; calibLCV.diag.sigma=sig_d;
calibLCV.diag.name='diagonal';

%% Analytic truth: the age recursion, in matrix form
%   mewz(:,1)     = Mew(:,1) + Rho(:,:,1)*z0,  z0=0
%   SigmaSqz(1)   = SigmaSq(:,:,1)
%   mewz(:,j)     = Mew(:,j) + Rho(:,:,j)*mewz(:,j-1)
%   SigmaSqz(j)   = Rho(:,:,j)*SigmaSqz(j-1)*Rho(:,:,j)' + SigmaSq(:,:,j)
% and the lag-one cross moment, which gives the per-variable autocorrelation:
%   Cov(Z_j,Z_{j-1}) = Rho(:,:,j)*SigmaSqz(j-1)
%
% Recomputed here and again in the subcodes that use it, NOT read from otheroutputs. The command
% returns otheroutputs.mewz and otheroutputs.SigmaSqz; comparing those against this is a real check,
% and using the command's own answer as the truth would make every accuracy number a tautology.
cnamesLCV={'vary','frozen','diag'};
for c_c=1:3
    cn=cnamesLCV{c_c};
    Mw=calibLCV.(cn).Mew_J; Rh=calibLCV.(cn).Rho_J; Sg=calibLCV.(cn).SigmaSq_J;
    mewz=zeros(M,J); SigmaSqz=zeros(M,M,J);
    mewz(:,1)=Mw(:,1); SigmaSqz(:,:,1)=Sg(:,:,1);
    for j_c=2:J
        mewz(:,j_c)=Mw(:,j_c)+Rh(:,:,j_c)*mewz(:,j_c-1);
        SigmaSqz(:,:,j_c)=Rh(:,:,j_c)*SigmaSqz(:,:,j_c-1)*Rh(:,:,j_c)'+Sg(:,:,j_c);
    end
    sigmaz=zeros(M,J); acz=nan(M,J); rhoz=nan(1,J);
    for j_c=1:J
        sigmaz(:,j_c)=sqrt(diag(SigmaSqz(:,:,j_c)));
        rhoz(j_c)=SigmaSqz(1,2,j_c)/(sigmaz(1,j_c)*sigmaz(2,j_c)); % contemporaneous correlation between the two variables
    end
    for j_c=2:J
        Clag=Rh(:,:,j_c)*SigmaSqz(:,:,j_c-1); % Cov(Z_j,Z_{j-1})
        acz(:,j_c)=diag(Clag)./(sigmaz(:,j_c).*sigmaz(:,j_c-1));
    end
    calibLCV.(cn).mewz=mewz; calibLCV.(cn).SigmaSqz=SigmaSqz;
    calibLCV.(cn).sigmaz=sigmaz; calibLCV.(cn).autocorr=acz; calibLCV.(cn).crosscorr=rhoz;
    % A calibration whose VAR is explosive would make every age profile diverge and every accuracy
    % number meaningless, so check it here rather than discovering it in the results.
    maxeig=0;
    for j_c=1:J
        maxeig=max(maxeig,max(abs(eig(Rh(:,:,j_c)))));
    end
    fprintf('\nP8 calibration %s: largest |eig(Rho)| over the ages is %2.4f, this should be below one: %i \n',cn,maxeig,maxeig<1)
    fprintf('P8 calibration %s: sd of variable 1 runs %2.4f to %2.4f, variable 2 %2.4f to %2.4f \n',cn,sigmaz(1,1),sigmaz(1,J),sigmaz(2,1),sigmaz(2,J))
    fprintf('P8 calibration %s: contemporaneous correlation between the variables runs %+2.4f to %+2.4f \n',cn,rhoz(1),rhoz(J))
end

% The frozen recursion must converge TOWARD the stationary covariance computed above - two
% independent routes to the same object, one iterative and one a linear solve. It does not reach it
% at finite J: starting from SigmaSqz(1)=SigmaSq, the deviation decays like the largest eigenvalue of
% kron(Rho,Rho), which is max(|eig(Rho)|)^2, so at age J it is still of order that to the J-1. The
% bar below is set from that prediction with an order of magnitude of slack, and the predicted value
% is printed next to the measured one so the two can be compared rather than just passed.
frozgap=max(abs(calibLCV.frozen.SigmaSqz(:,:,J)-calibLCV.frozen.statSigmaSqz),[],'all');
frozpred=max(abs(calibLCV.frozen.SigmaSqz(:,:,1)-calibLCV.frozen.statSigmaSqz),[],'all')*max(abs(eig(Rho_f)))^(2*(J-1));
fprintf('\nP8 frozen: age-J covariance against the stationary solve, gap %2.3e against the %2.3e the decay rate predicts \n',frozgap,frozpred)
fprintf('P8 frozen: that gap is below the bar the decay rate implies [T2], this should be below %g: %2.3e \n',10*frozpred,frozgap)

% P8 CANNOT USE THE SAME znums AS THE UNIVARIATE BLOCKS. Everywhere else the cost is linear in
% znum; here the joint state count is znum^M and the transition matrix is
% (znum^M)-by-(znum^M)-by-(J-1), so P2's top entry of 101 would mean a 10201-by-10201 matrix at each
% of 24 ages. This list tops out with the joint chain at 289 states - a real solve, and a
% manageable one. It is also why J is 25 here rather than P6's 41.
znums=[5,7,9,13,17];

calibLCV.timethreshold=5; calibLCV.nreps=3; calibLCV.entropytol=10^(-7);
