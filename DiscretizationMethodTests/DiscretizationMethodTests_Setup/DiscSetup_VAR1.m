% Setup for P5: VAR(1) with normally distributed innovations
%
%    z_t = Mew + Rho*z_{t-1} + e_t,   e_t ~ N(0,SigmaSq),   z_t an M-vector
%
% Two calibrations. The FULL one has off-diagonal terms in both Rho and SigmaSq, so the
% cross-covariance is a real target that a pair of independent AR(1)s cannot reach. The DIAGONAL
% one exists for the kron identity in the cross-tests, which is the only exact statement available
% about a multivariate chain.

calibVAR=struct();

%% The full M=2 calibration
calibVAR.full.Mew=[0.01;-0.02];
calibVAR.full.Rho=[0.9,0.05;0.1,0.8];
calibVAR.full.SigmaSq=[0.04,0.01;0.01,0.09];

%% The diagonal M=2 calibration (for the kron identity)
calibVAR.diag.Mew=[0.01;-0.02];
calibVAR.diag.Rho=diag([0.9,0.7]);
calibVAR.diag.SigmaSq=diag([0.04,0.09]);

%% An M=3 calibration, for the one point in the sweep that goes past two variables
calibVAR.three.Mew=[0.01;-0.02;0];
calibVAR.three.Rho=[0.8,0.05,0;0.1,0.7,0.05;0,0.1,0.6];
calibVAR.three.SigmaSq=[0.04,0.01,0;0.01,0.09,0.01;0,0.01,0.0625];

%% Analytic truth, for each calibration
% Mean:            zmean = (I-Rho)^(-1)*Mew
% Variance:        Var(z) = Rho*Var(z)*Rho' + SigmaSq, solved by vectorizing:
%                    vec(Var(z)) = (I - kron(Rho,Rho))^(-1) * vec(SigmaSq)
% Autocovariance:  Cov(z_t,z_{t-1}) = Rho*Var(z), so the autocorrelation of variable i is
%                    (Rho*Var(z))_{ii} / Var(z)_{ii}
%
% The variance solve is NOT taken on trust. B19 was a wrong unconditional variance that produced a
% grid five times too wide and printed nothing suspicious, so the residual of the defining equation
% is computed and reported here before anything downstream uses it.
cnames={'full','diag','three'};
for c_c=1:3
    cn=cnames{c_c};
    Rho=calibVAR.(cn).Rho; SigmaSq=calibVAR.(cn).SigmaSq; Mew=calibVAR.(cn).Mew;
    M=size(Rho,1);
    if max(abs(eig(Rho)))>=1
        error('P5 calibration %s is not a stationary VAR',cn)
    end
    [~,posdef]=chol(SigmaSq);
    if posdef
        error('P5 calibration %s has a SigmaSq that is not positive definite',cn)
    end
    zmean=(eye(M)-Rho)\Mew;
    SigmaSqz=reshape((eye(M^2)-kron(Rho,Rho))\SigmaSq(:),M,M);
    resid=max(abs(Rho*SigmaSqz*Rho'+SigmaSq-SigmaSqz),[],'all');
    acvec=diag(Rho*SigmaSqz)./diag(SigmaSqz);
    calibVAR.(cn).M=M;
    calibVAR.(cn).zmean=zmean;
    calibVAR.(cn).SigmaSqz=SigmaSqz;
    calibVAR.(cn).sigmaz=sqrt(diag(SigmaSqz));
    calibVAR.(cn).autocorr=acvec;
    fprintf('\nP5 calibration %s (M=%i): the Lyapunov solution satisfies its own equation, residual %2.3e \n',cn,M,resid)
    fprintf('P5 calibration %s: sigmaz = [',cn); fprintf('%2.4f ',calibVAR.(cn).sigmaz); fprintf('], autocorr = [')
    fprintf('%2.4f ',acvec); fprintf('] \n')
end

%% The cross-correlation is the moment that distinguishes a VAR method from a pair of AR(1)s
S=calibVAR.full.SigmaSqz;
calibVAR.full.crosscorr=S(1,2)/sqrt(S(1,1)*S(2,2));
fprintf('P5 calibration full: the cross-correlation of the two variables is %2.4f \n',calibVAR.full.crosscorr)
if abs(calibVAR.full.crosscorr)<0.1
    warning('P5 full calibration has almost no cross-correlation, so the comparison against independent AR(1)s will be close to vacuous')
end

%% Grid sizes. The state count is znum^M and the transition matrix is znum^(2M), so the ceiling
%% here bites much harder than in P4. At M=2 it allows znum up to 31 (961 states); at M=3 it
%% allows znum up to 9 (729 states).
znums=[5,9,15,31];
calibVAR.Ncap=1000;   % skip any (znum,M) with znum^M above this, and log the skip
calibVAR.nreps=3;
calibVAR.timethreshold=0.5;
calibVAR.entropytol=10^(-7);
