% Setup for P10: AR(p), with gaussian innovations and with gaussian-mixture innovations
%
%    z' = mew + Rho(1)*z + Rho(2)*zlag1 + ... + Rho(p)*zlag(p-1) + e
%
% with e~N(0,sigma^2) for discretizeARp_FarmerToda and e a gaussian mixture for
% discretizeARpwGM_FarmerToda. Both commands were written for this block; neither existed before.
%
% THE STATE IS p-DIMENSIONAL AND ONLY ONE DIMENSION IS STOCHASTIC. The chain lives on
% (z,zlag1,...,zlag(p-1)) with znum^p states, and from any state only znum of those are reachable -
% the new z is drawn and every lag shifts along one place. That sparsity is a structural property
% the commands must have, not an accident of the parameters, and it is asserted directly.
%
% TRUTH COMES FROM THE MA(inf) REPRESENTATION, which gives all four moments exactly. Write
% z_t = mewz + sum_j psi_j e_{t-j} with psi_0=1 and psi_j = sum_k Rho(k)*psi_{j-k}. Cumulants are
% additive over independent terms and scale by the power, so
%       kappa_n(z) = (sum_j psi_j^n) * kappa_n(e)
% for every n>=2. That gives the variance, and with a gaussian mixture it gives the skewness and
% excess kurtosis too - the same cumulant argument the life-cycle blocks use for their age profiles,
% applied along the MA lags instead of along age. Autocovariances come from the same weights,
% gamma_j = kappa_2(e) * sum_i psi_i*psi_{i+j}.
%
% The variance is computed a SECOND, independent way, by the Lyapunov equation on the companion
% form, which is what the commands themselves use. The two must agree; that they do is checked below
% rather than assumed, and it is the check that says the truth here is right before any command is
% measured against it.

calibARp=struct();

%% The calibrations
% p=2 and p=3, because a command that only worked at two lags would pass every p=2 check. The
% coefficients are chosen to be comfortably stationary and to have genuinely different dynamics:
% AR2 is smoothly persistent, AR3 has an oscillating component so the MA weights change sign.
calibARp.AR2.Rho=[0.70,0.10];
calibARp.AR3.Rho=[0.50,0.25,-0.20];
calibARp.AR2.name='AR(2)';
calibARp.AR3.name='AR(3)';

% mew is the INTERCEPT, not the mean, and it is deliberately NON-ZERO. This is the convention pin:
% discretizeAR1wGM_FarmerToda used to read mew as the unconditional mean and was changed to the
% intercept to match the other ten commands in the family; a port of Toda's code that did not make
% the same change would reintroduce exactly that, and it is silent whenever mew=0.
calibARp.mew=0.06;
calibARp.sigma=0.05; % gaussian innovation

% The mixture: two components, asymmetric, so skewness and excess kurtosis are both non-zero and a
% method that only matched two moments would be visibly worse than one that matched four.
% NOT [-0.02;0.03], which was the first choice and is WRONG for the purpose: 0.6*(-0.02)+0.4*0.03
% is exactly zero, so the mixture was mean zero, mewz and Ez coincided, and the check that they
% differ could not pass. The run of 2026-09-17 duly failed it ten times while this file printed
% "the mixture has mean +0.000000, which is NOT zero" - a label contradicting its own number.
calibARp.mixprobs_i=[0.6;0.4];
calibARp.mu_i=[-0.02;0.05];
calibARp.sigma_i=[0.03;0.07];
% The mixture must be checked for mean zero or not - it is NOT here, which is what makes the
% distinction between the grid centre mew/(1-sum(Rho)) and the true mean (mew+E(e))/(1-sum(Rho))
% observable. If E(e) were zero the two would coincide and the check would be vacuous.
calibARp.Ee=calibARp.mixprobs_i'*calibARp.mu_i;
% ASSERTED, not asserted in a comment. The previous version of this line said the mean was not zero
% in its text while printing zero, which is exactly the kind of claim that should be a check.
fprintf('\nP10: the mixture has mean %+2.6f \n',calibARp.Ee)
fprintf('P10: that mean is NOT zero, so the grid centre and the true mean of z differ [T0], this should be one: %i \n',abs(calibARp.Ee)>1e-10)

%% Analytic truth, both innovation types, both calibrations
cnamesARp={'AR2','AR3'};
for c_c=1:2
    cn=cnamesARp{c_c};
    Rho=calibARp.(cn).Rho; p=length(Rho);
    calibARp.(cn).p=p;

    % companion form, and the persistence measure the commands branch their defaults on
    F=zeros(p,p); F(1,:)=Rho;
    if p>1
        F(2:p,1:p-1)=eye(p-1);
    end
    maxabseig=max(abs(eig(F)));
    calibARp.(cn).F=F; calibARp.(cn).maxabseig=maxabseig;

    % MA(inf) weights. 4000 terms is far past where maxabseig^j is at machine precision for these
    % calibrations; the truncation is checked below rather than trusted.
    nMA=4000;
    psi=zeros(1,nMA); psi(1)=1;
    for j_c=2:nMA
        acc=0;
        for k_c=1:p
            if j_c-k_c>=1
                acc=acc+Rho(k_c)*psi(j_c-k_c);
            end
        end
        psi(j_c)=acc;
    end
    calibARp.(cn).psi=psi;
    calibARp.(cn).S2=sum(psi.^2); calibARp.(cn).S3=sum(psi.^3); calibARp.(cn).S4=sum(psi.^4);
    fprintf('P10 %s: largest |eig(companion)| %2.4f, and the MA weight at lag %i is %2.1e (truncation is harmless) \n',calibARp.(cn).name,maxabseig,nMA,abs(psi(end)))

    for t_c=1:2
        if t_c==1
            tn='gauss';
            % gaussian innovation: k2=sigma^2, k3=0, k4=0
            k1=0; k2=calibARp.sigma^2; k3=0; k4=0;
        else
            tn='gm';
            % gaussian mixture: central moments, then cumulants
            mp=calibARp.mixprobs_i; mu=calibARp.mu_i; sd=calibARp.sigma_i;
            k1=mp'*mu;
            m2=mp'*((mu-k1).^2+sd.^2);
            m3=mp'*((mu-k1).^3+3*(mu-k1).*sd.^2);
            m4=mp'*((mu-k1).^4+6*((mu-k1).^2).*sd.^2+3*sd.^4);
            k2=m2; k3=m3; k4=m4-3*m2^2;
        end
        % z's cumulants, by the MA sum
        K2=calibARp.(cn).S2*k2;
        K3=calibARp.(cn).S3*k3;
        K4=calibARp.(cn).S4*k4;
        Ez=(calibARp.mew+k1)/(1-sum(Rho));
        mewz=calibARp.mew/(1-sum(Rho)); % the grid centre the commands use

        % the same variance a second way, by Lyapunov on the companion form
        Q=zeros(p,p); Q(1,1)=k2;
        Sig=reshape((eye(p^2)-kron(F,F))\Q(:),[p,p]);

        % autocovariances and autocorrelations, lags 1..p
        gam=zeros(1,p+1);
        for j_c=0:p
            gam(j_c+1)=k2*sum(psi(1:nMA-j_c).*psi(1+j_c:nMA));
        end

        calibARp.(cn).(tn).Ez=Ez;
        calibARp.(cn).(tn).mewz=mewz;
        calibARp.(cn).(tn).var=K2;
        calibARp.(cn).(tn).sd=sqrt(K2);
        calibARp.(cn).(tn).skew=K3/K2^1.5;
        calibARp.(cn).(tn).exkurt=K4/K2^2;
        calibARp.(cn).(tn).gamma=gam;
        calibARp.(cn).(tn).autocorr=gam/gam(1);
        calibARp.(cn).(tn).k2=k2;

        % THE TRUTH CHECKS ITSELF. Two independent routes to Var(z) - the MA sum and the Lyapunov
        % solve - and if they disagree the truth is wrong and every accuracy number below it is
        % meaningless. This is the one place in the block where a failure means the SETUP is broken
        % rather than a command.
        fprintf('P10 %s %s: Var(z) by MA sum %2.8f against Lyapunov %2.8f [T1], this should be zero: %2.3e \n',calibARp.(cn).name,tn,K2,Sig(1,1),abs(K2-Sig(1,1)))
        fprintf('P10 %s %s: E(z)=%+2.6f, grid centre=%+2.6f, sd=%2.6f, skewness=%+2.4f, excess kurtosis=%+2.4f \n',calibARp.(cn).name,tn,Ez,mewz,sqrt(K2),K3/K2^1.5,K4/K2^2)
    end
end

%% Grid sizes
% pi_z is znum^p square, so cost grows like znum^(2p). At p=2 the list below tops out at 441 states
% and at p=3 at 1331, both of which are real solves and quick ones. P2's top entry of 101 would be
% a 10201-by-10201 matrix at p=2 and a 1030301-by-1030301 one at p=3.
calibARp.znums2=[5,7,9,13,17,21]; % used at p=2
calibARp.znums3=[5,7,9,11];       % used at p=3
calibARp.timethreshold=5; calibARp.nreps=3;
