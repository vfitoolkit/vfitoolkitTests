% Setup for P7: life-cycle (age-dependent) AR(1) with gaussian-mixture innovations
%
%    z(j) = mew(j) + rho(j)*z(j-1) + e(j),   e(j) ~ sum_i mixprobs_i(i,j)*N(mu_i(i,j),sigma_i(i,j)^2)
%
% THE TRUTH HERE IS A CUMULANT RECURSION, and it is what makes this block stronger than P6. Since
% e(j) is independent of z(j-1), cumulants add over the sum and scale as k_n(a*X)=a^n*k_n(X):
%    k1(z_j) = mew(j) + rho(j)*k1(z_{j-1}) + k1(e_j)
%    kn(z_j) =          rho(j)^n*kn(z_{j-1}) + kn(e_j)      for n>=2
% The mixture's own cumulants are closed form, so mean, variance, SKEWNESS and EXCESS KURTOSIS all
% get exact age profiles. P6 had truth for three moments; the two extra ones here are the whole
% reason gaussian mixtures are used, and the reason this block exists.

calibLCGM=struct();
calibLCGM.J=41;
J=calibLCGM.J;
agevec=(1:J);
nmix=2;

%% The age-varying calibration
% Deliberately NOT mean zero. P3 could not distinguish the two grid-centring conventions in this
% family because its mixture was mean zero by construction; this one can. See item 1 of
% DiscretizationMethods_todo.md - discretizeLifeCycleAR1wGM_KFTT centres its grid on a recursion
% with no E(e) term, discretizeLifeCycleAR1wGM_Tauchen includes it, and they coincide only when
% E(e)=0. The offset is E(e)/(1-rho), which persistence makes large.
calibLCGM.vary.rho=0.95-0.05*(agevec-1)/(J-1);
calibLCGM.vary.mew=0.02*ones(1,J);
mixp=[0.75;0.25];
calibLCGM.vary.mixprobs_i=repmat(mixp,1,J);
% a small common component and a rarer, more negative, more volatile one: a downside-risk mixture
calibLCGM.vary.mu_i=[0.05*ones(1,J); -0.10-0.05*(agevec-1)/(J-1)];
calibLCGM.vary.sigma_i=[0.10+0.04*sin(pi*(agevec-1)/(J-1)); 0.30*ones(1,J)];

%% The frozen calibration: constant parameters AND a mean-zero mixture
% Mean zero so that the two centring conventions coincide, which is what makes the frozen identity
% stateable for BOTH pairs at once - KFTT against discretizeAR1wGM_FarmerToda, and Tauchen against
% discretizeAR1wGM_Tauchen.
fp=[0.7;0.3]; fmu=[0.09; -0.7*0.09/0.3]; fsig=[0.12;0.28];
calibLCGM.frozen.rho=0.90*ones(1,J);
calibLCGM.frozen.mew=zeros(1,J);
calibLCGM.frozen.mixprobs_i=repmat(fp,1,J);
calibLCGM.frozen.mu_i=repmat(fmu,1,J);
calibLCGM.frozen.sigma_i=repmat(fsig,1,J);
fprintf('\nP7 frozen mixture has mean %2.3e (it must be zero for the frozen identity to be stateable) \n',sum(fp.*fmu))

%% The cumulant recursion, for both calibrations
cnamesLCGM={'vary','frozen'};
for c_c=1:2
    cn=cnamesLCGM{c_c};
    mew=calibLCGM.(cn).mew; rho=calibLCGM.(cn).rho;
    p=calibLCGM.(cn).mixprobs_i; m=calibLCGM.(cn).mu_i; s=calibLCGM.(cn).sigma_i;
    % cumulants of the mixture innovation at each age, from its raw moments
    r1=sum(p.*m,1);
    r2=sum(p.*(m.^2+s.^2),1);
    r3=sum(p.*(m.^3+3*m.*s.^2),1);
    r4=sum(p.*(m.^4+6*(m.^2).*(s.^2)+3*s.^4),1);
    c2=r2-r1.^2;
    c3=r3-3*r1.*r2+2*r1.^3;
    c4=r4-4*r1.*r3+6*(r1.^2).*r2-3*r1.^4;
    ek1=r1; ek2=c2; ek3=c3; ek4=c4-3*c2.^2; % cumulants of e(j)
    % the recursion for z, from z(0)=0 (all cumulants of a point mass at zero are zero)
    k1=zeros(1,J); k2=zeros(1,J); k3=zeros(1,J); k4=zeros(1,J);
    k1(1)=mew(1)+ek1(1); k2(1)=ek2(1); k3(1)=ek3(1); k4(1)=ek4(1);
    for j_c=2:J
        k1(j_c)=mew(j_c)+rho(j_c)*k1(j_c-1)+ek1(j_c);
        k2(j_c)=rho(j_c)^2*k2(j_c-1)+ek2(j_c);
        k3(j_c)=rho(j_c)^3*k3(j_c-1)+ek3(j_c);
        k4(j_c)=rho(j_c)^4*k4(j_c-1)+ek4(j_c);
    end
    calibLCGM.(cn).mewz=k1;
    calibLCGM.(cn).varz=k2;
    calibLCGM.(cn).sigmaz=sqrt(k2);
    calibLCGM.(cn).skewz=k3./k2.^1.5;
    calibLCGM.(cn).exkurtz=k4./k2.^2;
    calibLCGM.(cn).emean=ek1;
    calibLCGM.(cn).evar=ek2;
    acz=nan(1,J);
    for j_c=2:J
        acz(j_c)=rho(j_c)*sqrt(k2(j_c-1))/sqrt(k2(j_c));
    end
    calibLCGM.(cn).autocorr=acz;
    fprintf('P7 calibration %s: E(e) runs %+2.4f to %+2.4f; z has sd %2.4f to %2.4f, skewness %+2.4f to %+2.4f, excess kurtosis %2.4f to %2.4f \n',cn,ek1(1),ek1(J),sqrt(k2(1)),sqrt(k2(J)),calibLCGM.(cn).skewz(1),calibLCGM.(cn).skewz(J),calibLCGM.(cn).exkurtz(1),calibLCGM.(cn).exkurtz(J))
end

% Two limit checks on the recursion, rather than trusting it.
% (i) a one-component mixture is a normal, so its 3rd and 4th cumulants must be zero. This tests
% the central-moment-to-cumulant algebra above, which is where an error would be silent.
r1t=0.3; s1t=0.2;
k3test=1*(r1t^3+3*r1t*s1t^2)-3*r1t*(r1t^2+s1t^2)+2*r1t^3;
k4test=1*(r1t^4+6*r1t^2*s1t^2+3*s1t^4)-4*r1t*(r1t^3+3*r1t*s1t^2)+6*r1t^2*(r1t^2+s1t^2)-3*r1t^4-3*(s1t^2)^2;
fprintf('P7 cumulant formulas, limit check (a single normal has zero 3rd and 4th cumulants): %2.3e %2.3e \n',abs(k3test),abs(k4test))
% (ii) the recursion against an INDEPENDENT closed form. Any check of the form "set rho=0 and the
% recursion gives the innovation cumulants" is vacuous - that is the same line of algebra compared
% against itself. What is not vacuous is the stationary fixed point: with rho and the mixture held
% constant, iterating kn(j)=rho^n*kn(j-1)+kn(e) converges to kn(e)/(1-rho^n), which is derived
% separately rather than read off the loop. The frozen calibration is run out far enough to get
% there and compared against it.
rf=calibLCGM.frozen.rho(1);
pf=calibLCGM.frozen.mixprobs_i(:,1); mf=calibLCGM.frozen.mu_i(:,1); sf=calibLCGM.frozen.sigma_i(:,1);
q1=sum(pf.*mf); q2=sum(pf.*(mf.^2+sf.^2)); q3=sum(pf.*(mf.^3+3*mf.*sf.^2)); q4=sum(pf.*(mf.^4+6*(mf.^2).*(sf.^2)+3*sf.^4));
d2=q2-q1^2; d3=q3-3*q1*q2+2*q1^3; d4=q4-4*q1*q3+6*q1^2*q2-3*q1^4;
ef=[q1,d2,d3,d4-3*d2^2]; % the four cumulants of the frozen innovation
kk=zeros(1,4);
for i_c=1:2000
    kk=[calibLCGM.frozen.mew(1)+rf*kk(1)+ef(1), rf^2*kk(2)+ef(2), rf^3*kk(3)+ef(3), rf^4*kk(4)+ef(4)];
end
kfix=[(calibLCGM.frozen.mew(1)+ef(1))/(1-rf), ef(2)/(1-rf^2), ef(3)/(1-rf^3), ef(4)/(1-rf^4)];
fprintf('P7 cumulant recursion, checked against its stationary fixed point kn(e)/(1-rho^n): %2.3e \n',max(abs(kk-kfix)))

if abs(calibLCGM.vary.skewz(J))<0.1
    warning('P7 age-varying calibration has almost no skewness; the whole point of the block is the higher moments')
end

%% Grid sizes
znums=[5,9,15,31,51];
% No even-znum list here: B15 was a parity bug in the Fella-Gallipoli-Pan family, which P6 covers.
% Neither of this block's commands has a parity-dependent construction.
calibLCGM.nreps=3;
calibLCGM.timethreshold=0.5;
calibLCGM.entropytol=10^(-7);
calibLCGM.Jbig=200;
