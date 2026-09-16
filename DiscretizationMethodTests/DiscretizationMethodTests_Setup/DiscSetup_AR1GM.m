% Setup for P3: AR(1) with gaussian-mixture innovations
%
%    z' = mew + rho*z + e,   e ~ sum_i mixprobs_i * N(mu_i, sigma_i^2)
%
% (mew is the intercept, as everywhere else in the family - discretizeAR1wGM_FarmerToda was the
% last command to move onto that convention.)
%
% The point of a gaussian mixture is the moments beyond the second, so the calibration has to have
% some. A near-symmetric mixture would make every ordering in this block vacuous, because
% nMoments=2 and nMoments=4 would agree and the gaussian methods would not be misspecified in any
% way that shows.

calibGM=struct();
calibGM.mew=0;
calibGM.rho=0.9;
calibGM.mixprobs=[0.8;0.2];
calibGM.mu=[0.05;-0.2];
calibGM.sigma=[0.1;0.4];

%% Analytic moments of the innovation, from the mixture parameters
p=calibGM.mixprobs; mu=calibGM.mu; sd=calibGM.sigma;
m1=sum(p.*mu);                                   % mean
m2=sum(p.*(mu.^2+sd.^2));                        % uncentred 2nd
m3=sum(p.*(mu.^3+3*mu.*sd.^2));                  % uncentred 3rd
m4=sum(p.*(mu.^4+6*(mu.^2).*sd.^2+3*sd.^4));     % uncentred 4th
c2=m2-m1^2;
c3=m3-3*m1*m2+2*m1^3;
c4=m4-4*m1*m3+6*(m1^2)*m2-3*m1^4;
calibGM.e.mean=m1;
calibGM.e.var=c2;
calibGM.e.k3=c3;              % third cumulant = third central moment
calibGM.e.k4=c4-3*c2^2;       % fourth cumulant = fourth central moment - 3*variance^2
calibGM.e.skew=c3/c2^1.5;
calibGM.e.exkurt=calibGM.e.k4/c2^2;

%% Analytic moments of z itself
% z is a linear process, z = sum_k rho^k e_{t-k}, and cumulants are additive over independent
% terms and homogeneous of their own degree. So the n-th cumulant of z is the n-th cumulant of e
% divided by (1-rho^n). That gives exact skewness and excess kurtosis for the stationary
% distribution, which is what this block measures the methods against.
rho=calibGM.rho;
calibGM.z.mean=(calibGM.mew+m1)/(1-rho);
calibGM.z.k2=calibGM.e.var/(1-rho^2);
calibGM.z.k3=calibGM.e.k3/(1-rho^3);
calibGM.z.k4=calibGM.e.k4/(1-rho^4);
calibGM.z.var=calibGM.z.k2;
calibGM.z.skew=calibGM.z.k3/calibGM.z.k2^1.5;
calibGM.z.exkurt=calibGM.z.k4/calibGM.z.k2^2;

%% The anti-vacuity guard, printed so a future change to the calibration cannot quietly break it
fprintf('\nP3 calibration: innovation skewness %2.4f, excess kurtosis %2.4f \n',calibGM.e.skew,calibGM.e.exkurt)
fprintf('                z skewness %2.4f, excess kurtosis %2.4f \n',calibGM.z.skew,calibGM.z.exkurt)
if abs(calibGM.e.skew)<0.5 || abs(calibGM.e.exkurt)<1
    warning('P3 calibration is too close to gaussian; the orderings in this block will be vacuous')
end

%% A symmetric mixture, for the checks that need one
% The centro-symmetry check does not apply to a skewed process - the whole point of the main
% calibration is that it is not symmetric - so it would simply be lost for these two commands
% unless a symmetric-mixture config is carried alongside. Two components, equal weight, mirrored
% means: fat-tailed but symmetric.
% NOTE: mirrored means are NOT enough. [0.5;0.5] on N(0.15,0.1^2) and N(-0.15,0.3^2) has mirrored
% means but different spreads, so the mixture is skewed and the check it is meant to support fails
% by 0.44 rather than by a rounding error. Use a SCALE mixture instead: both components centred at
% zero, different variances, so it is symmetric by construction and fat-tailed, which also keeps it
% a meaningful test of a mixture method rather than a normal in disguise.
calibGM.sym.mixprobs=[0.8;0.2];
calibGM.sym.mu=[0;0];
calibGM.sym.sigma=[0.1;0.4];

%% Gaussian comparators, fitted to the mixture's first two moments
% This is the Kirkby (2025) comparison: what does a practitioner lose by discretizing a
% gaussian-mixture process with a gaussian method? The comparator is not "some gaussian AR(1)" but
% the gaussian AR(1) with the SAME mean and variance, which is what someone would actually fit.
calibGM.gaussian.mew=calibGM.mew+m1;   % so that E(z) matches
calibGM.gaussian.rho=rho;
calibGM.gaussian.sigma=sqrt(c2);       % innovation std dev matched

znums=[5,9,15,31,51,101];
calibGM.nreps=5;
calibGM.timethreshold=0.5;
calibGM.entropytol=10^(-7);
calibGM.entropytol4=10^(-4);  % nMoments=4 solves land further out; measured in P2
