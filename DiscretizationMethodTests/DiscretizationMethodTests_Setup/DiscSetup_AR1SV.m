% Setup for P4: AR(1) with log-AR(1) stochastic volatility
%
%    z_t = rho*z_{t-1} + u_t,               u_t ~ N(0,exp(x_t))
%    x_t = (1-phi)*xBar + phi*x_{t-1} + e,  e ~ N(0,sigmae^2)
%
% xBar is set inside the commands so that the unconditional standard deviation of u is sigmau.

calibSV=struct();
calibSV.rho=0.95;     % persistence of z
calibSV.phi=0.9;      % persistence of the log-volatility process
calibSV.sigmau=0.2;   % unconditional standard deviation of u
calibSV.sigmae=0.3;   % standard deviation of the innovations to log volatility

rho=calibSV.rho; phi=calibSV.phi; sigmau=calibSV.sigmau; sigmae=calibSV.sigmae;

%% The volatility process, which is an ordinary gaussian AR(1)
calibSV.sigmaX=(sigmae^2)/(1-phi^2);           % unconditional variance of x
calibSV.xBar=2*log(sigmau)-calibSV.sigmaX/2;   % unconditional mean of x
calibSV.x.mean=calibSV.xBar;
calibSV.x.var=calibSV.sigmaX;
calibSV.x.autocorr=phi;

%% z: mean, variance and autocorrelation are all closed form
calibSV.z.mean=0;
calibSV.z.var=sigmau^2/(1-rho^2);   % because E[exp(x)]=exp(xBar+sigmaX/2)=sigmau^2
calibSV.z.autocorr=rho;

%% z: the EXCESS KURTOSIS, which is the moment stochastic volatility exists to produce
% Conditional on the whole path of x, z_t is gaussian with variance V=sum_k rho^(2k)*exp(x_{t-k}),
% so E[z^4]=3*E[V^2] and the excess kurtosis is 3*E[V^2]/E[V]^2-3. Writing that out,
%    E[V]   = exp(xBar+sigmaX/2)/(1-rho^2)
%    E[V^2] = exp(2*xBar+sigmaX) * sum_{j,k>=0} rho^(2(j+k)) * exp(phi^|j-k| * sigmaX)
% and the exponential prefactors cancel in the ratio, leaving
%    excess kurtosis = 3*S*(1-rho^2)^2 - 3,   S = sum_{j,k>=0} rho^(2(j+k))*exp(phi^|j-k|*sigmaX)
% The double sum converges geometrically in rho^2, so truncating is exact to machine precision.
%
% Two limits check the formula, and both are asserted below rather than taken on trust:
%    sigmaX -> 0 (no stochastic volatility)  =>  S = 1/(1-rho^2)^2, so excess kurtosis = 0
%    rho = 0 (z is iid, a scale mixture)     =>  S = exp(sigmaX),   so excess kurtosis = 3*exp(sigmaX)-3
Kmax=2000; % rho^(2*2000) is far below machine precision for any rho used here
kk=(0:Kmax)';
S=0;
for j_c=0:Kmax
    S=S+sum(rho.^(2*(j_c+kk)).*exp(phi.^abs(j_c-kk)*calibSV.sigmaX));
end
calibSV.z.exkurt=3*S*(1-rho^2)^2-3;

% limit check 1: no stochastic volatility
S0=0;
for j_c=0:Kmax
    S0=S0+sum(rho.^(2*(j_c+kk)));
end
fprintf('\nP4 kurtosis formula, limit check (sigmaX=0 must give excess kurtosis zero): %2.3e \n',abs(3*S0*(1-rho^2)^2-3))
% limit check 2: rho=0 gives 3*exp(sigmaX)-3
% This must re-evaluate the SERIES at rho=0 and compare that against the closed form. Comparing
% the closed form against itself is identically zero and tests nothing.
S1=0;
for j_c=0:Kmax
    S1=S1+sum((0.^(2*(j_c+kk))).*exp(phi.^abs(j_c-kk)*calibSV.sigmaX));
end
fprintf('P4 kurtosis formula, limit check (rho=0 must give 3*exp(sigmaX)-3): %2.3e \n',abs((3*S1-3)-(3*exp(calibSV.sigmaX)-3)))
fprintf('P4 calibration: z has variance %2.4f, autocorrelation %2.2f, excess kurtosis %2.4f \n',calibSV.z.var,calibSV.z.autocorr,calibSV.z.exkurt)
if calibSV.z.exkurt<0.5
    warning('P4 calibration produces little excess kurtosis; the point of stochastic volatility is the fat tails, so the block will be close to vacuous')
end

%% Grid sizes. This is a two-dimensional command, so per the nested-sweep rule the xnum and znum
%% sweeps run independently over the same range, subject to a ceiling on prod(znum) because the
%% joint transition matrix is prod(znum)-by-prod(znum) and the cost scales as its square.
znums=[5,9,15,31,51,101];
calibSV.Ncap=1000;   % skip any (xnum,znum) with xnum*znum above this, and log the skip
calibSV.nreps=5;
calibSV.timethreshold=0.5;
calibSV.entropytol=10^(-7);
