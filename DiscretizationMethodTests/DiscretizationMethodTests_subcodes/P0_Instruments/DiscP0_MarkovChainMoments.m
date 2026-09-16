function output=DiscP0_MarkovChainMoments(calibP0,figure_c)
% P0: closed-form validation of MarkovChainMoments()
%
% Everything else in this bank measures a discretization by running it through
% MarkovChainMoments(). This checks the instrument itself against chains whose moments are known
% in closed form, so a failure here is unambiguously the instrument and not the discretization.
%
% Two chains:
%   (a) a hand-built two-state chain: stationary distribution, mean, variance and autocorrelation
%       all have closed forms
%   (b) a Rouwenhorst chain: matches the AR(1)'s mean, variance and autocorrelation EXACTLY for
%       any znum, by construction (Kopecky & Suen 2010), so it is a closed-form target too

fprintf('\n========== P0: MarkovChainMoments against closed form ========== \n')

output=struct();

%% (a) Two-state chain, everything in closed form
p=calibP0.two.p; q=calibP0.two.q;
z_grid=[calibP0.two.z1; calibP0.two.z2];
pi_z=[1-p, p; q, 1-q];

% Closed form
statdist_true=[q; p]/(p+q);
mean_true=statdist_true(1)*z_grid(1)+statdist_true(2)*z_grid(2);
var_true=statdist_true(1)*(z_grid(1)-mean_true)^2+statdist_true(2)*(z_grid(2)-mean_true)^2;
autocorr_true=1-p-q; % standard result for a two-state chain

[mcmean,mcvar,mcautocorr,mcstatdist]=MarkovChainMoments(z_grid,pi_z);

fprintf('two-state chain: stationary dist [T1], this should be zero: %2.8e \n',max(abs(mcstatdist(:)-statdist_true(:))))
fprintf('two-state chain: mean [T1], this should be zero: %2.8e \n',abs(mcmean-mean_true))
fprintf('two-state chain: variance [T1], this should be zero: %2.8e \n',abs(mcvar-var_true))
fprintf('two-state chain: autocorrelation [T1], this should be zero: %2.8e \n',abs(mcautocorr-autocorr_true))
output.two.mean=[mcmean,mean_true];
output.two.variance=[mcvar,var_true];
output.two.autocorrelation=[mcautocorr,autocorr_true];
output.two.statdist=[mcstatdist(:),statdist_true(:)];

% The autocorrelation row is the one that matters most, and it needs a probe to show that it
% discriminates something the other three do not.
%
% Scaling BOTH off-diagonals of a two-state chain by k leaves the stationary distribution exactly
% unchanged - it is proportional to [k*q, k*p], i.e. to [q, p] - while changing the persistence
% from 1-(p+q) to 1-k*(p+q). So this alternative chain has the SAME stationary distribution, mean
% and variance as the one above, and a different autocorrelation. If the autocorrelation check
% could not tell them apart, it would not be testing the dynamics at all.
%
% (Transposing pi_z would not work as a probe: the transpose of a transition matrix is not itself
% a transition matrix unless the chain is doubly stochastic, and MarkovChainMoments correctly
% rejects it at its own assert that the rows sum to one.)
kscale=0.5;
pi_z_alt=[1-kscale*p, kscale*p; kscale*q, 1-kscale*q];
autocorr_alt_true=1-kscale*(p+q);
[mcmean_alt,mcvar_alt,mcautocorr_alt,mcstatdist_alt]=MarkovChainMoments(z_grid,pi_z_alt);
fprintf('two-state chain, rescaled: stationary dist is UNCHANGED [T1], this should be zero: %2.8e \n',max(abs(mcstatdist_alt(:)-mcstatdist(:))))
fprintf('two-state chain, rescaled: mean is UNCHANGED [T1], this should be zero: %2.8e \n',abs(mcmean_alt-mcmean))
fprintf('two-state chain, rescaled: variance is UNCHANGED [T1], this should be zero: %2.8e \n',abs(mcvar_alt-mcvar))
fprintf('two-state chain, rescaled: autocorrelation matches its own closed form [T1], this should be zero: %2.8e \n',abs(mcautocorr_alt-autocorr_alt_true))
fprintf('two-state chain, rescaled: autocorrelation must DIFFER from the original, this should NOT be zero: %2.8e \n',abs(mcautocorr_alt-mcautocorr))
fprintf('   (if that last line prints zero, the autocorrelation check is not testing the dynamics) \n')
output.two.rescaled=[mcmean_alt,mcvar_alt,mcautocorr_alt,autocorr_alt_true];
output.two.statdist_rescaled=[mcstatdist(:),mcstatdist_alt(:)];

%% (b) Rouwenhorst chain: mean, variance and autocorrelation are exact by construction
mew=calibP0.rouw.mew; rho=calibP0.rouw.rho; sigma=calibP0.rouw.sigma; znum=calibP0.rouw.znum;
[rz_grid,rpi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum);

rmean_true=mew/(1-rho);              % note: mew is the intercept, so E(z)=mew/(1-rho)
rvar_true=sigma^2/(1-rho^2);
rautocorr_true=rho;

[rmcmean,rmcvar,rmcautocorr,rmcstatdist]=MarkovChainMoments(rz_grid,rpi_z);

fprintf('Rouwenhorst chain: mean [T1], this should be zero: %2.8e \n',abs(rmcmean-rmean_true))
fprintf('Rouwenhorst chain: variance [T1], this should be zero: %2.8e \n',abs(rmcvar-rvar_true))
fprintf('Rouwenhorst chain: autocorrelation [T1], this should be zero: %2.8e \n',abs(rmcautocorr-rautocorr_true))
output.rouw.mean=[rmcmean,rmean_true];
output.rouw.variance=[rmcvar,rvar_true];
output.rouw.autocorrelation=[rmcautocorr,rautocorr_true];

% The Rouwenhorst stationary distribution is binomial, which is another closed form
sbinom=zeros(znum,1);
for jj=1:znum
    sbinom(jj)=nchoosek(znum-1,jj-1);
end
sbinom=sbinom/2^(znum-1);
fprintf('Rouwenhorst chain: stationary dist is binomial [T1], this should be zero: %2.8e \n',max(abs(rmcstatdist(:)-sbinom(:))))
output.rouw.statdist=[rmcstatdist(:),sbinom(:)];

%% Option sweep on MarkovChainMoments
% eigenvector=1 (a direct eigen-solve) and eigenvector=0 (iterate the distribution to a fixed
% point) are two different solves for the same object. They agree only to the ITERATIVE one's
% tolerance, NOT to machine precision, so these are not T1 checks: mcmomentsoptions.Tolerance
% defaults to 1e-8 and that is the order the differences below should come out at. Asserting
% anything tighter would be asserting something the algorithm never promised.
mcopts1=struct(); mcopts1.eigenvector=1;
mcopts0=struct(); mcopts0.eigenvector=0; mcopts0.Tolerance=10^(-8);
[m1,v1,a1,s1]=MarkovChainMoments(rz_grid,rpi_z,mcopts1);
[m0,v0,a0,s0]=MarkovChainMoments(rz_grid,rpi_z,mcopts0);
tol0=10^(-8);
fprintf('eigenvector=1 vs 0 at Tolerance=%g: stationary dist [T2, expect order Tolerance], this should be below %g: %2.8e \n',tol0,100*tol0,max(abs(s1(:)-s0(:))))
fprintf('eigenvector=1 vs 0 at Tolerance=%g: mean [T2], this should be below %g: %2.8e \n',tol0,100*tol0,abs(m1-m0))
fprintf('eigenvector=1 vs 0 at Tolerance=%g: variance [T2], this should be below %g: %2.8e \n',tol0,100*tol0,abs(v1-v0))
fprintf('eigenvector=1 vs 0 at Tolerance=%g: autocorrelation [T2], this should be below %g: %2.8e \n',tol0,100*tol0,abs(a1-a0))
output.opt.eigenvector=[m1,v1,a1;m0,v0,a0];

% calcautocorrelation=0 skips the autocorrelation but must not disturb the other three
mcoptsNoAC=struct(); mcoptsNoAC.calcautocorrelation=0;
[m2,v2,~,s2]=MarkovChainMoments(rz_grid,rpi_z,mcoptsNoAC);
fprintf('calcautocorrelation=0: mean unchanged [T0], this should be zero: %2.8e \n',abs(m2-m1))
fprintf('calcautocorrelation=0: variance unchanged [T0], this should be zero: %2.8e \n',abs(v2-v1))
fprintf('calcautocorrelation=0: stationary dist unchanged [T0], this should be zero: %2.8e \n',max(abs(s2(:)-s1(:))))

% Tolerance controls the iterative path, so tightening it must move the answer TOWARDS the direct
% eigen-solve. That is a sharper statement than "the two paths roughly agree": it says the
% iterative path is converging to the right thing, and that Tolerance is what is binding.
tol3=10^(-12);
mcoptsTol=struct(); mcoptsTol.eigenvector=0; mcoptsTol.Tolerance=tol3;
[m3,v3,a3,s3]=MarkovChainMoments(rz_grid,rpi_z,mcoptsTol);
fprintf('Tolerance %g vs %g (eigenvector=0) changes the variance, this should NOT be zero: %2.8e \n',tol3,tol0,abs(v3-v0))
fprintf('   (if it is zero, Tolerance is not binding and the checks below prove nothing) \n')
fprintf('eigenvector=0 at Tolerance=%g vs the eigen-solve: stationary dist [T2, expect order Tolerance], this should be below %g: %2.8e \n',tol3,100*tol3,max(abs(s3(:)-s1(:))))
fprintf('eigenvector=0 at Tolerance=%g vs the eigen-solve: variance [T2], this should be below %g: %2.8e \n',tol3,100*tol3,abs(v3-v1))
fprintf('tightening Tolerance moves the iterative solve CLOSER to the eigen-solve [T2], this should be one: %i \n',abs(v3-v1)<abs(v0-v1))
output.opt.tolerance=[m0,v0,a0;m3,v3,a3];

% parallel=0 and parallel=1 are both cpu paths and must agree exactly
mcoptsP0=struct(); mcoptsP0.parallel=0;
mcoptsP1=struct(); mcoptsP1.parallel=1;
[m4,v4,a4,s4]=MarkovChainMoments(rz_grid,rpi_z,mcoptsP0);
[m5,v5,a5,s5]=MarkovChainMoments(rz_grid,rpi_z,mcoptsP1);
fprintf('parallel=0 vs 1: mean [T0], this should be zero: %2.8e \n',abs(m4-m5))
fprintf('parallel=0 vs 1: variance [T0], this should be zero: %2.8e \n',abs(v4-v5))
fprintf('parallel=0 vs 1: autocorrelation [T0], this should be zero: %2.8e \n',abs(a4-a5))
fprintf('parallel=0 vs 1: stationary dist [T0], this should be zero: %2.8e \n',max(abs(s4(:)-s5(:))))

%% Figure
figure(figure_c)
subplot(2,2,1); bar([mcstatdist(:),statdist_true(:)])
title('two-state: stationary dist'); legend('MarkovChainMoments','closed form','Location','best'); xlabel('z index')
subplot(2,2,2); bar([rmcstatdist(:),sbinom(:)])
title('Rouwenhorst: stationary dist'); legend('MarkovChainMoments','binomial','Location','best'); xlabel('z index')
subplot(2,2,3)
bar([abs(mcmean-mean_true),abs(mcvar-var_true),abs(mcautocorr-autocorr_true)])
set(gca,'YScale','log'); set(gca,'XTickLabel',{'mean','var','autocorr'})
title('two-state: |error| vs closed form')
subplot(2,2,4)
bar([abs(rmcmean-rmean_true),abs(rmcvar-rvar_true),abs(rmcautocorr-rautocorr_true)])
set(gca,'YScale','log'); set(gca,'XTickLabel',{'mean','var','autocorr'})
title('Rouwenhorst: |error| vs exact')

end
