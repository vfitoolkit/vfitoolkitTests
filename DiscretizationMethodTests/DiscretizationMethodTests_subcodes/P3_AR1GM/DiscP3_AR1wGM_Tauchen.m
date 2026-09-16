function output=DiscP3_AR1wGM_Tauchen(calib,znums,figure_c)
% P3: discretizeAR1wGM_Tauchen - TEST-FIRST
%
% This command does not exist yet. It is to be written as part of taking this bank green, and the
% checks below are its specification: they were authored before a line of it, so its signature and
% conventions are pinned here rather than being whatever the implementation happens to do.
%
% WHAT THE IMPLEMENTATION MUST SATISFY
%
%   signature   [z_grid,pi_z] = discretizeAR1wGM_Tauchen(mew,rho,mixprobs_i,mu_i,sigma_i,znum,Tauchen_q,tauchenoptions)
%
%   process     z' = mew + rho*z + e,  e ~ sum_i mixprobs_i*N(mu_i,sigma_i^2)
%               mew is the INTERCEPT, so E(z)=(mew+E(e))/(1-rho). Every command in the AR(1)
%               family now reads mew this way; discretizeAR1wGM_FarmerToda was the last to move.
%               Do not reintroduce the (1-rho)*mew reading here.
%
%   grid        evenly spaced, centred on E(z), half-width Tauchen_q*sd(z), where
%               sd(z)=sqrt(Var(e)/(1-rho^2)) and Var(e) is the MIXTURE's variance,
%               sum_i p_i*(mu_i^2+sigma_i^2) - (sum_i p_i*mu_i)^2.
%
%   pi_z        cdf mass in each bin under the mixture's conditional density, i.e. the weighted sum
%               over components of the normal cdf differences, with the two outer bins running to
%               -Inf and +Inf so that rows sum to one exactly:
%                   pi(i,j) = sum_k p_k*[Phi((z_j+w_j - condMean - mu_k)/sigma_k)
%                                        - Phi((z_j-w_j^- - condMean - mu_k)/sigma_k)]
%               with condMean = mew + rho*z_i and w the half-spacings.
%
%   options     tauchenoptions.parallel, as elsewhere in the family.
%
% Until it exists this subcode ERRORS, and takes the rest of the run with it. That is deliberate:
% a missing implementation should stop everything loudly rather than print a note and let the
% diary go on looking healthy. There is no graceful-degradation path here on purpose.

fprintf('\n========== P3: discretizeAR1wGM_Tauchen ========== \n')

output=struct();
mew=calib.mew; rho=calib.rho;
p=calib.mixprobs; mu=calib.mu; sd=calib.sigma;
zstar=(mew+calib.e.mean)/(1-rho);
sdz=sqrt(calib.z.var);
nz=length(znums);
Tauchen_q=3;
err_var=zeros(1,nz); err_skew=zeros(1,nz); err_exkurt=zeros(1,nz); runtime=zeros(1,nz);
err_condmean=zeros(1,nz); err_condvar=zeros(1,nz); err_condskew=zeros(1,nz); err_condkurt=zeros(1,nz); err_condsd_w=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    tauchenoptions=struct();

    % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
    tic;
    [z_grid,pi_z]=discretizeAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,Tauchen_q,tauchenoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid,pi_z]=discretizeAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,Tauchen_q,tauchenoptions);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);

    % --- invariants
    fprintf('znum=%i: size of z_grid [T0], this should be zero: %i \n',znum,any(size(z_grid)~=[znum,1]))
    fprintf('znum=%i: size of pi_z [T0], this should be zero: %i \n',znum,any(size(pi_z)~=[znum,znum]))
    fprintf('znum=%i: z_grid is strictly ascending [T0], this should be one: %i \n',znum,issorted(z_grid,'strictascend'))
    fprintf('znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z(:)<0)+any(pi_z(:)>1))
    % Exact, not approximate: the outer bins run to +-Inf, so the rows are complete by construction
    fprintf('znum=%i: rows of pi_z sum to one [T1], this should be zero: %2.8e \n',znum,max(abs(sum(pi_z,2)-1)))
    fprintf('znum=%i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
    % The grid specification, pinned
    fprintf('znum=%i: grid is centred on E(z)=(mew+E(e))/(1-rho) [T1], this should be zero: %2.8e \n',znum,abs(mean([z_grid(1),z_grid(end)])-zstar))
    fprintf('znum=%i: grid half-width is Tauchen_q*sd(z) [T1], this should be zero: %2.8e \n',znum,abs((z_grid(end)-z_grid(1))/2-Tauchen_q*sdz))

    % --- accuracy
    [mcmean,mcvar,~,statdist]=MarkovChainMoments(z_grid,pi_z);
    sk=sum(statdist.*(z_grid-mcmean).^3)/mcvar^1.5;
    ek=sum(statdist.*(z_grid-mcmean).^4)/mcvar^2-3;
    err_var(c_c)=abs(mcvar-calib.z.var);
    err_skew(c_c)=abs(sk-calib.z.skew);
    err_exkurt(c_c)=abs(ek-calib.z.exkurt);
    fprintf('znum=%i: variance error %2.3e, skewness error %2.3e, excess kurtosis error %2.3e, runtime %2.6f s \n',znum,err_var(c_c),err_skew(c_c),err_exkurt(c_c),runtime(c_c))

    % --- CONDITIONAL moments, to four orders. This is the most direct test there is of what this
    % command computes: conditional on z, the next period value is the mixture with every component
    % mean shifted by mew+rho*z, so the conditional moments ARE the mixture's moments. Every other
    % check in this subcode looks at the stationary distribution, which is two steps removed.
    % Truth, for every row i: the conditional distribution is the mixture shifted to mew+rho*z_i,
    % so its central moments are the mixture's own central moments, independent of i.
    condmean_true=mew+rho*z_grid;
    condmean_hat=pi_z*z_grid;
    condc2=zeros(znum,1); condc3=zeros(znum,1); condc4=zeros(znum,1);
    for z_c=1:znum
        dz=z_grid-condmean_hat(z_c);
        condc2(z_c)=sum(pi_z(z_c,:)'.*dz.^2);
        condc3(z_c)=sum(pi_z(z_c,:)'.*dz.^3);
        condc4(z_c)=sum(pi_z(z_c,:)'.*dz.^4);
    end
    ec2=calib.e.var; ec3=calib.e.k3; ec4=calib.e.k4+3*ec2^2; % central moments of the mixture
    [~,~,~,sd_c]=MarkovChainMoments(z_grid,pi_z);
    err_condmean(c_c)=max(abs(condmean_hat-condmean_true));
    err_condvar(c_c)=max(abs(condc2-ec2));
    err_condskew(c_c)=max(abs(condc3./condc2.^1.5-ec3/ec2^1.5));
    err_condkurt(c_c)=max(abs(condc4./condc2.^2-ec4/ec2^2));
    err_condsd_w(c_c)=abs(sum(sd_c.*sqrt(condc2))-sqrt(ec2));
    fprintf('znum=%i: conditional mean, worst row [T2] %2.3e; conditional variance, worst row [T2] %2.3e \n',znum,err_condmean(c_c),err_condvar(c_c))
    fprintf('znum=%i: conditional skewness, worst row [T2] %2.3e (truth %2.4f); conditional excess kurtosis, worst row [T2] %2.3e (truth %2.4f) \n',znum,err_condskew(c_c),ec3/ec2^1.5,err_condkurt(c_c),ec4/ec2^2-3)
    fprintf('znum=%i: weighted-average conditional sd error [T2] %2.3e \n',znum,err_condsd_w(c_c))

    output.sweep(c_c).znum=znum;
    output.sweep(c_c).z_grid=z_grid;
    output.sweep(c_c).pi_z=pi_z;
end
output.err_condmean=err_condmean; output.err_condvar=err_condvar;
output.err_condskew=err_condskew; output.err_condkurt=err_condkurt; output.err_condsd_w=err_condsd_w;
output.err_var=err_var; output.err_skew=err_skew; output.err_exkurt=err_exkurt;
output.runtime=runtime; output.znums=znums;

%% Reduction: a one-component mixture IS a normal, so this must reproduce discretizeAR1_Tauchen
znum=15;
[zgA,pzA]=discretizeAR1_Tauchen(mew,rho,sqrt(calib.e.var),znum,Tauchen_q,struct());
[zgB,pzB]=discretizeAR1wGM_Tauchen(mew,rho,1,0,sqrt(calib.e.var),znum,Tauchen_q,struct());
fprintf('nmix=1 vs discretizeAR1_Tauchen: z_grid [T1], this should be zero: %2.8e \n',max(abs(zgA-zgB)))
fprintf('nmix=1 vs discretizeAR1_Tauchen: pi_z [T1], this should be zero: %2.8e \n',max(abs(pzA-pzB),[],'all'))

%% Symmetry, on a symmetric mixture
[zgS,pzS]=discretizeAR1wGM_Tauchen(0,rho,calib.sym.mixprobs,calib.sym.mu,calib.sym.sigma,znum,Tauchen_q,struct());
fprintf('symmetric mixture: z_grid symmetric about zero [T1], this should be zero: %2.8e \n',max(abs((zgS+flipud(zgS))/2)))
fprintf('symmetric mixture: pi_z is centrosymmetric [T1], this should be zero: %2.8e \n',max(abs(pzS-rot90(pzS,2)),[],'all'))


%% REGRESSION BARS ON THE ACCURACY NUMBERS ABOVE
% The moment errors in the sweep above are printed with a value but no verdict, so a defect can sit
% in the output while the run summary reports a clean pass - which is what happened with B30 on
% 2026-08-28. These are regression bars, not accuracy claims: a single threshold cannot serve both
% a coarse grid, where a large error is legitimate, and a fine one. So each family gets a loose bar
% on the worst of the sweep and a tight one at the finest grid. Both come from the 2026-08-28 run
% with about an order of magnitude of headroom.
fprintf('worst variance error over the sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_var(:)))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',0.1,err_var(end))
fprintf('worst skewness error over the sweep [T2], this should be below %g: %2.3e \n',1,max(err_skew(:)))
fprintf('skewness error at the finest grid [T2], this should be below %g: %2.3e \n',1,err_skew(end))
fprintf('worst excess kurtosis error over the sweep [T2], this should be below %g: %2.3e \n',10,max(err_exkurt(:)))
fprintf('excess kurtosis error at the finest grid [T2], this should be below %g: %2.3e \n',10,err_exkurt(end))

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(znums,max(err_skew,10^(-18)),'o-',znums,max(err_exkurt,10^(-18)),'s-',znums,max(err_var,10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('AR1wGM\_Tauchen: |error|'); xlabel('znum'); legend('skewness','excess kurtosis','variance','Location','best')
subplot(1,2,2)
plot(znums,runtime,'o-'); set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('znum')

end
