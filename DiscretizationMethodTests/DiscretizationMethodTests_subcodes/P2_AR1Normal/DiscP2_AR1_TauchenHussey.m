function output=DiscP2_AR1_TauchenHussey(calib,znums,figure_c)
% P2: discretizeAR1_TauchenHussey
%
% P2 owns this command's option sweep. It has one option, baseSigma, and Floden (2008) is
% organised entirely around the three choices for it - so the sweep here IS his Table 1's three
% Tauchen-Hussey columns:
%
%   baseSigma = sigma_eps    Tauchen & Hussey's own suggestion (the conditional std dev)
%   baseSigma = sigma_z      the unconditional std dev; Klein (2007) suggests this
%   baseSigma = w*sigma_eps + (1-w)*sigma_z with w = 1/2 + rho/4    Floden's third variant,
%                            which is the toolkit's default
%
% Floden's finding is that the second is badly behaved at high persistence - it puts the nodes so
% far from the mean that the implied autocorrelation is too high and the conditional variance too
% low - and that the third is the robust compromise. The persistent calibration is where that
% shows, and this subcode asserts the ordering rather than just reporting it.
%
% One more thing this subcode is the first to exercise: znum=101 on this command. Its local
% gausshermite routine used to return unconverged nodes there, silently, giving a grid that was
% neither monotone nor symmetric. That was found by the invariant block below.

fprintf('\n========== P2: discretizeAR1_TauchenHussey ========== \n')

output=struct();
nz=length(znums);

for cal_c=1:3
    cname=calib.names{cal_c};
    mew=calib.(cname).mew; rho=calib.(cname).rho; sigma=calib.(cname).sigma;
    zstar=mew/(1-rho); varz=sigma^2/(1-rho^2);
    fprintf('\n--- calibration %s: mew=%g, rho=%g, sigma=%g (so E(z)=%g, sd(z)=%g) --- \n',cname,mew,rho,sigma,zstar,sqrt(varz))

    err_mean=zeros(1,nz); err_var=zeros(1,nz); err_ac=zeros(1,nz);
    err_condmean=zeros(1,nz); err_condvar=zeros(1,nz); err_condsd_w=zeros(1,nz);
    err_skew=zeros(1,nz); err_exkurt=zeros(1,nz); runtime=zeros(1,nz); nrepsused=zeros(1,nz);
    for c_c=1:nz
        znum=znums(c_c);
        tauchenhusseyoptions=struct();

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
        tic;
        [z_grid,pi_z]=discretizeAR1_TauchenHussey(mew,rho,sigma,znum,tauchenhusseyoptions);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z]=discretizeAR1_TauchenHussey(mew,rho,sigma,znum,tauchenhusseyoptions);
            treps(r_c)=toc;
        end
        runtime(c_c)=median(treps);
        nrepsused(c_c)=nreps;

        % --- invariants
        fprintf('%s znum=%i: size of z_grid [T0], this should be zero: %i \n',cname,znum,any(size(z_grid)~=[znum,1]))
        fprintf('%s znum=%i: size of pi_z [T0], this should be zero: %i \n',cname,znum,any(size(pi_z)~=[znum,znum]))
        % The check that found the unconverged-node bug at znum=101
        fprintf('%s znum=%i: z_grid is strictly ascending [T0], this should be one: %i \n',cname,znum,issorted(z_grid,'strictascend'))
        fprintf('%s znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',cname,znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        % pi_z is explicitly row-normalised by this command, so the row sums are exact by construction
        fprintf('%s znum=%i: rows of pi_z sum to one [T1], this should be zero: %2.8e \n',cname,znum,max(abs(sum(pi_z,2)-1)))
        fprintf('%s znum=%i: no NaN or Inf [T0], this should be zero: %i \n',cname,znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
        fprintf('%s znum=%i: z_grid symmetric about mew/(1-rho) [T1], this should be zero: %2.8e \n',cname,znum,max(abs((z_grid+flipud(z_grid))/2-zstar)))
        fprintf('%s znum=%i: pi_z is centrosymmetric [T1], this should be zero: %2.8e \n',cname,znum,max(abs(pi_z-rot90(pi_z,2)),[],'all'))

        % --- accuracy
        [mcmean,mcvar,mcautocorr]=MarkovChainMoments(z_grid,pi_z);
        err_mean(c_c)=abs(mcmean-zstar);
        err_var(c_c)=abs(mcvar-varz);
        err_ac(c_c)=abs(mcautocorr-rho);
        % The mean is exact: the quadrature nodes are symmetric about zstar and the row
        % normalisation preserves that, so the first moment comes out exactly right.
        fprintf('%s znum=%i: mean is exact by symmetry [T1], this should be zero: %2.8e \n',cname,znum,err_mean(c_c))
        fprintf('%s znum=%i: variance error [T2] %2.3e, autocorr error [T2] %2.3e, runtime %2.6f s \n',cname,znum,err_var(c_c),err_ac(c_c),runtime(c_c))

        % --- CONDITIONAL moments. The checks above are all about the stationary distribution, but
        % half of what this literature argues about is conditional: Floden (2008) reports sigma_eps
        % beside sigma_z precisely because methods buy one at the other's expense, and his finding
        % about the Tauchen-Hussey sigma_z variant is that it gets the unconditional variance right
        % WHILE getting the conditional variance and the autocorrelation wrong. A method can pass
        % every unconditional check here and still misrepresent every conditional distribution.
        %
        % Truth: E[z'|z_i] = mew + rho*z_i, and Var[z'|z_i] = sigma^2, for every i.
        % Reported two ways - the worst row, which is where a method fails first (the grid edges,
        % where the conditional distribution is truncated), and the stationary-weighted average,
        % which is what Floden's sigma_eps row is.
        condmean_true=mew+rho*z_grid;
        condmean_hat=pi_z*z_grid;
        condvar_hat=zeros(znum,1);
        for z_c=1:znum
            condvar_hat(z_c)=sum(pi_z(z_c,:)'.*((z_grid-condmean_hat(z_c)).^2));
        end
        [~,~,~,sd_c]=MarkovChainMoments(z_grid,pi_z);
        wcondsd=sum(sd_c.*sqrt(condvar_hat)); % the stationary-weighted average conditional std dev
        err_condmean(c_c)=max(abs(condmean_hat-condmean_true));
        err_condvar(c_c)=max(abs(condvar_hat-sigma^2));
        err_condsd_w(c_c)=abs(wcondsd-sigma);
        fprintf('%s znum=%i: conditional mean, worst row [T2] %2.3e; conditional variance, worst row [T2] %2.3e \n',cname,znum,err_condmean(c_c),err_condvar(c_c))
        fprintf('%s znum=%i: weighted-average conditional sd %2.6f against sigma %2.6f, error [T2] %2.3e \n',cname,znum,wcondsd,sigma,err_condsd_w(c_c))

        % --- unconditional skewness and excess kurtosis. Truth is zero and zero for a gaussian
        % AR(1), and nothing else in this bank checks that a discretization does not MANUFACTURE
        % higher moments. Not hypothetical: P3 found the mixture method returning a platykurtic
        % chain on a coarse grid, so a method inventing shape is a live failure mode.
        [~,~,~,sd_u]=MarkovChainMoments(z_grid,pi_z);
        m_u=sum(sd_u.*z_grid);
        v_u=sum(sd_u.*(z_grid-m_u).^2);
        err_skew(c_c)=abs(sum(sd_u.*(z_grid-m_u).^3)/v_u^1.5);
        err_exkurt(c_c)=abs(sum(sd_u.*(z_grid-m_u).^4)/v_u^2-3);
        fprintf('%s znum=%i: unconditional skewness is zero [T2], this should be small: %2.3e \n',cname,znum,err_skew(c_c))
        fprintf('%s znum=%i: unconditional excess kurtosis is zero [T2], this should be small: %2.3e \n',cname,znum,err_exkurt(c_c))

        output.(cname).sweep(c_c).znum=znum;
        output.(cname).sweep(c_c).z_grid=z_grid;
        output.(cname).sweep(c_c).pi_z=pi_z;
    end
    output.(cname).err_condmean=err_condmean;
    output.(cname).err_condvar=err_condvar;
    output.(cname).err_condsd_w=err_condsd_w;
    output.(cname).err_skew=err_skew;
    output.(cname).err_exkurt=err_exkurt;
    output.(cname).err_mean=err_mean;
    output.(cname).err_var=err_var;
    output.(cname).err_ac=err_ac;
    output.(cname).runtime=runtime;
    output.(cname).nrepsused=nrepsused;

    %% Regression bars on the accuracy numbers above, for this calibration
    % These are set INSIDE the calibration loop on purpose: err_* is re-declared each time round, so
    % a block after the loop would only ever see the last calibration. They are regression bars, not
    % accuracy claims - a single threshold cannot serve both a coarse grid and a fine one - and come
    % from the 2026-08-28 run with about an order of magnitude of headroom. Without them the moment
    % errors are printed with no verdict, which is how B30 hid inside a clean-looking run summary.
    fprintf('%s: worst variance error over the sweep [T2], this should be below %g: %2.3e \n',cname,1,max(err_var))
    fprintf('%s: variance error at the finest grid [T2], this should be below %g: %2.3e \n',cname,0.1,err_var(end))
    fprintf('%s: worst autocorrelation error over the sweep [T2], this should be below %g: %2.3e \n',cname,1,max(err_ac))
    fprintf('%s: autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',cname,1,err_ac(end))
end

%% The baseSigma sweep: Floden's three variants
% Reported per calibration, and the ordering asserted on the persistent one, which is the regime
% Floden's paper is about.
znum=9;
bsname={'sigma_eps (Tauchen-Hussey)','sigma_z (Klein)','weighted (Floden, the default)'};
err_ac_bs=zeros(3,3); err_var_bs=zeros(3,3); % (calibration, variant)
for cal_c=1:3
    cname=calib.names{cal_c};
    mew=calib.(cname).mew; rho=calib.(cname).rho; sigma=calib.(cname).sigma;
    zstar=mew/(1-rho); varz=sigma^2/(1-rho^2); sdz=sigma/sqrt(1-rho^2);
    w=0.5+rho/4;
    bslist=[sigma, sdz, w*sigma+(1-w)*sdz];
    for b_c=1:3
        opts=struct(); opts.baseSigma=bslist(b_c);
        [z_grid,pi_z]=discretizeAR1_TauchenHussey(mew,rho,sigma,znum,opts);
        [~,mcvar,mcac]=MarkovChainMoments(z_grid,pi_z);
        err_var_bs(cal_c,b_c)=abs(mcvar-varz);
        err_ac_bs(cal_c,b_c)=abs(mcac-rho);
        fprintf('%s, baseSigma=%s: variance error %2.3e, autocorr error %2.3e \n',cname,bsname{b_c},err_var_bs(cal_c,b_c),err_ac_bs(cal_c,b_c))
    end
end
output.err_var_bs=err_var_bs;
output.err_ac_bs=err_ac_bs;
output.bsname=bsname;

% Floden's finding, asserted on the persistent calibration: the sigma_z variant is the bad one at
% high persistence. This is the whole point of his paper and of the toolkit's default.
%
% It turns out to be a stronger finding than his table shows. At rho=0.99 the sigma_z variant does
% not merely lose accuracy - the nodes sit so far apart that the off-diagonal transition
% probabilities underflow, the chain becomes numerically ABSORBING, its stationary distribution
% collapses onto a point, and the autocorrelation comes back as 0/0 = NaN. Floden pg 518: it
% "chooses nodes far from the mean, and with a small number of nodes the implied transition
% probabilities from one node to another are miniscule".
%
% So the check has to be written in a way that a NaN cannot silently satisfy. Comparing
% err_weighted < err_sigmaz is FALSE when the second is NaN, which reads as the ordering failing
% when in fact the ordering held in the strongest possible way. State the degeneracy directly.
fprintf('at rho=%g, the weighted default gives a non-degenerate chain [T2], this should be one: %i \n',calib.persistent.rho,isfinite(err_ac_bs(2,3))&&isfinite(err_var_bs(2,3)))
fprintf('at rho=%g, the sigma_z variant DEGENERATES (autocorrelation not finite, or variance error \n',calib.persistent.rho)
fprintf('   larger than the variance itself), this should be one: %i \n',~isfinite(err_ac_bs(2,2))||err_var_bs(2,2)>calib.persistent.sigma^2/(1-calib.persistent.rho^2))
fprintf('at rho=%g, the weighted default beats the sigma_z variant on the variance [T2], this should be one: %i \n',calib.persistent.rho,err_var_bs(2,3)<err_var_bs(2,2))
% And at moderate persistence, where nothing degenerates, the ordering is an ordinary comparison
fprintf('at rho=%g (moderate), the weighted default beats both other variants on the autocorrelation [T2], this should be one: %i \n',calib.moderate.rho,err_ac_bs(1,3)<err_ac_bs(1,1)&&err_ac_bs(1,3)<err_ac_bs(1,2))
% And the default really is the weighted variant
mew=calib.moderate.mew; rho=calib.moderate.rho; sigma=calib.moderate.sigma;
w=0.5+rho/4; sdz=sigma/sqrt(1-rho^2);
[zgD,pzD]=discretizeAR1_TauchenHussey(mew,rho,sigma,9,struct());
optsW=struct(); optsW.baseSigma=w*sigma+(1-w)*sdz;
[zgW,pzW]=discretizeAR1_TauchenHussey(mew,rho,sigma,9,optsW);
fprintf('the default baseSigma is w*sigma_eps+(1-w)*sigma_z with w=1/2+rho/4 [T0], this should be zero: %2.8e \n',max(abs(zgD-zgW))+max(abs(pzD-pzW),[],'all'))


%% Figure
figure(figure_c)
subplot(1,3,1)
plot(znums,max(output.moderate.err_var,10^(-18)),'o-',znums,max(output.persistent.err_var,10^(-18)),'s-',znums,max(output.drift.err_var,10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('TauchenHussey: |variance error|'); xlabel('znum'); legend(calib.names,'Location','best')
subplot(1,3,2)
bar(err_ac_bs); set(gca,'YScale','log'); set(gca,'XTickLabel',calib.names)
title('|autocorrelation error| by baseSigma'); legend(bsname,'Location','best')
subplot(1,3,3)
plot(znums,output.moderate.runtime,'o-',znums,output.persistent.runtime,'s-',znums,output.drift.runtime,'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('znum')

end
