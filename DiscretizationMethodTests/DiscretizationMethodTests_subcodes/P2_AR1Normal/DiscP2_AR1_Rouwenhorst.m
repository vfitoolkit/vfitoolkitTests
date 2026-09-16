function output=DiscP2_AR1_Rouwenhorst(calib,znums,figure_c)
% P2: discretizeAR1_Rouwenhorst
%
% P2 owns this command's option sweep, which is small - the method has one option (parallel) and
% no hyperparameter, which is part of its appeal.
%
% What makes this subcode worth its length is that Rouwenhorst matches the AR(1)'s mean, variance
% AND first-order autocorrelation EXACTLY, for any znum, by construction (Kopecky & Suen 2010).
% So all three are T1 exactness checks rather than T2 accuracy ones - the only method in this block
% for which that is true of the autocorrelation.
%
% The drift calibration is the one that matters most here. Rouwenhorst's grid is built around
% mew/(1-rho); if it were built around mew instead, the mean check below would be out by a factor
% of ten on that calibration and by nothing at all on the other two.

fprintf('\n========== P2: discretizeAR1_Rouwenhorst ========== \n')

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
        rouwenhorstoptions=struct();

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
        tic;
        [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,rouwenhorstoptions);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,rouwenhorstoptions);
            treps(r_c)=toc;
        end
        runtime(c_c)=median(treps);
        nrepsused(c_c)=nreps;

        % --- invariants
        fprintf('%s znum=%i: size of z_grid [T0], this should be zero: %i \n',cname,znum,any(size(z_grid)~=[znum,1]))
        fprintf('%s znum=%i: size of pi_z [T0], this should be zero: %i \n',cname,znum,any(size(pi_z)~=[znum,znum]))
        fprintf('%s znum=%i: z_grid is strictly ascending [T0], this should be one: %i \n',cname,znum,issorted(z_grid,'strictascend'))
        fprintf('%s znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',cname,znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        fprintf('%s znum=%i: rows of pi_z sum to one [T1], this should be zero: %2.8e \n',cname,znum,max(abs(sum(pi_z,2)-1)))
        fprintf('%s znum=%i: no NaN or Inf [T0], this should be zero: %i \n',cname,znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
        fprintf('%s znum=%i: z_grid symmetric about mew/(1-rho) [T1], this should be zero: %2.8e \n',cname,znum,max(abs((z_grid+flipud(z_grid))/2-zstar)))
        fprintf('%s znum=%i: pi_z is centrosymmetric [T1], this should be zero: %2.8e \n',cname,znum,max(abs(pi_z-rot90(pi_z,2)),[],'all'))
        % Rouwenhorst's grid spans exactly +-sqrt(znum-1)*sd(z) about the mean
        fprintf('%s znum=%i: grid half-width is sqrt(znum-1)*sd(z) [T1], this should be zero: %2.8e \n',cname,znum,abs((z_grid(end)-z_grid(1))/2-sqrt(znum-1)*sqrt(varz)))

        % --- the exactness guarantee. All three moments, at every znum. T1, not T2.
        [mcmean,mcvar,mcautocorr]=MarkovChainMoments(z_grid,pi_z);
        err_mean(c_c)=abs(mcmean-zstar);
        err_var(c_c)=abs(mcvar-varz);
        err_ac(c_c)=abs(mcautocorr-rho);
        fprintf('%s znum=%i: mean matched exactly [T1], this should be zero: %2.8e \n',cname,znum,err_mean(c_c))
        fprintf('%s znum=%i: variance matched exactly [T1], this should be zero: %2.8e \n',cname,znum,err_var(c_c))
        fprintf('%s znum=%i: autocorrelation matched exactly [T1], this should be zero: %2.8e \n',cname,znum,err_ac(c_c))
        fprintf('%s znum=%i: runtime %2.6f s (nreps=%i) \n',cname,znum,runtime(c_c),nreps)

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
    fprintf('%s: worst mean error over the sweep [T2], this should be below %g: %2.3e \n',cname,1e-08,max(err_mean))
    fprintf('%s: mean error at the finest grid [T2], this should be below %g: %2.3e \n',cname,1e-08,err_mean(end))
    fprintf('%s: worst variance error over the sweep [T2], this should be below %g: %2.3e \n',cname,1e-08,max(err_var))
    fprintf('%s: variance error at the finest grid [T2], this should be below %g: %2.3e \n',cname,1e-08,err_var(end))
    fprintf('%s: worst autocorrelation error over the sweep [T2], this should be below %g: %2.3e \n',cname,1e-08,max(err_ac))
    fprintf('%s: autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',cname,1e-08,err_ac(end))
end

%% The stationary distribution is binomial, which is a second closed form
mew=calib.moderate.mew; rho=calib.moderate.rho; sigma=calib.moderate.sigma;
for znum=[5,9,15]
    [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,struct());
    [~,~,~,mcstatdist]=MarkovChainMoments(z_grid,pi_z);
    sbinom=zeros(znum,1);
    for jj=1:znum
        sbinom(jj)=nchoosek(znum-1,jj-1);
    end
    sbinom=sbinom/2^(znum-1);
    fprintf('znum=%i: stationary dist is binomial [T1], this should be zero: %2.8e \n',znum,max(abs(mcstatdist(:)-sbinom(:))))
end

%% The transition matrix does not depend on sigma, only the grid does
% p=q=(1+rho)/2 is a function of rho alone, so scaling sigma must rescale the grid about zstar and
% leave pi_z untouched. Cheap, and it isolates the grid construction from the transition.
znum=15;
[zgA,pzA]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,struct());
[zgB,pzB]=discretizeAR1_Rouwenhorst(mew,rho,2*sigma,znum,struct());
zstar=mew/(1-rho);
fprintf('doubling sigma leaves pi_z unchanged [T0], this should be zero: %2.8e \n',max(abs(pzA-pzB),[],'all'))
fprintf('doubling sigma doubles the grid about its mean [T1], this should be zero: %2.8e \n',max(abs((zgB-zstar)-2*(zgA-zstar))))

%% Option sweep
opts0=struct(); opts0.parallel=0;
opts1=struct(); opts1.parallel=1;
[zga,pza]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,opts0);
[zgb,pzb]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,opts1);
fprintf('parallel=0 vs 1: z_grid [T0], this should be zero: %2.8e \n',max(abs(zga-zgb)))
fprintf('parallel=0 vs 1: pi_z [T0], this should be zero: %2.8e \n',max(abs(pza-pzb),[],'all'))
if gpuDeviceCount>0
    opts2=struct(); opts2.parallel=2;
    [zgc,pzc]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,opts2);
    fprintf('parallel=2 returns gpuArrays [T0], this should be one: %i \n',isa(zgc,'gpuArray')&&isa(pzc,'gpuArray'))
    % Rouwenhorst builds on the cpu and moves the result, so unlike Tauchen this one IS exact
    fprintf('parallel=1 vs 2: z_grid [T0], this should be zero: %2.8e \n',max(abs(zgb-gather(zgc))))
    fprintf('parallel=1 vs 2: pi_z [T0], this should be zero: %2.8e \n',max(abs(pzb-gather(pzc)),[],'all'))
else
    fprintf('parallel=2 skipped: no gpu on this machine \n')
end


%% Figure
figure(figure_c)
subplot(1,3,1)
plot(znums,max(output.moderate.err_var,10^(-18)),'o-',znums,max(output.persistent.err_var,10^(-18)),'s-',znums,max(output.drift.err_var,10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('Rouwenhorst: |variance error|'); xlabel('znum'); legend(calib.names,'Location','best')
subplot(1,3,2)
plot(znums,max(output.moderate.err_ac,10^(-18)),'o-',znums,max(output.persistent.err_ac,10^(-18)),'s-',znums,max(output.drift.err_ac,10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('|autocorrelation error| (exact by construction)'); xlabel('znum')
subplot(1,3,3)
plot(znums,output.moderate.runtime,'o-',znums,output.persistent.runtime,'s-',znums,output.drift.runtime,'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('znum')

end
