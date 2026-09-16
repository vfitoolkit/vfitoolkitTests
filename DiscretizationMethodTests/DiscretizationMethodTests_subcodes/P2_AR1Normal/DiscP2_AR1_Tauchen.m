function output=DiscP2_AR1_Tauchen(calib,znums,figure_c)
% P2: discretizeAR1_Tauchen
%
% P2 owns this command's option sweep. Same template as P1's method subcodes: config sweep,
% inlined invariant block, inlined timing protocol, accuracy against analytic truth.
%
% Truth for a stationary AR(1) is closed form: E(z)=mew/(1-rho), Var(z)=sigma^2/(1-rho^2),
% autocorrelation=rho. Tauchen assigns probability by cdf mass between midpoints, so it
% approximates all three - these are T2 checks.
%
% Two things here are exercised for the first time by this bank:
%   - the mew~=0 calibration, which is what pins the grid-centring convention. At mew=0 a command
%     centred on mew rather than on mew/(1-rho) looks identical to a correct one.
%   - tauchenoptions.dshift, which was dead until recently (its guard used exist() on a dotted
%     name, which never resolves), so the code path below has never run before.

fprintf('\n========== P2: discretizeAR1_Tauchen ========== \n')

output=struct();
nz=length(znums);
Tauchen_q=3;

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
        tauchenoptions=struct();

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
        tic;
        [z_grid,pi_z]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,tauchenoptions);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,tauchenoptions);
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
        % Symmetry, stated about zstar rather than about zero. The process is mirror-symmetric
        % about its unconditional mean (put y=z-zstar and y'=rho*y+e, whose law is invariant under
        % y -> -y), and every method here builds a grid symmetric about zstar and assigns
        % probability by a rule depending only on grid points relative to conditional means. So the
        % chain inherits it. Stated this way the check applies at mew~=0 too, which is where it has
        % teeth: a grid centred on mew instead of on mew/(1-rho) fails it.
        fprintf('%s znum=%i: z_grid symmetric about mew/(1-rho) [T1], this should be zero: %2.8e \n',cname,znum,max(abs((z_grid+flipud(z_grid))/2-zstar)))
        fprintf('%s znum=%i: pi_z is centrosymmetric [T1], this should be zero: %2.8e \n',cname,znum,max(abs(pi_z-rot90(pi_z,2)),[],'all'))
        fprintf('%s znum=%i: grid half-width is Tauchen_q*sd(z) [T1], this should be zero: %2.8e \n',cname,znum,abs((z_grid(end)-z_grid(1))/2-Tauchen_q*sqrt(varz)))

        % --- accuracy against closed form
        [mcmean,mcvar,mcautocorr]=MarkovChainMoments(z_grid,pi_z);
        err_mean(c_c)=abs(mcmean-zstar);
        err_var(c_c)=abs(mcvar-varz);
        err_ac(c_c)=abs(mcautocorr-rho);
        fprintf('%s znum=%i: mean error [T2] %2.3e, variance error [T2] %2.3e, autocorr error [T2] %2.3e, runtime %2.6f s \n',cname,znum,err_mean(c_c),err_var(c_c),err_ac(c_c),runtime(c_c))

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
        output.(cname).sweep(c_c).moments=[mcmean,mcvar,mcautocorr];
    end
    % The mean is exact by symmetry OF THE DISCRETIZATION - the grid is symmetric about zstar and
    % the bin masses are symmetric with it, and both of those are asserted at T1 above and hold at
    % machine precision. What is measured here is different: the mean comes back through
    % MarkovChainMoments, whose stationary-distribution solve has its own error, and that error
    % grows as the chain becomes nearly reducible. At rho=0.99 with znum=5 it reaches 6e-09 while
    % every larger grid is at 1e-13 or better. So this is a check on the SOLVE, not on the
    % discretization, and it is tiered accordingly.
    fprintf('%s: mean through MarkovChainMoments, worst over the sweep, this should be below 1e-07: %2.8e \n',cname,max(err_mean))
    fprintf('%s: and at the largest znum it is back to machine precision [T1], this should be zero: %2.8e \n',cname,err_mean(end))
    % The variance converges to the TRUNCATION FLOOR at fixed Tauchen_q, not to zero (P1 measured
    % this: it falls, then flattens, then can tick back up). So assert the fall over the first half
    % of the sweep only.
    fprintf('%s: variance error falls over the first half of the sweep [T2], this should be one: %i \n',cname,err_var(4)<err_var(1))
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
    fprintf('%s: worst variance error over the sweep [T2], this should be below %g: %2.3e \n',cname,0.1,max(err_var))
    fprintf('%s: variance error at the finest grid [T2], this should be below %g: %2.3e \n',cname,0.01,err_var(end))
    fprintf('%s: worst autocorrelation error over the sweep [T2], this should be below %g: %2.3e \n',cname,0.1,max(err_ac))
    fprintf('%s: autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',cname,0.001,err_ac(end))
end

%% Option sweep, on the moderate calibration
mew=calib.moderate.mew; rho=calib.moderate.rho; sigma=calib.moderate.sigma;
zstar=mew/(1-rho); sdz=sigma/sqrt(1-rho^2);
znum=15;

% Tauchen_q: wider grid, less truncation, coarser spacing. A trade, not an ordering.
for qq=[2,3,4]
    [z_grid,pi_z]=discretizeAR1_Tauchen(mew,rho,sigma,znum,qq,struct());
    [~,mcvar,~]=MarkovChainMoments(z_grid,pi_z);
    fprintf('Tauchen_q=%i: grid half-width matches Tauchen_q*sd(z) [T1], this should be zero: %2.8e \n',qq,abs((z_grid(end)-z_grid(1))/2-qq*sdz))
    fprintf('Tauchen_q=%i: variance error [T2] %2.3e \n',qq,abs(mcvar-sigma^2/(1-rho^2)))
end

% z_grid option: passing back the grid the command just built must reproduce its pi_z EXACTLY.
% This is the newest code in the family and the least covered.
[zg0,pz0]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,struct());
optsug=struct(); optsug.z_grid=zg0;
[zg1,pz1]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,optsug);
% The grid itself round-trips exactly. pi_z does NOT, and should not be expected to: the usergrid
% path builds the bin edges as midpoints between adjacent grid points, (z(1:end-1)+z(2:end))/2,
% while the default path builds them as z +- omega/2 with omega=z(2)-z(1). Algebraically identical
% on an evenly spaced grid, but not the same floating-point operations, so this is T1 not T0.
fprintf('z_grid option round-trip: z_grid [T0], this should be zero: %2.8e \n',max(abs(zg0-zg1)))
fprintf('z_grid option round-trip: pi_z [T1], this should be zero: %2.8e \n',max(abs(pz0-pz1),[],'all'))
% and it must ignore Tauchen_q when a grid is supplied. Compared against the usergrid path itself,
% so that this isolates Tauchen_q rather than re-measuring the midpoint-vs-omega difference.
optsug2=struct(); optsug2.z_grid=zg0;
[~,pz2]=discretizeAR1_Tauchen(mew,rho,sigma,znum,7,optsug2);
fprintf('z_grid option makes Tauchen_q irrelevant [T0], this should be zero: %2.8e \n',max(abs(pz1-pz2),[],'all'))
% wrong length must error
try
    optsbad=struct(); optsbad.z_grid=zg0(1:end-1);
    discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,optsbad);
    fprintf('z_grid of the wrong length should have errored, this should be one: 0 \n')
catch
    fprintf('z_grid of the wrong length errors as documented, this should be one: 1 \n')
end

% dshift. This option was dead until its guard was fixed, so this is the first time the code path
% has run. It shifts the bin edges by a constant, so on a FIXED grid, shifting the edges by c is
% the same as moving the conditional mean by -c, i.e. the same as replacing mew by mew-c. That is
% a real oracle for an option that otherwise has none.
cshift=0.05;
optsA=struct(); optsA.z_grid=zg0; optsA.dshift=cshift;
[~,pzA]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,optsA);
optsB=struct(); optsB.z_grid=zg0; % same grid, no shift, but mew moved instead
[~,pzB]=discretizeAR1_Tauchen(mew-cshift,rho,sigma,znum,Tauchen_q,optsB);
fprintf('dshift=%g equals shifting mew by -%g on the same grid [T1], this should be zero: %2.8e \n',cshift,cshift,max(abs(pzA-pzB),[],'all'))
optsZ=struct(); optsZ.z_grid=zg0; optsZ.dshift=0;
[~,pzZ]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,optsZ);
% Compared against the usergrid path WITHOUT dshift, not against the default path: otherwise this
% would be re-measuring the midpoint-vs-omega difference above rather than isolating dshift.
fprintf('dshift=0 is the same as not setting dshift [T0], this should be zero: %2.8e \n',max(abs(pz1-pzZ),[],'all'))

% parallel
opts0=struct(); opts0.parallel=0;
opts1=struct(); opts1.parallel=1;
[zga,pza]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,opts0);
[zgb,pzb]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,opts1);
fprintf('parallel=0 vs 1: z_grid [T0], this should be zero: %2.8e \n',max(abs(zga-zgb)))
fprintf('parallel=0 vs 1: pi_z [T0], this should be zero: %2.8e \n',max(abs(pza-pzb),[],'all'))
if gpuDeviceCount>0
    opts2=struct(); opts2.parallel=2;
    [zgc,pzc]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,opts2);
    fprintf('parallel=2 returns gpuArrays [T0], this should be one: %i \n',isa(zgc,'gpuArray')&&isa(pzc,'gpuArray'))
    % z_grid is built on the cpu in both branches and then moved, so it is a strict zero. pi_z can
    % sit a couple of ulps off: erfc is not correctly rounded, so the cpu and gpu libraries may
    % disagree in the last bit. Was T1 until the gpu branch stopped building the cdf from 1+erf.
    fprintf('parallel=1 vs 2: z_grid [T0], this should be zero: %2.8e \n',max(abs(zgb-gather(zgc))))
    fprintf('parallel=1 vs 2: pi_z [T0], this should be zero: %2.8e \n',max(abs(pzb-gather(pzc)),[],'all'))
else
    fprintf('parallel=2 skipped: no gpu on this machine \n')
end


%% Figure
figure(figure_c)
for cal_c=1:3
    cname=calib.names{cal_c};
    subplot(2,3,cal_c)
    plot(znums,max(output.(cname).err_var,10^(-18)),'o-',znums,max(output.(cname).err_ac,10^(-18)),'s-')
    set(gca,'XScale','log'); set(gca,'YScale','log')
    title(['Tauchen: ',cname]); xlabel('znum'); legend('variance','autocorr','Location','best')
    subplot(2,3,3+cal_c)
    plot(znums,output.(cname).runtime,'o-')
    set(gca,'XScale','log'); set(gca,'YScale','log')
    title(['runtime (s): ',cname]); xlabel('znum')
end

end
