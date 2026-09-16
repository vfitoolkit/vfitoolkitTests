function output=DiscP2_AR1_FarmerToda(calib,znums,figure_c)
% P2: discretizeAR1_FarmerToda
%
% P2 owns this command's option sweep, and it is the largest in the bank: four grid methods, four
% moment targets, plus nSigmas, the user-grid option and parallel.
%
% This subcode also carries the EXACT CONDITIONAL MOMENT block, which is the sharpest check
% anywhere in the bank. Farmer-Toda does not merely approximate the conditional moments, it
% TARGETS them: for each row i of pi_z and each k=1..nMoments,
%
%     sum_j pi(i,j)*(z_j - condmean_i)^k  ==  T_k,     condmean_i = mew + rho*z_i
%     with T = [0, sigma^2, 0, 3*sigma^4]
%
% and the method reports, in otheroutputs.nMoments_grid, how many moments it actually achieved
% from each grid point (it falls back at the edges, where the target is not attainable). So the
% check is gated on that: a row that fell back to 2 moments is not asked about the 3rd and 4th.
%
% Calibration-independent, and it tests the entropy solver through its own contract rather than
% through an accuracy threshold. The count of fallback rows is reported per config, because a
% config that falls back everywhere is a config that is testing nothing.

fprintf('\n========== P2: discretizeAR1_FarmerToda ========== \n')

output=struct();
nz=length(znums);
entropytol=calib.entropytol;

for cal_c=1:3
    cname=calib.names{cal_c};
    mew=calib.(cname).mew; rho=calib.(cname).rho; sigma=calib.(cname).sigma;
    zstar=mew/(1-rho); varz=sigma^2/(1-rho^2);
    fprintf('\n--- calibration %s: mew=%g, rho=%g, sigma=%g (so E(z)=%g, sd(z)=%g) --- \n',cname,mew,rho,sigma,zstar,sqrt(varz))

    err_mean=zeros(1,nz); err_var=zeros(1,nz); err_ac=zeros(1,nz);
    err_condmean=zeros(1,nz); err_condvar=zeros(1,nz); err_condsd_w=zeros(1,nz);
    err_skew=zeros(1,nz); err_exkurt=zeros(1,nz); runtime=zeros(1,nz); nrepsused=zeros(1,nz); fallbackfrac=zeros(1,nz);
    for c_c=1:nz
        znum=znums(c_c);
        farmertodaoptions=struct();
        farmertodaoptions.verbose=0; % the per-row fallback warnings would swamp the diary

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
        tic;
        [z_grid,pi_z,otheroutputs]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,farmertodaoptions);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z,otheroutputs]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,farmertodaoptions);
            treps(r_c)=toc;
        end
        runtime(c_c)=median(treps);
        nrepsused(c_c)=nreps;

        % --- invariants
        fprintf('%s znum=%i: size of z_grid [T0], this should be zero: %i \n',cname,znum,any(size(z_grid)~=[znum,1]))
        fprintf('%s znum=%i: size of pi_z [T0], this should be zero: %i \n',cname,znum,any(size(pi_z)~=[znum,znum]))
        fprintf('%s znum=%i: z_grid is strictly ascending [T0], this should be one: %i \n',cname,znum,issorted(z_grid,'strictascend'))
        fprintf('%s znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',cname,znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        fprintf('%s znum=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',cname,znum,entropytol,max(abs(sum(pi_z,2)-1)))
        fprintf('%s znum=%i: no NaN or Inf [T0], this should be zero: %i \n',cname,znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
        fprintf('%s znum=%i: z_grid symmetric about mew/(1-rho) [T1], this should be zero: %2.8e \n',cname,znum,max(abs((z_grid+flipud(z_grid))/2-zstar)))
        fprintf('%s znum=%i: pi_z is centrosymmetric, this should be below %g: %2.8e \n',cname,znum,entropytol,max(abs(pi_z-rot90(pi_z,2)),[],'all'))

        % --- the exact conditional moment block, gated on what the method says it achieved
        nMg=otheroutputs.nMoments_grid(:);
        fallbackfrac(c_c)=mean(nMg<2); % the default target is 2 moments
        Tbar=[0; sigma^2; 0; 3*sigma^4];
        worstk=zeros(1,4);
        for z_c=1:znum
            condMean=mew+rho*z_grid(z_c);
            for kk=1:min(4,nMg(z_c))
                got=sum(pi_z(z_c,:)'.*((z_grid-condMean).^kk));
                % scale by sigma^kk so the tolerance means the same thing at every order
                worstk(kk)=max(worstk(kk),abs(got-Tbar(kk))/sigma^kk);
            end
        end
        fprintf('%s znum=%i: conditional moments hit their targets where the method says it achieved them \n',cname,znum)
        fprintf('   worst relative error, moment 1: %2.3e, moment 2: %2.3e (should be below %g) \n',worstk(1),worstk(2),entropytol)
        fprintf('%s znum=%i: rows that fell back below 2 moments: %2.1f%% (a config at 100%% is testing nothing) \n',cname,znum,100*fallbackfrac(c_c))

        % --- accuracy against closed form
        [mcmean,mcvar,mcautocorr]=MarkovChainMoments(z_grid,pi_z);
        err_mean(c_c)=abs(mcmean-zstar);
        err_var(c_c)=abs(mcvar-varz);
        err_ac(c_c)=abs(mcautocorr-rho);
        fprintf('%s znum=%i: mean error %2.3e, variance error %2.3e, autocorr error %2.3e, runtime %2.6f s \n',cname,znum,err_mean(c_c),err_var(c_c),err_ac(c_c),runtime(c_c))

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
        output.(cname).sweep(c_c).nMoments_grid=nMg;
        output.(cname).sweep(c_c).worstk=worstk;
    end
    % THE ACCURACY OF THIS METHOD IS NOT MONOTONE IN znum, and that is worth stating loudly,
    % because it is the opposite of what anyone would assume. As the grid grows the entropy solve
    % fails on more and more rows, those rows silently fall back to matching fewer moments, and the
    % unconditional variance degrades with them. Measured on the moderate calibration: the variance
    % error runs 2.7e-12 at znum=5 but 1.8e-03 at znum=101, while the fraction of rows falling back
    % below the 2-moment target goes 0%, 0%, 0%, 6.5%, 27.5%, 45.5%. At the largest grid, plain
    % Tauchen is the more accurate of the two.
    fprintf('%s: fallback rate rises with znum: ',cname)
    for c_c=1:nz
        fprintf('%i:%2.1f%%  ',znums(c_c),100*fallbackfrac(c_c));
    end
    fprintf(' \n')
    fprintf('%s: variance error at the smallest and largest znum: %2.3e and %2.3e \n',cname,err_var(1),err_var(end))
    fprintf('%s: accuracy degrades as the solve starts failing - the variance error at the highest \n',cname)
    fprintf('   fallback rate is worse than at the lowest [T2], this should be one: %i \n',err_var(fallbackfrac==max(fallbackfrac))>=err_var(fallbackfrac==min(fallbackfrac)))
    % That is a real ordering, not a threshold: it ties the accuracy loss to its cause.

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
    fprintf('%s: worst mean error over the sweep [T2], this should be below %g: %2.3e \n',cname,1e-06,max(err_mean))
    fprintf('%s: mean error at the finest grid [T2], this should be below %g: %2.3e \n',cname,1e-06,err_mean(end))
    fprintf('%s: worst variance error over the sweep [T2], this should be below %g: %2.3e \n',cname,0.1,max(err_var))
    fprintf('%s: variance error at the finest grid [T2], this should be below %g: %2.3e \n',cname,0.01,err_var(end))
    fprintf('%s: worst autocorrelation error over the sweep [T2], this should be below %g: %2.3e \n',cname,0.1,max(err_ac))
    fprintf('%s: autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',cname,0.01,err_ac(end))
    output.(cname).fallbackfrac=fallbackfrac;
end

%% Option sweep, on the moderate calibration
mew=calib.moderate.mew; rho=calib.moderate.rho; sigma=calib.moderate.sigma;
zstar=mew/(1-rho); varz=sigma^2/(1-rho^2);
znum=15;
methodlist={'even','gauss-legendre','clenshaw-curtis','gauss-hermite'};
nMomentslist=[1,2,4];
errgrid=zeros(length(methodlist),length(nMomentslist));
for m_c=1:length(methodlist)
    for n_c=1:length(nMomentslist)
        opts=struct(); opts.method=methodlist{m_c}; opts.nMoments=nMomentslist(n_c); opts.verbose=0;
        [z_grid,pi_z,oo]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,opts);
        % The entropy solve lands further from the exact answer the harder the problem is, and
        % nMoments is the biggest lever on that: measured here, the asymmetry of pi_z runs 1e-13 to
        % 1e-9 at nMoments 1 and 2, but 3e-6 to 6e-6 at nMoments=4. So the tolerance is
        % nMoments-dependent rather than one number. (gauss-hermite is the exception at 2e-10,
        % because its prior q=W is exactly symmetric where the others use W.*normpdf, which is not.)
        if nMomentslist(n_c)>=4
            soltol=10^(-4);
        else
            soltol=10^(-7);
        end
        % Invariants must hold for every grid method, including clenshaw-curtis, which is the only
        % one that builds its nodes in descending order and flips them.
        fprintf('method=%s, nMoments=%i: z_grid is strictly ascending [T0], this should be one: %i \n',methodlist{m_c},nMomentslist(n_c),issorted(z_grid,'strictascend'))
        fprintf('method=%s, nMoments=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',methodlist{m_c},nMomentslist(n_c),soltol,max(abs(sum(pi_z,2)-1)))
        fprintf('method=%s, nMoments=%i: pi_z is centrosymmetric, this should be below %g: %2.8e \n',methodlist{m_c},nMomentslist(n_c),soltol,max(abs(pi_z-rot90(pi_z,2)),[],'all'))
        [~,mcvar,~]=MarkovChainMoments(z_grid,pi_z);
        errgrid(m_c,n_c)=abs(mcvar-varz);
        fprintf('method=%s, nMoments=%i: variance error %2.3e, rows below 2 moments %2.1f%% \n',methodlist{m_c},nMomentslist(n_c),errgrid(m_c,n_c),100*mean(oo.nMoments_grid(:)<2))
    end
end
output.errgrid=errgrid;
output.methodlist=methodlist;
output.nMomentslist=nMomentslist;
% Asking for 2 moments must not do worse on the variance than asking for 1
fprintf('nMoments=2 beats nMoments=1 on the variance, for every grid method [T2], this should be one: %i \n',all(errgrid(:,2)<=errgrid(:,1)))

% nSigmas
for ns=[2,3]
    opts=struct(); opts.method='even'; opts.nSigmas=ns; opts.verbose=0;
    [z_grid,~,~]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,opts);
    fprintf('nSigmas=%i: grid half-width is nSigmas*sd(z) [T1], this should be zero: %2.8e \n',ns,abs((z_grid(end)-z_grid(1))/2-ns*sqrt(varz)))
end

% z_grid option round-trip
optsbase=struct(); optsbase.method='even'; optsbase.verbose=0;
[zg0,pz0,~]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,optsbase);
optsug=struct(); optsug.z_grid=zg0; optsug.verbose=0;
[zg1,pz1,~]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,optsug);
fprintf('z_grid option round-trip: z_grid [T0], this should be zero: %2.8e \n',max(abs(zg0(:)-zg1(:))))
% Not T0 on pi_z: with a user grid the prior q is built as for 'even' regardless of method, and
% the entropy solve starts from a different initial guess, so it converges to the same answer to
% solver tolerance rather than bit-for-bit.
fprintf('z_grid option round-trip: pi_z, this should be below %g: %2.8e \n',entropytol,max(abs(pz0-pz1),[],'all'))

% error paths
try
    discretizeAR1_FarmerToda(mew,rho,sigma,2,struct());
    fprintf('znum<3 should have errored, this should be one: 0 \n')
catch
    fprintf('znum<3 errors as documented, this should be one: 1 \n')
end
try
    optsbad=struct(); optsbad.nMoments=5;
    discretizeAR1_FarmerToda(mew,rho,sigma,15,optsbad);
    fprintf('nMoments=5 should have errored, this should be one: 0 \n')
catch
    fprintf('nMoments=5 errors as documented, this should be one: 1 \n')
end
try
    discretizeAR1_FarmerToda(mew,1.0,sigma,15,struct());
    fprintf('rho=1 should have errored, this should be one: 0 \n')
catch
    fprintf('rho=1 errors as documented, this should be one: 1 \n')
end


%% Figure
figure(figure_c)
subplot(1,3,1)
plot(znums,max(output.moderate.err_var,10^(-18)),'o-',znums,max(output.persistent.err_var,10^(-18)),'s-',znums,max(output.drift.err_var,10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('FarmerToda: |variance error|'); xlabel('znum'); legend(calib.names,'Location','best')
subplot(1,3,2)
bar(errgrid); set(gca,'YScale','log'); set(gca,'XTickLabel',methodlist)
title('|variance error| by grid method'); legend('nMoments=1','nMoments=2','nMoments=4','Location','best')
subplot(1,3,3)
plot(znums,output.moderate.runtime,'o-',znums,output.persistent.runtime,'s-',znums,output.drift.runtime,'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('znum'); legend(calib.names,'Location','best')

end
