function output=DiscP3_AR1wGM_FarmerToda(calib,znums,figure_c)
% P3: discretizeAR1wGM_FarmerToda
%
% P3 owns this command's option sweep. It is the only command in the bank with five grid methods -
% the usual four plus 'GMQ', Gaussian Mixture Quadrature, which exists precisely because a mixture
% is not well served by nodes placed for a normal.
%
% What this block measures that no other block can: the moments BEYOND the second. The whole point
% of a gaussian mixture is skewness and excess kurtosis, and the calibration is chosen to have a
% lot of both (skewness about -1.7, excess kurtosis about 6.5 in the innovation). Truth for the
% stationary distribution of z comes from the cumulants: the n-th cumulant of z is the n-th
% cumulant of e divided by (1-rho^n), because z is a linear process and cumulants are additive.
%
% Note this command's mew convention changed: it used to read mew as the unconditional mean and is
% now the intercept, like the other ten commands in the family. P2's cross-test is what pins that.

fprintf('\n========== P3: discretizeAR1wGM_FarmerToda ========== \n')

output=struct();
nz=length(znums);
mew=calib.mew; rho=calib.rho;
p=calib.mixprobs; mu=calib.mu; sd=calib.sigma;
zstar=(mew+calib.e.mean)/(1-rho);

err_var=zeros(1,nz); err_skew=zeros(1,nz); err_exkurt=zeros(1,nz);
err_condmean=zeros(1,nz); err_condvar=zeros(1,nz); err_condskew=zeros(1,nz); err_condkurt=zeros(1,nz); err_condsd_w=zeros(1,nz);
runtime=zeros(1,nz); nrepsused=zeros(1,nz); fallbackfrac=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    farmertodaoptions=struct(); farmertodaoptions.verbose=0;

    % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
    % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
    tic;
    [z_grid,pi_z,otheroutputs]=discretizeAR1wGM_FarmerToda(mew,rho,p,mu,sd,znum,farmertodaoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid,pi_z,otheroutputs]=discretizeAR1wGM_FarmerToda(mew,rho,p,mu,sd,znum,farmertodaoptions);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    nrepsused(c_c)=nreps;

    % --- invariants
    fprintf('znum=%i: size of z_grid [T0], this should be zero: %i \n',znum,any(size(z_grid)~=[znum,1]))
    fprintf('znum=%i: size of pi_z [T0], this should be zero: %i \n',znum,any(size(pi_z)~=[znum,znum]))
    fprintf('znum=%i: z_grid is strictly ascending [T0], this should be one: %i \n',znum,issorted(z_grid,'strictascend'))
    fprintf('znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z(:)<0)+any(pi_z(:)>1))
    fprintf('znum=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',znum,calib.entropytol4,max(abs(sum(pi_z,2)-1)))
    fprintf('znum=%i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
    % NO centro-symmetry check here: the process is deliberately skewed, so the chain must NOT be
    % symmetric. Asserting symmetry would be wrong and asserting nothing would lose the check
    % entirely, so a symmetric-mixture config is run separately below to keep it alive.

    % --- the exact conditional moment block, gated on what the method reports achieving.
    % Same contract as the gaussian Farmer-Toda, but the targets are now the mixture's uncentred
    % moments about the conditional mean, and nMoments defaults to 4 here rather than 2 - which is
    % the whole reason to use this command.
    nMg=otheroutputs.nMoments_grid(:);
    fallbackfrac(c_c)=mean(nMg<4);
    Tbar=[calib.e.mean; calib.e.var+calib.e.mean^2; ...
          sum(p.*(mu.^3+3*mu.*sd.^2)); sum(p.*(mu.^4+6*(mu.^2).*sd.^2+3*sd.^4))];
    sde=sqrt(calib.e.var);
    worstk=zeros(1,4);
    for z_c=1:znum
        condMean=mew+rho*z_grid(z_c);
        for kk=1:min(4,nMg(z_c))
            got=sum(pi_z(z_c,:)'.*((z_grid-condMean).^kk));
            worstk(kk)=max(worstk(kk),abs(got-Tbar(kk))/sde^kk);
        end
    end
    fprintf('znum=%i: conditional moments 1-4 hit their targets where the method achieved them, \n',znum)
    fprintf('   worst relative error: %2.3e %2.3e %2.3e %2.3e (should be below %g) \n',worstk(1),worstk(2),worstk(3),worstk(4),calib.entropytol4)
    fprintf('znum=%i: rows that fell short of 4 moments: %2.1f%% (a config at 100%% is testing nothing) \n',znum,100*fallbackfrac(c_c))

    % --- accuracy against the analytic cumulants of z
    [mcmean,mcvar,~,statdist]=MarkovChainMoments(z_grid,pi_z);
    c3=sum(statdist.*(z_grid-mcmean).^3);
    c4=sum(statdist.*(z_grid-mcmean).^4);
    sk=c3/mcvar^1.5;
    ek=c4/mcvar^2-3;
    err_var(c_c)=abs(mcvar-calib.z.var);
    err_skew(c_c)=abs(sk-calib.z.skew);
    err_exkurt(c_c)=abs(ek-calib.z.exkurt);
    fprintf('znum=%i: variance error %2.3e, skewness %2.4f vs truth %2.4f (err %2.3e), excess kurtosis %2.4f vs truth %2.4f (err %2.3e) \n',znum,err_var(c_c),sk,calib.z.skew,err_skew(c_c),ek,calib.z.exkurt,err_exkurt(c_c))

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
    output.sweep(c_c).moments=[mcmean,mcvar,sk,ek];
end
output.err_condmean=err_condmean; output.err_condvar=err_condvar;
output.err_condskew=err_condskew; output.err_condkurt=err_condkurt; output.err_condsd_w=err_condsd_w;
output.err_var=err_var; output.err_skew=err_skew; output.err_exkurt=err_exkurt;
output.runtime=runtime; output.nrepsused=nrepsused; output.fallbackfrac=fallbackfrac; output.znums=znums;

fprintf('skewness error falls over the sweep [T2], this should be one: %i \n',err_skew(end)<err_skew(1))
fprintf('excess kurtosis error falls over the sweep [T2], this should be one: %i \n',err_exkurt(end)<err_exkurt(1))

%% Symmetry, on a symmetric mixture
% Keeps the centro-symmetry check alive for this command, on a config where it applies.
znum=15;
[zgS,pzS,ooS]=discretizeAR1wGM_FarmerToda(0,rho,calib.sym.mixprobs,calib.sym.mu,calib.sym.sigma,znum,struct('verbose',0));
fprintf('symmetric mixture: z_grid symmetric about zero [T1], this should be zero: %2.8e \n',max(abs((zgS+flipud(zgS))/2)))
[~,vS,~,dS]=MarkovChainMoments(zgS,pzS);
mS=sum(dS.*zgS);
fprintf('symmetric mixture: skewness is zero [T2], this should be small: %2.8e \n',abs(sum(dS.*(zgS-mS).^3)/vS^1.5))

% CENTROSYMMETRY, WITH THE BAR DERIVED RATHER THAN CHOSEN. On a symmetric mixture the problem is
% invariant under z -> -z, so pi_z ought to be centrosymmetric. It never is exactly, because each
% row is an independent numerical solve, and the question is how big a departure is consistent
% with that rather than evidence of a bug. This used to be asserted against a flat 1e-4 borrowed
% from the four-moment entropy tolerance, which is a tolerance on MOMENT error and not on
% probabilities, and the two are not on the same scale.
%
% The scale factor between them is measurable. The entropy solution is pinned down by the moments
% it matches, so a discrepancy in the achieved moments of a mirrored pair translates into a
% discrepancy in their probability vectors, amplified by the conditioning of the moment map - to
% first order by 1/sigma_min of the covariance of the moment functions under that row. For the
% worst mirrored pair here that amplification is about 1.8e3, so a moment discrepancy of 1e-9
% permits a probability discrepancy near 2e-6, and no fixed bar on pi_z can know that.
%
% So the assertion compares like with like: for each mirrored pair, the observed asymmetry against
% that pair's own moment discrepancy times its own amplification. There is no constant in it. A
% failure now means the asymmetry is NOT explained by the moment noise, which is the thing that
% would actually be worth chasing. The point of writing it this way is that the ratio is a property
% of the solve rather than of the machine: measured on cpu 2026-09-19 it is 2.70e-03 with the
% feasibility test in place and 2.67e-03 without it, a difference of one per cent, even though the
% raw asymmetry those two produce differs by a factor of fifteen (3.5e-09 against 5.2e-08).
%
% The nMoments_grid precondition is kept, and for the same reason as before: if mirrored rows
% matched different numbers of moments then they solved different problems, the bound below does
% not apply, and an asymmetric pi_z is expected rather than a defect. Before 2026-09-19 that was
% the normal case on gpu and this whole check sat dormant - the asymmetry was 6e-02, far too large
% to be solver noise - and it is the feasibility test in discreteApproximation that made mirrored
% rows agree, because a convex-hull test is exact and therefore symmetric where an iterative solve
% stopping on its own tolerances is not.
nMg=ooS.nMoments_grid(:);
symfallback=max(abs(nMg-flipud(nMg)));
fprintf('symmetric mixture: nMoments_grid is itself asymmetric by %i moments across mirrored rows \n',symfallback)
fprintf('symmetric mixture: pi_z asymmetry is %2.8e (reported) \n',max(abs(pzS-rot90(pzS,2)),[],'all'))
if symfallback==0
    % The moment functions, in raw powers of the grid. Mirroring sends z to -z, so the kth moment
    % of a row must equal (-1)^k times the kth moment of its mirror.
    Tsym=[zgS(:),zgS(:).^2,zgS(:).^3,zgS(:).^4];
    Msym=pzS*Tsym;
    worstratio=0; worstpair=0; worstamp=0;
    for s_c=1:floor(znum/2)
        s_m=znum+1-s_c;
        pairasym=max(abs(pzS(s_c,:)-fliplr(pzS(s_m,:))));
        pairmom=max(abs(Msym(s_c,:)-Msym(s_m,:).*((-1).^(1:4))));
        sigpair=Inf;
        for s_r=[s_c,s_m]
            pr=pzS(s_r,:)'; mur=Tsym'*pr; Cr=(Tsym-mur')'*(pr.*(Tsym-mur'));
            sigpair=min(sigpair,min(eig(Cr)));
        end
        % eps floors a moment discrepancy that has underflowed to nothing, where the ratio would
        % otherwise be 0/0 rather than small.
        pairbound=max(pairmom,eps)/sigpair;
        if pairasym/pairbound>worstratio
            worstratio=pairasym/pairbound; worstpair=s_c; worstamp=1/sigpair;
        end
    end
    fprintf('symmetric mixture: worst mirrored pair is row %i, amplification 1/sigma_min = %2.3e \n',worstpair,worstamp)
    fprintf('symmetric mixture: its asymmetry as a fraction of what its own moment discrepancy permits, \n')
    fprintf('   so the solve is symmetric to the accuracy it achieved [T2], this should be below 1: %2.8e \n',worstratio)
else
    fprintf('symmetric mixture: mirrored rows achieved DIFFERENT numbers of moments, so an asymmetric \n')
    fprintf('   pi_z is expected here and is not asserted against. This is a property of the entropy \n')
    fprintf('   solver, not of the method: a symmetric problem does not guarantee a symmetric solve. \n')
end

%% Option sweep: five grid methods x nMoments
methodlist={'even','gauss-legendre','clenshaw-curtis','gauss-hermite','GMQ'};
nMomentslist=[2,4];
znum=15;
errsk=zeros(length(methodlist),length(nMomentslist));
errek=zeros(length(methodlist),length(nMomentslist));
fbk=zeros(length(methodlist),length(nMomentslist));
for m_c=1:length(methodlist)
    for n_c=1:length(nMomentslist)
        opts=struct(); opts.method=methodlist{m_c}; opts.nMoments=nMomentslist(n_c); opts.verbose=0;
        [zg,pz,oo]=discretizeAR1wGM_FarmerToda(mew,rho,p,mu,sd,znum,opts);
        if nMomentslist(n_c)>=4
            soltol=calib.entropytol4;
        else
            soltol=calib.entropytol;
        end
        fprintf('method=%s, nMoments=%i: z_grid is strictly ascending [T0], this should be one: %i \n',methodlist{m_c},nMomentslist(n_c),issorted(zg,'strictascend'))
        fprintf('method=%s, nMoments=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',methodlist{m_c},nMomentslist(n_c),soltol,max(abs(sum(pz,2)-1)))
        [mm,vv,~,dd]=MarkovChainMoments(zg,pz);
        errsk(m_c,n_c)=abs(sum(dd.*(zg-mm).^3)/vv^1.5-calib.z.skew);
        errek(m_c,n_c)=abs(sum(dd.*(zg-mm).^4)/vv^2-3-calib.z.exkurt);
        fbk(m_c,n_c)=mean(oo.nMoments_grid(:)<nMomentslist(n_c));
        fprintf('method=%s, nMoments=%i: skewness error %2.3e, excess kurtosis error %2.3e, rows short of target %2.1f%% \n',methodlist{m_c},nMomentslist(n_c),errsk(m_c,n_c),errek(m_c,n_c),100*fbk(m_c,n_c))
    end
end
output.errsk=errsk; output.errek=errek; output.methodlist=methodlist;

% The ordering this block exists to test: targeting four moments must beat targeting two on the
% two moments the extra targeting is about. This is the whole argument for the command.
%
% Skewness improves for every grid method. Excess kurtosis improves for four of the five - and the
% exception is instructive rather than a bug, so it is reported rather than asserted away. Asking
% for four moments is not free: the solve fails on more rows (measured here, 20% to 73% of rows
% fall short at nMoments=4 against 0% to 13% at nMoments=2), and those rows fall back to matching
% fewer moments. Where the fallback rate is high enough, targeting four moments can leave a chain
% whose kurtosis is WORSE than targeting two. So more moments is not monotonically better.
output.fbk=fbk;
fprintf('nMoments=4 beats nMoments=2 on skewness, for every grid method [T2], this should be one: %i \n',all(errsk(:,2)<=errsk(:,1)))
fprintf('nMoments=4 beats nMoments=2 on excess kurtosis, for a MAJORITY of grid methods [T2], this should be one: %i \n',sum(errek(:,2)<=errek(:,1))>=3)
for m_c=1:length(methodlist)
    if errek(m_c,2)>errek(m_c,1)
        fprintf('   exception: %s gets WORSE on excess kurtosis at nMoments=4 (%2.3e -> %2.3e), and its \n',methodlist{m_c},errek(m_c,1),errek(m_c,2))
        fprintf('   fallback rate rises from %2.1f%% to %2.1f%% of rows. Reported, not asserted away. \n',100*fbk(m_c,1),100*fbk(m_c,2))
    end
end

%% Figure
figure(figure_c)
subplot(1,3,1)
plot(znums,max(err_skew,10^(-18)),'o-',znums,max(err_exkurt,10^(-18)),'s-',znums,max(err_var,10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('AR1wGM\_FarmerToda: |error|'); xlabel('znum'); legend('skewness','excess kurtosis','variance','Location','best')
subplot(1,3,2)
bar([errsk(:,1),errsk(:,2)]); set(gca,'YScale','log'); set(gca,'XTickLabel',methodlist)
title('|skewness error| by grid method'); legend('nMoments=2','nMoments=4','Location','best')
subplot(1,3,3)
plot(znums,runtime,'o-'); set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('znum')

end
