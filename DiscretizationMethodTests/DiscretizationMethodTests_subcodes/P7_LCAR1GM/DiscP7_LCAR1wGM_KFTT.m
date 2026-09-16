function output=DiscP7_LCAR1wGM_KFTT(calib,znums,figure_c)
% P7: discretizeLifeCycleAR1wGM_KFTT
%
% The most complex command in the discretization surface, and until this block the least tested
% thing in it: before P7 it appeared in exactly one subcode of the whole bank, and only for the B27
% GMQ regression, so its own behaviour had never been measured.
%
% ALL FOUR MOMENTS GET AN EXACT AGE PROFILE HERE, from the cumulant recursion in DiscSetup_LCAR1GM.
% That is the difference from P6, which had truth for mean, variance and autocorrelation only. The
% two extra moments - skewness and excess kurtosis - are the reason gaussian mixtures are used at
% all, so they are what this block is really measuring; the first two are hygiene.
%
% MarkovChainMoments_FHorz returns mean, variance and autocorrelation but not the higher moments,
% so skewness and excess kurtosis are computed here from the age-conditional distribution, which
% the same command returns as its fourth output.

fprintf('\n========== P7: discretizeLifeCycleAR1wGM_KFTT ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho;
p=calib.vary.mixprobs_i; mu=calib.vary.mu_i; sd=calib.vary.sigma_i;
mewzT=calib.vary.mewz; varzT=calib.vary.varz; skewT=calib.vary.skewz; exkurtT=calib.vary.exkurtz;
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_skew=zeros(1,nz); err_exkurt=zeros(1,nz);
runtime=zeros(1,nz); fallback=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    ko=struct(); ko.verbose=0;

    tic;
    [z_grid_J,pi_z_J,jequaloneDistz,oo]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,ko);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,ko);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);

    % --- invariants
    fprintf('znum=%3i: size of z_grid_J is znum-by-J [T0], this should be zero: %i \n',znum,any(size(z_grid_J)~=[znum,J]))
    fprintf('znum=%3i: size of pi_z_J is znum-by-znum-by-(J-1) [T0], this should be zero: %i \n',znum,any(size(pi_z_J)~=[znum,znum,J-1]))
    asc=1;
    for j_c=1:J
        asc=asc*issorted(z_grid_J(:,j_c),'strictascend');
    end
    fprintf('znum=%3i: the grid at every age is strictly ascending [T0], this should be one: %i \n',znum,asc)
    fprintf('znum=%3i: pi_z_J is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z_J(:)<0)+any(pi_z_J(:)>1))
    fprintf('znum=%3i: rows of pi_z_J sum to one at every age, this should be below %g: %2.8e \n',znum,calib.entropytol,max(abs(sum(pi_z_J,2)-1),[],'all'))
    fprintf('znum=%3i: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',znum,abs(sum(jequaloneDistz)-1))
    fprintf('znum=%3i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid_J(:)))+any(~isfinite(pi_z_J(:))))
    nMg=oo.nMoments_grid;
    fallback(c_c)=mean(nMg(:)<4);
    fprintf('znum=%3i: conditional distributions matching fewer than 4 moments: %2.1f%% \n',znum,100*fallback(c_c))

    % --- the four age profiles
    [mcmean,mcvar,~,mcdist]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    skj=zeros(1,J); ekj=zeros(1,J);
    for j_c=1:J
        d=gather(mcdist(:,j_c)); d=d/sum(d);
        g=z_grid_J(:,j_c);
        m1=sum(d.*g); v=sum(d.*(g-m1).^2);
        skj(j_c)=sum(d.*(g-m1).^3)/v^1.5;
        ekj(j_c)=sum(d.*(g-m1).^4)/v^2-3;
    end
    [err_mean(c_c),jm]=max(abs(mcmean(:)'-mewzT));
    [err_var(c_c),jv]=max(abs(mcvar(:)'-varzT));
    [err_skew(c_c),js]=max(abs(skj-skewT));
    [err_exkurt(c_c),jk]=max(abs(ekj-exkurtT));
    fprintf('znum=%3i: worst age for the mean is j=%i, error [T2] %2.3e \n',znum,jm,err_mean(c_c))
    fprintf('znum=%3i: worst age for the variance is j=%i, error [T2] %2.3e \n',znum,jv,err_var(c_c))
    fprintf('znum=%3i: worst age for the SKEWNESS is j=%i, error [T2] %2.3e \n',znum,js,err_skew(c_c))
    fprintf('znum=%3i: worst age for the EXCESS KURTOSIS is j=%i, error [T2] %2.3e \n',znum,jk,err_exkurt(c_c))
    fprintf('znum=%3i: runtime %2.6f s \n',znum,runtime(c_c))
    output.sweep(c_c).var=mcvar(:)'; output.sweep(c_c).skew=skj; output.sweep(c_c).exkurt=ekj;
end

output.err_mean=err_mean; output.err_var=err_var; output.err_skew=err_skew; output.err_exkurt=err_exkurt;
output.runtime=runtime; output.fallback=fallback; output.znums=znums;

% Convergence, on the moment the block exists for. Reported as a floor rather than as monotonicity,
% for the reason P6 established: these methods are solver-limited at fine grids, and demanding a
% monotone sequence of a quantity that has stopped depending on the grid asks the wrong question.
fprintf('\nthe worst-age skewness error falls over the sweep [T2], this should be one: %i \n',err_skew(end)<err_skew(1))
fprintf('the worst-age excess kurtosis error falls over the sweep [T2], this should be one: %i \n',err_exkurt(end)<err_exkurt(1))
fprintf('   skewness       :');
for c_c=1:nz
    fprintf(' %2.1e',err_skew(c_c));
end
fprintf(' \n   excess kurtosis:');
for c_c=1:nz
    fprintf(' %2.1e',err_exkurt(c_c));
end
fprintf(' \n')


%% REGRESSION BARS ON THE ACCURACY NUMBERS ABOVE
% Every moment error printed in the sweep is reported with a value but no verdict, which means a
% defect can sit in the output while the run summary reports a clean pass. That is exactly what
% happened on 2026-08-28: B30 put a family's mean error at 0.234, flat across a tenfold grid
% refinement, and DiscSummary reported "1 failed" because nothing asserted on it.
%
% These bars are NOT accuracy claims - they are regression bars. A single accuracy threshold cannot
% work across a sweep, because a coarse grid is legitimately bad and a fine one is not; a bar loose
% enough for znum=5 catches nothing at znum=51. So each family gets two: the worst over the whole
% sweep against a loose bar, which catches gross regressions anywhere, and the finest grid against a
% tight one, which is where an accuracy claim can actually be made. Both are set from the run of
% 2026-08-28 with roughly an order of magnitude of headroom, so ordinary drift does not fire them.
fprintf('worst mean error over the whole sweep [T2], this should be below %g: %2.3e \n',0.001,max(err_mean))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',1e-05,err_mean(end))
fprintf('worst variance error over the whole sweep [T2], this should be below %g: %2.3e \n',0.01,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',0.0001,err_var(end))
fprintf('worst skewness error over the whole sweep [T2], this should be below %g: %2.3e \n',1,max(err_skew))
fprintf('skewness error at the finest grid [T2], this should be below %g: %2.3e \n',0.05,err_skew(end))
fprintf('worst excess kurtosis error over the whole sweep [T2], this should be below %g: %2.3e \n',5,max(err_exkurt))
fprintf('excess kurtosis error at the finest grid [T2], this should be below %g: %2.3e \n',0.2,err_exkurt(end))

%% kfttoptions.method='gauss-hermite'
% COVERAGE. The sweep above runs on the default method, which is 'even', and DiscP7_crosstests runs
% 'GMQ' for the B27 regression. That leaves the 'gauss-hermite' branch of the prior in
% discretizeLifeCycleAR1wGM_KFTT as the one arm of its switch that nothing in the bank ever entered.
% It matters now because that line was one of the sites rewritten when the discretization commands
% were taken off the Statistics Toolbox: every other rewritten line has been run, and this one had
% not, so a typo in it would have shipped unnoticed.
%
% WHAT CAN BE ASSERTED, and what cannot. The method does two things here, and only one of them is
% the rewritten line. It picks the GRID - 'even' spaces it evenly over +-nSigmas*sigma_z while
% 'gauss-hermite' puts it on mewz+sqrt(2)*sigmaz*(Gauss-Hermite nodes), so the two grids are
% genuinely different and nothing may assert they agree - and it picks the PRIOR that the maximum
% entropy step starts from, which is the rewritten line. So the structural invariants are the only
% things asserted here. The accuracy numbers below are printed without bars, because on this
% calibration they are not close and are not meant to be: see the finding recorded below.
%
% THE FINDING, from the run of 2026-09-16, which is the first time anything in the bank entered this
% branch. On P7's calibration gauss-hermite is not a viable choice and the numbers are not marginal:
%
%    matching fewer than 4 moments    81.0%    against 18.2% for the default 'even'
%    worst-age excess kurtosis error  8.917    against 0.136 for 'even', so 65 times worse
%    worst-age variance error         2.2e-01  against 1.6e-06 for 'even'
%
% That is the command's own warning being cashed out. It emits "Model is persistent; even-spaced
% grid is recommended" whenever rho(jj)>0.8, and P7's rho runs 0.95 down to 0.90 across the ages
% (DiscSetup_LCAR1GM), so the warning applies at EVERY age, not just some of them. The
% mechanism is the grid, not the prior: Gauss-Hermite nodes crowd the centre and reach far out in
% the tails, which is right for integrating against a single normal and wrong as a support for a
% fat-tailed mixture whose conditional distributions then cannot match four moments on it. The
% entropy solve falls back, and a chain assembled mostly from fallback rows has an excess kurtosis
% that is an artefact rather than a property.
%
% So do NOT read the comparison below as "which method is better" - that is settled, and by a wide
% margin. It is here so the branch is executed and its output shown to be well formed, and so the
% size of the gap is on the record rather than being rediscovered.
znum_gh=znums(min(3,nz));
kg=struct(); kg.verbose=0; kg.method='gauss-hermite';
[zg_gh,pz_gh,jd_gh,oo_gh]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum_gh,J,kg);
ke=struct(); ke.verbose=0; ke.method='even';
[zg_ev,pz_ev,jd_ev,oo_ev]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum_gh,J,ke);
zg_gh=gather(zg_gh); pz_gh=gather(pz_gh); jd_gh=gather(jd_gh);
zg_ev=gather(zg_ev); pz_ev=gather(pz_ev); jd_ev=gather(jd_ev);

fprintf('\nmethod=gauss-hermite, znum=%i: it runs [T0], this should be zero: %i \n',znum_gh,any(size(pz_gh)~=[znum_gh,znum_gh,J-1]))
fprintf('method=gauss-hermite: no NaN or Inf [T0], this should be zero: %i \n',any(~isfinite(zg_gh(:)))+any(~isfinite(pz_gh(:))))
fprintf('method=gauss-hermite: pi_z_J is in [0,1] [T0], this should be zero: %i \n',any(pz_gh(:)<0)+any(pz_gh(:)>1))
fprintf('method=gauss-hermite: rows of pi_z_J sum to one, this should be below %g: %2.8e \n',calib.entropytol,max(abs(sum(pz_gh,2)-1),[],'all'))
fprintf('method=gauss-hermite: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',abs(sum(jd_gh)-1))
fprintf('method=gauss-hermite: the grid is strictly ascending at every age [T0], this should be one: %i \n',all(arrayfun(@(j) issorted(zg_gh(:,j),'strictascend'),1:J)))
fprintf('method=gauss-hermite: conditional distributions matching fewer than 4 moments: %2.1f%% (method=even gives %2.1f%%) \n',100*mean(oo_gh.nMoments_grid(:)<4),100*mean(oo_ev.nMoments_grid(:)<4))

% The option has to have CHANGED something, or it was ignored and every check here is vacuous. The
% grid is the visible half: gauss-hermite nodes are not evenly spaced, so the spacing that is
% constant under method=even must not be constant here.
fprintf('method=gauss-hermite: the grid is NOT evenly spaced, i.e. the option reached the switch [T0], this should be one: %i \n',max(abs(diff(diff(zg_gh(:,1)))))>1e-10)
fprintf('method=gauss-hermite: method=even by contrast IS evenly spaced [T0], this should be one: %i \n',max(abs(diff(diff(zg_ev(:,1)))))<1e-10)

% ...and yet the moments it solved for should be the same, wherever both solves got all four.
[m_gh,v_gh,~,d_gh]=MarkovChainMoments_FHorz(zg_gh,pz_gh,jd_gh);
[m_ev,v_ev,~,d_ev]=MarkovChainMoments_FHorz(zg_ev,pz_ev,jd_ev);
sk_gh=zeros(1,J); ek_gh=zeros(1,J); sk_ev=zeros(1,J); ek_ev=zeros(1,J);
for j_c=1:J
    d=gather(d_gh(:,j_c)); d=d/sum(d); g=zg_gh(:,j_c);
    m1=sum(d.*g); vv=sum(d.*(g-m1).^2);
    sk_gh(j_c)=sum(d.*(g-m1).^3)/vv^1.5; ek_gh(j_c)=sum(d.*(g-m1).^4)/vv^2-3;
    d=gather(d_ev(:,j_c)); d=d/sum(d); g=zg_ev(:,j_c);
    m1=sum(d.*g); vv=sum(d.*(g-m1).^2);
    sk_ev(j_c)=sum(d.*(g-m1).^3)/vv^1.5; ek_ev(j_c)=sum(d.*(g-m1).^4)/vv^2-3;
end
fprintf('method=gauss-hermite vs even: worst age-conditional mean gap [T2] %2.3e, variance gap %2.3e \n',max(abs(m_gh(:)-m_ev(:))),max(abs(v_gh(:)-v_ev(:))))
fprintf('method=gauss-hermite vs even: worst skewness gap [T2] %2.3e, excess kurtosis gap %2.3e \n',max(abs(sk_gh-sk_ev)),max(abs(ek_gh-ek_ev)))
fprintf('   (gaps between two discretizations on DIFFERENT grids, not errors against truth. They are \n')
fprintf('   small only to the extent both solves matched their four moments, which the fallback line \n')
fprintf('   above reports; where either fell back, the gap is the fallback talking.) \n')

% Against the exact age profiles from the cumulant recursion, both methods side by side. No bar:
% the gap is expected to be large for the reason set out at the top of this block, and the sweep
% above is where the accuracy claim for this command is actually made. What these lines are for is
% to keep the size of the gap visible, so that if it ever narrows sharply - which would mean the
% fallback rate had collapsed - somebody notices.
fprintf('method=gauss-hermite vs truth: worst-age mean error [T2] %2.3e (method=even gives %2.3e) \n',max(abs(m_gh(:)'-mewzT)),max(abs(m_ev(:)'-mewzT)))
fprintf('method=gauss-hermite vs truth: worst-age variance error [T2] %2.3e (method=even gives %2.3e) \n',max(abs(v_gh(:)'-varzT)),max(abs(v_ev(:)'-varzT)))
fprintf('method=gauss-hermite vs truth: worst-age skewness error [T2] %2.3e (method=even gives %2.3e) \n',max(abs(sk_gh-skewT)),max(abs(sk_ev-skewT)))
fprintf('method=gauss-hermite vs truth: worst-age excess kurtosis error [T2] %2.3e (method=even gives %2.3e) \n',max(abs(ek_gh-exkurtT)),max(abs(ek_ev-exkurtT)))

%% Figure
figure(figure_c)
subplot(1,3,1)
plot(1:J,skewT,'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,output.sweep(c_c).skew,'--')
end
hold off
xlabel('age j'); ylabel('skewness'); title('KFTT: skewness profile')
legend([{'truth'},arrayfun(@(x) ['znum=',num2str(x)],znums,'UniformOutput',false)],'Location','best')
subplot(1,3,2)
plot(1:J,exkurtT,'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,output.sweep(c_c).exkurt,'--')
end
hold off
xlabel('age j'); ylabel('excess kurtosis'); title('KFTT: excess kurtosis profile')
subplot(1,3,3)
semilogy(znums,max(err_var,10^(-18)),'o-',znums,max(err_skew,10^(-18)),'s-',znums,max(err_exkurt,10^(-18)),'d-')
xlabel('znum'); ylabel('worst-age error'); title('KFTT: convergence')
legend({'variance','skewness','excess kurtosis'},'Location','best')

end
