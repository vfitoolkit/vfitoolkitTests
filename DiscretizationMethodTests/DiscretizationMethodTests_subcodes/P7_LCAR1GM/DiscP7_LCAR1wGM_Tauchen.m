function output=DiscP7_LCAR1wGM_Tauchen(calib,znums,figure_c)
% P7: discretizeLifeCycleAR1wGM_Tauchen
%
% Written test-first alongside this block, as discretizeAR1wGM_Tauchen was for P3 and
% discretizeAR1wSV_Tauchen for P4. It is the Tauchen counterpart of the KFTT command: the
% age-dependent extension of discretizeAR1wGM_Tauchen.
%
% TWO DIFFERENCES FROM ITS KFTT SIBLING that this subcode has to respect.
%
% It takes Tauchen_q rather than nSigmas, in the argument position before the options struct, and
% it has no maximum entropy step - so there is no nMoments, no nMoments_grid and no fallback rate.
% What a bad grid costs it shows up directly in the moments instead of partly in a solver.
%
% It centres its grid on E(z_j) INCLUDING the mean of the mixture, where the KFTT command uses a
% recursion with no E(e) term. This block's calibration is deliberately NOT mean zero, so the two
% conventions genuinely differ here and that difference is measured in DiscP7_crosstests rather
% than hidden by a calibration that cannot see it.

fprintf('\n========== P7: discretizeLifeCycleAR1wGM_Tauchen ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho;
p=calib.vary.mixprobs_i; mu=calib.vary.mu_i; sd=calib.vary.sigma_i;
mewzT=calib.vary.mewz; varzT=calib.vary.varz; skewT=calib.vary.skewz; exkurtT=calib.vary.exkurtz;
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_skew=zeros(1,nz); err_exkurt=zeros(1,nz);
runtime=zeros(1,nz);
Tauchen_q=3;

for c_c=1:nz
    znum=znums(c_c);
    to=struct(); to.verbose=0;

    tic;
    [z_grid_J,pi_z_J,jequaloneDistz,oo]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,Tauchen_q,to);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,Tauchen_q,to);
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
    fprintf('znum=%3i: rows of pi_z_J sum to one at every age [T1], this should be zero: %2.8e \n',znum,max(abs(sum(pi_z_J,2)-1),[],'all'))
    fprintf('znum=%3i: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',znum,abs(sum(jequaloneDistz)-1))
    fprintf('znum=%3i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid_J(:)))+any(~isfinite(pi_z_J(:))))

    % --- the grid is centred on E(z_j) and scaled by sd(z_j). For a Tauchen grid that is an
    % identity, so it catches a wrong mewz or sigmaz before any moment is computed - the same check
    % that caught B19 in P5 and that pins this command's centring convention.
    gmid=(z_grid_J(end,:)+z_grid_J(1,:))/2;
    ghw=(z_grid_J(end,:)-z_grid_J(1,:))/2;
    fprintf('znum=%3i: the grid is centred on E(z_j), mixture mean included [T1], this should be zero: %2.8e \n',znum,max(abs(gmid-mewzT)))
    fprintf('znum=%3i: the grid half-width is Tauchen_q*sd(z_j) [T1], this should be zero: %2.8e \n',znum,max(abs(ghw-Tauchen_q*sqrt(varzT))))
    fprintf('znum=%3i: otheroutputs.mew_z and sigma_z match the recursion [T1], these should be zero: %2.8e %2.8e \n',znum,max(abs(gather(oo.mew_z)-mewzT)),max(abs(gather(oo.sigma_z)-sqrt(varzT))))

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
output.runtime=runtime; output.znums=znums; output.Tauchen_q=Tauchen_q;

fprintf('\nthe worst-age skewness error falls over the sweep [T2], this should be one: %i \n',err_skew(end)<err_skew(1))
fprintf('the worst-age excess kurtosis error falls over the sweep [T2], this should be one: %i \n',err_exkurt(end)<err_exkurt(1))
fprintf('   (at a FIXED Tauchen_q of %g the grid never widens, so these converge to a truncation \n',Tauchen_q)
fprintf('    floor rather than to zero - the same behaviour P2 documents for the gaussian Tauchen \n')
fprintf('    method, and DiscP7_gridwidth is where the width is varied) \n')


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
fprintf('worst mean error over the whole sweep [T2], this should be below %g: %2.3e \n',1,max(err_mean))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',0.02,err_mean(end))
fprintf('worst variance error over the whole sweep [T2], this should be below %g: %2.3e \n',0.5,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',0.05,err_var(end))
fprintf('worst skewness error over the whole sweep [T2], this should be below %g: %2.3e \n',1,max(err_skew))
fprintf('skewness error at the finest grid [T2], this should be below %g: %2.3e \n',0.4,err_skew(end))
fprintf('worst excess kurtosis error over the whole sweep [T2], this should be below %g: %2.3e \n',3,max(err_exkurt))
fprintf('excess kurtosis error at the finest grid [T2], this should be below %g: %2.3e \n',2.5,err_exkurt(end))

%% Figure
figure(figure_c)
subplot(1,3,1)
plot(1:J,skewT,'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,output.sweep(c_c).skew,'--')
end
hold off
xlabel('age j'); ylabel('skewness'); title('wGM\_Tauchen: skewness profile')
legend([{'truth'},arrayfun(@(x) ['znum=',num2str(x)],znums,'UniformOutput',false)],'Location','best')
subplot(1,3,2)
plot(1:J,exkurtT,'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,output.sweep(c_c).exkurt,'--')
end
hold off
xlabel('age j'); ylabel('excess kurtosis'); title('wGM\_Tauchen: excess kurtosis profile')
subplot(1,3,3)
semilogy(znums,max(err_var,10^(-18)),'o-',znums,max(err_skew,10^(-18)),'s-',znums,max(err_exkurt,10^(-18)),'d-')
xlabel('znum'); ylabel('worst-age error'); title('wGM\_Tauchen: convergence')
legend({'variance','skewness','excess kurtosis'},'Location','best')

end
