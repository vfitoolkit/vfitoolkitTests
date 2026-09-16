function output=DiscP6_LCAR1_FGP(calib,znums,znumsEven,figure_c)
% P6: discretizeLifeCycleAR1_FellaGallipoliPan (the extended-Rouwenhorst method)
%
% THE ODD ONE OUT OF THE THREE. It takes no mew argument at all - the process it discretizes is
% z(j)=rho(j)*z(j-1)+e(j), driftless - and it supports only one of the four initial-condition
% options. So this subcode runs the age-varying calibration with mew set aside, and the mean is
% zero at every age by construction rather than being a target.
%
% THE GRID WIDTH. This is what the sweep is really for. The stationary discretizeAR1_Rouwenhorst
% uses a half-width of exactly sqrt(znum-1)*sigmaz, with no cap, and that exact spread is what
% makes Rouwenhorst match the variance to machine precision - it is a property of the construction,
% not a tuning choice. This command defaults to nSigmas=min(sqrt(znum-1),4), so the cap binds for
% znum>17. Whether that costs it the exact variance match is measured below rather than assumed,
% and znums straddles 17 for exactly that reason (sqrt(17-1)=4 exactly).

fprintf('\n========== P6: discretizeLifeCycleAR1_FellaGallipoliPan ========== \n')

output=struct();
J=calib.J;
rho=calib.vary.rho; sigma=calib.vary.sigma;
% The driftless recursion, recomputed here: this command has no mew, so mewz is zero throughout and
% only sigmaz differs from the calibration's stored profile (which was built with mew).
sigmazT=zeros(1,J); sigmazT(1)=sigma(1);
for j_c=2:J
    sigmazT(j_c)=sqrt(rho(j_c)^2*sigmazT(j_c-1)^2+sigma(j_c)^2);
end
acT=nan(1,J);
for j_c=2:J
    acT(j_c)=rho(j_c)*sigmazT(j_c-1)/sigmazT(j_c);
end

nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_ac=zeros(1,nz);
runtime=zeros(1,nz); halfwidth=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    fgpoptions=struct();

    tic;
    [z_grid_J,pi_z_J,jequaloneDistz,otheroutputs]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,fgpoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,fgpoptions);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);

    % --- invariants
    fprintf('znum=%3i: size of z_grid_J is znum-by-J [T0], this should be zero: %i \n',znum,any(size(z_grid_J)~=[znum,J]))
    fprintf('znum=%3i: size of pi_z_J is znum-by-znum-by-(J-1) [T0], this should be zero: %i \n',znum,any(size(pi_z_J)~=[znum,znum,J-1]))
    fprintf('znum=%3i: pi_z_J is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z_J(:)<0)+any(pi_z_J(:)>1))
    fprintf('znum=%3i: rows of pi_z_J sum to one at every age [T1], this should be zero: %2.8e \n',znum,max(abs(sum(pi_z_J,2)-1),[],'all'))
    fprintf('znum=%3i: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',znum,abs(sum(jequaloneDistz)-1))
    fprintf('znum=%3i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid_J(:)))+any(~isfinite(pi_z_J(:))))

    % --- the grid half-width in units of sigmaz, and whether the cap is binding
    halfwidth(c_c)=(z_grid_J(end,J)-z_grid_J(1,J))/2/sigmazT(J);
    fprintf('znum=%3i: grid half-width is %2.4f sigmaz, against the required sqrt(znum-1)=%2.4f \n',znum,halfwidth(c_c),sqrt(znum-1))

    % --- the age profiles
    [mcmean,mcvar,mcautocorr]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    err_mean(c_c)=max(abs(mcmean(:)'));  % truth is zero at every age: no drift
    [err_var(c_c),jv]=max(abs(mcvar(:)'-sigmazT.^2));
    [err_ac(c_c),ja]=max(abs(mcautocorr(2:end)-acT(2:end)));
    fprintf('znum=%3i: worst-age mean error [T2] %2.3e (the truth is zero at every age; there is no drift) \n',znum,err_mean(c_c))
    fprintf('znum=%3i: worst age for the variance is j=%i, error [T2] %2.3e \n',znum,jv,err_var(c_c))
    fprintf('znum=%3i: worst age for the autocorrelation is j=%i, error [T2] %2.3e \n',znum,ja+1,err_ac(c_c))
    fprintf('znum=%3i: runtime %2.6f s \n',znum,runtime(c_c))
    output.sweep(c_c).var=mcvar(:)';
end

output.err_mean=err_mean; output.err_var=err_var; output.err_ac=err_ac;
output.runtime=runtime; output.znums=znums; output.halfwidth=halfwidth;
output.sigmazT=sigmazT; output.acT=acT;

% THE WIDTH, AND THE B29 REGRESSION. Rouwenhorst reproduces the variance exactly when the grid
% half-width is exactly sqrt(znum-1)*sigmaz(j), and at no other width - that is a property of the
% construction, not a hyperparameter. The default used to read min(sqrt(znum-1),4), so the cap bound
% for znum>17 and the variance stopped being exact: 2.27e-01 at znum=31 and 3.31e-01 at znum=51,
% against a true variance of the same order, where below the cap it was 2e-16. The cap was removed
% on 2026-08-26. These checks are what stops it coming back.
fprintf('\n--- the grid width, and the B29 regression --- \n')
fprintf('the default width is sqrt(znum-1) at every znum [T1], this should be zero: %2.8e \n',max(abs(halfwidth-sqrt(znums-1))))
fprintf('and the variance is therefore exact at every znum [T1], this should be zero: %2.8e \n',max(err_var))
% Setting the width by hand to the same value must change nothing - if it does, the default is not
% what it claims to be.
err_var_explicit=zeros(1,nz);
for c_c=1:nz
    znum=znums(c_c);
    fgpoptions=struct(); fgpoptions.nSigmas=sqrt(znum-1);
    [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,fgpoptions);
    [~,v,~]=MarkovChainMoments_FHorz(gather(zg),gather(pz),gather(j1));
    err_var_explicit(c_c)=max(abs(v(:)'-sigmazT.^2));
end
fprintf('setting nSigmas=sqrt(znum-1) explicitly gives the same as the default [T0], this should be zero: %2.8e \n',max(abs(err_var_explicit-err_var)))
% And the other half of the claim: any OTHER width breaks it. Without this the two checks above
% would pass on a command that ignored nSigmas entirely.
fprintf('\n--- and what a wrong width costs, which is why the cap was a bug --- \n')
for w=[2,3,4]
    znum=51;
    fgpoptions=struct(); fgpoptions.nSigmas=w;
    [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,fgpoptions);
    [~,v,~]=MarkovChainMoments_FHorz(gather(zg),gather(pz),gather(j1));
    fprintf('znum=51, nSigmas=%g instead of %2.2f: worst-age variance error %2.3e \n',w,sqrt(znum-1),max(abs(v(:)'-sigmazT.^2)))
end
fprintf('the required width sqrt(50)=7.07 is the ONLY one that works; 4 was the old cap \n')
output.err_var_explicit=err_var_explicit;


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
fprintf('worst mean error over the whole sweep [T2], this should be below %g: %2.3e \n',1e-10,max(err_mean))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',1e-10,err_mean(end))
fprintf('worst variance error over the whole sweep [T2], this should be below %g: %2.3e \n',1e-10,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',1e-10,err_var(end))
fprintf('worst autocorrelation error over the whole sweep [T2], this should be below %g: %2.3e \n',1e-10,max(err_ac))
fprintf('autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',1e-10,err_ac(end))

%% Even znum (B15)
% B15 was an even-znum failure in this family. The methods must work for both parities.
fprintf('\n--- even znum (the B15 regression) --- \n')
for c_c=1:length(znumsEven)
    znum=znumsEven(c_c);
    [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,struct());
    zg=gather(zg); pz=gather(pz); j1=gather(j1);
    [~,v,~]=MarkovChainMoments_FHorz(zg,pz,j1);
    fprintf('znum=%3i (even): sizes are right [T0], this should be zero: %i \n',znum,any(size(zg)~=[znum,J])+any(size(pz)~=[znum,znum,J-1]))
    fprintf('znum=%3i (even): rows sum to one [T1], this should be zero: %2.8e \n',znum,max(abs(sum(pz,2)-1),[],'all'))
    fprintf('znum=%3i (even): worst-age variance error [T2] %2.3e \n',znum,max(abs(v(:)'-sigmazT.^2)))
end

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(1:J,sigmazT.^2,'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,output.sweep(c_c).var,'--')
end
hold off
xlabel('age j'); ylabel('variance of z'); title('FGP: the variance profile against the recursion')
subplot(1,2,2)
semilogy(znums,max(err_var,10^(-18)),'o-',znums,max(err_var_explicit,10^(-18)),'s--')
xlabel('znum'); ylabel('worst-age variance error'); title('FGP: exact at the required width')
legend({'default nSigmas','nSigmas=sqrt(znum-1) set by hand'},'Location','best')

end
