function output=DiscP6_LCAR1_FGPTauchen(calib,znums,znumsEven,figure_c)
% P6: discretizeLifeCycleAR1_FellaGallipoliPanTauchen (the extended-Tauchen method)
%
% Unlike its Rouwenhorst sibling this one does take mew, and supports all four initial-condition
% options, so it runs the full age-varying calibration.
%
% Its nSigmas defaults to min(sqrt(znum-1),3), and that cap is NOT the same issue as the one in
% DiscP6_LCAR1_FGP. For Tauchen, nSigmas is a free hyperparameter - the method has no exactness
% property that a particular width would deliver - so capping it is a choice about where to
% truncate, not a broken construction. What it does mean is the familiar Tauchen truncation floor:
% at fixed width, refining znum stops helping once the truncated tail dominates. P2 measured that
% for the stationary method and it is measured here too, by sweeping nSigmas at fixed znum.

fprintf('\n========== P6: discretizeLifeCycleAR1_FellaGallipoliPanTauchen ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho; sigma=calib.vary.sigma;
mewzT=calib.vary.mewz; sigmazT=calib.vary.sigmaz; acT=calib.vary.autocorr;
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_ac=zeros(1,nz); runtime=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    fgpoptions=struct();

    tic;
    [z_grid_J,pi_z_J,jequaloneDistz,otheroutputs]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,fgpoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,fgpoptions);
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

    % --- the grid is centred on the recursion's mean and scaled by its standard deviation. That
    % is an identity for a Tauchen grid, so it catches a wrong mewz or sigmaz before any moment is
    % computed - the same check that caught B19 in P5.
    gmid=(z_grid_J(end,:)+z_grid_J(1,:))/2;
    ghw=(z_grid_J(end,:)-z_grid_J(1,:))/2;
    nS=min(sqrt(znum-1),3);
    fprintf('znum=%3i: the grid is centred on mewz(j) at every age [T1], this should be zero: %2.8e \n',znum,max(abs(gmid-mewzT)))
    fprintf('znum=%3i: the grid half-width is nSigmas*sigmaz(j) at every age [T1], this should be zero: %2.8e \n',znum,max(abs(ghw-nS*sigmazT)))

    % --- the age profiles
    [mcmean,mcvar,mcautocorr]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    [err_mean(c_c),jm]=max(abs(mcmean(:)'-mewzT));
    [err_var(c_c),jv]=max(abs(mcvar(:)'-sigmazT.^2));
    [err_ac(c_c),ja]=max(abs(mcautocorr(2:end)-acT(2:end)));
    fprintf('znum=%3i: worst age for the mean is j=%i, error [T2] %2.3e \n',znum,jm,err_mean(c_c))
    fprintf('znum=%3i: worst age for the variance is j=%i, error [T2] %2.3e \n',znum,jv,err_var(c_c))
    fprintf('znum=%3i: worst age for the autocorrelation is j=%i, error [T2] %2.3e \n',znum,ja+1,err_ac(c_c))
    fprintf('znum=%3i: runtime %2.6f s \n',znum,runtime(c_c))
    output.sweep(c_c).var=mcvar(:)';
end

output.err_mean=err_mean; output.err_var=err_var; output.err_ac=err_ac;
output.runtime=runtime; output.znums=znums;

%% The truncation floor: sweep nSigmas at fixed znum
% At a FIXED nSigmas, refining znum converges to a truncation floor rather than to zero - P2
% established that for the stationary Tauchen method and pinned the convergence ORDER instead. The
% same thing shows here as a U shape in nSigmas at fixed znum: too narrow truncates the tails, too
% wide coarsens the spacing, and the error is minimised somewhere in between.
fprintf('\n--- nSigmas at fixed znum=31 --- \n')
nSlist=[2,2.5,3,3.5,4,5];
err_nS=zeros(1,length(nSlist));
for s_c=1:length(nSlist)
    fgpoptions=struct(); fgpoptions.nSigmas=nSlist(s_c);
    [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,31,J,fgpoptions);
    zg=gather(zg); pz=gather(pz); j1=gather(j1);
    [~,v,~]=MarkovChainMoments_FHorz(zg,pz,j1);
    err_nS(s_c)=max(abs(v(:)'-sigmazT.^2));
    fprintf('nSigmas=%3.1f: worst-age variance error %2.3e \n',nSlist(s_c),err_nS(s_c))
end
[~,bs]=min(err_nS);
fprintf('the error is minimised at nSigmas=%3.1f, and rises on both sides of it [T2], this should be one: %i \n',nSlist(bs),(bs>1)&&(bs<length(nSlist)))
fprintf('   (that U shape is the point: it is why refining znum at a fixed nSigmas stops helping, \n')
fprintf('    and why the default cap of 3 is a truncation choice rather than an error) \n')
output.err_nS=err_nS; output.nSlist=nSlist;


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
fprintf('worst variance error over the whole sweep [T2], this should be below %g: %2.3e \n',0.5,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',0.01,err_var(end))
fprintf('worst autocorrelation error over the whole sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_ac))
fprintf('autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',0.005,err_ac(end))

%% Scalar versus age-vector nSigmas (the B17 regression)
% B17 was age-dependent nSigmas. A scalar and a constant vector of the same value must give
% bit-identical output; anything else means the scalar is not being expanded the way it claims.
fprintf('\n--- scalar versus age-vector nSigmas (the B17 regression) --- \n')
o1=struct(); o1.nSigmas=3;
o2=struct(); o2.nSigmas=3*ones(J,1);
[zA,pA,jA]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,15,J,o1);
[zB,pB,jB]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,15,J,o2);
fprintf('z_grid_J [T0], this should be zero: %2.8e \n',max(abs(gather(zA)-gather(zB)),[],'all'))
fprintf('pi_z_J [T0], this should be zero: %2.8e \n',max(abs(gather(pA)-gather(pB)),[],'all'))
fprintf('jequaloneDistz [T0], this should be zero: %2.8e \n',max(abs(gather(jA)-gather(jB))))
% and a genuinely age-varying nSigmas, which only the vector form can express
o3=struct(); o3.nSigmas=linspace(2,4,J)';
[zC,pC,jC]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,15,J,o3);
zC=gather(zC);
ghwC=(zC(end,:)-zC(1,:))/2;
fprintf('age-varying nSigmas: the half-width tracks nSigmas(j)*sigmaz(j) [T1], this should be zero: %2.8e \n',max(abs(ghwC-linspace(2,4,J).*sigmazT)))

%% Even znum (B15)
fprintf('\n--- even znum (the B15 regression) --- \n')
for c_c=1:length(znumsEven)
    znum=znumsEven(c_c);
    [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,struct());
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
xlabel('age j'); ylabel('variance of z'); title('FGP-Tauchen: the variance profile')
subplot(1,2,2)
semilogy(nSlist,max(err_nS,10^(-18)),'o-')
xlabel('nSigmas'); ylabel('worst-age variance error'); title('FGP-Tauchen: the truncation floor, znum=31')

end
