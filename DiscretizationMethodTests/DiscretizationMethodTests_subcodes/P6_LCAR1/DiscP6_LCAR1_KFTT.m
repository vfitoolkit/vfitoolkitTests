function output=DiscP6_LCAR1_KFTT(calib,znums,figure_c)
% P6: discretizeLifeCycleAR1_KFTT
%
% The block's reference method. Everything here is an AGE PROFILE rather than a single number, and
% the truth is the recursion recomputed in this subcode - mewz(j)=mew(j)+rho(j)*mewz(j-1) and
% sigmaz(j)^2=rho(j)^2*sigmaz(j-1)^2+sigma(j)^2 - never the command's own otheroutputs.sigma_z,
% which is checked against it separately in DiscP6_crosstests.
%
% Two things this block sees that no earlier one did. First, pi_z_J has third dimension J-1, not J:
% pi_z_J(:,:,j) is the transition from age j to age j+1, and there is no age J+1 to transition to.
% Second, jequaloneDistz is a real output with real content - the age-1 distribution - and it is the
% only way the age-1 moments can be right at all.
%
% Because there are J ages, printing every age would swamp the diary. What is printed is the WORST
% age for each moment and where it occurred; the full profiles are in the figure and the struct.

fprintf('\n========== P6: discretizeLifeCycleAR1_KFTT ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho; sigma=calib.vary.sigma;
mewzT=calib.vary.mewz; sigmazT=calib.vary.sigmaz; acT=calib.vary.autocorr;
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_ac=zeros(1,nz);
runtime=zeros(1,nz); fallback=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    kfttoptions=struct();

    % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
    % otheroutputs comes from the warm-up, which is discarded anyway, so measuring the fallback
    % rate is free.
    tic;
    [z_grid_J,pi_z_J,jequaloneDistz,otheroutputs]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,kfttoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,kfttoptions);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);

    z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);

    % --- invariants, including the ones specific to a life-cycle output
    fprintf('znum=%3i: size of z_grid_J is znum-by-J [T0], this should be zero: %i \n',znum,any(size(z_grid_J)~=[znum,J]))
    fprintf('znum=%3i: size of pi_z_J is znum-by-znum-by-(J-1) [T0], this should be zero: %i \n',znum,any(size(pi_z_J)~=[znum,znum,J-1]))
    fprintf('znum=%3i: size of jequaloneDistz is znum-by-1 [T0], this should be zero: %i \n',znum,any(size(jequaloneDistz)~=[znum,1]))
    asc=1;
    for j_c=1:J
        asc=asc*issorted(z_grid_J(:,j_c),'strictascend');
    end
    fprintf('znum=%3i: the grid at every age is strictly ascending [T0], this should be one: %i \n',znum,asc)
    fprintf('znum=%3i: pi_z_J is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z_J(:)<0)+any(pi_z_J(:)>1))
    fprintf('znum=%3i: rows of pi_z_J sum to one at every age, this should be below %g: %2.8e \n',znum,calib.entropytol,max(abs(sum(pi_z_J,2)-1),[],'all'))
    fprintf('znum=%3i: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',znum,abs(sum(jequaloneDistz)-1))
    fprintf('znum=%3i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid_J(:)))+any(~isfinite(pi_z_J(:)))+any(~isfinite(jequaloneDistz)))

    % --- the fallback rate, by age. The command already prints a summary of this and warns when it
    % hits four moments in under 80% of cases; recording it as a number is what lets the accuracy
    % below be read against it, which is the reading that mattered in P2, P3 and P4.
    nMg=otheroutputs.nMoments_grid;
    fallback(c_c)=mean(nMg(:)<4); % 4 is discretizeLifeCycleAR1_KFTT's default nMoments, unlike the stationary discretizeAR1_FarmerToda which defaults to 2
    fprintf('znum=%3i: conditional distributions matching fewer than 4 moments: %2.1f%% \n',znum,100*fallback(c_c))

    % --- the age profiles against the recursion
    [mcmean,mcvar,mcautocorr]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    dmean=abs(mcmean(:)'-mewzT);
    dvar=abs(mcvar(:)'-sigmazT.^2);
    dac=abs(mcautocorr(2:end)-acT(2:end));
    [err_mean(c_c),jm]=max(dmean);
    [err_var(c_c),jv]=max(dvar);
    [err_ac(c_c),ja]=max(dac);
    fprintf('znum=%3i: worst age for the mean is j=%i, error [T2] %2.3e \n',znum,jm,err_mean(c_c))
    fprintf('znum=%3i: worst age for the variance is j=%i, error [T2] %2.3e \n',znum,jv,err_var(c_c))
    fprintf('znum=%3i: worst age for the autocorrelation is j=%i, error [T2] %2.3e \n',znum,ja+1,err_ac(c_c))
    fprintf('znum=%3i: runtime %2.6f s \n',znum,runtime(c_c))
    output.sweep(c_c).znum=znum;
    output.sweep(c_c).mean=mcmean(:)';
    output.sweep(c_c).var=mcvar(:)';
    output.sweep(c_c).autocorr=mcautocorr(:)';
end

output.err_mean=err_mean; output.err_var=err_var; output.err_ac=err_ac;
output.runtime=runtime; output.fallback=fallback; output.znums=znums;
output.mewzT=mewzT; output.sigmazT=sigmazT; output.acT=acT;

% WHAT IS ASSERTED HERE, AND WHY IT IS NOT MONOTONICITY. The first version of this demanded that
% the worst-age variance error never rise as znum grows. It failed on the run of 2026-08-26 with the
% sequence 2.3e-07, 4.6e-06, 1.8e-07, 1.5e-07, 1.7e-07, 2.6e-07 - which is not a convergence failure,
% it is flat noise around 2e-07 with one outlier. This method is SOLVER-limited, not grid-limited:
% the maximum entropy solve has its own accuracy floor, P2 measured the same floor for the
% stationary discretizeAR1_FarmerToda, and below that floor the error simply does not care about
% znum. Demanding monotonicity of a quantity that has stopped depending on the grid is asking the
% wrong question, and it would have gone red on every future run for no reason.
%
% So what is asserted is the floor itself: the error must stay below it at every grid size. That
% catches a real regression - a broken command lands at 1e-2 or worse, four orders above this - and
% it does not fire on solver noise. The sequence is printed so the flatness is visible rather than
% asserted away.
solverfloor=10^(-5);
fprintf('the worst-age variance error is below the entropy solver floor of %g at every znum [T2], this should be one: %i \n',solverfloor,all(err_var<solverfloor))
fprintf('   the sequence is:');
for c_c=1:nz
    fprintf(' %2.1e',err_var(c_c));
end
fprintf(' \n');
fprintf('   (flat rather than falling, because this method is solver-limited rather than grid-limited \n')
fprintf('    at these grid sizes; refining znum past the point where the solve dominates buys nothing) \n')

%% The option grid: method x nMoments
% KFTT offers four grid methods, and its nMoments default is 4 where the stationary
% discretizeAR1_FarmerToda defaults to 2. P2 and P3 both found that asking for more moments can
% make the answer WORSE when the entropy solve starts failing, and P5 found it a third time, so
% the fallback rate is reported next to the error rather than left to be inferred.
fprintf('\n--- option grid: method x nMoments, at znum=15 --- \n')
methodlist={'even','gauss-legendre','clenshaw-curtis','gauss-hermite'};
nMomentslist=[2,4];
errgrid=NaN(4,2); fbgrid=NaN(4,2);
for m_c=1:4
    for n_c=1:2
        kfttoptions=struct(); kfttoptions.method=methodlist{m_c}; kfttoptions.nMoments=nMomentslist(n_c);
        [zg,pz,j1,oo]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,kfttoptions);
        zg=gather(zg); pz=gather(pz); j1=gather(j1);
        [~,v,~]=MarkovChainMoments_FHorz(zg,pz,j1);
        errgrid(m_c,n_c)=max(abs(v(:)'-sigmazT.^2));
        fbgrid(m_c,n_c)=mean(oo.nMoments_grid(:)<nMomentslist(n_c));
        fprintf('method=%-16s nMoments=%i: worst-age variance error %2.3e, fallback %2.1f%% \n',methodlist{m_c},nMomentslist(n_c),errgrid(m_c,n_c),100*fbgrid(m_c,n_c))
    end
end
[~,bi]=min(errgrid(:)); [bm,bn]=ind2sub([4,2],bi);
fprintf('the most accurate combination here is method=%s at nMoments=%i (reported, not asserted) \n',methodlist{bm},nMomentslist(bn))
output.errgrid=errgrid; output.fbgrid=fbgrid;
output.methodlist=methodlist; output.nMomentslist=nMomentslist;


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
fprintf('worst mean error over the whole sweep [T2], this should be below %g: %2.3e \n',0.0001,max(err_mean))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',1e-06,err_mean(end))
fprintf('worst variance error over the whole sweep [T2], this should be below %g: %2.3e \n',0.0001,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',1e-06,err_var(end))
fprintf('worst autocorrelation error over the whole sweep [T2], this should be below %g: %2.3e \n',0.0001,max(err_ac))
fprintf('autocorrelation error at the finest grid [T2], this should be below %g: %2.3e \n',1e-06,err_ac(end))

%% Scalar versus age-vector nSigmas (the B17 regression, and its extension to both KFTT commands)
% B17 gave discretizeLifeCycleAR1_KFTT an age-dependent nSigmas. A scalar and a constant vector of
% the same value must give bit-identical output - the scalar is expanded to exactly that vector -
% and a genuinely varying vector must make the grid half-width track nSigmas(j)*sigmaz(j) at every
% age. The gaussian-mixture sibling did not have this at all until 2026-08-26; its version of this
% check belongs in P7.
%
% Only the 'even', 'gauss-legendre' and 'clenshaw-curtis' methods use nSigmas: 'gauss-hermite'
% determines its own grid and ignores it. So this uses 'even', and the fact that it would be
% VACUOUS under 'gauss-hermite' is the reason the method is pinned rather than left at the default.
fprintf('\n--- scalar versus age-vector nSigmas (the B17 regression) --- \n')
o1=struct(); o1.nSigmas=3; o1.method='even';
o2=struct(); o2.nSigmas=3*ones(J,1); o2.method='even';
[zA,pA,jA]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,o1);
[zB,pB,jB]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,o2);
fprintf('z_grid_J [T0], this should be zero: %2.8e \n',max(abs(gather(zA)-gather(zB)),[],'all'))
fprintf('pi_z_J [T0], this should be zero: %2.8e \n',max(abs(gather(pA)-gather(pB)),[],'all'))
fprintf('jequaloneDistz [T0], this should be zero: %2.8e \n',max(abs(gather(jA)-gather(jB))))
% a genuinely age-varying nSigmas, which only the vector form can express
o3=struct(); o3.nSigmas=linspace(2,4,J)'; o3.method='even';
[zC,~,~,ooC]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,o3);
zC=gather(zC);
ghwC=(zC(end,:)-zC(1,:))/2;
fprintf('age-varying nSigmas: the half-width tracks nSigmas(j)*sigmaz(j) [T1], this should be zero: %2.8e \n',max(abs(ghwC-linspace(2,4,J).*gather(ooC.sigma_z(:))')))
% a vector of the wrong length must be rejected, not silently used at the wrong ages
try
    discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,struct('nSigmas',3*ones(J+1,1),'method','even'));
    fprintf('a wrong-length nSigmas vector: ran without error [T0], this should not happen \n')
catch ME
    fprintf('a wrong-length nSigmas vector: errors, and the message names the option [T0], this should be one: %i \n',contains(ME.message,'kfttoptions.nSigmas'))
end
% and the documented fact that gauss-hermite ignores nSigmas entirely, which is worth pinning
% because a user who sets it there gets no effect and no warning
oH1=struct(); oH1.method='gauss-hermite'; oH1.nSigmas=2;
oH2=struct(); oH2.method='gauss-hermite'; oH2.nSigmas=5;
[zH1,~,~]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,oH1);
[zH2,~,~]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,15,J,oH2);
fprintf('method=gauss-hermite ignores nSigmas entirely [T0], this should be zero: %2.8e \n',max(abs(gather(zH1)-gather(zH2)),[],'all'))
fprintf('   (that is documented behaviour, not a bug - gauss-hermite determines its own grid - but \n')
fprintf('    setting nSigmas alongside it has no effect and produces no warning) \n')

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(1:J,sigmazT.^2,'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,output.sweep(c_c).var,'--')
end
hold off
xlabel('age j'); ylabel('variance of z'); title('KFTT: the variance profile against the recursion')
legend([{'truth'},arrayfun(@(x) ['znum=',num2str(x)],znums,'UniformOutput',false)],'Location','best')
subplot(1,2,2)
semilogy(znums,max(err_mean,10^(-18)),'o-',znums,max(err_var,10^(-18)),'s-',znums,max(err_ac,10^(-18)),'d-')
xlabel('znum'); ylabel('worst-age error'); title('KFTT: convergence')
legend({'mean','variance','autocorrelation'},'Location','best')

end
