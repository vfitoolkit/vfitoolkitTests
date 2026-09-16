function output=DiscP1_IID_TanakaToda(calib,znums,figure_c)
% P1: discretizeIID_TanakaToda, run at DEFAULTS only
%
% This command is not swept here. It is currently an exact copy of discretizeIIDNormal_TanakaToda,
% which owns the option sweep, and the plan is for discretizeIID_TanakaToda to be generalized to
% discretize a wide range of distributions while the IIDNormal one stays normal-only.
%
% So what this subcode does is run it at defaults across the sweep and record the outputs; the
% actual comparison against its twin is the T0 check in DiscP1_crosstests. That check is worth
% keeping after the generalization lands: it becomes the statement that the general command still
% reproduces the normal-only one on the normal case, which is exactly the regression anyone
% generalizing them would want.

fprintf('\n========== P1: discretizeIID_TanakaToda (defaults only) ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;
% (this command has no Tauchen_q input)
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); runtime=zeros(1,nz); nrepsused=zeros(1,nz); nMomentsused=zeros(1,nz);

for c_c=1:nz
    enum=znums(c_c);

    % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
    % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
    tic;
    [e_grid,pi_e,otheroutputs]=discretizeIID_TanakaToda(mew,sigma,enum,struct());
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [e_grid,pi_e]=discretizeIID_TanakaToda(mew,sigma,enum,struct());
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    nrepsused(c_c)=nreps;
    nMomentsused(c_c)=otheroutputs.nMoments;
    fprintf('enum=%i: the entropy solve matched %i moments of the 2 requested [T2], this should be one: %i \n',enum,otheroutputs.nMoments,otheroutputs.nMoments==2)

    % --- invariants
    fprintf('enum=%i: size of e_grid [T0], this should be zero: %i \n',enum,any(size(e_grid)~=[enum,1]))
    fprintf('enum=%i: size of pi_e [T0], this should be zero: %i \n',enum,any(size(pi_e)~=[enum,1]))
    fprintf('enum=%i: e_grid is strictly ascending [T0], this should be one: %i \n',enum,issorted(e_grid,'strictascend'))
    fprintf('enum=%i: pi_e is in [0,1] [T0], this should be zero: %i \n',enum,any(pi_e<0)+any(pi_e>1))
    fprintf('enum=%i: pi_e sums to one [T1], this should be zero: %2.8e \n',enum,abs(sum(pi_e)-1))
    fprintf('enum=%i: no NaN or Inf [T0], this should be zero: %i \n',enum,any(~isfinite(e_grid))+any(~isfinite(pi_e)))
    fprintf('enum=%i: e_grid symmetric about mew [T1], this should be zero: %2.8e \n',enum,max(abs((e_grid+flipud(e_grid))/2-mew)))
    fprintf('enum=%i: pi_e equals its own reversal [T1], this should be zero: %2.8e \n',enum,max(abs(pi_e-flipud(pi_e))))

    m=sum(pi_e.*e_grid);
    v=sum(pi_e.*(e_grid-m).^2);
    err_mean(c_c)=abs(m-mew);
    err_var(c_c)=abs(v-sigma^2);
    fprintf('enum=%i: mean error %2.3e, variance error %2.3e, runtime %2.6f s (nreps=%i) \n',enum,err_mean(c_c),err_var(c_c),runtime(c_c),nreps)

    output.sweep(c_c).enum=enum;
    output.sweep(c_c).e_grid=e_grid;
    output.sweep(c_c).pi_e=pi_e;
end
output.err_mean=err_mean;
output.err_var=err_var;
output.runtime=runtime;
output.nrepsused=nrepsused;
output.nMomentsused=nMomentsused;
output.znums=znums;

%% Figure
figure(figure_c)
subplot(1,2,1); plot(znums,max(err_mean,10^(-18)),'o-',znums,max(err_var,10^(-18)),'s-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('discretizeIID_TanakaToda: |error|','Interpreter','none'); xlabel('enum'); legend('mean','variance','Location','best')
subplot(1,2,2); plot(znums,runtime,'o-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('enum')

end
