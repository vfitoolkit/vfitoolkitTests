function output=DiscP1_IIDNormal_TanakaToda(calib,znums,figure_c)
% P1: discretizeIIDNormal_TanakaToda
%
% P1 owns the option sweep for this command. Same template as DiscP1_IIDNormal_Tauchen.
%
% The important difference from Tauchen: Tanaka-Toda is a maximum-entropy method that TARGETS the
% moments, so with nMoments>=2 it matches the mean and variance of the iid normal EXACTLY rather
% than approximately. Those checks are therefore T1, not T2 - and that sharper statement is the
% single most valuable thing this subcode asserts. With nMoments=4 it also matches the third and
% fourth central moments (0 and 3*sigma^4) exactly.
%
% Caveat on those checks: the method falls back to fewer moments when the entropy solve fails. It
% used to give no record of that, so a failure here could only be read together with whatever
% warning happened to print above it. As of B25 it returns otheroutputs.nMoments - a scalar, not a
% grid, since there is no loop here - and the sweep below asserts on it directly.
%
% "Exactly" means "to the entropy solver's tolerance", and that tolerance DEGRADES as the problem
% gets harder - with enum, and with nMoments. Measured: the variance match runs 3e-17 at enum=5 but
% 1.3e-09 at enum=101, and at nMoments=4 it is around 5e-11 by enum=15. So the threshold here is
% 1e-7, not the 1e-10 of a T1 check, and the degradation itself is reported: it is a real property
% of the method and worth seeing in the diary rather than hidden behind a pass.
entropytol=10^(-7);

fprintf('\n========== P1: discretizeIIDNormal_TanakaToda ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;

%% Main sweep: enum, at the default options
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); runtime=zeros(1,nz); nrepsused=zeros(1,nz); nMomentsused=zeros(1,nz);

for c_c=1:nz
    enum=znums(c_c);
    tanakatodaoptions=struct();

    % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
    % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
    tic;
    [e_grid,pi_e,otheroutputs]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,tanakatodaoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [e_grid,pi_e]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,tanakatodaoptions);
        if isfield(tanakatodaoptions,'parallel')
            if tanakatodaoptions.parallel==2
                wait(gpuDevice) % otherwise the timing measures kernel launch, not kernel execution
            end
        end
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
    fprintf('enum=%i: pi_e equals its own reversal, this should be below %g: %2.8e \n',enum,entropytol,max(abs(pi_e-flipud(pi_e))))

    % --- the moment-matching guarantee. T1, not T2: with nMoments=2 (the default) this method
    % targets these two moments and should hit them to solver tolerance at every enum.
    m=sum(pi_e.*e_grid);
    v=sum(pi_e.*(e_grid-m).^2);
    err_mean(c_c)=abs(m-mew);
    err_var(c_c)=abs(v-sigma^2);
    fprintf('enum=%i: mean matched to the entropy solver tolerance, this should be below %g: %2.8e \n',enum,entropytol,err_mean(c_c))
    fprintf('enum=%i: variance matched to the entropy solver tolerance, this should be below %g: %2.8e \n',enum,entropytol,err_var(c_c))
    fprintf('enum=%i: runtime %2.6f s (nreps=%i) \n',enum,runtime(c_c),nreps)

    output.sweep(c_c).enum=enum;
    output.sweep(c_c).e_grid=e_grid;
    output.sweep(c_c).pi_e=pi_e;
    output.sweep(c_c).moments=[m,v];
end
output.err_mean=err_mean;
output.err_var=err_var;
output.runtime=runtime;
output.nrepsused=nrepsused;
output.nMomentsused=nMomentsused;
output.znums=znums;
% Report the degradation explicitly rather than letting it hide inside a pass
fprintf('the entropy solve gets harder as the grid grows: variance match is %2.2e at enum=%i and %2.2e at enum=%i \n',err_var(1),znums(1),err_var(end),znums(end))

%% Option sweep at a fixed enum: method x nMoments
enum=15;
methodlist={'even','gauss-legendre','clenshaw-curtis','gauss-hermite'};
nMomentslist=[1,2,4];
err4=zeros(length(methodlist),length(nMomentslist)); % error on the 4th central moment
for m_c=1:length(methodlist)
    for n_c=1:length(nMomentslist)
        opts=struct();
        opts.method=methodlist{m_c};
        opts.nMoments=nMomentslist(n_c);
        opts.verbose=1; % so that any fallback to fewer moments is visible in the diary
        [e_grid,pi_e]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,opts);

        % invariants must hold for every method, including the two that build their grids by
        % quadrature and the one that flips them (clenshaw-curtis)
        fprintf('method=%s, nMoments=%i: e_grid is strictly ascending [T0], this should be one: %i \n',methodlist{m_c},nMomentslist(n_c),issorted(e_grid,'strictascend'))
        fprintf('method=%s, nMoments=%i: pi_e sums to one [T1], this should be zero: %2.8e \n',methodlist{m_c},nMomentslist(n_c),abs(sum(pi_e)-1))
        fprintf('method=%s, nMoments=%i: pi_e equals its own reversal, this should be below %g: %2.8e \n',methodlist{m_c},nMomentslist(n_c),entropytol,max(abs(pi_e-flipud(pi_e))))

        m1=sum(pi_e.*e_grid);
        m2=sum(pi_e.*(e_grid-m1).^2);
        m3=sum(pi_e.*(e_grid-m1).^3);
        m4=sum(pi_e.*(e_grid-m1).^4);
        if nMomentslist(n_c)>=2
            fprintf('method=%s, nMoments=%i: mean and variance matched, these should be below %g: %2.8e %2.8e \n',methodlist{m_c},nMomentslist(n_c),entropytol,abs(m1-mew),abs(m2-sigma^2))
        else
            fprintf('method=%s, nMoments=%i: mean matched, this should be below %g: %2.8e \n',methodlist{m_c},nMomentslist(n_c),entropytol,abs(m1-mew))
        end
        if nMomentslist(n_c)==4
            fprintf('method=%s, nMoments=4: 3rd and 4th central moments matched, these should be below %g: %2.8e %2.8e \n',methodlist{m_c},entropytol,abs(m3-0),abs(m4-3*sigma^4))
        end
        err4(m_c,n_c)=abs(m4-3*sigma^4);
    end
end
output.err4=err4;
output.methodlist=methodlist;
output.nMomentslist=nMomentslist;

% Ordering: asking for 4 moments must not do WORSE on the 4th moment than asking for 2.
fprintf('nMoments=4 beats nMoments=2 on the 4th central moment, for every method [T2], this should be one: %i \n',all(err4(:,3)<=err4(:,2)))

% parallel: 0 and 1 are both cpu paths
opts0=struct(); opts0.parallel=0;
opts1=struct(); opts1.parallel=1;
[eg0,pe0]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,opts0);
[eg1,pe1]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,opts1);
fprintf('parallel=0 vs 1: e_grid [T0], this should be zero: %2.8e \n',max(abs(eg0-eg1)))
fprintf('parallel=0 vs 1: pi_e [T0], this should be zero: %2.8e \n',max(abs(pe0-pe1)))
if gpuDeviceCount>0
    opts2=struct(); opts2.parallel=2;
    [eg2,pe2]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,opts2);
    fprintf('parallel=2 returns gpuArrays [T0], this should be one: %i \n',isa(eg2,'gpuArray')&&isa(pe2,'gpuArray'))
    fprintf('parallel=1 vs 2: e_grid [T1], this should be zero: %2.8e \n',max(abs(eg1-gather(eg2))))
    fprintf('parallel=1 vs 2: pi_e [T1], this should be zero: %2.8e \n',max(abs(pe1-gather(pe2))))
else
    fprintf('parallel=2 skipped: no gpu on this machine \n')
end

%% Error paths: the documented errors must actually fire
try
    discretizeIIDNormal_TanakaToda(mew,sigma,2,struct());
    fprintf('enum<3 should have errored, this should be one: 0 \n')
catch
    fprintf('enum<3 errors as documented, this should be one: 1 \n')
end
try
    optsbad=struct(); optsbad.nMoments=5;
    discretizeIIDNormal_TanakaToda(mew,sigma,15,optsbad);
    fprintf('nMoments=5 should have errored, this should be one: 0 \n')
catch
    fprintf('nMoments=5 errors as documented, this should be one: 1 \n')
end

%% Figure
figure(figure_c)
subplot(1,3,1); plot(znums,max(err_mean,10^(-18)),'o-',znums,max(err_var,10^(-18)),'s-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('discretizeIIDNormal\_TanakaToda: |error|'); xlabel('enum'); legend('mean','variance','Location','best')
subplot(1,3,2); plot(znums,runtime,'o-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('enum')
subplot(1,3,3); bar(err4)
set(gca,'YScale','log'); set(gca,'XTickLabel',methodlist)
title('|error| on 4th central moment'); legend('nMoments=1','nMoments=2','nMoments=4','Location','best')

end
