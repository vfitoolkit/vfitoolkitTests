function output=DiscP1_IIDNormal_Tauchen(calib,znums,figure_c)
% P1: discretizeIIDNormal_Tauchen
%
% P1 owns the option sweep for this command (it is native here). The structure of this subcode is
% the TEMPLATE that every other method subcode in this bank copies: config sweep, inline invariant
% block, inline timing, accuracy against analytic truth. The invariant block and the timing
% protocol are written out in full here rather than shared, per the toolkit's no-helper-functions
% rule - which means the copies must be kept identical by hand. Every assertion below has a
% distinctive message string so that a grep across the subcodes finds all of them when one changes.
%
% Truth for an iid normal is closed form: mean=mew, variance=sigma^2. Tauchen assigns probability
% by cdf mass between midpoints, so it does NOT match those exactly - the moment checks are T2 and
% are paired with a convergence assertion in enum.

fprintf('\n========== P1: discretizeIIDNormal_Tauchen ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;

%% Main sweep: enum, at the default Tauchen_q
% One sweep serves both accuracy and runtime: the timed calls are the calls the accuracy analysis
% reads. Discretization is deterministic, so the repeats all return identical (grid,pi).
Tauchen_q=3;
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); runtime=zeros(1,nz); nrepsused=zeros(1,nz);

for c_c=1:nz
    enum=znums(c_c);
    tauchenoptions=struct();

    % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
    % nreps drops to 1 for any config whose warm-up exceeded calib.timethreshold seconds.
    tic;
    [e_grid,pi_e]=discretizeIIDNormal_Tauchen(mew,sigma,enum,Tauchen_q,tauchenoptions);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [e_grid,pi_e]=discretizeIIDNormal_Tauchen(mew,sigma,enum,Tauchen_q,tauchenoptions);
        if isfield(tauchenoptions,'parallel')
            if tauchenoptions.parallel==2
                wait(gpuDevice) % otherwise the timing measures kernel launch, not kernel execution
            end
        end
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    nrepsused(c_c)=nreps;

    % --- invariants
    fprintf('enum=%i: size of e_grid [T0], this should be zero: %i \n',enum,any(size(e_grid)~=[enum,1]))
    fprintf('enum=%i: size of pi_e [T0], this should be zero: %i \n',enum,any(size(pi_e)~=[enum,1]))
    fprintf('enum=%i: e_grid is strictly ascending [T0], this should be one: %i \n',enum,issorted(e_grid,'strictascend'))
    fprintf('enum=%i: pi_e is in [0,1] [T0], this should be zero: %i \n',enum,any(pi_e<0)+any(pi_e>1))
    fprintf('enum=%i: pi_e sums to one [T1], this should be zero: %2.8e \n',enum,abs(sum(pi_e)-1))
    fprintf('enum=%i: no NaN or Inf [T0], this should be zero: %i \n',enum,any(~isfinite(e_grid))+any(~isfinite(pi_e)))
    % symmetry. The iid normal is symmetric about mew, and this command builds a grid symmetric
    % about mew, so pi_e must equal its own reversal. This is the cheap check that catches an
    % off-by-one or a one-tail-only correction: pi_e could still sum to one with either of those.
    fprintf('enum=%i: e_grid symmetric about mew [T1], this should be zero: %2.8e \n',enum,max(abs((e_grid+flipud(e_grid))/2-mew)))
    fprintf('enum=%i: pi_e equals its own reversal [T1], this should be zero: %2.8e \n',enum,max(abs(pi_e-flipud(pi_e))))

    % --- accuracy against closed form
    % The MEAN is exact for any enum, not approximate: the grid is symmetric about mew and the bin
    % masses are symmetric with it, so the mean is mew by symmetry. That makes it a T1 check, and
    % asserting that it "converges" would be asserting something already exact.
    m=sum(pi_e.*e_grid);
    v=sum(pi_e.*(e_grid-m).^2);
    err_mean(c_c)=abs(m-mew);
    err_var(c_c)=abs(v-sigma^2);
    fprintf('enum=%i: mean is exact by symmetry [T1], this should be zero: %2.8e \n',enum,err_mean(c_c))
    fprintf('enum=%i: variance error [T2] %2.3e, runtime %2.6f s (nreps=%i) \n',enum,err_var(c_c),runtime(c_c),nreps)

    output.sweep(c_c).enum=enum;
    output.sweep(c_c).e_grid=e_grid;
    output.sweep(c_c).pi_e=pi_e;
    output.sweep(c_c).moments=[m,v];
end

% --- convergence, and what it converges TO.
% With Tauchen_q held fixed the grid never reaches beyond mew +- Tauchen_q*sigma no matter how many
% points are added, and all the tail mass is lumped onto the two end points. So the variance error
% does NOT go to zero: it falls while the interior spacing is the binding error, then flattens out
% at the TRUNCATION FLOOR set by Tauchen_q. Asserting monotone convergence to zero here would be
% asserting something the method does not do at fixed q - and would fail, because the error can
% tick back up once the floor is reached.
fprintf('variance error falls over the first half of the sweep [T2], this should be one: %i \n',err_var(4)<err_var(1))
fprintf('variance error has flattened by the end of the sweep (within a factor of 10) [T2], this should be one: %i \n',err_var(end)/err_var(end-1)<10 && err_var(end-1)/err_var(end)<10)
fprintf('the floor it flattens to is the tail mass beyond +-Tauchen_q*sigma, which for q=%i is %2.3e \n',Tauchen_q,2*(normcdf(-Tauchen_q)*(Tauchen_q*sigma)^2))
% To get genuine convergence, Tauchen_q has to GROW with enum, so that the truncation floor keeps
% receding while the interior spacing shrinks. With q growing, the error that remains is the
% interior discretization error, which is O(spacing^2). So the sharp statement is not "the error
% falls" but "error/spacing^2 tends to a constant", which pins the convergence ORDER rather than
% just its direction.
err_var_growq=zeros(1,nz); spacing_growq=zeros(1,nz);
for c_c=1:nz
    enum=znums(c_c);
    qgrow=1.2*log(enum); % Floden's spacing rule, which grows with enum
    [eg,pe]=discretizeIIDNormal_Tauchen(mew,sigma,enum,qgrow,struct());
    mg=sum(pe.*eg);
    err_var_growq(c_c)=abs(sum(pe.*(eg-mg).^2)-sigma^2);
    spacing_growq(c_c)=(eg(2)-eg(1))/sigma;
end
ratio_growq=err_var_growq./(spacing_growq.^2);
for c_c=1:nz
    fprintf('growing q: enum=%i, q=%2.3f, spacing/sigma=%2.4f, variance error %2.3e, error/spacing^2 %2.3e \n',znums(c_c),1.2*log(znums(c_c)),spacing_growq(c_c),err_var_growq(c_c),ratio_growq(c_c))
end
% The enum=5 entry is NOT a usable baseline and the assertions below start at enum=9. At
% q=1.2*log(5)=1.93 the grid is so narrow that the truncation error, which understates the
% variance, very nearly cancels the coarse-spacing error, which overstates it. That gives an
% accidentally tiny error (~6e-05 on this calibration) and makes the sequence RISE from enum=5 to
% enum=9 before it starts converging. Using a fluke cancellation as the reference point would make
% this check unpassable, which is exactly what it did on the first run.
fprintf('with q growing, the variance error falls monotonically from enum=%i onward [T2], this should be one: %i \n',znums(2),all(diff(err_var_growq(2:end))<0))
fprintf('with q growing, the variance error is O(spacing^2): error/spacing^2 is constant to within 15%% over the last three enum [T2], this should be one: %i \n',(max(ratio_growq(end-2:end))/min(ratio_growq(end-2:end))-1)<0.15)
fprintf('   (that constant is %2.3e; it pins the convergence ORDER, not just the direction) \n',ratio_growq(end))
output.err_var_growq=err_var_growq;
output.spacing_growq=spacing_growq;
output.ratio_growq=ratio_growq;
output.err_mean=err_mean;
output.err_var=err_var;
output.runtime=runtime;
output.nrepsused=nrepsused;
output.znums=znums;

%% Option sweep at a fixed enum
enum=15;

% Tauchen_q sets how far the grid reaches. A wider grid captures more of the tail (better
% variance) at the cost of coarser spacing, so this is a trade rather than an ordering.
for qq=[2,3,4]
    [e_grid,pi_e]=discretizeIIDNormal_Tauchen(mew,sigma,enum,qq,struct());
    m=sum(pi_e.*e_grid); v=sum(pi_e.*(e_grid-m).^2);
    fprintf('Tauchen_q=%i: grid half-width %2.4f (should be %2.4f), mean error [T2] %2.3e, variance error [T2] %2.3e \n',qq,(e_grid(end)-e_grid(1))/2,qq*sigma,abs(m-mew),abs(v-sigma^2))
    fprintf('Tauchen_q=%i: grid half-width matches Tauchen_q*sigma [T1], this should be zero: %2.8e \n',qq,abs((e_grid(end)-e_grid(1))/2-qq*sigma))
end

% Tauchen_q=[] takes the default, min(sqrt(enum-1),4). This command had no default until now - the
% argument was simply required - while discretizeAR1_Tauchen has had one since the width question
% was settled. Checked across enum so BOTH limbs of the min are exercised: below enum=17 the
% sqrt(enum-1) limb binds, at and above it the cap of 4 does, and a rule that only ever tested one
% side would not notice if the other were wrong.
for ee=[5,9,17,31]
    [egd,ped]=discretizeIIDNormal_Tauchen(mew,sigma,ee,[],struct());
    egd=gather(egd); ped=gather(ped);
    qdef=min(sqrt(ee-1),4);
    [egx,pex]=discretizeIIDNormal_Tauchen(mew,sigma,ee,qdef,struct());
    egx=gather(egx); pex=gather(pex);
    fprintf('enum=%2i: Tauchen_q=[] gives the half-width min(sqrt(enum-1),4)=%4.2f [T1], this should be zero: %2.8e \n',ee,qdef,abs((egd(end)-egd(1))/2-qdef*sigma))
    fprintf('enum=%2i: and is identical to passing that number explicitly: e_grid [T0], this should be zero: %2.8e \n',ee,max(abs(egd-egx)))
    fprintf('enum=%2i: and pi_e likewise [T0], this should be zero: %2.8e \n',ee,max(abs(ped-pex)))
end
fprintf('the cap binds from enum=17 upward, so the sweep above exercises both limbs of the min [T0], this should be one: %i \n',(min(sqrt(17-1),4)==4)&&(min(sqrt(9-1),4)<4))

% parallel: 0 and 1 are both cpu paths, and 2 (gpu) now runs the same erfc expression, so all
% three are T0.
opts0=struct(); opts0.parallel=0;
opts1=struct(); opts1.parallel=1;
[eg0,pe0]=discretizeIIDNormal_Tauchen(mew,sigma,enum,Tauchen_q,opts0);
[eg1,pe1]=discretizeIIDNormal_Tauchen(mew,sigma,enum,Tauchen_q,opts1);
fprintf('parallel=0 vs 1: e_grid [T0], this should be zero: %2.8e \n',max(abs(eg0-eg1)))
fprintf('parallel=0 vs 1: pi_e [T0], this should be zero: %2.8e \n',max(abs(pe0-pe1)))
if gpuDeviceCount>0
    opts2=struct(); opts2.parallel=2;
    [eg2,pe2]=discretizeIIDNormal_Tauchen(mew,sigma,enum,Tauchen_q,opts2);
    fprintf('parallel=2 returns gpuArrays [T0], this should be one: %i \n',isa(eg2,'gpuArray')&&isa(pe2,'gpuArray'))
    fprintf('parallel=1 vs 2: e_grid [T0], this should be zero: %2.8e \n',max(abs(eg1-gather(eg2))))
    fprintf('parallel=1 vs 2: pi_e [T0], this should be zero: %2.8e \n',max(abs(pe1-gather(pe2))))
    % e_grid is built on the cpu in both branches and then moved, so it is a strict zero. pi_e can
    % sit a couple of ulps off: erfc is not correctly rounded, so the cpu and gpu libraries may
    % disagree in the last bit. Was T1 until the gpu branch stopped building the cdf from 1+erf.
else
    fprintf('parallel=2 skipped: no gpu on this machine \n')
end

%% Figure
figure(figure_c)
subplot(1,3,1); plot(znums,err_mean,'o-',znums,err_var,'s-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('discretizeIIDNormal\_Tauchen: |error|'); xlabel('enum'); legend('mean','variance','Location','best')
subplot(1,3,2); plot(znums,runtime,'o-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('enum')
subplot(1,3,3); bar(output.sweep(3).e_grid,output.sweep(3).pi_e)
title(['pi\_e at enum=',num2str(output.sweep(3).enum)]); xlabel('e')

end
