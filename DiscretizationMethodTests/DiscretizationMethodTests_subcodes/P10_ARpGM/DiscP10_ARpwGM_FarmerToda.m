function output=DiscP10_ARpwGM_FarmerToda(calib,cn,figure_c)
% P10: discretizeARpwGM_FarmerToda
%
% A command written for this block, so nothing here is a regression test - it is the first exercise.
%
% WHAT THE MIXTURE ADDS over the gaussian sibling, and why the block runs both. With gaussian
% innovations z is gaussian, so its skewness and excess kurtosis are zero and matching more than two
% conditional moments buys nothing. With a mixture they are not zero, and they are the reason to use
% a mixture at all - so they are what this subcode is really measuring, and mean and variance are
% hygiene. The truth for both comes from the MA(inf) cumulant sum kappa_n(z)=(sum_j psi_j^n)*kappa_n(e).
%
% AND THE MIXTURE IS NOT MEAN ZERO, deliberately. That makes the grid centre mew/(1-sum(Rho)) differ
% from the true mean (mew+E(e))/(1-sum(Rho)), so a command that confused the two would show it. The
% gaussian sibling cannot see that distinction, since a gaussian innovation here has mean zero.
%
% WHAT IS DIFFERENT FROM EVERY EARLIER BLOCK. The state is p-dimensional but only ONE dimension is
% stochastic: the new z is drawn and the p-1 lags shift along one place. So pi_z is structurally
% sparse - exactly znum non-zeros per row out of znum^p, at positions fixed by the shift - and that
% sparsity is a property the command must have rather than something the parameters happen to give.
% It is asserted directly below, index by index, because it is the one check that would catch a
% command that built a correct-looking stochastic matrix with the lags wired up wrongly.
%
% AND THE AUTOCORRELATIONS COME FOR FREE. Dimension k of the state IS z lagged k-1 periods, so under
% the stationary distribution the autocorrelation of z at lag j is just the cross-correlation
% between dimensions 1 and j+1 of the same state vector. No propagation needed, and it is a sharper
% check than it looks: it ties the lag structure to the moments.

fprintf('\n========== P10: discretizeARpwGM_FarmerToda, %s ========== \n',calib.(cn).name)

output=struct();
Rho=calib.(cn).Rho; p=calib.(cn).p;
mew=calib.mew;
mixprobs_i=calib.mixprobs_i; mu_i=calib.mu_i; sigma_i=calib.sigma_i;
T=calib.(cn).gm; % the analytic truth for this calibration with gaussian-mixture innovations
if p==2
    znums=calib.znums2;
else
    znums=calib.znums3;
end
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_ac=zeros(1,nz); runtime=zeros(1,nz); fallback=zeros(1,nz);
err_skew=zeros(1,nz); err_exkurt=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    fo=struct(); fo.verbose=0;

    tic;
    [z_grid,pi_z,oo]=discretizeARpwGM_FarmerToda(mew,Rho,mixprobs_i,mu_i,sigma_i,znum,fo);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid,pi_z]=discretizeARpwGM_FarmerToda(mew,Rho,mixprobs_i,mu_i,sigma_i,znum,fo);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    z_grid=gather(z_grid); pi_z=gather(pi_z);
    Nz=znum^p;

    % --- invariants
    fprintf('znum=%3i: z_grid is STACKED, p*znum-by-1 [T0], this should be zero: %i \n',znum,any(size(z_grid)~=[p*znum,1]))
    fprintf('znum=%3i: size of pi_z is znum^p-by-znum^p [T0], this should be zero: %i \n',znum,any(size(pi_z)~=[Nz,Nz]))
    % The p stacked blocks must be IDENTICAL: they are the same variable at different lags, so they
    % share one grid. A command that built p separate grids would still return the right shape.
    blockgap=0;
    for m_c=2:p
        blockgap=max(blockgap,max(abs(z_grid((m_c-1)*znum+1:m_c*znum)-z_grid(1:znum))));
    end
    fprintf('znum=%3i: the p stacked blocks are identical copies of one grid [T1], this should be zero: %2.8e \n',znum,blockgap)
    fprintf('znum=%3i: the shared grid is strictly ascending [T0], this should be one: %i \n',znum,issorted(z_grid(1:znum),'strictascend'))
    fprintf('znum=%3i: pi_z is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z(:)<0)+any(pi_z(:)>1))
    fprintf('znum=%3i: rows of pi_z sum to one, this should be below 1e-07: %2.8e \n',znum,max(abs(sum(pi_z,2)-1)))
    fprintf('znum=%3i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
    fallback(c_c)=mean(gather(oo.nMoments_grid)<4); % 4 is this command's default nMoments
    fprintf('znum=%3i: states matching fewer than the requested 4 moments: %2.1f%% \n',znum,100*fallback(c_c))

    % --- THE SPARSITY, index by index
    % From joint state ii the reachable states are exactly znum*mod(ii-1,znum^(p-1))+(1:znum).
    % Checking the pattern rather than just the count is what makes this catch a wired-up-wrongly
    % lag structure: a command that shifted the wrong way would have znum non-zeros in each row too.
    nnzrow=sum(pi_z~=0,2);
    patternok=1;
    for ii=1:Nz
        destbase=znum*mod(ii-1,znum^(p-1));
        should=false(1,Nz); should(destbase+(1:znum))=true;
        if any((pi_z(ii,:)~=0) & ~should)
            patternok=0; break
        end
    end
    fprintf('znum=%3i: every row has at most znum non-zeros [T0], this should be one: %i \n',znum,max(nnzrow)<=znum)
    fprintf('znum=%3i: and they sit exactly where the lag shift puts them [T0], this should be one: %i \n',znum,patternok)
    fprintf('znum=%3i: which is %i non-zeros out of %i, a density of %2.3f%% \n',znum,sum(nnzrow),Nz*Nz,100*sum(nnzrow)/(Nz*Nz))

    % --- otheroutputs against the setup's independently computed truth
    fprintf('znum=%3i: otheroutputs.mewz is the grid centre mew/(1-sum(Rho)) [T1], this should be zero: %2.8e \n',znum,abs(gather(oo.mewz)-T.mewz))
    fprintf('znum=%3i: otheroutputs.Ez is the true mean (mew+E(e))/(1-sum(Rho)) [T1], this should be zero: %2.8e \n',znum,abs(gather(oo.Ez)-T.Ez))
    fprintf('znum=%3i: and the two DIFFER, so the mixture being off-centre is visible [T0], this should be one: %i \n',znum,abs(gather(oo.Ez)-gather(oo.mewz))>1e-10)
    fprintf('znum=%3i: otheroutputs.sigmaz matches the MA(inf) truth [T1], this should be zero: %2.8e \n',znum,abs(gather(oo.sigmaz)-T.sd))
    fprintf('znum=%3i: otheroutputs.maxabseig matches the companion eigenvalue [T1], this should be zero: %2.8e \n',znum,abs(gather(oo.maxabseig)-calib.(cn).maxabseig))

    % --- the stationary distribution, and the moments it carries
    d=ones(Nz,1)/Nz;
    for it_c=1:100000
        dn=pi_z'*d;
        if max(abs(dn-d))<1e-15
            d=dn; break
        end
        d=dn;
    end
    d=d/sum(d);
    zv=CreateGridvals(znum*ones(p,1),z_grid,1); % (Nz)-by-p, dimension 1 fastest
    m1=d'*zv(:,1);
    v1=d'*((zv(:,1)-m1).^2);
    sk=(d'*((zv(:,1)-m1).^3))/v1^1.5;
    ek=(d'*((zv(:,1)-m1).^4))/v1^2-3;
    err_mean(c_c)=abs(m1-T.Ez);
    err_var(c_c)=abs(v1-T.var);
    fprintf('znum=%3i: E(z) %+2.6f against truth %+2.6f, error [T2] %2.3e \n',znum,m1,T.Ez,err_mean(c_c))
    fprintf('znum=%3i: Var(z) %2.6f against truth %2.6f, error [T2] %2.3e \n',znum,v1,T.var,err_var(c_c))
    err_skew(c_c)=abs(sk-T.skew); err_exkurt(c_c)=abs(ek-T.exkurt);
    fprintf('znum=%3i: skewness %+2.4f against truth %+2.4f, error [T2] %2.3e \n',znum,sk,T.skew,err_skew(c_c))
    fprintf('znum=%3i: excess kurtosis %+2.4f against truth %+2.4f, error [T2] %2.3e \n',znum,ek,T.exkurt,err_exkurt(c_c))

    % --- the lag structure, read straight off the state vector
    % Dimension k is z lagged k-1, so the stationary marginal of every dimension must be the SAME
    % distribution. If it is not, the shift is wrong in a way the row sums would never reveal.
    margapm=0; margapv=0;
    for m_c=2:p
        mk=d'*zv(:,m_c);
        vk=d'*((zv(:,m_c)-mk).^2);
        margapm=max(margapm,abs(mk-m1)); margapv=max(margapv,abs(vk-v1));
    end
    fprintf('znum=%3i: every lag dimension has the same stationary mean as dimension 1 [T1], this should be zero: %2.8e \n',znum,margapm)
    fprintf('znum=%3i: and the same variance [T1], this should be zero: %2.8e \n',znum,margapv)
    acj=zeros(1,p-1);
    for j_c=1:p-1
        mk=d'*zv(:,j_c+1); vk=d'*((zv(:,j_c+1)-mk).^2);
        acj(j_c)=(d'*((zv(:,1)-m1).*(zv(:,j_c+1)-mk)))/sqrt(v1*vk);
    end
    err_ac(c_c)=max(abs(acj-T.autocorr(2:p)));
    fprintf('znum=%3i: autocorrelations at lags 1 to %i:',znum,p-1);
    for j_c=1:p-1
        fprintf(' %+2.5f (truth %+2.5f)',acj(j_c),T.autocorr(j_c+1));
    end
    fprintf(', worst error [T2] %2.3e \n',err_ac(c_c))
    fprintf('znum=%3i: runtime %2.6f s (nreps=%i) \n',znum,runtime(c_c),nreps)
end
output.err_mean=err_mean; output.err_var=err_var; output.err_ac=err_ac;
output.err_skew=err_skew; output.err_exkurt=err_exkurt;
output.runtime=runtime; output.znums=znums; output.fallback=fallback;

%% Convergence
% THE CONVERGENCE ASSERTIONS ARE FLOOR-AWARE, and they have to be. The maximum entropy solve
% ACCEPTS a solution once norm(momentError)<=1e-5, so nothing downstream of it is exact to better
% than that - on the run of 2026-09-17 the autocorrelation errors across the sweep were 7.8e-08,
% 1.7e-06, 7.5e-07, 2.8e-07, 2.2e-07, 4.4e-07: all at the solver's noise floor, in no particular
% order, and "the last is smaller than the first" duly failed while nothing was wrong. So each
% assertion below is: the error either FELL, or it is already at the floor and there is nothing
% left to fall. This is the third time in this bank that ordering two numbers at a floor has
% produced a red with no content; the pattern is worth recognising on sight.
solverfloor=1e-05; % what the entropy solve accepts, so the finest error this sweep can show
fprintf('\nvariance error falls over the sweep, or is already at the solver floor [T2], this should be one: %i \n',(err_var(end)<err_var(1))||(err_var(end)<solverfloor))
fprintf('autocorrelation error falls over the sweep, or is already at the solver floor [T2], this should be one: %i \n',(err_ac(end)<err_ac(1))||(err_ac(end)<solverfloor))
fprintf('excess kurtosis error falls over the sweep, or is already at the solver floor [T2], this should be one: %i \n',(err_exkurt(end)<err_exkurt(1))||(err_exkurt(end)<solverfloor))
fprintf('   mean    :');
for c_c=1:nz
    fprintf(' %2.1e',err_mean(c_c));
end
fprintf(' \n   variance:');
for c_c=1:nz
    fprintf(' %2.1e',err_var(c_c));
end
fprintf(' \n   autocorr:');
for c_c=1:nz
    fprintf(' %2.1e',err_ac(c_c));
end
fprintf(' \n   skewness:');
for c_c=1:nz
    fprintf(' %2.1e',err_skew(c_c));
end
fprintf(' \n   ex kurt :');
for c_c=1:nz
    fprintf(' %2.1e',err_exkurt(c_c));
end
fprintf(' \n')

%% REGRESSION BARS
% First run of a new command, so these are set loose and will be tightened once there is a
% measurement to set them from. They exist so a gross regression cannot pass silently, which is what
% the verdict-less accuracy lines above would otherwise allow.
fprintf('\nworst mean error over the sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_mean))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',0.01,err_mean(end))
fprintf('worst variance error over the sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',0.01,err_var(end))
fprintf('worst autocorrelation error over the sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_ac))
fprintf('worst skewness error over the sweep [T2], this should be below %g: %2.3e \n',1,max(err_skew))
fprintf('worst excess kurtosis error over the sweep [T2], this should be below %g: %2.3e \n',2,max(err_exkurt))

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(znums,max(err_mean,1e-18),'o-',znums,max(err_var,1e-18),'s-',znums,max(err_ac,1e-18),'^-',znums,max(err_exkurt,1e-18),'d-')
set(gca,'XScale','log')
xlabel('znum per dimension'); ylabel('|error|')
title(['ARpwGM\_FarmerToda: ',calib.(cn).name]); legend('mean','variance','autocorr','ex kurtosis','Location','best')
subplot(1,2,2)
loglog(znums,runtime,'o-')
xlabel('znum per dimension'); ylabel('seconds'); title('runtime')
sgtitle(['P10: discretizeARpwGM\_FarmerToda, ',calib.(cn).name])

end
