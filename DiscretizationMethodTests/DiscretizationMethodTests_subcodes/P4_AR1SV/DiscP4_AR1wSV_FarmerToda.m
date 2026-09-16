function output=DiscP4_AR1wSV_FarmerToda(calib,znums,figure_c)
% P4: discretizeAR1wSV_FarmerToda
%
% This block is the bank's first two-dimensional command, so as much of it is about the STACKED
% OUTPUT CONVENTION as about accuracy: z_grid is [x_grid; z_grid] of length xnum+znum, pi_z is the
% joint transition on (x,z) of size (xnum*znum)^2, and x varies FASTEST in the joint index. With
% xnum==znum an error in that stacking is invisible - every wrong slice still has the right size -
% so the sweep runs xnum and znum INDEPENDENTLY over the same range.
%
% Truth: E[z]=0, Var(z)=sigmau^2/(1-rho^2), autocorr(z)=rho, and the excess kurtosis of z from the
% convergent series in DiscSetup_AR1SV. The x block is an ordinary gaussian AR(1) with mean xBar,
% variance sigmae^2/(1-phi^2) and autocorrelation phi, so it has its own exact truth.

fprintf('\n========== P4: discretizeAR1wSV_FarmerToda ========== \n')

output=struct();
rho=calib.rho; phi=calib.phi; sigmau=calib.sigmau; sigmae=calib.sigmae;
Tauchen_q=3;
nz=length(znums);

err_var=zeros(nz,nz); err_ac=zeros(nz,nz); err_exkurt=zeros(nz,nz);
err_xvar=zeros(nz,nz); err_xac=zeros(nz,nz);
runtime=zeros(nz,nz); skipped=zeros(nz,nz);
fallback=NaN(nz,nz); % fraction of rows that fell back to one conditional moment

for x_c=1:nz
    for z_c=1:nz
        xnum=znums(x_c); znum=znums(z_c);
        if xnum*znum>calib.Ncap
            skipped(x_c,z_c)=1;
            fprintf('xnum=%i, znum=%i: SKIPPED, prod=%i exceeds the ceiling of %i (the joint matrix would be %i-by-%i) \n',xnum,znum,xnum*znum,calib.Ncap,xnum*znum,xnum*znum)
            continue
        end
        opts=struct(); opts.nSigmas=Tauchen_q;

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        tic;
        [z_grid,pi_z,otheroutputs]=discretizeAR1wSV_FarmerToda(rho,phi,sigmau,sigmae,xnum,znum,opts);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z]=discretizeAR1wSV_FarmerToda(rho,phi,sigmau,sigmae,xnum,znum,opts);
            treps(r_c)=toc;
        end
        runtime(x_c,z_c)=median(treps);

        % --- how many rows actually got two conditional moments. Before B25 this command kept no
        % record of its fallbacks and warned once per row instead: the run of 2026-08-25 emitted
        % 1348 warnings and 5676 lines, over half the diary, and the information in them was not
        % recoverable from the outputs. It now returns otheroutputs.nMoments_grid, so the fallback
        % rate is a measured quantity like any other. It matters because a row that fell back to
        % one moment has an unconstrained conditional variance, which is the suspected mechanism
        % behind the variance error stalling rather than converging. Taken from the warm-up call
        % above, which is discarded for timing anyway, so measuring it costs nothing.
        nMg=otheroutputs.nMoments_grid(:);
        fallback(x_c,z_c)=mean(nMg<2);
        fprintf('xnum=%i znum=%i: rows matching fewer than 2 conditional moments: %2.1f%% (%i of %i) \n',xnum,znum,100*fallback(x_c,z_c),sum(nMg<2),numel(nMg))

        % --- invariants, including the ones specific to a stacked two-dimensional output
        fprintf('xnum=%i znum=%i: size of z_grid is xnum+znum [T0], this should be zero: %i \n',xnum,znum,any(size(z_grid)~=[xnum+znum,1]))
        fprintf('xnum=%i znum=%i: size of pi_z is (xnum*znum)^2 [T0], this should be zero: %i \n',xnum,znum,any(size(pi_z)~=[xnum*znum,xnum*znum]))
        xg=z_grid(1:xnum); zg=z_grid(xnum+1:end);
        fprintf('xnum=%i znum=%i: the x block of the grid is strictly ascending [T0], this should be one: %i \n',xnum,znum,issorted(xg,'strictascend'))
        fprintf('xnum=%i znum=%i: the z block of the grid is strictly ascending [T0], this should be one: %i \n',xnum,znum,issorted(zg,'strictascend'))
        fprintf('xnum=%i znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',xnum,znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        fprintf('xnum=%i znum=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',xnum,znum,calib.entropytol,max(abs(sum(pi_z,2)-1)))
        fprintf('xnum=%i znum=%i: no NaN or Inf [T0], this should be zero: %i \n',xnum,znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))
        % z is mean zero, so its grid must be symmetric about zero
        fprintf('xnum=%i znum=%i: the z grid is symmetric about zero [T1], this should be zero: %2.8e \n',xnum,znum,max(abs((zg+flipud(zg))/2)))
        % and the x grid symmetric about xBar, since x is a gaussian AR(1)
        fprintf('xnum=%i znum=%i: the x grid is symmetric about xBar [T1], this should be zero: %2.8e \n',xnum,znum,max(abs((xg+flipud(xg))/2-calib.xBar)))

        % --- the joint stationary distribution, and the two marginals read off it.
        % This is what pins the ordering convention: if x did not vary fastest, the marginals below
        % would come out transposed and their moments would be wrong.
        % The joint stationary distribution, by iteration. MarkovChainMoments is deliberately NOT
        % used here: it computes moments against a single z_grid, and this command's grid is a
        % STACKED PAIR of grids rather than the values of one variable, so its mean and variance
        % would be meaningless. Only the distribution is wanted, and iterating for it is inlined
        % rather than shared, per the no-helper-functions rule.
        statdist=ones(xnum*znum,1)/(xnum*znum);
        for i_c=1:10000
            statdistnew=pi_z'*statdist;
            if max(abs(statdistnew-statdist))<10^(-14)
                statdist=statdistnew;
                break
            end
            statdist=statdistnew;
        end
        statdist=statdist/sum(statdist);
        Pxz=reshape(statdist,[xnum,znum]); % x fastest, so this reshape IS the ordering assumption
        xmarg=sum(Pxz,2); zmarg=sum(Pxz,1)';
        fprintf('xnum=%i znum=%i: the joint stationary distribution sums to one [T1], this should be zero: %2.8e \n',xnum,znum,abs(sum(statdist)-1))

        mx=sum(xmarg.*xg); vx=sum(xmarg.*(xg-mx).^2);
        mz=sum(zmarg.*zg); vz=sum(zmarg.*(zg-mz).^2);
        ek=sum(zmarg.*(zg-mz).^4)/vz^2-3;
        err_xvar(x_c,z_c)=abs(vx-calib.x.var);
        err_var(x_c,z_c)=abs(vz-calib.z.var);
        err_exkurt(x_c,z_c)=abs(ek-calib.z.exkurt);
        fprintf('xnum=%i znum=%i: E[x] %2.4f vs xBar %2.4f; Var(x) %2.4f vs truth %2.4f (err [T2] %2.3e) \n',xnum,znum,mx,calib.xBar,vx,calib.x.var,err_xvar(x_c,z_c))
        fprintf('xnum=%i znum=%i: E[z] %2.3e (truth zero); Var(z) %2.4f vs truth %2.4f (err [T2] %2.3e) \n',xnum,znum,mz,vz,calib.z.var,err_var(x_c,z_c))
        fprintf('xnum=%i znum=%i: excess kurtosis of z %2.4f vs truth %2.4f (err [T2] %2.3e) \n',xnum,znum,ek,calib.z.exkurt,err_exkurt(x_c,z_c))
        fprintf('xnum=%i znum=%i: runtime %2.6f s \n',xnum,znum,runtime(x_c,z_c))

        output.sweep(x_c,z_c).xnum=xnum;
        output.sweep(x_c,z_c).znum=znum;
        output.sweep(x_c,z_c).moments=[mz,vz,ek,mx,vx];
    end
end
output.err_var=err_var; output.err_ac=err_ac; output.err_exkurt=err_exkurt;
output.err_xvar=err_xvar; output.runtime=runtime; output.skipped=skipped; output.znums=znums;
output.fallback=fallback;

% Excess kurtosis is the moment stochastic volatility exists to produce, so it gets the
% convergence assertion. The diagonal (xnum==znum) is used, since convergence in one dimension
% while the other is held coarse is not what "refining the grid" means.
diagek=zeros(1,nz); ok=false(1,nz);
for c_c=1:nz
    if ~skipped(c_c,c_c)
        diagek(c_c)=err_exkurt(c_c,c_c); ok(c_c)=true;
    end
end
% NOTE ON THE FORM OF THIS ASSERTION. It used to compare the first diagonal cell against the
% last, and that was too weak to be worth having: on the run of 2026-08-25 the FarmerToda errors
% were 1.759, 1.268, 1.305, 1.349 - rising monotonically over the last three cells - and a
% first-versus-last test passed on that. What is wanted is that refining the grid never makes the
% answer worse, so every consecutive pair is checked.
dek=diagek(ok);
worst=0; worstat=0;
for c_c=2:length(dek)
    if dek(c_c)-dek(c_c-1)>worst
        worst=dek(c_c)-dek(c_c-1); worstat=c_c;
    end
end
fprintf('excess kurtosis error never rises between consecutive diagonal cells [T2], this should be one: %i \n',worst<=10^(-12))
if worst>10^(-12)
    zd=znums(ok);
    fprintf('   the worst increase is from xnum=znum=%i to %i, %2.4e to %2.4e \n',zd(worstat-1),zd(worstat),dek(worstat-1),dek(worstat))
end
fprintf('   (diagonal means xnum==znum; refining one dimension while the other stays coarse is a \n')
fprintf('    different question, and the off-diagonal cells are here to test the stacking, not \n')
fprintf('    convergence) \n')

%% Figure
figure(figure_c)
subplot(1,2,1)
imagesc(log10(max(err_exkurt,10^(-18)))); colorbar
set(gca,'XTick',1:nz,'XTickLabel',znums,'YTick',1:nz,'YTickLabel',znums)
xlabel('znum'); ylabel('xnum'); title('discretizeAR1wSV_FarmerToda: log10 |excess kurtosis error|','Interpreter','none')
subplot(1,2,2)
imagesc(log10(max(runtime,10^(-18)))); colorbar
set(gca,'XTick',1:nz,'XTickLabel',znums,'YTick',1:nz,'YTickLabel',znums)
xlabel('znum'); ylabel('xnum'); title('log10 runtime (s)')

end

