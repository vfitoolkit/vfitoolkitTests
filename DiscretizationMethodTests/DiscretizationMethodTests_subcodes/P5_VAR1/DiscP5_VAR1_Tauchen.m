function output=DiscP5_VAR1_Tauchen(calib,znums,figure_c)
% P5: discretizeVAR1_Tauchen
%
% OUTPUT CONVENTION. This command returns z_grid STACKED, sum(znum)-by-1, one block per variable,
% and builds its joint index with CreateGridvals, so variable 1 varies FASTEST. Its sibling
% discretizeVAR1_FarmerToda returns the grid JOINT, (znum^M)-by-M, and uses allcomb2, whose default
% has the LAST variable varying fastest. Two commands for the same process, two shapes and two
% orderings. Neither subcode assumes the other's; both are pinned against kron in DiscP5_crosstests.
%
% Note this is not literally Tauchen's method: Tauchen forms the product of the marginal cdfs,
% this evaluates the multivariate cdf directly (mvncdf). The two coincide exactly when SigmaSq is
% diagonal, which is what makes the kron identity an exact statement here rather than an approximate
% one, and makes it a test of the code rather than of the method.
%
% Truth: zmean=(I-Rho)^(-1)*Mew, the Lyapunov solution for Var(z), and autocorrelation
% (Rho*Var(z))_ii/Var(z)_ii - all computed in DiscSetup_VAR1 and checked there against their own
% defining equations.

fprintf('\n========== P5: discretizeVAR1_Tauchen ========== \n')

output=struct();
nz=length(znums);
Mlist=[1,2,3];
err_mean=NaN(nz,3); err_var=NaN(nz,3); err_ac=NaN(nz,3); err_cross=NaN(nz,3);
err_sigmaz=NaN(nz,3); runtime=NaN(nz,3); skipped=ones(nz,3);
Tauchen_q=3;

for m_c=1:3
    M=Mlist(m_c);
    if M==1
        Mew=calib.full.Mew(1); Rho=calib.full.Rho(1,1); SigmaSq=calib.full.SigmaSq(1,1);
        zmeanT=Mew/(1-Rho); SigmaSqzT=SigmaSq/(1-Rho^2); acT=Rho;
    elseif M==2
        Mew=calib.full.Mew; Rho=calib.full.Rho; SigmaSq=calib.full.SigmaSq;
        zmeanT=calib.full.zmean; SigmaSqzT=calib.full.SigmaSqz; acT=calib.full.autocorr;
    else
        Mew=calib.three.Mew; Rho=calib.three.Rho; SigmaSq=calib.three.SigmaSq;
        zmeanT=calib.three.zmean; SigmaSqzT=calib.three.SigmaSqz; acT=calib.three.autocorr;
    end
    sigmazT=sqrt(diag(SigmaSqzT));

    for c_c=1:nz
        znum=znums(c_c);
        if znum^M>calib.Ncap
            fprintf('M=%i, znum=%i: SKIPPED, znum^M=%i exceeds the ceiling of %i (the joint matrix would be %i-by-%i) \n',M,znum,znum^M,calib.Ncap,znum^M,znum^M)
            continue
        end
        skipped(c_c,m_c)=0;
        tauchenoptions=struct(); tauchenoptions.verbose=0;

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        tic;
        [z_grid,pi_z]=discretizeVAR1_Tauchen(Mew,Rho,SigmaSq,znum,Tauchen_q,tauchenoptions);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z]=discretizeVAR1_Tauchen(Mew,Rho,SigmaSq,znum,Tauchen_q,tauchenoptions);
            treps(r_c)=toc;
        end
        runtime(c_c,m_c)=median(treps);

        z_grid=gather(z_grid); pi_z=gather(pi_z);

        % --- invariants, stacked-output form
        fprintf('M=%i znum=%i: size of z_grid is sum(znum) [T0], this should be zero: %i \n',M,znum,any(size(z_grid)~=[M*znum,1]))
        fprintf('M=%i znum=%i: size of pi_z is (znum^M)^2 [T0], this should be zero: %i \n',M,znum,any(size(pi_z)~=[znum^M,znum^M]))
        asc=1;
        for v_c=1:M
            asc=asc*issorted(z_grid((v_c-1)*znum+1:v_c*znum),'strictascend');
        end
        fprintf('M=%i znum=%i: every variable block of the grid is strictly ascending [T0], this should be one: %i \n',M,znum,asc)
        fprintf('M=%i znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',M,znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        fprintf('M=%i znum=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',M,znum,calib.entropytol,max(abs(sum(pi_z,2)-1)))
        fprintf('M=%i znum=%i: no NaN or Inf [T0], this should be zero: %i \n',M,znum,any(~isfinite(z_grid))+any(~isfinite(pi_z(:))))

        % --- THE GRID WIDTH ITSELF. This is the B19 regression, and it is checked before any
        % moment is computed, because B19 was a wrong sigmaz that made the grid five times too
        % wide while every probability still summed to one. The grid is
        % linspace(-1,1,znum)*Tauchen_q*sigmaz + zmean, so its half-width divided by Tauchen_q must
        % be sigmaz exactly. That is an identity, not an approximation.
        halfwidth=zeros(M,1); gridmid=zeros(M,1);
        for v_c=1:M
            blk=z_grid((v_c-1)*znum+1:v_c*znum);
            halfwidth(v_c)=(blk(end)-blk(1))/2;
            gridmid(v_c)=(blk(end)+blk(1))/2;
        end
        err_sigmaz(c_c,m_c)=max(abs(halfwidth/Tauchen_q-sigmazT));
        fprintf('M=%i znum=%i: the grid half-width is Tauchen_q*sigmaz [T1], this should be zero: %2.8e \n',M,znum,err_sigmaz(c_c,m_c))
        fprintf('M=%i znum=%i: the grid is centred on the unconditional mean [T1], this should be zero: %2.8e \n',M,znum,max(abs(gridmid-zmeanT)))

        % --- moments of the chain
        statdist=ones(znum^M,1)/(znum^M);
        for i_c=1:10000
            sn=pi_z'*statdist;
            if max(abs(sn-statdist))<10^(-14), statdist=sn; break; end
            statdist=sn;
        end
        statdist=statdist/sum(statdist);
        zvals=CreateGridvals(znum*ones(M,1),z_grid,1); % (znum^M)-by-M, variable 1 fastest
        mz=(statdist'*zvals)';
        dz=zvals-mz';
        Vz=(dz.*statdist)'*dz;
        % one-step-ahead cross moment, for the autocorrelation
        Ez1=pi_z*zvals; % E[z' | z], (znum^M)-by-M
        Cz=((zvals-mz').*statdist)'*(Ez1-mz');
        acz=diag(Cz)./diag(Vz);
        err_mean(c_c,m_c)=max(abs(mz-zmeanT));
        err_var(c_c,m_c)=max(abs(diag(Vz)-diag(SigmaSqzT)));
        err_ac(c_c,m_c)=max(abs(acz-acT));
        fprintf('M=%i znum=%i: mean error [T2] %2.3e, variance error [T2] %2.3e, autocorrelation error [T2] %2.3e \n',M,znum,err_mean(c_c,m_c),err_var(c_c,m_c),err_ac(c_c,m_c))
        if M>1
            err_cross(c_c,m_c)=abs(Vz(1,2)-SigmaSqzT(1,2));
            fprintf('M=%i znum=%i: the CROSS-covariance error [T2] %2.3e (this is the moment a pair of AR(1)s cannot reach) \n',M,znum,err_cross(c_c,m_c))
        end
        fprintf('M=%i znum=%i: runtime %2.6f s \n',M,znum,runtime(c_c,m_c))
        output.sweep(c_c,m_c).moments=[mz(:)',diag(Vz)',acz(:)'];
    end
end

output.err_mean=err_mean; output.err_var=err_var; output.err_ac=err_ac;
output.err_cross=err_cross; output.err_sigmaz=err_sigmaz;
output.runtime=runtime; output.skipped=skipped; output.znums=znums; output.Mlist=Mlist;

% Convergence, at M=2, on the cross-covariance. Refining the grid must not make it worse.
col=2; d=err_cross(:,col); ok=~skipped(:,col) & ~isnan(d);
dd=d(ok); worst=0; worstat=0;
for c_c=2:length(dd)
    if dd(c_c)-dd(c_c-1)>worst
        worst=dd(c_c)-dd(c_c-1); worstat=c_c;
    end
end
fprintf('M=2: the cross-covariance error never rises as znum grows [T2], this should be one: %i \n',worst<=10^(-12))
if worst>10^(-12)
    zd=znums(ok);
    fprintf('   the worst increase is from znum=%i to %i, %2.4e to %2.4e \n',zd(worstat-1),zd(worstat),dd(worstat-1),dd(worstat))
end

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(znums,max(err_var(:,2),10^(-18)),'o-',znums,max(err_ac(:,2),10^(-18)),'s-',znums,max(err_cross(:,2),10^(-18)),'d-')
xlabel('znum'); ylabel('error'); title('discretizeVAR1\_Tauchen, M=2')
legend({'variance','autocorrelation','cross-covariance'},'Location','best')
subplot(1,2,2)
semilogy(znums,max(err_sigmaz(:,1),10^(-18)),'o-',znums,max(err_sigmaz(:,2),10^(-18)),'s-',znums,max(err_sigmaz(:,3),10^(-18)),'d-')
xlabel('znum'); ylabel('|grid half-width/Tauchen\_q - sigmaz|'); title('the B19 regression')
legend({'M=1','M=2','M=3'},'Location','best')

end
