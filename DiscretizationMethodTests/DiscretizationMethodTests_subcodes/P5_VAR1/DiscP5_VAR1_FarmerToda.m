function output=DiscP5_VAR1_FarmerToda(calib,znums,figure_c)
% P5: discretizeVAR1_FarmerToda
%
% OUTPUT CONVENTION. This command returns z_grid JOINT, (znum^M)-by-M, one row per state tuple -
% the toolkit's [prod(n_z), length(n_z)] convention - where discretizeVAR1_Tauchen returns it
% stacked, sum(znum)-by-1. Because the grid lists every tuple explicitly, the MOMENTS below are
% ordering-agnostic: they read values out of the rows and never assume which variable moves fastest.
% The ordering itself is a separate question and is settled in DiscP5_crosstests against kron.
%
% The grid is built in a transformed space (x = chol(SigmaSq)^(-1) z, so the innovations are
% standard normal and independent) and mapped back at the end. That is exactly what
% discretizeVAR1_Tauchen fails to do, which is B19: the same Lyapunov solve returns Var(x), not
% Var(z), and this command is correct only because it never leaves the transformed space until the
% mapping back.

fprintf('\n========== P5: discretizeVAR1_FarmerToda ========== \n')

output=struct();
nz=length(znums);
Mlist=[1,2,3];
err_mean=NaN(nz,3); err_var=NaN(nz,3); err_ac=NaN(nz,3); err_cross=NaN(nz,3);
fallback=NaN(nz,3); runtime=NaN(nz,3); skipped=ones(nz,3);

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

    for c_c=1:nz
        znum=znums(c_c);
        if znum^M>calib.Ncap
            fprintf('M=%i, znum=%i: SKIPPED, znum^M=%i exceeds the ceiling of %i (the joint matrix would be %i-by-%i) \n',M,znum,znum^M,calib.Ncap,znum^M,znum^M)
            continue
        end
        skipped(c_c,m_c)=0;
        opts=struct(); opts.parallel=1; opts.verbose=0;

        % --- timing: one warm-up call, discarded, then nreps timed calls, report the median.
        % otheroutputs is taken from the warm-up, which is discarded anyway, so the fallback
        % measurement is free (see B25).
        tic;
        [z_grid,pi_z,otheroutputs]=discretizeVAR1_FarmerToda(Mew,Rho,SigmaSq,znum,opts);
        twarm=toc;
        if twarm>calib.timethreshold
            nreps=1;
        else
            nreps=calib.nreps;
        end
        treps=zeros(1,nreps);
        for r_c=1:nreps
            tic;
            [z_grid,pi_z]=discretizeVAR1_FarmerToda(Mew,Rho,SigmaSq,znum,opts);
            treps(r_c)=toc;
        end
        runtime(c_c,m_c)=median(treps);

        z_grid=gather(z_grid); pi_z=gather(pi_z);

        % --- invariants, joint-output form
        fprintf('M=%i znum=%i: size of z_grid is (znum^M)-by-M [T0], this should be zero: %i \n',M,znum,any(size(z_grid)~=[znum^M,M]))
        fprintf('M=%i znum=%i: size of pi_z is (znum^M)^2 [T0], this should be zero: %i \n',M,znum,any(size(pi_z)~=[znum^M,znum^M]))
        fprintf('M=%i znum=%i: pi_z is in [0,1] [T0], this should be zero: %i \n',M,znum,any(pi_z(:)<0)+any(pi_z(:)>1))
        fprintf('M=%i znum=%i: rows of pi_z sum to one, this should be below %g: %2.8e \n',M,znum,calib.entropytol,max(abs(sum(pi_z,2)-1)))
        fprintf('M=%i znum=%i: no NaN or Inf [T0], this should be zero: %i \n',M,znum,any(~isfinite(z_grid(:)))+any(~isfinite(pi_z(:))))
        % The grid is NOT a tensor grid in z-space, and it is worth being precise about why, because
        % the obvious invariant - "each variable takes exactly znum distinct values" - is FALSE here
        % and would have been a wrong test. The command builds a tensor grid in a transformed space
        % and maps it back with z = C*y + mu where C = chol(SigmaSq,'lower')*U and U is an orthogonal
        % rotation from minVarTrace. A general linear map of a tensor grid is not a tensor grid; the
        % command's own header says as much ("It is NOT a kronecker-product grid"). Variable v takes
        % up to znum^M distinct values.
        %
        % What IS exactly true is that z_grid(:,v) = mu_v + sum_k C(v,k)*y1D(k,i_k), which is
        % ADDITIVELY SEPARABLE in the M grid indices. Additive separability is equivalent to every
        % mixed second difference vanishing, and that is an identity, so it is a T1 check. It fails
        % if the grid is assembled in the wrong order, or if the mapping back is applied to the
        % wrong thing.
        sepresid=0;
        if M==2
            for v_c=1:M
                A=reshape(z_grid(:,v_c),[znum,znum]);
                sepresid=max(sepresid,max(abs(A-A(:,1)-A(1,:)+A(1,1)),[],'all'));
            end
        elseif M==3
            for v_c=1:M
                A=reshape(z_grid(:,v_c),[znum,znum,znum]);
                sepresid=max(sepresid,max(abs(A-A(:,1,1)-reshape(A(1,:,1),[1,znum,1])-reshape(A(1,1,:),[1,1,znum])+2*A(1,1,1)),[],'all'));
            end
        end
        if M>1
            fprintf('M=%i znum=%i: the grid is additively separable in the grid indices [T1], this should be zero: %2.8e \n',M,znum,sepresid)
        end
        % Every row of pi_z is a product distribution: the maximum entropy problem is solved one
        % variable at a time in the transformed space, where the innovations are independent, and
        % the row is assembled with prod(allcomb2(temp),2). So each row, reshaped to the grid shape,
        % is rank one. This holds whatever SigmaSq is, and it is the structural fact that makes the
        % kron identity in DiscP5_crosstests exact rather than approximate.
        if M==2
            rk=0;
            for r_c=1:znum^M
                sv=svd(reshape(pi_z(r_c,:),[znum,znum]));
                rk=max(rk,sv(2)/sv(1));
            end
            fprintf('M=%i znum=%i: every row of pi_z is a product distribution (rank one) [T1], this should be zero: %2.8e \n',M,znum,rk)
        end

        % --- the maximum entropy fallback rate (available since B25)
        nMg=otheroutputs.nMoments_grid(:);
        fallback(c_c,m_c)=mean(nMg<2);
        fprintf('M=%i znum=%i: (row,variable) pairs matching fewer than 2 moments: %2.1f%% (%i of %i) \n',M,znum,100*fallback(c_c,m_c),sum(nMg<2),numel(nMg))

        % --- moments of the chain, read from the joint grid rows (ordering-agnostic)
        statdist=ones(znum^M,1)/(znum^M);
        for i_c=1:10000
            sn=pi_z'*statdist;
            if max(abs(sn-statdist))<10^(-14), statdist=sn; break; end
            statdist=sn;
        end
        statdist=statdist/sum(statdist);
        zvals=z_grid;
        mz=(statdist'*zvals)';
        dz=zvals-mz';
        Vz=(dz.*statdist)'*dz;
        Ez1=pi_z*zvals;
        Cz=((zvals-mz').*statdist)'*(Ez1-mz');
        acz=diag(Cz)./diag(Vz);
        err_mean(c_c,m_c)=max(abs(mz-zmeanT));
        err_var(c_c,m_c)=max(abs(diag(Vz)-diag(SigmaSqzT)));
        err_ac(c_c,m_c)=max(abs(acz-acT));
        fprintf('M=%i znum=%i: mean error [T2] %2.3e, variance error [T2] %2.3e, autocorrelation error [T2] %2.3e \n',M,znum,err_mean(c_c,m_c),err_var(c_c,m_c),err_ac(c_c,m_c))
        if M>1
            err_cross(c_c,m_c)=abs(Vz(1,2)-SigmaSqzT(1,2));
            fprintf('M=%i znum=%i: the CROSS-covariance error [T2] %2.3e \n',M,znum,err_cross(c_c,m_c))
        end
        fprintf('M=%i znum=%i: runtime %2.6f s \n',M,znum,runtime(c_c,m_c))
    end
end

output.err_mean=err_mean; output.err_var=err_var; output.err_ac=err_ac;
output.err_cross=err_cross; output.fallback=fallback;
output.runtime=runtime; output.skipped=skipped; output.znums=znums; output.Mlist=Mlist;

%% The option grid: method x nMoments, at a fixed grid size
% P2 and P3 both found the entropy solve degrading as more moments are demanded, and doing WORSE at
% nMoments=4 than at 2 when the fallback rate is high. This is the same measurement for the VAR
% case, where the solve runs once per variable per row and so has more chances to fail.
fprintf('\n--- option grid: method x nMoments, at M=2, znum=9 --- \n')
methodlist={'even','quantile','gauss-hermite'};
nMomentslist=[1,2,4];
Mew=calib.full.Mew; Rho=calib.full.Rho; SigmaSq=calib.full.SigmaSq;
SigmaSqzT=calib.full.SigmaSqz;
errgrid=NaN(3,3); crossgrid=NaN(3,3); fbgrid=NaN(3,3);
for m_c=1:3
    for n_c=1:3
        opts=struct(); opts.method=methodlist{m_c}; opts.nMoments=nMomentslist(n_c);
        opts.parallel=1; opts.verbose=0;
        [zg,pz,oo]=discretizeVAR1_FarmerToda(Mew,Rho,SigmaSq,9,opts);
        zg=gather(zg); pz=gather(pz);
        sd=ones(9^2,1)/(9^2);
        for i_c=1:10000
            sn=pz'*sd;
            if max(abs(sn-sd))<10^(-14), sd=sn; break; end
            sd=sn;
        end
        sd=sd/sum(sd);
        mz=(sd'*zg)'; dz=zg-mz'; Vz=(dz.*sd)'*dz;
        errgrid(m_c,n_c)=max(abs(diag(Vz)-diag(SigmaSqzT)));
        crossgrid(m_c,n_c)=abs(Vz(1,2)-SigmaSqzT(1,2));
        fbgrid(m_c,n_c)=mean(oo.nMoments_grid(:)<nMomentslist(n_c));
        fprintf('method=%s, nMoments=%i: variance error %2.3e, cross-covariance error %2.3e, fallback %2.1f%% \n',methodlist{m_c},nMomentslist(n_c),errgrid(m_c,n_c),crossgrid(m_c,n_c),100*fbgrid(m_c,n_c))
    end
end
% Reported, not asserted. Farmer & Toda recommend 'even' for persistent processes and
% 'gauss-hermite' otherwise, and the command's own default follows that rule; the paper makes no
% claim that more moments is always better, and P2 and P3 have both measured cases where it is not.
[~,bi]=min(errgrid(:)); [bm,bn]=ind2sub([3,3],bi);
fprintf('the most accurate combination here is method=%s at nMoments=%i (reported, not asserted) \n',methodlist{bm},nMomentslist(bn))
if any(eig(Rho)>0.8)
    defaultmethod='even';
else
    defaultmethod='gauss-hermite';
end
fprintf('the command defaults to ''even'' when any eigenvalue of Rho exceeds 0.8, which for this \n')
fprintf('   calibration (max |eig(Rho)| = %2.4f) means it defaults to ''%s'' \n',max(abs(eig(Rho))),defaultmethod)
output.errgrid=errgrid; output.crossgrid=crossgrid; output.fbgrid=fbgrid;
output.methodlist=methodlist; output.nMomentslist=nMomentslist;


%% REGRESSION BARS ON THE ACCURACY NUMBERS ABOVE
% The moment errors in the sweep above are printed with a value but no verdict, so a defect can sit
% in the output while the run summary reports a clean pass - which is what happened with B30 on
% 2026-08-28. These are regression bars, not accuracy claims: a single threshold cannot serve both
% a coarse grid, where a large error is legitimate, and a fine one. So each family gets a loose bar
% on the worst of the sweep and a tight one at the finest grid. Both come from the 2026-08-28 run
% with about an order of magnitude of headroom.
fprintf('worst mean error over the sweep [T2], this should be below %g: %2.3e \n',1e-06,max(err_mean(:)))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',1e-06,max(err_mean(end,:)))
fprintf('worst variance error over the sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_var(:)))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',1e-06,max(err_var(end,:)))
fprintf('worst cross-covariance error over the sweep [T2], this should be below %g: %2.3e \n',0.1,max(err_cross(:)))
fprintf('cross-covariance error at the finest grid [T2], this should be below %g: %2.3e \n',1e-06,max(err_cross(end,:)))

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(znums,max(err_var(:,2),10^(-18)),'o-',znums,max(err_ac(:,2),10^(-18)),'s-',znums,max(err_cross(:,2),10^(-18)),'d-')
xlabel('znum'); ylabel('error'); title('discretizeVAR1\_FarmerToda, M=2')
legend({'variance','autocorrelation','cross-covariance'},'Location','best')
subplot(1,2,2)
imagesc(log10(max(errgrid,10^(-18)))); colorbar
set(gca,'XTick',1:3,'XTickLabel',nMomentslist,'YTick',1:3,'YTickLabel',methodlist)
xlabel('nMoments'); title('log10 variance error, M=2, znum=9')

end
