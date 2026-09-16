function output=DiscP5_comparator(calib,znums,figure_c)
% P5: two independent AR(1) discretizations, kron'd together (misspecified on purpose)
%
% This is what a user gets by discretizing each variable of a VAR separately and taking the product
% chain: the diagonal of Rho and the diagonal of SigmaSq are used, the off-diagonals are discarded.
% It should hit each variable's own variance and autocorrelation about as well as any AR(1) method
% does, and miss the cross-covariance entirely - it is zero by construction. That gap IS the reason
% the VAR commands exist, so it is measured rather than asserted from theory.
%
% Note the comparator is misspecified in TWO ways at once, and they are not the same thing:
%   (i) it throws away the off-diagonal of Rho, so each variable's own dynamics are wrong too - the
%       truth for variable 1 is not an AR(1) with coefficient Rho(1,1);
%  (ii) it throws away the off-diagonal of SigmaSq, so the innovations are independent.
% Only (ii) forces the cross-covariance to zero. (i) means even the marginal moments are off, which
% is why the variance errors below do not go to zero as the grid is refined.

fprintf('\n========== P5: two independent AR(1)s, kron''d (misspecified) ========== \n')

output=struct();
Rho=calib.full.Rho; SigmaSq=calib.full.SigmaSq; Mew=calib.full.Mew;
SigmaSqzT=calib.full.SigmaSqz; zmeanT=calib.full.zmean; acT=calib.full.autocorr;
nz=length(znums);
nmlist={'AR1_Tauchen','AR1_Rouwenhorst','AR1_FarmerToda'};
err_var=NaN(3,nz); err_ac=NaN(3,nz); err_cross=NaN(3,nz);

fprintf('the truth has cross-covariance %2.6f (cross-correlation %2.4f); this comparator gives zero \n',SigmaSqzT(1,2),calib.full.crosscorr)

for c_c=1:nz
    znum=znums(c_c);
    if znum^2>calib.Ncap
        fprintf('znum=%i: SKIPPED, znum^2=%i exceeds the ceiling of %i \n',znum,znum^2,calib.Ncap)
        continue
    end
    for m_c=1:3
        P=cell(1,2); g=cell(1,2);
        for v_c=1:2
            % each variable on its own, as an AR(1) with its own diagonal coefficients
            mew_v=Mew(v_c); rho_v=Rho(v_c,v_c); sigma_v=sqrt(SigmaSq(v_c,v_c));
            switch m_c
                case 1
                    [g{v_c},P{v_c}]=discretizeAR1_Tauchen(mew_v,rho_v,sigma_v,znum,3,struct());
                case 2
                    [g{v_c},P{v_c}]=discretizeAR1_Rouwenhorst(mew_v,rho_v,sigma_v,znum,struct());
                case 3
                    [g{v_c},P{v_c}]=discretizeAR1_FarmerToda(mew_v,rho_v,sigma_v,znum,struct('verbose',0));
            end
            g{v_c}=gather(g{v_c}); P{v_c}=gather(P{v_c});
        end
        % Variable 1 fastest, matching CreateGridvals, so the product chain is kron(P2,P1)
        pi_z=kron(P{2},P{1});
        zvals=CreateGridvals([znum;znum],[g{1};g{2}],1);
        sd=ones(znum^2,1)/(znum^2);
        for i_c=1:10000
            sn=pi_z'*sd;
            if max(abs(sn-sd))<10^(-14), sd=sn; break; end
            sd=sn;
        end
        sd=sd/sum(sd);
        mz=(sd'*zvals)'; dz=zvals-mz'; Vz=(dz.*sd)'*dz;
        Ez1=pi_z*zvals; Cz=((zvals-mz').*sd)'*(Ez1-mz');
        acz=diag(Cz)./diag(Vz);
        err_var(m_c,c_c)=max(abs(diag(Vz)-diag(SigmaSqzT)));
        err_ac(m_c,c_c)=max(abs(acz-acT));
        err_cross(m_c,c_c)=abs(Vz(1,2)-SigmaSqzT(1,2));
        fprintf('%s, znum=%i: variance error %2.3e, autocorrelation error %2.3e, cross-covariance %2.3e (truth %2.6f) \n',nmlist{m_c},znum,err_var(m_c,c_c),err_ac(m_c,c_c),Vz(1,2),SigmaSqzT(1,2))
    end
end

% The one assertion. Refining the grid cannot fix a misspecification, and the cross-covariance is
% zero by construction here, so the error must be exactly |truth| at every grid size for all three.
last=find(~isnan(err_cross(1,:)),1,'last');
fprintf('the cross-covariance error is |truth| at every grid size, for all three [T1], this should be zero: %2.8e \n',max(abs(err_cross(:,1:last)-abs(SigmaSqzT(1,2))),[],'all'))
fprintf('refining the grid does not reduce the cross-covariance error for any of the three [T2], this should be one: %i \n',all(err_cross(:,last)>=err_cross(:,1)-10^(-12)))

output.err_var=err_var; output.err_ac=err_ac; output.err_cross=err_cross;
output.nmlist=nmlist; output.znums=znums;

%% Figure
figure(figure_c)
semilogy(znums,max(err_var(1,:),10^(-18)),'o-',znums,max(err_var(2,:),10^(-18)),'s-',znums,max(err_var(3,:),10^(-18)),'d-')
xlabel('znum'); ylabel('variance error'); title('independent AR(1)s on a VAR(1) (misspecified)')
legend(nmlist,'Location','best','Interpreter','none')

end
