function output=DiscP3_gaussianmethods(calib,znums,figure_c)
% P3: the four gaussian AR(1) methods applied to the gaussian-mixture process
%
% Not swept here - P2 owns their option sweeps - and run at defaults as MISSPECIFIED comparators.
%
% This is the Kirkby (2025) comparison, and the comparator has to be chosen honestly. It is not
% "some gaussian AR(1)", it is the gaussian AR(1) a practitioner would actually fit: same
% unconditional mean, same innovation variance. Anything the mixture methods gain over this is a
% gain from modelling the higher moments, not from being handed a better-calibrated process.
%
% What these methods cannot do is produce ANY skewness or excess kurtosis: their grids are
% symmetric and their transition matrices centrosymmetric, so the stationary distribution is
% symmetric and its third central moment is zero by construction. Their skewness error is
% therefore the truth itself. That is not a defect of the implementations - it is the point of the
% block, and it is why the ordering below is asserted rather than reported.

fprintf('\n========== P3: the gaussian methods on a gaussian-mixture process ========== \n')
fprintf('comparator is the gaussian AR(1) with the same mean and innovation variance: \n')
fprintf('   mew=%g, rho=%g, sigma=%g (against the mixture''s sd of %g) \n',calib.gaussian.mew,calib.gaussian.rho,calib.gaussian.sigma,sqrt(calib.e.var))
fprintf('truth for z: variance %2.4f, skewness %2.4f, excess kurtosis %2.4f \n',calib.z.var,calib.z.skew,calib.z.exkurt)

output=struct();
mewG=calib.gaussian.mew; rhoG=calib.gaussian.rho; sigG=calib.gaussian.sigma;
nz=length(znums);
nmlist={'AR1_Tauchen','AR1_Rouwenhorst','AR1_TauchenHussey','AR1_FarmerToda'};
err_var=zeros(4,nz); err_skew=zeros(4,nz); err_exkurt=zeros(4,nz); runtime=zeros(4,nz);

for c_c=1:nz
    znum=znums(c_c);
    for m_c=1:4
        tic;
        if m_c==1
            [z_grid,pi_z]=discretizeAR1_Tauchen(mewG,rhoG,sigG,znum,3,struct());
        elseif m_c==2
            [z_grid,pi_z]=discretizeAR1_Rouwenhorst(mewG,rhoG,sigG,znum,struct());
        elseif m_c==3
            [z_grid,pi_z]=discretizeAR1_TauchenHussey(mewG,rhoG,sigG,znum,struct());
        else
            [z_grid,pi_z,~]=discretizeAR1_FarmerToda(mewG,rhoG,sigG,znum,struct('verbose',0));
        end
        runtime(m_c,c_c)=toc;

        [mcmean,mcvar,~,statdist]=MarkovChainMoments(z_grid,pi_z);
        sk=sum(statdist.*(z_grid-mcmean).^3)/mcvar^1.5;
        ek=sum(statdist.*(z_grid-mcmean).^4)/mcvar^2-3;
        err_var(m_c,c_c)=abs(mcvar-calib.z.var);
        err_skew(m_c,c_c)=abs(sk-calib.z.skew);
        err_exkurt(m_c,c_c)=abs(ek-calib.z.exkurt);

        % These methods produce a symmetric chain, so their skewness is zero by construction, not
        % by accident. Assert that, because if it were ever non-zero the comparison below would be
        % measuring something other than what it claims.
        fprintf('%s, znum=%i: skewness of the fitted gaussian chain is zero [T2], this should be small: %2.3e \n',nmlist{m_c},znum,abs(sk))
        fprintf('%s, znum=%i: so its skewness error IS the truth, %2.4f, and its excess kurtosis error is %2.3e \n',nmlist{m_c},znum,err_skew(m_c,c_c),err_exkurt(m_c,c_c))
    end
end
output.err_var=err_var; output.err_skew=err_skew; output.err_exkurt=err_exkurt;
output.runtime=runtime; output.nmlist=nmlist; output.znums=znums;

% Adding grid points cannot help: the deficiency is structural, not a resolution problem. This is
% the sharpest way to say what the mixture methods are for.
fprintf('refining the grid does not reduce the skewness error for any gaussian method [T2], this should be one: %i \n',all(err_skew(:,end)>0.5*err_skew(:,1)))
fprintf('   (skewness error at znum=%i vs znum=%i, by method: \n',znums(1),znums(end))
for m_c=1:4
    fprintf('    %-18s %2.4f -> %2.4f \n',nmlist{m_c},err_skew(m_c,1),err_skew(m_c,end))
end

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(znums,err_skew(1,:),'o-',znums,err_skew(2,:),'s-',znums,err_skew(3,:),'^-',znums,err_skew(4,:),'v-')
set(gca,'XScale','log')
title('gaussian methods: |skewness error| (flat by construction)'); xlabel('znum'); legend(nmlist,'Location','best','Interpreter','none')
subplot(1,2,2)
plot(znums,max(err_var(1,:),10^(-18)),'o-',znums,max(err_var(2,:),10^(-18)),'s-',znums,max(err_var(3,:),10^(-18)),'^-',znums,max(err_var(4,:),10^(-18)),'v-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('|variance error| (they do target this)'); xlabel('znum'); legend(nmlist,'Location','best','Interpreter','none')

end
