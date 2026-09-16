function output=DiscP4_comparator(calib,znums,figure_c)
% P4: the plain AR(1) at the average volatility, as a misspecified comparator
%
% What does a practitioner lose by ignoring the stochastic volatility? The honest comparator is not
% "some AR(1)" but the one they would actually fit: same persistence, and innovation standard
% deviation set to sigmau, the unconditional standard deviation of u. That process has exactly the
% right mean, variance and autocorrelation - and ZERO excess kurtosis, where the truth is positive.
%
% So this is the same shape of comparison as P3's gaussian-versus-mixture one, and it isolates the
% same thing: what the extra machinery buys is the tail behaviour, and nothing else.

fprintf('\n========== P4: the plain AR(1) at the average volatility (misspecified) ========== \n')
fprintf('comparator: rho=%g, sigma=%g. It matches Var(z)=%2.4f and autocorr=%2.2f by construction, \n',calib.rho,calib.sigmau,calib.z.var,calib.z.autocorr)
fprintf('and has zero excess kurtosis where the truth is %2.4f. \n',calib.z.exkurt)

output=struct();
rho=calib.rho; sigmau=calib.sigmau;
nz=length(znums);
nmlist={'AR1_Tauchen','AR1_Rouwenhorst','AR1_FarmerToda'};
err_var=zeros(3,nz); err_ac=zeros(3,nz); err_exkurt=zeros(3,nz);

for c_c=1:nz
    znum=znums(c_c);
    for m_c=1:3
        if m_c==1
            [z_grid,pi_z]=discretizeAR1_Tauchen(0,rho,sigmau,znum,3,struct());
        elseif m_c==2
            [z_grid,pi_z]=discretizeAR1_Rouwenhorst(0,rho,sigmau,znum,struct());
        else
            [z_grid,pi_z,~]=discretizeAR1_FarmerToda(0,rho,sigmau,znum,struct('verbose',0));
        end
        [mcmean,mcvar,mcac,statdist]=MarkovChainMoments(z_grid,pi_z);
        ek=sum(statdist.*(z_grid-mcmean).^4)/mcvar^2-3;
        err_var(m_c,c_c)=abs(mcvar-calib.z.var);
        err_ac(m_c,c_c)=abs(mcac-calib.z.autocorr);
        err_exkurt(m_c,c_c)=abs(ek-calib.z.exkurt);
        fprintf('%s, znum=%i: variance error %2.3e, autocorr error %2.3e, excess kurtosis %2.4f vs truth %2.4f (err %2.3e) \n',nmlist{m_c},znum,err_var(m_c,c_c),err_ac(m_c,c_c),ek,calib.z.exkurt,err_exkurt(m_c,c_c))
    end
end
output.err_var=err_var; output.err_ac=err_ac; output.err_exkurt=err_exkurt;
output.nmlist=nmlist; output.znums=znums;

% Refining the grid cannot help: the deficiency is structural. A homoskedastic gaussian AR(1) has
% no excess kurtosis at any grid size, so its error stays at roughly the truth however many points
% are added. Same statement as P3 makes about the gaussian methods on a mixture.
fprintf('refining the grid does not reduce the kurtosis error for any of the three [T2], this should be one: %i \n',all(err_exkurt(:,end)>0.5*err_exkurt(:,1)))
for m_c=1:3
    fprintf('   %-18s excess kurtosis error %2.4f at znum=%i, %2.4f at znum=%i \n',nmlist{m_c},err_exkurt(m_c,1),znums(1),err_exkurt(m_c,end),znums(end))
end

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(znums,err_exkurt(1,:),'o-',znums,err_exkurt(2,:),'s-',znums,err_exkurt(3,:),'^-')
set(gca,'XScale','log'); title('misspecified comparator: |excess kurtosis error|'); xlabel('znum')
legend(nmlist,'Location','best','Interpreter','none')
subplot(1,2,2)
plot(znums,max(err_var(1,:),10^(-18)),'o-',znums,max(err_var(2,:),10^(-18)),'s-',znums,max(err_var(3,:),10^(-18)),'^-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('|variance error| (which it does match)'); xlabel('znum'); legend(nmlist,'Location','best','Interpreter','none')

end
