function output=DiscP7_comparator(calib,znums,figure_c)
% P7: the gaussian life-cycle method on a mixture process (misspecified on purpose)
%
% THIS IS THE BLOCK'S HEADLINE. discretizeLifeCycleAR1_KFTT is given the mixture's first two moments
% age by age - the right mean and the right variance of the innovation, as a normal - and asked to
% discretize. It should reproduce the mean and variance profiles about as well as the mixture
% methods do, and miss the skewness and excess kurtosis profiles at EVERY age, because a normal has
% neither. That gap is Kirkby (2025)'s point in life-cycle form, and P3 measured the stationary
% version of it.
%
% What refining the grid does here is the thing to watch: it converges, to a process whose third and
% fourth cumulants are zero. More grid points cannot manufacture a skewness the innovation does not
% have, so the skewness error converges to |truth| rather than to zero.

fprintf('\n========== P7: the gaussian life-cycle method on a mixture process (misspecified) ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho;
mewzT=calib.vary.mewz; varzT=calib.vary.varz; skewT=calib.vary.skewz; exkurtT=calib.vary.exkurtz;
% the gaussian innovation that matches the mixture's first two moments at each age
emean=calib.vary.emean; evar=calib.vary.evar;
sigma_matched=sqrt(evar);
mew_matched=mew+emean; % fold the mixture's mean into the intercept, since a normal here is mean zero
fprintf('the mixture is matched by a normal with the same mean and variance at every age \n')
fprintf('truth: skewness runs %+2.4f to %+2.4f and excess kurtosis %2.4f to %2.4f; a normal has zero of both \n',skewT(1),skewT(J),exkurtT(1),exkurtT(J))

nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_skew=zeros(1,nz); err_exkurt=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    ko=struct(); ko.verbose=0;
    [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(mew_matched,rho,sigma_matched,znum,J,ko);
    z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);
    [mcmean,mcvar,~,mcdist]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    skj=zeros(1,J); ekj=zeros(1,J);
    for j_c=1:J
        d=gather(mcdist(:,j_c)); d=d/sum(d);
        g=z_grid_J(:,j_c);
        m1=sum(d.*g); v=sum(d.*(g-m1).^2);
        skj(j_c)=sum(d.*(g-m1).^3)/v^1.5;
        ekj(j_c)=sum(d.*(g-m1).^4)/v^2-3;
    end
    err_mean(c_c)=max(abs(mcmean(:)'-mewzT));
    err_var(c_c)=max(abs(mcvar(:)'-varzT));
    err_skew(c_c)=max(abs(skj-skewT));
    err_exkurt(c_c)=max(abs(ekj-exkurtT));
    fprintf('znum=%3i: worst-age errors - mean %2.3e, variance %2.3e, SKEWNESS %2.3e, EXCESS KURTOSIS %2.3e \n',znum,err_mean(c_c),err_var(c_c),err_skew(c_c),err_exkurt(c_c))
    if c_c==nz
        output.skew=skj; output.exkurt=ekj;
    end
end

% The two claims, and they are different. Matching two moments is enough for two moments, so the
% variance should be good. It is not enough for the next two, and no amount of grid is.
fprintf('\nthe variance profile is still reproduced well [T2], this should be one: %i \n',err_var(end)<0.05*max(varzT))
fprintf('the skewness error is close to |truth| - the matched normal has none [T2], this should be one: %i \n',err_skew(end)>0.5*max(abs(skewT)))
fprintf('refining the grid does not reduce the skewness error [T2], this should be one: %i \n',err_skew(end)>0.5*err_skew(1))
fprintf('   skewness error over the sweep:');
for c_c=1:nz
    fprintf(' %2.2e',err_skew(c_c));
end
fprintf(' \n   against a true skewness of up to %2.4f in absolute value \n',max(abs(skewT)))

output.err_mean=err_mean; output.err_var=err_var; output.err_skew=err_skew; output.err_exkurt=err_exkurt;
output.znums=znums;

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(1:J,skewT,'k-','LineWidth',1.5); hold on
plot(1:J,output.skew,'r--','LineWidth',1.2); hold off
xlabel('age j'); ylabel('skewness'); title('a matched normal cannot reach the skewness')
legend({'truth (mixture)','gaussian KFTT, moments matched'},'Location','best')
subplot(1,2,2)
plot(1:J,exkurtT,'k-','LineWidth',1.5); hold on
plot(1:J,output.exkurt,'r--','LineWidth',1.2); hold off
xlabel('age j'); ylabel('excess kurtosis'); title('nor the excess kurtosis')

end
