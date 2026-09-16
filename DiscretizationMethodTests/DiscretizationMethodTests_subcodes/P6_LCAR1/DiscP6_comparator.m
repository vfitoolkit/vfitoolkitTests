function output=DiscP6_comparator(calib,znums,figure_c)
% P6: a stationary AR(1) held fixed across ages (misspecified on purpose)
%
% This is what a user gets by ignoring age-dependence: discretize once, at the average parameters,
% and use the same grid and transition matrix at every age. It should be about right in the middle
% of the life cycle, where the average parameters are closest to the actual ones, and wrong at both
% ends - and it can never reproduce the RISING variance profile, because a stationary chain
% started from its own stationary distribution has a flat one. That gap is the reason the
% age-dependent methods exist.
%
% Note what is NOT claimed here. The comparator is not a bad discretization of anything - it is an
% exact discretization of the wrong process. Refining znum makes it converge, just to the wrong
% answer, and the checks below say so rather than treating it as an accuracy failure.

fprintf('\n========== P6: a stationary AR(1) at the average parameters (misspecified) ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho; sigma=calib.vary.sigma;
sigmazT=calib.vary.sigmaz; mewzT=calib.vary.mewz;
mbar=mean(mew); rbar=mean(rho); sbar=mean(sigma);
statsigmaz=sbar/sqrt(1-rbar^2); statmewz=mbar/(1-rbar);
fprintf('average parameters: mew=%+2.4f, rho=%2.4f, sigma=%2.4f, giving a stationary sd of %2.4f and mean %+2.4f \n',mbar,rbar,sbar,statsigmaz,statmewz)
fprintf('the truth has a sd rising from %2.4f at age 1 to %2.4f at age %i, so a flat profile cannot match both ends \n',sigmazT(1),sigmazT(J),J)

nz=length(znums);
nmlist={'AR1_Tauchen','AR1_Rouwenhorst','AR1_FarmerToda'};
err_var=zeros(3,nz); err_mean=zeros(3,nz); err_var_mid=zeros(3,nz);
jmid=round(J/2);

for c_c=1:nz
    znum=znums(c_c);
    for m_c=1:3
        switch m_c
            case 1
                [zg,pz]=discretizeAR1_Tauchen(mbar,rbar,sbar,znum,3,struct());
            case 2
                [zg,pz]=discretizeAR1_Rouwenhorst(mbar,rbar,sbar,znum,struct());
            case 3
                [zg,pz]=discretizeAR1_FarmerToda(mbar,rbar,sbar,znum,struct('verbose',0));
        end
        zg=gather(zg); pz=gather(pz);
        % use it as a life-cycle chain: same grid and transition at every age, started from its own
        % stationary distribution
        z_grid_J=repmat(zg(:),1,J);
        pi_z_J=repmat(pz,1,1,J-1);
        [~,~,~,statdist]=MarkovChainMoments(zg,pz);
        jequaloneDistz=statdist(:);
        [mcmean,mcvar,~]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
        err_mean(m_c,c_c)=max(abs(mcmean(:)'-mewzT));
        err_var(m_c,c_c)=max(abs(mcvar(:)'-sigmazT.^2));
        err_var_mid(m_c,c_c)=abs(mcvar(jmid)-sigmazT(jmid)^2);
        if c_c==nz
            output.profile(m_c).var=mcvar(:)';
            output.profile(m_c).mean=mcmean(:)';
        end
    end
    fprintf('znum=%3i: worst-age variance error - Tauchen %2.3e, Rouwenhorst %2.3e, FarmerToda %2.3e \n',znum,err_var(1,c_c),err_var(2,c_c),err_var(3,c_c))
end

% THE ASSERTION, AND WHY IT IS NOT "REFINING DOES NOT HELP". The first version of this demanded
% that the worst-age variance error not fall as znum grows. It failed on the run of 2026-08-26,
% because AR1_Tauchen runs 6.34e-01 down to 3.62e-01 - a real reduction. But that reduction is
% Tauchen converging to ITS OWN stationary answer, from a coarse-grid error towards the pure
% misspecification error, not towards the truth. Rouwenhorst and FarmerToda are already exact at
% znum=5, so they sit flat on that floor from the start.
%
% The right claim is therefore about the FLOOR, not the direction: refining converges to a positive
% error and stays there. Two things are asserted. The error at the largest grid is bounded away
% from zero - it must be, because the process being discretized is the wrong one - and the last two
% grid sizes agree, which is what "converged, to the wrong answer" means.
lastcol=err_var(:,end); prevcol=err_var(:,end-1);
fprintf('\nthe error converges to a POSITIVE floor rather than to zero [T2], this should be one: %i \n',all(lastcol>0.1*max(sigmazT.^2)))
fprintf('and it has converged: the last two grid sizes agree to within 5%% [T2], this should be one: %i \n',all(abs(lastcol-prevcol)<0.05*lastcol))
for m_c=1:3
    fprintf('   %-16s worst-age variance error %2.4e at znum=%i, %2.4e at znum=%i \n',nmlist{m_c},err_var(m_c,1),znums(1),err_var(m_c,end),znums(end))
end
fprintf('   (the floor is about %2.3f, against a true variance peaking at %2.3f - so the misspecified \n',mean(lastcol),max(sigmazT.^2))
fprintf('    chain is wrong by most of the quantity it is trying to reproduce, at every grid size) \n')
% and the mid-life check, which is the honest version of "it is not useless"
fprintf('at mid-life (age %i) the best of the three is off by %2.3e, against %2.3e at the worst age \n',jmid,min(err_var_mid(:,end)),max(err_var(:,end)))

output.err_var=err_var; output.err_mean=err_mean; output.err_var_mid=err_var_mid;
output.nmlist=nmlist; output.znums=znums; output.jmid=jmid;

%% Figure
figure(figure_c)
plot(1:J,sigmazT.^2,'k-','LineWidth',1.5); hold on
for m_c=1:3
    plot(1:J,output.profile(m_c).var,'--')
end
hold off
xlabel('age j'); ylabel('variance of z'); title('a stationary AR(1) used across ages (misspecified)')
legend([{'truth'},nmlist],'Location','best','Interpreter','none')

end
