function output=DiscP2_floden(calib,figure_c)
% P2: reproduce Floden (2008) Table 1 - an EXTERNAL oracle
%
% Every other accuracy check in this bank compares the toolkit against a formula this bank
% computes for itself. This one compares it against numbers that an independent implementation
% published seventeen years ago, which is a different and stronger kind of evidence: it tests the
% toolkit against the literature rather than against itself.
%
% Floden, Martin (2008), "A note on the accuracy of Markov-chain approximations to highly
% persistent AR(1) processes", Economics Letters 99(3), 516-520.
%
% Two things have to line up before the numbers can even be compared:
%   1. Floden spaces his Tauchen nodes over +-1.2*sigma_z*log(n), not over +-Tauchen_q*sigma_z. So
%      the toolkit is called with Tauchen_q=1.2*log(znum). His table's z^n/sigma_z row is the
%      built-in check that this is right, and it is asserted below before anything else.
%   2. His sigma_e is the innovation std dev IMPLIED by the fitted chain, i.e.
%      sigma_z*sqrt(1-rho^2) using the chain's own moments, not an average of conditional std devs.
%
% TRANSCRIPTION STATUS (see DiscSetup_AR1Normal for the detail): the sigma_e and sigma_z rows were
% transcribed from the paper AND independently recomputed from Floden's description of his method,
% and the two agree to all four printed decimals in all nine cells - so those are ASSERTED. The rho
% row is NOT confirmed: the independent recomputation is systematically lower than the
% transcription, by 0.0004 to 0.0248, shrinking as n grows. So rho is REPORTED, with the
% discrepancy printed, until someone checks the paper. Do not promote it to an assertion without
% doing that.

fprintf('\n========== P2: Floden (2008) Table 1 as an external oracle ========== \n')
fprintf('Floden spaces his Tauchen nodes over +-1.2*sigma_z*log(n), so the toolkit is called with \n')
fprintf('Tauchen_q=1.2*log(znum). Table entries are printed to four decimals, so the tolerance is %g. \n',calib.floden.tol)

output=struct();
tol=calib.floden.tol;
nlist=calib.floden.nlist;

%% First: the spacing rule itself
for n_c=1:length(nlist)
    n=nlist(n_c);
    % Floden prints this row to four decimals, so the exact 1.2*log(n) can only agree to half a
    % unit in the last place. This is a rounding comparison, not a T1 identity.
    fprintf('spacing rule: 1.2*log(%i)=%2.6f, Floden Table 1 prints %2.4f, this should be below %g: %2.8e \n',n,1.2*log(n),calib.floden.znoversigmaz(n_c),calib.floden.roundtol,abs(1.2*log(n)-calib.floden.znoversigmaz(n_c)))
end
% If those three fail, nothing below means anything: either the transcription is wrong or the
% spacing rule has been misread, and the comparison would be against the wrong grid.

%% Now the table itself
err_se=zeros(3,3); err_se_direct=zeros(3,3); err_sz=zeros(3,3); dev_rho=zeros(3,3);
for p_c=1:3
    rho=calib.floden.rho(p_c);
    sigma=sqrt(calib.floden.sigmasq_e(p_c));
    fprintf('\n--- %s: rho=%g, sigma_eps^2=%g (so sigma_eps=%2.4f, sigma_z=%2.4f) --- \n',calib.floden.name{p_c},rho,calib.floden.sigmasq_e(p_c),sigma,sigma/sqrt(1-rho^2))
    for n_c=1:length(nlist)
        n=nlist(n_c);
        q=1.2*log(n);
        [z_grid,pi_z]=discretizeAR1_Tauchen(0,rho,sigma,n,q,struct());
        [~,mcvar,mcac,statdist]=MarkovChainMoments(z_grid,pi_z);
        sz_hat=sqrt(mcvar);
        rho_hat=mcac;
        se_hat=sz_hat*sqrt(1-rho_hat^2); % the innovation std dev implied by the fitted chain

        % Floden does not say which of two definitions his sigma_eps row uses, and the choice is
        % testable. The one above backs it out of the fitted AR(1); the alternative is the direct
        % one, the stationary-weighted average of the chain's own conditional standard deviations,
        % which is also what his Fig 2 plots node by node. Both are computed and both compared,
        % because whichever reproduces the table is the definition - and that also bears on the rho
        % row, which could not be confirmed from the PDF: if the direct definition matches while the
        % backed-out one does not, then his rho is not cov/var under the stationary distribution.
        condmean_hat=pi_z*z_grid;
        condvar_hat=zeros(n,1);
        for z_c=1:n
            condvar_hat(z_c)=sum(pi_z(z_c,:)'.*((z_grid-condmean_hat(z_c)).^2));
        end
        se_direct=sum(statdist.*sqrt(condvar_hat)); % weighted-average conditional std dev

        err_se(p_c,n_c)=abs(se_hat-calib.floden.tauchen_sigma_e(p_c,n_c));
        err_sz(p_c,n_c)=abs(sz_hat-calib.floden.tauchen_sigma_z(p_c,n_c));
        dev_rho(p_c,n_c)=abs(rho_hat-calib.floden.tauchen_rho_unconfirmed(p_c,n_c));

        fprintf('n=%2i: sigma_e toolkit %2.4f vs Floden %2.4f, this should be below %g: %2.8e \n',n,se_hat,calib.floden.tauchen_sigma_e(p_c,n_c),tol,err_se(p_c,n_c))
        fprintf('n=%2i: sigma_e by the DIRECT definition (weighted-average conditional sd) is %2.4f, difference from Floden %2.3e \n',n,se_direct,abs(se_direct-calib.floden.tauchen_sigma_e(p_c,n_c)))
        err_se_direct(p_c,n_c)=abs(se_direct-calib.floden.tauchen_sigma_e(p_c,n_c));
        fprintf('n=%2i: sigma_z toolkit %2.4f vs Floden %2.4f, this should be below %g: %2.8e \n',n,sz_hat,calib.floden.tauchen_sigma_z(p_c,n_c),tol,err_sz(p_c,n_c))
        fprintf('n=%2i: rho      toolkit %2.4f vs Floden %2.4f (REPORTED not asserted, transcription unconfirmed), difference %2.4f \n',n,rho_hat,calib.floden.tauchen_rho_unconfirmed(p_c,n_c),dev_rho(p_c,n_c))
        % The grid Floden's rule is supposed to produce. Divided by the TRUE sigma_z, not by the
        % chain's implied one: his z^n/sigma_z row is a statement about the grid construction, and
        % dividing by the implied sd instead turns an exact identity into a number that drifts with
        % the discretization error (it read 2.78 against 3.25 for STY at n=15 when computed that way,
        % which looked alarming and meant nothing).
        sz_true=sigma/sqrt(1-rho^2);
        fprintf('n=%2i: grid half-width over the TRUE sd(z) is %2.6f, and 1.2*log(n)=%2.6f [T1], this should be zero: %2.8e \n',n,(z_grid(end)-z_grid(1))/2/sz_true,1.2*log(n),abs((z_grid(end)-z_grid(1))/2/sz_true-1.2*log(n)))
    end
end

fprintf('\nsigma_e reproduced across all nine cells, this should be one: %i \n',all(err_se(:)<tol))
fprintf('   which of the two definitions matches Floden better, over the nine cells: \n')
fprintf('   backed out of the fitted AR(1), sigma_z*sqrt(1-rho^2): worst %2.3e, %i of 9 cells within %g \n',max(err_se(:)),sum(err_se(:)<tol),tol)
fprintf('   direct, weighted-average conditional sd:                worst %2.3e, %i of 9 cells within %g \n',max(err_se_direct(:)),sum(err_se_direct(:)<tol),tol)
fprintf('   (if the direct definition wins, Floden''s rho is not cov/var under the stationary \n')
fprintf('    distribution, which would explain the one row of his table that could not be confirmed) \n')
output.err_se_direct=err_se_direct;
fprintf('sigma_z reproduced across all nine cells, this should be one: %i \n',all(err_sz(:)<tol))
fprintf('rho: worst deviation from the (unconfirmed) transcription is %2.4f, at which the transcription \n',max(dev_rho(:)))
fprintf('     or Floden''s definition of the implied autocorrelation needs checking against the paper. \n')
fprintf('     Not a failure of the toolkit: sigma_e above is computed FROM rho, and it matches to 1e-4, \n')
fprintf('     which it could not do if the chain''s autocorrelation were wrong. \n')

output.err_se=err_se;
output.err_sz=err_sz;
output.dev_rho=dev_rho;

%% Floden's conclusion, as an assertion
% Read his pg 518 carefully, because the obvious paraphrase gets it backwards. What he says about
% the sigma_z variant is: "the resulting approximation is relatively accurate for the unconditional
% variance, but typically at the price of a much too high autocorrelation and a too low conditional
% variance." So sigma_z is GOOD on the unconditional variance and BAD on the autocorrelation - the
% opposite way round from what "the sigma_z variant is the bad one" would suggest.
%
% The Storesletten-Telmer-Yaron process, rho=0.98, is his hardest case.
rho=calib.floden.rho(3); sigma=sqrt(calib.floden.sigmasq_e(3)); n=9;
[zgT,pzT]=discretizeAR1_Tauchen(0,rho,sigma,n,1.2*log(n),struct());
[~,vT,aT]=MarkovChainMoments(zgT,pzT);
optsz=struct(); optsz.baseSigma=sigma/sqrt(1-rho^2);
[zgH,pzH]=discretizeAR1_TauchenHussey(0,rho,sigma,n,optsz);
[~,vH,aH]=MarkovChainMoments(zgH,pzH);
truevar=sigma^2/(1-rho^2);
fprintf('\nat rho=%g (Storesletten-Telmer-Yaron), n=%i: \n',rho,n)
fprintf('   unconditional variance: Tauchen %2.4f, TauchenHussey(sigma_z) %2.4f, truth %2.4f \n',vT,vH,truevar)
fprintf('   autocorrelation:        Tauchen %2.4f, TauchenHussey(sigma_z) %2.4f, truth %2.4f \n',aT,aH,rho)
fprintf('TauchenHussey(sigma_z) beats Tauchen on the UNCONDITIONAL VARIANCE, as Floden says [T2], this should be one: %i \n',abs(vH-truevar)<abs(vT-truevar))
fprintf('TauchenHussey(sigma_z) has a TOO HIGH autocorrelation, as Floden says [T2], this should be one: %i \n',aH>rho)
fprintf('   (Floden pg 518: it "chooses nodes far from the mean", which buys the unconditional \n')
fprintf('    variance at the cost of the dynamics) \n')
% And the collapse threshold, which his table does not show because he stops at rho=0.98: at
% rho=0.99 the same variant degenerates completely, its transition probabilities underflowing to
% an absorbing chain (see DiscP2_AR1_TauchenHussey). So the useful range of this variant ends
% somewhere between 0.98 and 0.99.
fprintf('at rho=0.98 the sigma_z variant still works; at rho=0.99 it degenerates entirely, so its \n')
fprintf('   useful range ends between the two. Floden stops at 0.98 and so does not show this. \n')

%% Figure
figure(figure_c)
subplot(1,2,1); bar(err_se'); set(gca,'YScale','log'); set(gca,'XTickLabel',{'n=5','n=9','n=15'})
title('|toolkit - Floden Table 1|: sigma_e'); legend(calib.floden.name,'Location','best','Interpreter','none')
hold on; plot([0.5,3.5],[tol,tol],'k--'); hold off
subplot(1,2,2); bar(err_sz'); set(gca,'YScale','log'); set(gca,'XTickLabel',{'n=5','n=9','n=15'})
title('|toolkit - Floden Table 1|: sigma_z'); legend(calib.floden.name,'Location','best','Interpreter','none')
hold on; plot([0.5,3.5],[tol,tol],'k--'); hold off

end
