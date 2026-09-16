function output=DiscP6_compare(calib,znums,outputP6,figure_c)
% P6: method comparison, life-cycle AR(1) with normal innovations
%
% The comparison is on the VARIANCE PROFILE, which is what age-dependence is about: the truth rises
% with age, and a method that got the average level right while missing the shape would be no use.
% The number reported is the worst age, because an average over ages would hide exactly the
% end-of-life-cycle failure the comparator has.

fprintf('\n========== P6: method comparison, life-cycle AR(1) with normal innovations ========== \n')

output=struct();
J=calib.J;
sigmazT=calib.vary.sigmaz;
fprintf('truth: the variance of z rises from %2.6f at age 1 to %2.6f at age %i \n',sigmazT(1)^2,sigmazT(J)^2,J)

nz=length(znums);
eK=outputP6.LCAR1_KFTT.err_var;
eF=outputP6.LCAR1_FGP.err_var;
eT=outputP6.LCAR1_FGPTauchen.err_var;
eC=outputP6.comparator.err_var;   % 3-by-nz

for c_c=1:nz
    fprintf('znum=%3i: worst-age variance error - KFTT %2.4e, FGP %2.4e, FGP-Tauchen %2.4e, best stationary %2.4e \n',znums(c_c),eK(c_c),eF(c_c),eT(c_c),min(eC(:,c_c)))
end

% Every age-dependent method must beat the misspecified stationary one. This is not a claim about
% which of the three is best - it is the claim that modelling age-dependence at all beats not
% modelling it, which is the only thing the comparator is here to establish.
bestC=min(eC,[],1);
fprintf('\nKFTT beats the best stationary comparator at every znum [T2], this should be one: %i \n',all(eK<bestC))
fprintf('FGP beats the best stationary comparator at every znum [T2], this should be one: %i \n',all(eF<bestC))
fprintf('FGP-Tauchen beats the best stationary comparator at every znum [T2], this should be one: %i \n',all(eT<bestC))

% The two Fella-Gallipoli-Pan commands both print, on every call, that they are "typically inferior
% to the KFTT method" and that KFTT is "strongly recommended". Unlike the comment in
% discretizeAR1wGM_Tauchen, that one is a claim this bank can check, so it is checked - and
% reported rather than asserted, because the calibration here is one calibration and the claim is
% about the general case.
fprintf('\nthe two FGP commands both recommend KFTT instead. On this calibration, at znum=%i: \n',znums(end))
fprintf('   KFTT %2.4e, FGP %2.4e, FGP-Tauchen %2.4e (worst-age variance error) \n',eK(end),eF(end),eT(end))
fprintf('   KFTT is the most accurate of the three: %i (reported, not asserted - one calibration) \n',(eK(end)<=eF(end))&&(eK(end)<=eT(end)))

output.eK=eK; output.eF=eF; output.eT=eT; output.bestC=bestC;

%% Figure
figure(figure_c)
semilogy(znums,max(eK,10^(-18)),'o-',znums,max(eF,10^(-18)),'s-',znums,max(eT,10^(-18)),'d-',znums,max(bestC,10^(-18)),'k--')
xlabel('znum'); ylabel('worst-age variance error'); title('life-cycle AR(1): the three methods and a stationary comparator')
legend({'KFTT','FGP','FGP-Tauchen','best stationary'},'Location','best')

end
