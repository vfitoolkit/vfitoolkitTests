function output=DiscP7_compare(calib,znums,outputP7,figure_c)
% P7: method comparison, life-cycle AR(1) with gaussian-mixture innovations
%
% On skewness and excess kurtosis, not variance. Every method here gets the variance about right -
% including the misspecified one, which is handed the right variance by construction - so comparing
% on it would say nothing. The higher moments are what separate them.

fprintf('\n========== P7: method comparison, life-cycle AR(1) with gaussian-mixture innovations ========== \n')

output=struct();
J=calib.J;
fprintf('truth: skewness runs %+2.4f to %+2.4f, excess kurtosis %2.4f to %2.4f \n',calib.vary.skewz(1),calib.vary.skewz(J),calib.vary.exkurtz(1),calib.vary.exkurtz(J))

nz=length(znums);
sK=outputP7.LCAR1wGM_KFTT.err_skew;   kK=outputP7.LCAR1wGM_KFTT.err_exkurt;
sT=outputP7.LCAR1wGM_Tauchen.err_skew; kT=outputP7.LCAR1wGM_Tauchen.err_exkurt;
sC=outputP7.comparator.err_skew;      kC=outputP7.comparator.err_exkurt;

for c_c=1:nz
    fprintf('znum=%3i: worst-age skewness error        - KFTT %2.3e, Tauchen %2.3e, matched normal %2.3e \n',znums(c_c),sK(c_c),sT(c_c),sC(c_c))
end
fprintf('\n')
for c_c=1:nz
    fprintf('znum=%3i: worst-age excess kurtosis error - KFTT %2.3e, Tauchen %2.3e, matched normal %2.3e \n',znums(c_c),kK(c_c),kT(c_c),kC(c_c))
end

% Both mixture methods must beat the matched normal on the higher moments at the finest grid. This
% is the claim the block exists to make; it is not a claim about which mixture method is better.
fprintf('\nKFTT beats the matched normal on skewness at the finest grid [T2], this should be one: %i \n',sK(end)<sC(end))
fprintf('wGM_Tauchen beats the matched normal on skewness at the finest grid [T2], this should be one: %i \n',sT(end)<sC(end))
fprintf('KFTT beats the matched normal on excess kurtosis at the finest grid [T2], this should be one: %i \n',kK(end)<kC(end))
fprintf('wGM_Tauchen beats the matched normal on excess kurtosis at the finest grid [T2], this should be one: %i \n',kT(end)<kC(end))

% The two mixture methods are NOT ordered against each other here. Both commands print a
% recommendation of KFTT on every call, and unlike most such comments this bank can check it - but
% one calibration is not enough to assert a general ordering, so it is reported.
fprintf('\nthe wGM_Tauchen command recommends KFTT instead. On this calibration, at znum=%i: \n',znums(end))
fprintf('   skewness        KFTT %2.3e vs Tauchen %2.3e; KFTT better: %i \n',sK(end),sT(end),sK(end)<=sT(end))
fprintf('   excess kurtosis KFTT %2.3e vs Tauchen %2.3e; KFTT better: %i \n',kK(end),kT(end),kK(end)<=kT(end))
fprintf('   (reported, not asserted - one calibration cannot establish a general ordering) \n')

output.sK=sK; output.sT=sT; output.sC=sC; output.kK=kK; output.kT=kT; output.kC=kC;

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(znums,max(sK,10^(-18)),'o-',znums,max(sT,10^(-18)),'s-',znums,max(sC,10^(-18)),'k--')
xlabel('znum'); ylabel('worst-age skewness error'); title('skewness')
legend({'wGM\_KFTT','wGM\_Tauchen','matched normal'},'Location','best')
subplot(1,2,2)
semilogy(znums,max(kK,10^(-18)),'o-',znums,max(kT,10^(-18)),'s-',znums,max(kC,10^(-18)),'k--')
xlabel('znum'); ylabel('worst-age excess kurtosis error'); title('excess kurtosis')

end
