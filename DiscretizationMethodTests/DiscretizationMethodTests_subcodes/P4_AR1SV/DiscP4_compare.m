function output=DiscP4_compare(calib,znums,outputP4,figure_c)
% P4: compare the two stochastic-volatility methods against each other and against the
% misspecified plain AR(1)

fprintf('\n========== P4: method comparison, AR(1) with stochastic volatility ========== \n')
fprintf('truth for z: variance %2.6f, autocorrelation %2.4f, excess kurtosis %2.6f \n',calib.z.var,calib.z.autocorr,calib.z.exkurt)

output=struct();
nz=length(znums);
eF=outputP4.AR1wSV_FarmerToda.err_exkurt;
eT=outputP4.AR1wSV_Tauchen.err_exkurt;
sF=outputP4.AR1wSV_FarmerToda.skipped;
eG=outputP4.comparator.err_exkurt;

% Compare along the diagonal, where both dimensions are refined together
dF=nan(1,nz); dT=nan(1,nz);
for c_c=1:nz
    if ~sF(c_c,c_c)
        dF(c_c)=eF(c_c,c_c); dT(c_c)=eT(c_c,c_c);
    end
end
for c_c=1:nz
    if ~isnan(dF(c_c))
        fprintf('xnum=znum=%3i: wSV_FarmerToda kurtosis error %2.4e, wSV_Tauchen %2.4e, best plain AR(1) %2.4e \n',znums(c_c),dF(c_c),dT(c_c),min(eG(:,c_c)))
    end
end

%% The ordering this block exists to establish
% Both SV methods must beat the misspecified comparator on excess kurtosis, because the comparator
% has none at all. Asserted only where the grid is fine enough, following P3: a coarse grid can
% overshoot and land further from the truth than zero does.
enough=(znums>=15) & ~isnan(dF);
fprintf('\nwSV_FarmerToda beats every plain AR(1) on excess kurtosis once znum>=15 [T2], this should be one: %i \n',all(dF(enough)<=min(eG(:,enough),[],1)))
fprintf('wSV_Tauchen beats every plain AR(1) on excess kurtosis once znum>=15 [T2], this should be one: %i \n',all(dT(enough)<=min(eG(:,enough),[],1)))
fprintf('   (stochastic volatility exists to generate fat tails from conditionally normal shocks, \n')
fprintf('    so a homoskedastic approximation cannot reach them at any grid size) \n')

%% NOT asserted: which of the two SV methods is better
% They are not discretizing quite the same object. discretizeAR1wSV_FarmerToda conditions the z
% transition on the EXPECTED next-period volatility E[exp(x_t)|x_{t-1}]; discretizeAR1wSV_Tauchen
% conditions on the REALIZED exp(x_t), which is the exact conditional law. So a difference between
% them is expected, and which is more accurate overall is a question for the numbers rather than a
% claim to assert in advance. Reported.
fprintf('\nNOT asserted: an ordering between the two SV methods. They condition the z transition on \n')
fprintf('   different things - FarmerToda on the expected next-period volatility, Tauchen on the \n')
fprintf('   realized one - so they are not two implementations of one object. At the largest \n')
fprintf('   diagonal grid: FarmerToda %2.4e, Tauchen %2.4e \n',dF(find(~isnan(dF),1,'last')),dT(find(~isnan(dT),1,'last')))

output.dF=dF; output.dT=dT; output.eG=eG;

%% Figure
figure(figure_c)
plot(znums,max(dF,10^(-18)),'ko-','LineWidth',2); hold on
plot(znums,max(dT,10^(-18)),'ks-','LineWidth',2)
for m_c=1:3
    plot(znums,max(eG(m_c,:),10^(-18)),'--')
end
hold off; set(gca,'XScale','log'); set(gca,'YScale','log')
title('|excess kurtosis error| (diagonal of the sweep)'); xlabel('xnum=znum')
legend([{'wSV\_FarmerToda','wSV\_Tauchen'},outputP4.comparator.nmlist],'Location','best','Interpreter','none')

end
