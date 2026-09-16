function output=DiscP3_compare(calib,znums,outputP3,figure_c)
% P3: compare the mixture methods against the gaussian ones, and against analytic truth
%
% Truth comes from the cumulants of z, which are exact: the n-th cumulant of z is the n-th cumulant
% of the innovation divided by (1-rho^n). No simulation, no threshold read off a previous run.

fprintf('\n========== P3: method comparison, AR(1) with gaussian-mixture innovations ========== \n')
fprintf('truth for z: variance %2.6f, skewness %2.6f, excess kurtosis %2.6f \n',calib.z.var,calib.z.skew,calib.z.exkurt)

output=struct();
eGM_sk=outputP3.AR1wGM_FarmerToda.err_skew;
eGM_ek=outputP3.AR1wGM_FarmerToda.err_exkurt;
eG_sk=outputP3.gaussianmethods.err_skew;
eG_ek=outputP3.gaussianmethods.err_exkurt;
nmlist=outputP3.gaussianmethods.nmlist;

for c_c=1:length(znums)
    fprintf('znum=%3i: wGM_FarmerToda skewness error %2.4e, best gaussian method %2.4e \n',znums(c_c),eGM_sk(c_c),min(eG_sk(:,c_c)))
end

%% The ordering this whole block exists to establish
fprintf('\nwGM_FarmerToda beats EVERY gaussian method on skewness, at EVERY znum [T2], this should be one: %i \n',all(all(eGM_sk<=eG_sk)))
% Kurtosis needs a fine enough grid before the mixture machinery pays off, and that is a finding
% rather than a caveat. At znum=5 the mixture method cannot represent a fat-tailed conditional
% distribution at all: measured, it returns excess kurtosis of -0.38 against a truth of +0.68, so
% its error (1.06) is WORSE than a gaussian method's (0.68-0.79), which simply returns about zero.
% From znum=15 it wins, and by znum=31 it wins by four orders of magnitude.
% Skewness has no such threshold - the mixture method wins there even at znum=5 - because the
% gaussian chains are exactly symmetric and so get the skewness maximally wrong at any grid size.
enough=(znums>=15);
fprintf('wGM_FarmerToda beats EVERY gaussian method on excess kurtosis, once znum>=15 [T2], this should be one: %i \n',all(all(eGM_ek(enough)<=eG_ek(:,enough))))
fprintf('   and at znum=%i it does NOT: wGM_FarmerToda %2.3e against the best gaussian method %2.3e \n',znums(1),eGM_ek(1),min(eG_ek(:,1)))
fprintf('   (a coarse grid overshoots: the mixture method returns a platykurtic chain where the \n')
fprintf('    gaussian ones return an almost mesokurtic one, and the truth is leptokurtic) \n')
fprintf('   (Kirkby 2025 is about exactly this: the gaussian methods cannot represent the higher \n')
fprintf('    moments at any grid size, because their chains are symmetric by construction, so the \n')
fprintf('    gap is structural rather than a matter of resolution) \n')

%% And where the gaussian methods are NOT at a disadvantage
% Both target the variance, so no ordering is predicted there and none is asserted. Reporting it
% is the honest thing, and it is also the useful thing: it says what the mixture machinery costs
% you on the moments the simpler methods already handle.
eGM_v=outputP3.AR1wGM_FarmerToda.err_var;
eG_v=outputP3.gaussianmethods.err_var;
fprintf('\nNOT asserted: any ordering on the variance. Both target it, so neither should dominate. \n')
fprintf('   at znum=%i: wGM_FarmerToda %2.3e, and by gaussian method: \n',znums(end),eGM_v(end))
for m_c=1:4
    fprintf('      %-18s %2.3e \n',nmlist{m_c},eG_v(m_c,end))
end

output.eGM_sk=eGM_sk; output.eG_sk=eG_sk; output.eGM_ek=eGM_ek; output.eG_ek=eG_ek;

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(znums,max(eGM_sk,10^(-18)),'ko-','LineWidth',2); hold on
for m_c=1:4
    plot(znums,max(eG_sk(m_c,:),10^(-18)),'--')
end
hold off; set(gca,'XScale','log'); set(gca,'YScale','log')
title('|skewness error|'); xlabel('znum'); legend([{'wGM\_FarmerToda'},nmlist],'Location','best','Interpreter','none')
subplot(1,2,2)
plot(znums,max(eGM_ek,10^(-18)),'ko-','LineWidth',2); hold on
for m_c=1:4
    plot(znums,max(eG_ek(m_c,:),10^(-18)),'--')
end
hold off; set(gca,'XScale','log'); set(gca,'YScale','log')
title('|excess kurtosis error|'); xlabel('znum')

end
