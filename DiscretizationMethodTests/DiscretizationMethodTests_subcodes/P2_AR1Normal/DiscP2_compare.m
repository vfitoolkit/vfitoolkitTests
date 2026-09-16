function output=DiscP2_compare(calib,znums,outputP2,figure_c)
% P2: compare the four methods against each other and against analytic truth
%
% Consumes what the method subcodes already computed.
%
% The orderings asserted here are the ones the literature actually predicts, in the regime it
% predicts them for. P1 taught this the hard way: asserting "Farmer-Toda beats Tauchen-Hussey" at
% rho=0 was asserting something neither paper claims, and it failed - correctly. So each ordering
% below names its source and its regime.

fprintf('\n========== P2: method comparison, AR(1) with normal innovations ========== \n')

output=struct();
nmlist={'Tauchen','Rouwenhorst','TauchenHussey','FarmerToda'};

for cal_c=1:3
    cname=calib.names{cal_c};
    mew=calib.(cname).mew; rho=calib.(cname).rho; sigma=calib.(cname).sigma;
    fprintf('\n--- calibration %s: rho=%g, truth: E(z)=%g, Var(z)=%g, autocorr=%g --- \n',cname,rho,mew/(1-rho),sigma^2/(1-rho^2),rho)

    ev=[outputP2.Tauchen.(cname).err_var; outputP2.Rouwenhorst.(cname).err_var; outputP2.TauchenHussey.(cname).err_var; outputP2.FarmerToda.(cname).err_var];
    ea=[outputP2.Tauchen.(cname).err_ac;  outputP2.Rouwenhorst.(cname).err_ac;  outputP2.TauchenHussey.(cname).err_ac;  outputP2.FarmerToda.(cname).err_ac];
    output.(cname).err_var=ev;
    output.(cname).err_ac=ea;

    for m_c=1:4
        fprintf('%-14s variance error at znum=%i: %2.3e, at znum=%i: %2.3e; autocorr error at znum=%i: %2.3e \n',nmlist{m_c},znums(1),ev(m_c,1),znums(end),ev(m_c,end),znums(end),ea(m_c,end))
    end

    % Rouwenhorst matches all three moments algebraically, at every znum and every calibration.
    % So it beats everything on the autocorrelation, everywhere. This is Kopecky & Suen (2010).
    fprintf('%s: Rouwenhorst beats all three others on the autocorrelation, at every znum [T2], this should be one: %i \n',cname,all(all(ea(2,:)<=ea([1,3,4],:))))
end

%% Orderings that are regime-specific, asserted only in their regime
% Farmer & Toda (2017), last paragraph of pg 678: their method outperforms Tauchen for almost all
% discretizations of a gaussian AR(1), but Rouwenhorst takes over above rho=0.99 - which is why
% discretizeAR1_FarmerToda prints a comment recommending Rouwenhorst at rho>=0.99.
evmod=output.moderate.err_var;
fbk=outputP2.FarmerToda.moderate.fallbackfrac;
solveok=(fbk<0.05);
fprintf('\nat rho=%g (moderate), FarmerToda beats Tauchen on the variance WHEREVER THE ENTROPY SOLVE \n',calib.moderate.rho)
fprintf('   SUCCEEDS (fallback under 5%% of rows, which is znum up to %i here) [T2], this should be one: %i \n',max(znums(solveok)),all(evmod(4,solveok)<=evmod(1,solveok)))
fprintf('   and at the largest znum, where %2.1f%% of rows fall back, Tauchen is the more accurate: \n',100*fbk(end))
fprintf('   FarmerToda %2.3e against Tauchen %2.3e \n',evmod(4,end),evmod(1,end))
fprintf('   (Farmer & Toda 2017 pg 678 claim their method beats Tauchen for almost all gaussian \n')
fprintf('    AR(1) discretizations, and it does - until the moment-matching problem becomes too \n')
fprintf('    hard to solve. The claim is about the method, and this is about its solver.) \n')

evper=output.persistent.err_var;
eaper=output.persistent.err_ac;
fprintf('at rho=%g (persistent), Rouwenhorst beats FarmerToda on the autocorrelation [T2], this should be one: %i \n',calib.persistent.rho,all(eaper(2,:)<=eaper(4,:)))
fprintf('   (Farmer & Toda 2017 pg 678: above rho=0.99 Rouwenhorst takes over, which is what the \n')
fprintf('    comment printed by discretizeAR1_FarmerToda at rho>=0.99 is telling the user) \n')

% Floden (2008): the Tauchen-Hussey sigma_z variant is the badly-behaved one at high persistence.
% Asserted in DiscP2_AR1_TauchenHussey, where the baseSigma sweep lives.

%% NOT asserted, and why
fprintf('\nNOT asserted: any ordering between Tauchen and TauchenHussey. Floden (2008) finds Tauchen \n')
fprintf('   "relatively robust" and the TauchenHussey variants better or worse depending on which \n')
fprintf('   baseSigma is used, so the comparison is against a single TauchenHussey column rather \n')
fprintf('   than against the method, and no ordering holds across all three baseSigma choices. \n')
fprintf('NOT asserted: convergence to zero for Tauchen. At fixed Tauchen_q the grid never reaches \n')
fprintf('   past +-q*sd(z), so its variance error converges to a truncation floor rather than to \n')
fprintf('   zero. That is measured in DiscP2_AR1_Tauchen, and P1 established the same thing for \n')
fprintf('   the iid case. \n')

%% The mew convention, seen from the comparison
% All four methods put E(z)=mew/(1-rho). On the drift calibration that is 1.0, ten times mew. If
% any method had centred its grid on mew instead, its mean error on this calibration would be
% about 0.9 and about zero on the other two.
fprintf('\nall four methods agree on E(z) at the drift calibration, where mew=%g but E(z)=%g: \n',calib.drift.mew,calib.drift.mew/(1-calib.drift.rho))
em=[outputP2.Tauchen.drift.err_mean; outputP2.Rouwenhorst.drift.err_mean; outputP2.TauchenHussey.drift.err_mean; outputP2.FarmerToda.drift.err_mean];
for m_c=1:4
    fprintf('   %-14s worst mean error over the sweep: %2.3e \n',nmlist{m_c},max(em(m_c,:)))
end
fprintf('no method is out by anything like mew*rho/(1-rho)=%g, which is what centring on mew would cost [T2], this should be one: %i \n',calib.drift.mew*calib.drift.rho/(1-calib.drift.rho),max(em(:))<0.01)

%% Figures
figure(figure_c)
for cal_c=1:3
    cname=calib.names{cal_c};
    ev=output.(cname).err_var; ea=output.(cname).err_ac;
    subplot(2,3,cal_c)
    plot(znums,max(ev(1,:),10^(-18)),'o-',znums,max(ev(2,:),10^(-18)),'s-',znums,max(ev(3,:),10^(-18)),'^-',znums,max(ev(4,:),10^(-18)),'v-')
    set(gca,'XScale','log'); set(gca,'YScale','log')
    title(['|variance error|: ',cname]); xlabel('znum')
    if cal_c==1
        legend(nmlist,'Location','best')
    end
    subplot(2,3,3+cal_c)
    plot(znums,max(ea(1,:),10^(-18)),'o-',znums,max(ea(2,:),10^(-18)),'s-',znums,max(ea(3,:),10^(-18)),'^-',znums,max(ea(4,:),10^(-18)),'v-')
    set(gca,'XScale','log'); set(gca,'YScale','log')
    title(['|autocorrelation error|: ',cname]); xlabel('znum')
end

end
