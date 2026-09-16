function output=DiscP1_compare(calib,znums,outputP1,figure_c)
% P1: compare the methods against each other and against analytic truth
%
% Consumes what the method subcodes already computed rather than recomputing it.
%
% Truth for an iid normal is closed form: mean=mew, variance=sigma^2.
%
% A bare threshold on an error only catches regressions - it says nothing about whether the method
% is accurate today. So the T2 numbers here are always paired with either an ORDERING that theory
% predicts, or a CONVERGENCE assertion in enum. Where theory predicts no ordering, this subcode
% says so and reports rather than asserting.

fprintf('\n========== P1: method comparison, iid normal ========== \n')
fprintf('truth: mean=%2.6f, variance=%2.6f \n',calib.mew,calib.sigma^2)

output=struct();

e_TA=outputP1.IIDNormal_Tauchen.err_var;
e_TT=outputP1.IIDNormal_TanakaToda.err_var;
e_AR=outputP1.AR1methods_rho0.err_var; % 4-by-nz: Tauchen, Rouwenhorst, FarmerToda, TauchenHussey

%% Orderings that theory predicts
% Tanaka-Toda TARGETS the mean and variance (maximum entropy moment matching), so it should match
% them to solver tolerance at every enum, where Tauchen only approximates them. This is the one
% ordering in this block that is a statement about the methods rather than about a threshold.
fprintf('TanakaToda beats Tauchen on the variance, at every enum [T2], this should be one: %i \n',all(e_TT<=e_TA))
fprintf('   worst enum for that comparison: %i \n',znums(find(e_TT>e_TA,1)))

% Rouwenhorst and Farmer-Toda also match exactly at rho=0 (Rouwenhorst by construction,
% Farmer-Toda by targeting), so they should beat Tauchen and Tauchen-Hussey too.
fprintf('Rouwenhorst beats AR1_Tauchen on the variance, at every znum [T2], this should be one: %i \n',all(e_AR(2,:)<=e_AR(1,:)))
fprintf('FarmerToda beats AR1_Tauchen on the variance, at every znum [T2], this should be one: %i \n',all(e_AR(3,:)<=e_AR(1,:)))

% NOT asserted: "Farmer-Toda beats Tauchen-Hussey". That is Farmer & Toda (2017)'s claim, and the
% reason discretizeAR1_TauchenHussey prints a comment recommending Farmer-Toda - but it is a claim
% about PERSISTENT processes, and this block is at rho=0. At rho=0 Tauchen-Hussey's Gauss-Hermite
% nodes integrate the mean and variance exactly, while Farmer-Toda reaches them only to its entropy
% solver's tolerance, so Tauchen-Hussey actually wins here. Asserting the ordering in this block
% would be asserting something neither paper claims. It belongs in P2, at rho>0.
fprintf('reported, not asserted: at rho=0 the variance error is %2.2e for TauchenHussey and %2.2e for FarmerToda (at the largest znum) \n',e_AR(4,end),e_AR(3,end))

%% Convergence, for the methods that approximate rather than match
fprintf('Tauchen variance error falls over the first half of the sweep [T2], this should be one: %i \n',e_TA(4)<e_TA(1))
% NOT asserted: monotonicity all the way to the largest enum. At a FIXED Tauchen_q the grid never
% reaches past mew +- q*sigma however many points are added, so the variance error converges to the
% truncation floor rather than to zero, and can tick back up once it gets there. Measured on this
% calibration: 1.7e-02, 3.9e-03, 9.6e-04, 1.4e-04, 3.4e-04, 4.2e-04 - falling to enum=31, then
% flat-to-rising. The genuine convergence statement needs Tauchen_q to grow with enum, and that is
% asserted in DiscP1_IIDNormal_Tauchen rather than here.

%% Where theory predicts nothing, say so
fprintf('NOT asserted: any ordering among the four TanakaToda grid methods (even, gauss-legendre, \n')
fprintf('   clenshaw-curtis, gauss-hermite). They all match the targeted moments exactly, so their \n')
fprintf('   differences show up only in moments beyond nMoments, and no ordering is predicted. \n')

output.err_var_Tauchen=e_TA;
output.err_var_TanakaToda=e_TT;
output.err_var_AR1methods=e_AR;
output.znums=znums;

%% Figure
% Accuracy on one axis, cost on the other: that is the actual method-selection question, so the
% two panels use the same colours and the same x-axis.
figure(figure_c)
subplot(1,2,1)
plot(znums,max(e_TA,10^(-18)),'o-',znums,max(e_TT,10^(-18)),'s-',znums,max(e_AR(2,:),10^(-18)),'^-',znums,max(e_AR(4,:),10^(-18)),'v-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('iid normal: |variance error|'); xlabel('enum')
legend('IIDNormal\_Tauchen','IIDNormal\_TanakaToda','AR1\_Rouwenhorst(rho=0)','AR1\_TauchenHussey(rho=0)','Location','best')
subplot(1,2,2)
plot(znums,outputP1.IIDNormal_Tauchen.runtime,'o-',znums,outputP1.IIDNormal_TanakaToda.runtime,'s-',znums,outputP1.AR1methods_rho0.runtime(2,:),'^-',znums,outputP1.AR1methods_rho0.runtime(4,:),'v-')
set(gca,'XScale','log'); set(gca,'YScale','log')
title('runtime (s)'); xlabel('enum')

end
