function output=DiscP5_compare(calib,znums,outputP5,figure_c)
% P5: method comparison, VAR(1) with normal innovations
%
% The comparison that matters here is the CROSS-covariance, not the variances. Both VAR methods and
% the misspecified comparator can get the marginal variances roughly right; only a genuine VAR
% method can get the covariance between the variables, and that is the whole reason for the block.

fprintf('\n========== P5: method comparison, VAR(1) with normal innovations ========== \n')

output=struct();
S=calib.full.SigmaSqz;
fprintf('truth at M=2: variances [%2.6f %2.6f], cross-covariance %2.6f, autocorrelations [%2.4f %2.4f] \n',S(1,1),S(2,2),S(1,2),calib.full.autocorr(1),calib.full.autocorr(2))

nz=length(znums);
eT=outputP5.VAR1_Tauchen.err_cross(:,2);      % M=2 column
eF=outputP5.VAR1_FarmerToda.err_cross(:,2);
vT=outputP5.VAR1_Tauchen.err_var(:,2);
vF=outputP5.VAR1_FarmerToda.err_var(:,2);
eC=outputP5.comparator.err_cross;             % 3-by-nz

for c_c=1:nz
    if isnan(eT(c_c)) && isnan(eF(c_c))
        continue
    end
    fprintf('znum=%3i: cross-covariance error - VAR1_Tauchen %2.4e, VAR1_FarmerToda %2.4e, best independent-AR(1) %2.4e \n',znums(c_c),eT(c_c),eF(c_c),min(eC(:,c_c)))
end
fprintf('\n')
for c_c=1:nz
    if isnan(vT(c_c)) && isnan(vF(c_c))
        continue
    end
    fprintf('znum=%3i: variance error          - VAR1_Tauchen %2.4e, VAR1_FarmerToda %2.4e \n',znums(c_c),vT(c_c),vF(c_c))
end

% Both VAR methods must beat the misspecified comparator on the cross-covariance at every grid size
% that ran. This is not a claim about which VAR method is better - it is the claim that modelling
% the covariance at all beats not modelling it, which is the only thing the block asserts about the
% comparator.
ok=~isnan(eT) & ~isnan(eF);
bestC=min(eC,[],1)';
fprintf('\nVAR1_Tauchen beats the best independent-AR(1) on the cross-covariance at every znum [T2], this should be one: %i \n',all(eT(ok)<bestC(ok)))
fprintf('VAR1_FarmerToda beats the best independent-AR(1) on the cross-covariance at every znum [T2], this should be one: %i \n',all(eF(ok)<bestC(ok)))

% NOT asserted: an ordering between the two VAR methods. They differ in more than accuracy - one
% evaluates the multivariate normal cdf on rectangles, the other solves a maximum entropy problem
% in a transformed basis - and they do not even use the same state ordering (see DiscP5_crosstests).
last=find(ok,1,'last');
if ~isempty(last)
    fprintf('\nNOT asserted: an ordering between the two VAR methods. At the largest grid that ran \n')
    fprintf('   (znum=%i): Tauchen %2.4e, FarmerToda %2.4e on the cross-covariance. \n',znums(last),eT(last),eF(last))
end

output.eT=eT; output.eF=eF; output.bestC=bestC;

%% Figure
figure(figure_c)
semilogy(znums,max(eT,10^(-18)),'o-',znums,max(eF,10^(-18)),'s-',znums,max(bestC,10^(-18)),'d--')
xlabel('znum'); ylabel('cross-covariance error'); title('VAR(1): the moment only a VAR method can reach')
legend({'VAR1\_Tauchen','VAR1\_FarmerToda','best independent AR(1)'},'Location','best')

end
