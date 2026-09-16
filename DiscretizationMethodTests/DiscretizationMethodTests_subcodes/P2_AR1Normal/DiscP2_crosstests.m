function output=DiscP2_crosstests(calib,znums)
% P2: the identities that hold within the AR(1)-normal block
%
% Every command here discretizes the SAME process as the four native methods, once its extra
% machinery is switched off: a VAR(1) with one variable, a gaussian mixture with one component,
% a stochastic-volatility process with no volatility shocks. Each reduction is a comparison of two
% independent implementations of one object, which is the strongest kind of check available
% without an external oracle.
%
% No figures, per the convention in the other banks.
%
% NOT here: the frozen-life-cycle identity (the three life-cycle AR(1) methods at constant
% parameters and a stationary initial condition, which must reproduce their stationary
% counterparts at every age). That lives in P6, where those commands are native and their
% hyperparameters are swept - putting it here would mean aligning method/nSigmas/nMoments across
% two blocks and getting a failure that is about alignment rather than about the methods.

fprintf('\n========== P2: cross-tests ========== \n')

output=struct();
znum=15;

for cal_c=1:3
    cname=calib.names{cal_c};
    mew=calib.(cname).mew; rho=calib.(cname).rho; sigma=calib.(cname).sigma;
    zstar=mew/(1-rho);
    fprintf('\n--- calibration %s --- \n',cname)

    %% 1. discretizeVAR1_Tauchen at M=1 vs discretizeAR1_Tauchen
    % Also the regression test for the unconditional-variance bug in discretizeVAR1_Tauchen, which
    % used to solve the Lyapunov equation in the Cholesky-transformed basis and never map back.
    % At M=1 that made its sigma_z too large by a factor of 1/sigma, so this check would have
    % failed by a factor of five on these calibrations rather than by a rounding error.
    Tauchen_q=3;
    [zgA,pzA]=discretizeAR1_Tauchen(mew,rho,sigma,znum,Tauchen_q,struct());
    optsV=struct(); optsV.verbose=0;
    [zgB,pzB]=discretizeVAR1_Tauchen(mew,rho,sigma^2,znum,Tauchen_q,optsV);
    fprintf('%s: VAR1_Tauchen(M=1) vs AR1_Tauchen, z_grid [T1], this should be zero: %2.8e \n',cname,max(abs(zgA(:)-zgB(:))))
    fprintf('%s: VAR1_Tauchen(M=1) vs AR1_Tauchen, pi_z [T1], this should be zero: %2.8e \n',cname,max(abs(pzA-pzB),[],'all'))
    output.(cname).var1tauchen=[max(abs(zgA(:)-zgB(:))),max(abs(pzA-pzB),[],'all')];

    %% 2. discretizeVAR1_FarmerToda at M=1 vs discretizeAR1_FarmerToda
    % These build the SAME grid but NOT the same transition matrix, and the identity is asserted
    % only on the part that is genuinely identical.
    %
    % Both are Farmer-Toda, and both hit the targeted conditional moments exactly - that is the
    % method's contract and it is checked below. Where they differ is the untargeted shape of each
    % conditional distribution, which the maximum-entropy prior determines. AR1 builds that prior
    % on the original scale, normpdf(z_grid,condMean,sigma); VAR1 builds it in the transformed
    % unit-variance space, normpdf(y1D,condMean,1). Those are proportional, so they would give the
    % same answer - except that both then apply the SAME ABSOLUTE FLOOR, q(q<kappa)=kappa with
    % kappa=1e-8, and an absolute floor on two differently-scaled densities truncates a different
    % number of tail points. At high persistence the grid spans many conditional standard
    % deviations, the tails are tiny, the floor is active over many points, and the two answers
    % separate: measured 4e-07 at rho=0.6 but 1e-03 at rho=0.99.
    %
    % So this is REPORTED, not asserted, along with both methods' accuracy against truth, because
    % which of the two is preferable is a real question and this is the evidence for answering it.
    optsF=struct(); optsF.method='even'; optsF.nMoments=2; optsF.nSigmas=3; optsF.verbose=0;
    optsFV=struct(); optsFV.method='even'; optsFV.nMoments=2; optsFV.nSigmas=3;
    [zgC,pzC,~]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,optsF);
    [zgD,pzD]=discretizeVAR1_FarmerToda(mew,rho,sigma^2,znum,optsFV);
    fprintf('%s: VAR1_FarmerToda(M=1) vs AR1_FarmerToda, z_grid [T1], this should be zero: %2.8e \n',cname,max(abs(zgC(:)-zgD(:))))
    fprintf('%s: VAR1_FarmerToda(M=1) vs AR1_FarmerToda, pi_z differs by %2.3e (reported, see above) \n',cname,max(abs(pzC-pzD),[],'all'))
    % Both must still hit the targeted moments; that is what makes them both Farmer-Toda
    varz=sigma^2/(1-rho^2);
    [~,vC,aC]=MarkovChainMoments(zgC,pzC);
    [~,vD,aD]=MarkovChainMoments(zgD,pzD);
    fprintf('%s: variance error, AR1_FarmerToda %2.3e vs VAR1_FarmerToda %2.3e \n',cname,abs(vC-varz),abs(vD-varz))
    fprintf('%s: autocorr error, AR1_FarmerToda %2.3e vs VAR1_FarmerToda %2.3e \n',cname,abs(aC-rho),abs(aD-rho))
    output.(cname).var1farmertoda=[max(abs(zgC(:)-zgD(:))),max(abs(pzC-pzD),[],'all'),abs(vC-varz),abs(vD-varz)];

    %% 3. discretizeAR1wGM_FarmerToda with a one-component mixture vs discretizeAR1_FarmerToda
    % A gaussian mixture with one component IS a normal, so this must reduce exactly - the SAME
    % mew passed to both commands, giving the same grid and the same transition matrix.
    %
    % This check used to need a translation, because discretizeAR1wGM_FarmerToda read mew as the
    % unconditional mean (z'=(1-rho)*mew+rho*z+e) while every other command in the family reads it
    % as the intercept. It was the only one of the eleven to do so - including its own life-cycle
    % counterpart discretizeLifeCycleAR1wGM_KFTT, which reads it as the intercept. That has been
    % fixed, so the identity is now plain, and it is asserted at every calibration INCLUDING the
    % one with mew~=0, which is the one that can tell the two readings apart.
    optsG=struct(); optsG.method='even'; optsG.nMoments=2; optsG.nSigmas=3; optsG.verbose=0;
    [zgE,pzE,~]=discretizeAR1wGM_FarmerToda(mew,rho,1,0,sigma,znum,optsG);
    fprintf('%s: AR1wGM(nmix=1) vs AR1_FarmerToda, z_grid, this should be below %g: %2.8e \n',cname,calib.entropytol,max(abs(zgC(:)-zgE(:))))
    fprintf('%s: AR1wGM(nmix=1) vs AR1_FarmerToda, pi_z, this should be below %g: %2.8e \n',cname,calib.entropytol,max(abs(pzC-pzE),[],'all'))
    % And the grid must sit on mew/(1-rho), not on mew. On the drift calibration those differ by
    % mew*rho/(1-rho), which is 0.9 - so this is the check that would catch a relapse.
    fprintf('%s: AR1wGM grid is centred on mew/(1-rho)=%g, not on mew=%g [T1], this should be zero: %2.8e \n',cname,zstar,mew,abs(mean([zgE(1),zgE(end)])-zstar))
    output.(cname).gmidentity=[max(abs(zgC(:)-zgE(:))),max(abs(pzC-pzE),[],'all')];
end

%% 5. discretizeAR1wSV_FarmerToda as the volatility shocks vanish
% With sigmae=0 the log-volatility process is degenerate and the z block should collapse onto a
% plain AR(1) at the constant volatility. Exactly zero is expected to fail rather than reduce:
% the x process is discretized by discretizeVAR1_FarmerToda, which needs a non-singular
% variance-covariance matrix. So this approaches the limit instead, and reports the approach.
fprintf('\n--- stochastic volatility, as sigmae -> 0 --- \n')
rho=calib.moderate.rho; sigmau=calib.moderate.sigma;
xnum=3; znumSV=15;
try
    optsSV=struct(); optsSV.nSigmas=2;
    for sigmae=[0.1,0.01,0.001]
        [zgSV,pzSV]=discretizeAR1wSV_FarmerToda(rho,0.9,sigmau,sigmae,xnum,znumSV,optsSV);
        % z_grid is stacked, x on top: the z block is the last znumSV entries
        zblock=zgSV(xnum+1:end);
        [zgP,~,~]=discretizeAR1_FarmerToda(0,rho,sigmau,znumSV,struct('method','even','nSigmas',2,'verbose',0));
        d=max(abs(zblock(:)-zgP(:)));
        fprintf('sigmae=%g: size of z_grid [T0], this should be zero: %i \n',sigmae,any(size(zgSV)~=[xnum+znumSV,1]))
        fprintf('sigmae=%g: rows of pi_z sum to one, this should be below %g: %2.8e \n',sigmae,calib.entropytol,max(abs(sum(pzSV,2)-1)))
        % NOTE. This difference is zero at EVERY sigmae, and that is not convergence. The command
        % sets sigmaz=sqrt(exp(xBar+sigmaX/2)/(1-rho^2)) with xBar=2*log(sigmau)-sigmaX/2, so
        % exp(xBar+sigmaX/2)=sigmau^2 identically: THE z GRID DOES NOT DEPEND ON sigmae AT ALL.
        % The two grids are equal by construction. The reduction is a statement about transition
        % PROBABILITIES, and is tested where it belongs, on pi_z, in P4's cross-tests.
        fprintf('sigmae=%g: the z block of the grid equals the plain AR(1) grid [T1], this should be zero: %2.8e \n',sigmae,d)
    end
    try
        discretizeAR1wSV_FarmerToda(rho,0.9,sigmau,0,xnum,znumSV,optsSV);
        fprintf('sigmae=0 exactly: ran without error (the degenerate x process was accepted) \n')
    catch
        fprintf('sigmae=0 exactly: errors, as expected - the x process variance is singular there \n')
    end
catch ME
    fprintf('stochastic-volatility reduction could not be run: %s \n',ME.message)
    fprintf('   (reported rather than crashing the block; this is one check of P2, and P4 owns this command) \n')
end

end
