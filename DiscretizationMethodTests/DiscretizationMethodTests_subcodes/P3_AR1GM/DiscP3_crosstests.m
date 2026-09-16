function output=DiscP3_crosstests(calib,znums)
% P3: the identities that hold within the gaussian-mixture block
%
% No figures, per the convention in the other banks.

fprintf('\n========== P3: cross-tests ========== \n')

output=struct();
mew=calib.mew; rho=calib.rho;
p=calib.mixprobs; mu=calib.mu; sd=calib.sigma;
znum=15;
optsG=struct(); optsG.method='even'; optsG.nMoments=2; optsG.nSigmas=3; optsG.verbose=0;

%% 1. Two identical components vs one component
% A mixture of two identical normals with weights summing to one IS that normal. This is the
% cheapest possible check on the mixture machinery: if the weighting or the moment accumulation
% were wrong, this would not reduce.
m0=0.05; s0=0.2;
[zg1,pz1,~]=discretizeAR1wGM_FarmerToda(mew,rho,1,m0,s0,znum,optsG);
[zg2,pz2,~]=discretizeAR1wGM_FarmerToda(mew,rho,[0.5;0.5],[m0;m0],[s0;s0],znum,optsG);
fprintf('two identical components vs one: z_grid [T1], this should be zero: %2.8e \n',max(abs(zg1-zg2)))
fprintf('two identical components vs one: pi_z, this should be below %g: %2.8e \n',calib.entropytol,max(abs(pz1-pz2),[],'all'))

%% 2. Unequal weights on identical components
% Same idea, but with weights that are not 0.5/0.5, so a bug that happened to cancel under equal
% weighting has nowhere to hide.
[zg3,pz3,~]=discretizeAR1wGM_FarmerToda(mew,rho,[0.7;0.3],[m0;m0],[s0;s0],znum,optsG);
fprintf('unequal weights on identical components vs one: z_grid [T1], this should be zero: %2.8e \n',max(abs(zg1-zg3)))
fprintf('unequal weights on identical components vs one: pi_z, this should be below %g: %2.8e \n',calib.entropytol,max(abs(pz1-pz3),[],'all'))

%% 3. A one-component mixture vs discretizeAR1_FarmerToda
% The reduction across commands. P2 asserts this too, at three calibrations, as the pin on the mew
% convention; it is repeated here at the mixture block's own rho so that a failure is attributable
% to this block rather than to P2's calibrations.
optsF=struct(); optsF.method='even'; optsF.nMoments=2; optsF.nSigmas=3; optsF.verbose=0;
[zg4,pz4,~]=discretizeAR1_FarmerToda(mew,rho,s0,znum,optsF);
[zg5,pz5,~]=discretizeAR1wGM_FarmerToda(mew,rho,1,0,s0,znum,optsG);
fprintf('nmix=1 (mean-zero) vs AR1_FarmerToda: z_grid [T1], this should be zero: %2.8e \n',max(abs(zg4-zg5)))
fprintf('nmix=1 (mean-zero) vs AR1_FarmerToda: pi_z, this should be below %g: %2.8e \n',calib.entropytol,max(abs(pz4-pz5),[],'all'))

%% 4. Component ordering must not matter
% A mixture is a set of weighted components, so permuting them is the same distribution. If any
% part of the implementation depended on the order they arrive in, this would catch it.
[zg6,pz6,~]=discretizeAR1wGM_FarmerToda(mew,rho,p,mu,sd,znum,optsG);
[zg7,pz7,~]=discretizeAR1wGM_FarmerToda(mew,rho,flipud(p),flipud(mu),flipud(sd),znum,optsG);
fprintf('permuting the mixture components: z_grid [T1], this should be zero: %2.8e \n',max(abs(zg6-zg7)))
fprintf('permuting the mixture components: pi_z, this should be below %g: %2.8e \n',calib.entropytol,max(abs(pz6-pz7),[],'all'))

%% 5. A zero-weight component
% Mathematically a component with probability zero cannot change the distribution, so adding one
% with a deliberately absurd mu and sigma should be a no-op. THE TWO COMMANDS DIFFER ON WHETHER
% THEY ACCEPT IT, and that difference is worth recording rather than working around:
%
%   discretizeAR1wGM_FarmerToda ERRORS. It builds a gmdistribution object to evaluate the mixture
%   density, and MATLAB's gmdistribution refuses a zero mixing proportion ("The mixing proportions
%   must be positive"). So the limitation comes from the Statistics Toolbox, not from the
%   discretization code.
%
%   discretizeAR1wGM_Tauchen ACCEPTS it. It sums component cdfs weighted by mixprobs directly, so
%   a zero weight simply contributes nothing.
%
% A user building a mixture programmatically - with a variable number of active components, say -
% will hit this, so it is checked both ways rather than assumed.
try
    discretizeAR1wGM_FarmerToda(mew,rho,[p;0],[mu;5],[sd;3],znum,optsG);
    fprintf('wGM_FarmerToda with a zero-weight component: ran without error (gmdistribution now allows it) \n')
catch
    fprintf('wGM_FarmerToda with a zero-weight component: errors, via gmdistribution, this should be one: 1 \n')
end
[zgT0,pzT0]=discretizeAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,3,struct());
[zgT1,pzT1]=discretizeAR1wGM_Tauchen(mew,rho,[p;0],[mu;5],[sd;3],znum,3,struct());
fprintf('wGM_Tauchen with a zero-weight component: z_grid [T1], this should be zero: %2.8e \n',max(abs(zgT0-zgT1)))
fprintf('wGM_Tauchen with a zero-weight component: pi_z [T1], this should be zero: %2.8e \n',max(abs(pzT0-pzT1),[],'all'))

%% 6. The two mixture commands against each other
% Both discretize the same process, by different methods, so they do NOT have to agree - what they
% must agree on is the moments each targets. Reported, with the Tauchen one's accuracy alongside.
[zgT,pzT]=discretizeAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,3,struct());
[~,vT,~,dT]=MarkovChainMoments(zgT,pzT);
mT=sum(dT.*zgT);
[~,vF,~,dF]=MarkovChainMoments(zg6,pz6);
mF=sum(dF.*zg6);
skT=sum(dT.*(zgT-mT).^3)/vT^1.5;
skF=sum(dF.*(zg6-mF).^3)/vF^1.5;
fprintf('wGM_Tauchen vs wGM_FarmerToda, both on the same process (reported, they need not agree): \n')
fprintf('   variance error %2.3e vs %2.3e; skewness error %2.3e vs %2.3e \n',abs(vT-calib.z.var),abs(vF-calib.z.var),abs(skT-calib.z.skew),abs(skF-calib.z.skew))

output.done=1;

end
