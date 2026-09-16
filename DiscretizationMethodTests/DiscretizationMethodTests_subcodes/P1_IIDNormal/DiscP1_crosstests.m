function output=DiscP1_crosstests(calib,znums)
% P1: the identities that hold within the iid normal block
%
% Because this bank is organised by process rather than by command, a reduction like "an AR(1) at
% rho=0 is an iid normal" sits here, next to the iid methods it reduces to, rather than in a
% separate pile of cross-tests.
%
% No figures, per the convention in the other banks.

fprintf('\n========== P1: cross-tests ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;
enum=15;
Tauchen_q=3;

%% 1. discretizeIID_Tauchen vs discretizeIIDNormal_Tauchen
% These are currently exact copies, so this is T0. It is worth keeping after discretizeIID_Tauchen
% is generalized to other distributions: it then becomes the statement that the general command
% still reproduces the normal-only one on the normal case.
[eg_a,pe_a]=discretizeIID_Tauchen(mew,sigma,enum,Tauchen_q,struct());
[eg_b,pe_b]=discretizeIIDNormal_Tauchen(mew,sigma,enum,Tauchen_q,struct());
fprintf('IID_Tauchen vs IIDNormal_Tauchen: e_grid [T0], this should be zero: %2.8e \n',max(abs(eg_a-eg_b)))
fprintf('IID_Tauchen vs IIDNormal_Tauchen: pi_e [T0], this should be zero: %2.8e \n',max(abs(pe_a-pe_b)))

%% 2. discretizeIID_TanakaToda vs discretizeIIDNormal_TanakaToda
[eg_c,pe_c]=discretizeIID_TanakaToda(mew,sigma,enum,struct());
[eg_d,pe_d]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,struct());
fprintf('IID_TanakaToda vs IIDNormal_TanakaToda: e_grid [T0], this should be zero: %2.8e \n',max(abs(eg_c-eg_d)))
fprintf('IID_TanakaToda vs IIDNormal_TanakaToda: pi_e [T0], this should be zero: %2.8e \n',max(abs(pe_c-pe_d)))

%% 3. discretizeAR1_Tauchen at rho=0 vs discretizeIIDNormal_Tauchen
% Same grid (at rho=0 the unconditional std dev is sigma, and both use mew +- Tauchen_q*sigma),
% and every row of the AR(1) transition matrix must equal the iid probability vector.
[zg,pz]=discretizeAR1_Tauchen(mew,0,sigma,enum,Tauchen_q,struct());
fprintf('AR1_Tauchen(rho=0) vs IIDNormal_Tauchen: grid [T1], this should be zero: %2.8e \n',max(abs(zg-eg_b)))
fprintf('AR1_Tauchen(rho=0) vs IIDNormal_Tauchen: every row of pi_z equals pi_e [T1], this should be zero: %2.8e \n',max(abs(pz-repmat(pe_b',enum,1)),[],'all'))

%% 4. discretizeAR1_FarmerToda at rho=0 vs discretizeIIDNormal_TanakaToda
% Tanaka-Toda is the iid analogue of Farmer-Toda, so this is the same method reached two ways.
% The defaults already line up at rho=0 (both pick gauss-hermite, both use nMoments=2 and the same
% nSigmas rule), but they are set explicitly here rather than relied on: a default that changes
% later should not silently turn this identity into a comparison of two different methods.
ftopts=struct(); ftopts.method='even'; ftopts.nMoments=2; ftopts.nSigmas=3;
ttopts=struct(); ttopts.method='even'; ttopts.nMoments=2; ttopts.nSigmas=3;
[zg2,pz2]=discretizeAR1_FarmerToda(mew,0,sigma,enum,ftopts);
[eg2,pe2]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,ttopts);
fprintf('AR1_FarmerToda(rho=0) vs IIDNormal_TanakaToda: grid [T1], this should be zero: %2.8e \n',max(abs(zg2-eg2)))
fprintf('AR1_FarmerToda(rho=0) vs IIDNormal_TanakaToda: every row of pi_z equals pi_e [T1], this should be zero: %2.8e \n',max(abs(pz2-repmat(pe2',enum,1)),[],'all'))

%% 5. MVNormal_ProbabilitiesOnGrid at M=1 vs discretizeIIDNormal_Tauchen
% Both allocate the cdf mass between adjacent midpoints, with the two outer bins running to
% +-Inf, so on the same grid they are two implementations of one formula and must agree.
% This is also the regression test for the M=1 branch of MVNormal_ProbabilitiesOnGrid, which used
% to pass the VARIANCE to normcdf where normcdf wants the standard deviation.
mvnopts=struct(); mvnopts.verbose=0;
P=MVNormal_ProbabilitiesOnGrid(eg_b,mew,sigma^2,enum,mvnopts);
fprintf('MVNormal_ProbabilitiesOnGrid(M=1) vs IIDNormal_Tauchen: pi_e [T1], this should be zero: %2.8e \n',max(abs(P(:)-pe_b(:))))
fprintf('MVNormal_ProbabilitiesOnGrid(M=1): sums to one [T1], this should be zero: %2.8e \n',abs(sum(P(:))-1))

%% 6. The reductions must hold across the whole sweep, not just at one enum
worst3=0; worst5=0;
for c_c=1:length(znums)
    en=znums(c_c);
    [egx,pex]=discretizeIIDNormal_Tauchen(mew,sigma,en,Tauchen_q,struct());
    [~,pzx]=discretizeAR1_Tauchen(mew,0,sigma,en,Tauchen_q,struct());
    worst3=max(worst3,max(abs(pzx-repmat(pex',en,1)),[],'all'));
    Px=MVNormal_ProbabilitiesOnGrid(egx,mew,sigma^2,en,mvnopts);
    worst5=max(worst5,max(abs(Px(:)-pex(:))));
end
fprintf('AR1_Tauchen(rho=0) identity, worst over the whole enum sweep [T1], this should be zero: %2.8e \n',worst3)
fprintf('MVNormal(M=1) identity, worst over the whole enum sweep [T1], this should be zero: %2.8e \n',worst5)

output.worst_AR1Tauchen_rho0=worst3;
output.worst_MVNormal_M1=worst5;

end
