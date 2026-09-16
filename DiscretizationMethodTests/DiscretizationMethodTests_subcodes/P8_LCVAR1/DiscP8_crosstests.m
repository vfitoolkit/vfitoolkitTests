function DiscP8_crosstests(calib,znum)
% P8: cross-tests for discretizeLifeCycleVAR1_Tauchen
%
% WITH ONE COMMAND IN THE BLOCK, THIS IS WHERE THE REAL ORACLES LIVE. Three independent
% implementations stand in for special cases the life-cycle VAR command must reproduce:
%   discretizeVAR1_Tauchen                        - the stationary VAR it extends
%   discretizeLifeCycleAR1_FellaGallipoliPanTauchen - the life-cycle AR(1) it generalises
%   discretizeLifeCycleAR1_FellaGallipoliPanTauchen, twice and kron'd - the independent-variables case
% None of these is an approximation of the life-cycle VAR command; each is separately written code
% that has to agree with it exactly where the special case bites.

fprintf('\n========== P8: cross-tests ========== \n')
J=calib.J; M=calib.M;

%% 1. The frozen-life-cycle identity, stationary-initial-condition form
% Constant parameters plus the STATIONARY initial condition means the age profile is flat, so the
% life-cycle command must equal the stationary one at EVERY age, not just in the limit. This is the
% form that surfaced B18. Both the grid and the transition matrix are asserted, because a command
% could get the grid right and assemble the transitions wrongly.
Rho=calib.frozen.Rho; SigmaSq=calib.frozen.SigmaSq;
tf=struct(); tf.verbose=0;
tf.initialj1mewz=calib.frozen.statmewz;
tf.initialj1SigmaSqz=calib.frozen.statSigmaSqz;
[zgF,pzF]=discretizeLifeCycleVAR1_Tauchen(calib.frozen.Mew_J,calib.frozen.Rho_J,calib.frozen.SigmaSq_J,znum,J,tf);
zgF=gather(zgF); pzF=gather(pzF);
% The stationary command takes Tauchen_q as a REQUIRED argument - it has no default at all - so it
% has to be handed the width the life-cycle command chose for itself, or the two build different
% grids for reasons that have nothing to do with the identity.
nS=min(sqrt(znum-1),3);
ts=struct(); ts.verbose=0;
[zgS,pzS]=discretizeVAR1_Tauchen(zeros(M,1),Rho,SigmaSq,znum,nS*ones(M,1),ts);
zgS=gather(zgS); pzS=gather(pzS);
gapg=0; gapp=0;
for j_c=1:J
    gapg=max(gapg,max(abs(zgF(:,j_c)-zgS)));
end
for j_c=1:J-1
    gapp=max(gapp,max(abs(pzF(:,:,j_c)-pzS),[],'all'));
end
fprintf('\n--- 1. frozen parameters + stationary initial condition vs discretizeVAR1_Tauchen --- \n')
fprintf('the grid is the same at EVERY age [T1], this should be zero: %2.8e \n',gapg)
fprintf('the transition matrix is the same at every age [T1], this should be zero: %2.8e \n',gapp)
fprintf('   (the life-cycle command is handed no width, so it takes its own default min(sqrt(znum-1),3); \n')
fprintf('   the stationary command has no default and is handed that same number explicitly.) \n')

%% 2. The frozen-life-cycle identity, transient form
% From the default z0=0 the age profile is a genuine transient climbing toward the stationary
% covariance, so only the LAST age approaches the stationary answer, and only asymptotically. The
% rate is known - the covariance deviation decays like max(|eig(Rho)|)^(2j) - so the check is against
% a bar computed from that rate rather than against zero. A large J is used to make the gap small.
Jbig=120;
Mew_big=zeros(M,Jbig); Rho_big=repmat(Rho,[1,1,Jbig]); Sig_big=repmat(SigmaSq,[1,1,Jbig]);
tt=struct(); tt.verbose=0;
[zgT,pzT]=discretizeLifeCycleVAR1_Tauchen(Mew_big,Rho_big,Sig_big,znum,Jbig,tt);
zgT=gather(zgT); pzT=gather(pzT);
lam=max(abs(eig(Rho)));
SigmaSq1=SigmaSq;
dev0=max(abs(SigmaSq1-calib.frozen.statSigmaSqz),[],'all');
predgap=dev0*lam^(2*(Jbig-1));
fprintf('\n--- 2. frozen parameters from z0=0, large J, vs discretizeVAR1_Tauchen at the last age --- \n')
fprintf('J=%i, and the covariance deviation decays like |eig(Rho)|^(2j) with |eig(Rho)|=%2.4f, so the \n',Jbig,lam)
fprintf('predicted covariance gap at the last age is %2.3e; the grid gap should be of the same order. \n',predgap)
lastgapg=max(abs(zgT(:,Jbig)-zgS));
fprintf('the grid at the last age matches the stationary grid [T2], this should be below %g: %2.3e \n',max(1e-6,100*predgap),lastgapg)
fprintf('the transition at the last age matches [T2], this should be below %g: %2.3e \n',max(1e-6,100*predgap),max(abs(pzT(:,:,Jbig-1)-pzS),[],'all'))
fprintf('and the gap SHRINKS with age: last age vs age 2 [T2], this should be one: %i \n',lastgapg<max(abs(zgT(:,2)-zgS)))

%% 3. M=1 against the life-cycle AR(1) command
% At one variable the VAR is an AR(1), and discretizeLifeCycleAR1_FellaGallipoliPanTauchen is a
% separately written extension of the SAME stationary method to age-dependent parameters. The two
% share a default width, min(sqrt(znum-1),3), which is the one place in the family that cap is 3
% rather than 4 - so this identity also pins that the two agree about the cap.
rho1=calib.diag.rho(1); sig1=calib.diag.sigma(1);
Mew1=zeros(1,J); Rho1=zeros(1,1,J); Sig1=zeros(1,1,J);
Rho1(1,1,:)=rho1; Sig1(1,1,:)=sig1^2;
t1=struct(); t1.verbose=0;
[zg1,pz1]=discretizeLifeCycleVAR1_Tauchen(Mew1,Rho1,Sig1,znum,J,t1);
zg1=gather(zg1); pz1=gather(pz1);
f1=struct(); f1.verbose=0;
[zgA,pzA]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(zeros(1,J),rho1*ones(1,J),sig1*ones(1,J),znum,J,f1);
zgA=gather(zgA); pzA=gather(pzA);
fprintf('\n--- 3. M=1 vs discretizeLifeCycleAR1_FellaGallipoliPanTauchen --- \n')
fprintf('z_grid_J agrees at every age [T1], this should be zero: %2.8e \n',max(abs(zg1-zgA),[],'all'))
fprintf('pi_z_J agrees at every age [T1], this should be zero: %2.8e \n',max(abs(pz1-pzA),[],'all'))
fprintf('both return pi_z_J with third dimension J-1 [T0], this should be one: %i \n',(size(pz1,3)==J-1)&&(size(pzA,3)==J-1))

%% 4. Diagonal Rho and SigmaSq: the joint chain is the kron of two univariate chains
% THIS IS THE STATE-ORDERING TEST, and it is not a formality. P5 found that the two stationary VAR
% commands order their joint index OPPOSITELY - discretizeVAR1_Tauchen matches CreateGridvals, with
% variable 1 varying fastest, while discretizeVAR1_FarmerToda uses allcomb2, with the last variable
% fastest. A life-cycle command that picked the wrong one would still produce a valid stochastic
% matrix with correct marginals, and only a test that reconstructs the joint from its parts catches
% it. With independent variables the joint transition must be a kron of the two univariate ones, and
% which kron ordering holds is exactly the convention question.
rd=calib.diag.rho; sd=calib.diag.sigma;
td=struct(); td.verbose=0;
[zgD,pzD]=discretizeLifeCycleVAR1_Tauchen(calib.diag.Mew_J,calib.diag.Rho_J,calib.diag.SigmaSq_J,znum,J,td);
zgD=gather(zgD); pzD=gather(pzD);
fa=struct(); fa.verbose=0;
[zga,pza]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(zeros(1,J),rd(1)*ones(1,J),sd(1)*ones(1,J),znum,J,fa);
[zgb,pzb]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(zeros(1,J),rd(2)*ones(1,J),sd(2)*ones(1,J),znum,J,fa);
zga=gather(zga); pza=gather(pza); zgb=gather(zgb); pzb=gather(pzb);
fprintf('\n--- 4. diagonal Rho and SigmaSq: the joint chain against a kron of two univariate chains --- \n')
fprintf('the stacked grid is the two univariate grids stacked [T1], this should be zero: %2.8e \n',max(abs(zgD-[zga;zgb]),[],'all'))
g12=0; g21=0;
for j_c=1:J-1
    g12=max(g12,max(abs(pzD(:,:,j_c)-kron(pzb(:,:,j_c),pza(:,:,j_c))),[],'all')); % variable 1 fastest
    g21=max(g21,max(abs(pzD(:,:,j_c)-kron(pza(:,:,j_c),pzb(:,:,j_c))),[],'all')); % variable 2 fastest
end
fprintf('joint equals kron(P2,P1), i.e. VARIABLE 1 VARIES FASTEST [T1], this should be zero: %2.8e \n',g12)
fprintf('joint equals kron(P1,P2), i.e. variable 2 varies fastest, the other convention: %2.8e \n',g21)
fprintf('exactly one of those two orderings holds [T0], this should be one: %i \n',(g12<1e-10)~=(g21<1e-10))
fprintf('   The first is CreateGridvals'' ordering, which is what the rest of the toolkit uses to unpack \n')
fprintf('   a joint z grid, and it is the one discretizeVAR1_Tauchen uses. If the second holds instead, \n')
fprintf('   this command disagrees with its own stationary counterpart about the joint index. \n')

%% 5. rho=0 at every age makes every row of the transition identical
% With no persistence the conditional distribution does not depend on the current state, so every row
% of pi_z_J(:,:,j) is the same. It is a weak check but it is independent of everything above, and it
% would catch a transition matrix built from the wrong conditional mean.
Rho0=zeros(M,M,J);
t0=struct(); t0.verbose=0;
[~,pz0]=discretizeLifeCycleVAR1_Tauchen(calib.diag.Mew_J,Rho0,calib.diag.SigmaSq_J,znum,J,t0);
pz0=gather(pz0);
rowgap=0;
for j_c=1:J-1
    Pj=pz0(:,:,j_c);
    rowgap=max(rowgap,max(abs(Pj-Pj(1,:)),[],'all'));
end
fprintf('\n--- 5. rho=0 at every age --- \n')
fprintf('every row of every transition slice is identical [T1], this should be zero: %2.8e \n',rowgap)

end
