function output=DiscP7_zeromean(calib)
% P7: kfttoptions.setmixturemutoenforcezeromean
%
% This option has never been exercised anywhere - not by the bank, not by any other test. It takes a
% mu_i with one FEWER row than mixprobs_i and sigma_i, appends the missing row, and solves it so the
% mixture has mean zero at every age:
%     mu_i(end,j) = -sum(mu_i(1:end-1,j).*mixprobs_i(1:end-1,j)) / mixprobs_i(end,j)
% The command then warns if the resulting mean is not zero.
%
% WHY IT MATTERS BEYOND ITS OWN CORRECTNESS. Its existence is the best evidence available about a
% design question the bank has been carrying since P3: the three gaussian-mixture commands disagree
% about whether the grid centres on mew/(1-rho) or on (mew+E(e))/(1-rho), and they coincide only for
% a mean-zero mixture. An option whose whole purpose is to force mean zero suggests the intended
% reading is that mixtures should be mean zero - which would make the KFTT centring right and the
% Tauchen one merely more general. That is item 1 of DiscretizationMethods_todo.md, and this subcode
% is the evidence for it rather than the decision.

fprintf('\n========== P7: setmixturemutoenforcezeromean ========== \n')

output=struct();
J=calib.J; znum=15;
mew=calib.vary.mew; rho=calib.vary.rho;
p=calib.vary.mixprobs_i; sd=calib.vary.sigma_i;
mu_full=calib.vary.mu_i;
nmix=size(p,1);

% Give it all but the last row of mu_i and let it solve for the last.
mu_short=mu_full(1:nmix-1,:);
ko=struct(); ko.verbose=0; ko.setmixturemutoenforcezeromean=1;
[zgA,pzA,j1A,ooA]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu_short,sd,znum,J,ko);
zgA=gather(zgA); pzA=gather(pzA); j1A=gather(j1A);

% What the solved-for row must be, computed here rather than read from the command
mu_solved=-(sum(mu_short.*p(1:nmix-1,:),1))./p(nmix,:);
mu_recon=[mu_short;mu_solved];
emean_recon=sum(p.*mu_recon,1);
fprintf('the reconstructed mixture has mean zero at every age [T1], this should be zero: %2.8e \n',max(abs(emean_recon)))

% The command should now be equivalent to being handed the full mu_i directly. That is the real
% check: the option is a convenience, so it must produce exactly what the explicit call produces.
ko2=struct(); ko2.verbose=0;
[zgB,pzB,j1B]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu_recon,sd,znum,J,ko2);
zgB=gather(zgB); pzB=gather(pzB); j1B=gather(j1B);
fprintf('setting the option equals passing the solved mu_i explicitly: z_grid_J [T0], this should be zero: %2.8e \n',max(abs(zgA-zgB),[],'all'))
fprintf('setting the option equals passing the solved mu_i explicitly: pi_z_J [T0], this should be zero: %2.8e \n',max(abs(pzA-pzB),[],'all'))
fprintf('setting the option equals passing the solved mu_i explicitly: jequaloneDistz [T0], this should be zero: %2.8e \n',max(abs(j1A-j1B)))

% And the process it produces really is mean zero in the innovation, so the two grid-centring
% conventions coincide: with E(e)=0 the KFTT grid centre equals E(z_j).
mewz_nodrift=zeros(1,J);
mewz_nodrift(1)=mew(1);
for j_c=2:J
    mewz_nodrift(j_c)=mew(j_c)+rho(j_c)*mewz_nodrift(j_c-1);
end
gmid=(zgA(end,:)+zgA(1,:))/2;
fprintf('with the mixture forced to mean zero the grid centre equals E(z_j) [T1], this should be zero: %2.8e \n',max(abs(gmid-mewz_nodrift)))
fprintf('   (that is the point: the two centring conventions in this family agree exactly when the \n')
fprintf('    mixture is mean zero, and this option exists to make that so - see the header) \n')

% The solve divides by mixprobs_i(end,:), so a tiny last weight is where it should break. Reported
% rather than asserted: what SHOULD happen for a near-zero weight is a design question, and the
% honest thing is to record what does happen.
fprintf('\n--- what happens when the last mixture weight is small, since the solve divides by it --- \n')
for w=[0.25,0.05,0.01]
    pw=[1-w;w]*ones(1,J);
    muw=mu_full(1,:); % one row short, as the option requires
    kow=struct(); kow.verbose=0; kow.setmixturemutoenforcezeromean=1;
    [zgw,~,~,~]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,pw,muw,sd,znum,J,kow);
    zgw=gather(zgw);
    musolved=-(muw.*pw(1,:))./pw(2,:);
    fprintf('last weight %4.2f: the solved component mean reaches %+2.3f, and the age-1 grid spans %2.3f \n',w,max(abs(musolved)),zgw(end,1)-zgw(1,1))
end
fprintf('   (a small last weight forces a large component mean to hold the average at zero, which \n')
fprintf('    widens the grid; nothing is asserted here because what should happen is a design call) \n')

output.emean_recon=emean_recon;
output.done=1;

end
