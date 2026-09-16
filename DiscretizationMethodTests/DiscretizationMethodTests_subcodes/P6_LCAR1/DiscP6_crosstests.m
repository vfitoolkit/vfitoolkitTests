function output=DiscP6_crosstests(calib,znums)
% P6: the identities that hold within the life-cycle AR(1) block
%
% The centrepiece is the FROZEN-LIFE-CYCLE IDENTITY. Give an age-dependent method constant
% parameters and it is discretizing a stationary process, so it must agree with the stationary
% method it extends. That gives P6 an oracle which is not an approximation of anything: the
% stationary commands are separately tested in P2, and the comparison is exact arithmetic against
% exact arithmetic rather than a discretization against a truth.
%
% Two forms, and they are different claims:
%   STATIONARY-INITIAL-CONDITION form. Constant parameters PLUS an initial condition equal to the
%      stationary distribution makes sigmaz(j) constant, so every age's grid is the stationary grid
%      and the identity holds at EVERY age. This is the strong form.
%   TRANSIENT form. Constant parameters from the default z(0)=0 makes sigmaz(j) rise towards the
%      stationary value, so the identity holds only at the LAST age, and only as J grows. This is
%      the form that surfaced B18 in the VAR(1) life-cycle command.

fprintf('\n========== P6: cross-tests ========== \n')

output=struct();
J=calib.J;
rho=calib.frozen.rho; sigma=calib.frozen.sigma; mew=calib.frozen.mew;
r=rho(1); sg=sigma(1); mw=mew(1);
statsigmaz=calib.frozen.statsigmaz;

%% 1. The frozen identity, stationary-initial-condition form
% Setting initialj0sigmaz to the stationary standard deviation freezes sigmaz(j) at that value for
% every j, so every age's grid coincides with the stationary grid and every transition matrix
% coincides with the stationary transition matrix.
fprintf('\n--- the frozen identity, stationary initial condition (holds at every age) --- \n')
fprintf('Every comparison below pins nSigmas and nMoments on BOTH sides. The life-cycle commands and \n')
fprintf('the stationary ones they extend have different defaults for both, so an unpinned comparison \n')
fprintf('would be measuring the defaults rather than the extension. \n')
for c_c=1:length(znums)
    znum=znums(c_c);

    % KFTT against discretizeAR1_FarmerToda. The two commands have DIFFERENT nMoments defaults -
    % KFTT uses 4, the stationary command uses 2 - so the identity is only stateable once both are
    % pinned to the same value, and it is checked at both settings. If they agree at one and not
    % the other, that is worth knowing rather than hiding behind a single choice.
    % nSigmas has to be pinned as well as nMoments, and for a less obvious reason. KFTT's default
    % is min(sqrt(2*(znum-1)),3) or min(sqrt(znum-1),3) depending on whether rho <= 1-2/(znum-1) -
    % and rho is a J-vector there, so the branch depends on the whole age path. The stationary
    % discretizeAR1_FarmerToda defaults to min(sqrt((znum-1)/2),2), a third formula. Left to their
    % defaults the two commands would build different grids and the identity would fail for a
    % reason that has nothing to do with the life-cycle extension being right or wrong.
    nSid=min(sqrt(znum-1),3);
    for nM=[2,4]
        ko=struct(); ko.initialj0sigmaz=statsigmaz; ko.nMoments=nM; ko.method='even'; ko.nSigmas=nSid;
        [zK,pK]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,ko);
        zK=gather(zK); pK=gather(pK);
        fo=struct(); fo.nMoments=nM; fo.method='even'; fo.nSigmas=nSid; fo.verbose=0;
        [zS,pS]=discretizeAR1_FarmerToda(mw,r,sg,znum,fo);
        zS=gather(zS); pS=gather(pS);
        dgrid=max(abs(zK-repmat(zS(:),1,J)),[],'all');
        dpi=max(abs(pK-repmat(pS,1,1,J-1)),[],'all');
        fprintf('znum=%3i, nMoments=%i: KFTT vs AR1_FarmerToda - z_grid at every age [T1] %2.8e, pi_z at every age [T1] %2.8e \n',znum,nM,dgrid,dpi)
        output.KFTT(c_c).nM(nM).dgrid=dgrid;
        output.KFTT(c_c).nM(nM).dpi=dpi;
    end

    % FellaGallipoliPan against discretizeAR1_Rouwenhorst. This one takes no mew, so the frozen
    % calibration is driftless. The stationary Rouwenhorst half-width is sqrt(znum-1)*sigmaz with
    % no cap, so nSigmas must be set to that here for the identity to be stateable at all above
    % znum=17 - which is itself the finding measured in DiscP6_LCAR1_FGP.
    go=struct(); go.initialj0sigmaz=statsigmaz; go.nSigmas=sqrt(znum-1);
    [zG,pG]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,go);
    zG=gather(zG); pG=gather(pG);
    [zR,pR]=discretizeAR1_Rouwenhorst(0,r,sg,znum,struct());
    zR=gather(zR); pR=gather(pR);
    dgrid=max(abs(zG-repmat(zR(:),1,J)),[],'all');
    dpi=max(abs(pG-repmat(pR,1,1,J-1)),[],'all');
    fprintf('znum=%3i: FGP vs AR1_Rouwenhorst - z_grid at every age [T1] %2.8e, pi_z at every age [T1] %2.8e \n',znum,dgrid,dpi)
    output.FGP(c_c).dgrid=dgrid; output.FGP(c_c).dpi=dpi;
    % and the same with the command's DEFAULT nSigmas, which is where the cap bites
    [zG2,pG2]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,struct('initialj0sigmaz',statsigmaz));
    d2=max(abs(gather(zG2)-repmat(zR(:),1,J)),[],'all');
    fprintf('znum=%3i: FGP vs AR1_Rouwenhorst at FGP''s DEFAULT nSigmas - z_grid %2.3e (zero only while sqrt(znum-1)<=4) \n',znum,d2)
    output.FGP(c_c).dgrid_default=d2;

    % FellaGallipoliPanTauchen against discretizeAR1_Tauchen. Tauchen_q must be set to the same
    % nSigmas the life-cycle command uses, which is min(sqrt(znum-1),3), not the bank's usual 3.
    nS=min(sqrt(znum-1),3);
    to=struct(); to.initialj0sigmaz=statsigmaz; to.nSigmas=nS;
    [zT,pT]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,to);
    zT=gather(zT); pT=gather(pT);
    [zA,pA]=discretizeAR1_Tauchen(mw,r,sg,znum,nS,struct());
    zA=gather(zA); pA=gather(pA);
    dgrid=max(abs(zT-repmat(zA(:),1,J)),[],'all');
    dpi=max(abs(pT-repmat(pA,1,1,J-1)),[],'all');
    fprintf('znum=%3i: FGP-Tauchen vs AR1_Tauchen - z_grid at every age [T1] %2.8e, pi_z at every age [T1] %2.8e \n',znum,dgrid,dpi)
    output.FGPT(c_c).dgrid=dgrid; output.FGPT(c_c).dpi=dpi;
end

%% 2. The frozen identity, transient form
% Default z(0)=0, constant parameters, large J. sigmaz(j) rises towards the stationary value, so
% the grids differ at early ages and converge at late ones. The identity is therefore a statement
% about the LAST age only, and it must get tighter as J grows - which is the check, since a single
% large J proves nothing about convergence.
fprintf('\n--- the frozen identity, transient form (holds at the last age, as J grows) --- \n')
znumT=15; nS=min(sqrt(znumT-1),3);
[zA,pA]=discretizeAR1_Tauchen(mw,r,sg,znumT,nS,struct());
zA=gather(zA); pA=gather(pA);
Jlist=[10,25,50,calib.Jbig];
dlast=zeros(1,length(Jlist));
for j_c=1:length(Jlist)
    Jt=Jlist(j_c);
    to=struct(); to.nSigmas=nS;
    [zT,pT]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew(1)*ones(1,Jt),r*ones(1,Jt),sg*ones(1,Jt),znumT,Jt,to);
    zT=gather(zT); pT=gather(pT);
    dlast(j_c)=max(abs(zT(:,Jt)-zA(:)));
    fprintf('J=%4i: the last-age grid against the stationary grid, difference %2.3e \n',Jt,dlast(j_c))
end
fprintf('the difference shrinks monotonically as J grows [T2], this should be one: %i \n',all(diff(dlast)<=10^(-14)))
fprintf('   (the early ages are SUPPOSED to differ; a process started from a point has not yet \n')
fprintf('    reached its stationary spread, and a method that matched at age 1 would be wrong) \n')
output.dlast=dlast; output.Jlist=Jlist;

%% 3. rho(j)=0 at every age: z is iid within each age, so every transition matrix has identical rows
fprintf('\n--- rho(j)=0 at every age --- \n')
znum0=11;
[z0,p0,j0]=discretizeLifeCycleAR1_KFTT(mew,zeros(1,J),sigma,znum0,J,struct());
z0=gather(z0); p0=gather(p0); j0=gather(j0);
maxrowdiff=0;
for j_c=1:J-1
    slice=p0(:,:,j_c);
    maxrowdiff=max(maxrowdiff,max(abs(slice-slice(1,:)),[],'all'));
end
fprintf('every transition matrix has identical rows [T1], this should be zero: %2.8e \n',maxrowdiff)
% and with no persistence the standard deviation at every age is just sigma(j)
[~,v0,~]=MarkovChainMoments_FHorz(z0,p0,j0);
fprintf('the standard deviation at every age is sigma(j) [T2], this should be small: %2.3e \n',max(abs(sqrt(v0(:)')-sigma)))

%% 4. The N_j-1 shape, and otheroutputs.sigma_z against the recursion
% All three commands emit pi_z_J with third dimension J-1, not J: pi_z_J(:,:,j) goes from age j to
% age j+1 and there is no age J+1. That was a deliberate change, and P6 is the first block to check
% it. Separately, each command reports otheroutputs.sigma_z, and comparing it against the recursion
% recomputed here is a real check - reading the truth off the command's own output would make every
% accuracy number in this block a tautology.
fprintf('\n--- the J-1 shape, and otheroutputs.sigma_z against the recursion --- \n')
mewv=calib.vary.mew; rhov=calib.vary.rho; sigmav=calib.vary.sigma;
sigmazT=calib.vary.sigmaz;
sigmazT_nodrift=zeros(1,J); sigmazT_nodrift(1)=sigmav(1);
for j_c=2:J
    sigmazT_nodrift(j_c)=sqrt(rhov(j_c)^2*sigmazT_nodrift(j_c-1)^2+sigmav(j_c)^2);
end
znumS=15;
[~,pKa,~,ooK]=discretizeLifeCycleAR1_KFTT(mewv,rhov,sigmav,znumS,J,struct());
[~,pGa,~,ooG]=discretizeLifeCycleAR1_FellaGallipoliPan(rhov,sigmav,znumS,J,struct());
[~,pTa,~,ooT]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mewv,rhov,sigmav,znumS,J,struct());
fprintf('KFTT: size(pi_z_J,3) is J-1 [T0], this should be zero: %i \n',size(pKa,3)~=J-1)
fprintf('FGP: size(pi_z_J,3) is J-1 [T0], this should be zero: %i \n',size(pGa,3)~=J-1)
fprintf('FGP-Tauchen: size(pi_z_J,3) is J-1 [T0], this should be zero: %i \n',size(pTa,3)~=J-1)
fprintf('KFTT: otheroutputs.sigma_z against the recursion [T1], this should be zero: %2.8e \n',max(abs(gather(ooK.sigma_z(:))'-sigmazT)))
fprintf('FGP: otheroutputs.sigma_z against the driftless recursion [T1], this should be zero: %2.8e \n',max(abs(gather(ooG.sigma_z(:))'-sigmazT_nodrift)))
fprintf('FGP-Tauchen: otheroutputs.sigma_z against the recursion [T1], this should be zero: %2.8e \n',max(abs(gather(ooT.sigma_z(:))'-sigmazT)))
fprintf('   (the sigma_z recursion does not involve mew at all, so KFTT and FGP-Tauchen must match \n')
fprintf('    the SAME profile that FGP does, drift or no drift) \n')
fprintf('KFTT and FGP report the same sigma_z profile [T1], this should be zero: %2.8e \n',max(abs(gather(ooK.sigma_z(:))-gather(ooG.sigma_z(:)))))
% KFTT also reports the MEAN profile, which the other two do not. It has its own recursion,
% mewz(j)=mew(j)+rho(j)*mewz(j-1), and unlike sigma_z it does depend on mew - so checking both
% separates a wrong drift from a wrong variance, which a single combined check could not.
mewzT=calib.vary.mewz;
fprintf('KFTT: otheroutputs.mew_z against the mean recursion [T1], this should be zero: %2.8e \n',max(abs(gather(ooK.mew_z(:))'-mewzT)))

%% 5. Error paths: an unsupported method name (B8), and the GMQ forwarding bug
% B8 was never a correctness problem. The switch on method has no `otherwise`, but the variables it
% assigns are assigned nowhere else, so an unsupported name always errored - just with
% "Unrecognized function or variable 'X1'", naming an internal variable in a file the caller did not
% write, with no mention of method. Fixed by validating the name up front, in the style
% discretizeVAR1_FarmerToda already used. The check here is that the message names the OPTION.
fprintf('\n--- error paths --- \n')
try
    discretizeLifeCycleAR1_KFTT(calib.vary.mew,calib.vary.rho,calib.vary.sigma,9,J,struct('method','not-a-real-method'));
    fprintf('an unsupported method name: ran without error [T0], this should not happen \n')
catch ME
    fprintf('an unsupported method name: errors, and the message names the option [T0], this should be one: %i \n',contains(ME.message,'kfttoptions.method'))
    fprintf('   the message is: %s \n',ME.message)
end

% THE GMQ FORWARDING BUG, found 2026-08-26 while scoping B8, and a genuine defect where B8 was not.
% discretizeLifeCycleAR1wGM_KFTT supports method='GMQ' (gaussian mixture quadrature) and supports
% initialj0sigmaz. Both are documented, both are legal, and together they used to crash: the command
% forwards its method to discretizeAR1_FarmerToda to build the period 0 distribution, and that
% command has no GMQ case, so the grid variable came back undefined. Two legal options on one
% command, failing inside a different command the caller never invoked.
%
% Fixed by translating GMQ to 'even' for that one forwarded call - period 0 is a plain normal, so
% there is no mixture for a mixture quadrature to be tuned to. This checks the combination runs and
% that the period 0 distribution it produces is sane.
fprintf('\n--- method=GMQ together with initialj0sigmaz, in the gaussian-mixture life-cycle command --- \n')
znumG=11;
% A mean-zero mixture, computed exactly rather than rounded: the KFTT commands centre their grid on
% mew/(1-rho), which is the unconditional mean only when the mixture has mean zero.
mixprobs=[0.7;0.3]; mu_i=[0.1;-0.7*0.1/0.3]; sigma_i=[0.15;0.35];
mixprobs_J=repmat(mixprobs,1,J); mu_J=repmat(mu_i,1,J); sigma_J=repmat(sigma_i,1,J);
ko=struct(); ko.method='GMQ'; ko.initialj0sigmaz=0.3;
[zG,pG,jG,ooG2]=discretizeLifeCycleAR1wGM_KFTT(calib.vary.mew,calib.vary.rho,mixprobs_J,mu_J,sigma_J,znumG,J,ko);
zG=gather(zG); pG=gather(pG); jG=gather(jG);
fprintf('it runs [T0], this should be zero: %i \n',any(size(zG)~=[znumG,J])+any(size(pG)~=[znumG,znumG,J-1]))
fprintf('rows of pi_z_J sum to one [T1], this should be zero: %2.8e \n',max(abs(sum(pG,2)-1),[],'all'))
fprintf('jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',abs(sum(jG)-1))
fprintf('otheroutputs.jequalzeroDistz sums to one [T1], this should be zero: %2.8e \n',abs(sum(gather(ooG2.jequalzeroDistz))-1))
% and the same call with method='even', which never had the bug: the period 0 distribution is built
% the same way in both, so these must now agree exactly.
ko2=ko; ko2.method='even';
[~,~,~,ooE]=discretizeLifeCycleAR1wGM_KFTT(calib.vary.mew,calib.vary.rho,mixprobs_J,mu_J,sigma_J,znumG,J,ko2);
fprintf('the period 0 distribution matches the method=''even'' call [T0], this should be zero: %2.8e \n',max(abs(gather(ooG2.jequalzeroDistz)-gather(ooE.jequalzeroDistz))))
fprintf('   (they must, because GMQ is translated to ''even'' for that one forwarded call - period 0 \n')
fprintf('    is a plain normal and there is no mixture for a mixture quadrature to be tuned to) \n')

output.done=1;

end
