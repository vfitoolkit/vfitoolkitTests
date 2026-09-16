function output=DiscP7_crosstests(calib,znums)
% P7: the identities that hold within the life-cycle gaussian-mixture block
%
% Two oracles here, and they are different in kind. The REDUCTION - a one-component mixture is a
% normal, so each command must reproduce its own gaussian counterpart exactly - and the FROZEN
% identity, where constant parameters plus a stationary initial condition mean the age-dependent
% method must reproduce the stationary method it extends, at every age.
%
% Every comparison pins the options on both sides. The life-cycle commands and the stationary ones
% they extend have different defaults for method, nMoments and nSigmas, so an unpinned comparison
% would be measuring the defaults rather than the extension. P6 learned that the hard way.

fprintf('\n========== P7: cross-tests ========== \n')

output=struct();
J=calib.J;

%% 1. nmix=1: a one-component mixture is a normal
% discretizeLifeCycleAR1wGM_KFTT must reproduce discretizeLifeCycleAR1_KFTT, and the Tauchen one
% must reproduce discretizeLifeCycleAR1_FellaGallipoliPanTauchen - but only after the mixture's
% mean is folded into the intercept, because a one-component mixture with a non-zero mu is a normal
% with a shifted mean, and the gaussian commands take their innovation to be mean zero.
fprintf('\n--- nmix=1 reduces to the gaussian life-cycle commands --- \n')
mew=calib.vary.mew; rho=calib.vary.rho;
mu1=calib.vary.mu_i(1,:); sd1=calib.vary.sigma_i(1,:);
for c_c=1:length(znums)
    znum=znums(c_c);
    nS=min(sqrt(znum-1),4);
    ko=struct(); ko.verbose=0; ko.method='even'; ko.nMoments=4; ko.nSigmas=nS;
    [zM,pM]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,ones(1,J),mu1,sd1,znum,J,ko);
    go=struct(); go.verbose=0; go.method='even'; go.nMoments=4; go.nSigmas=nS;
    [zG,pG]=discretizeLifeCycleAR1_KFTT(mew+mu1,rho,sd1,znum,J,go);
    fprintf('znum=%3i: wGM_KFTT(nmix=1) vs LifeCycleAR1_KFTT - z_grid_J [T1] %2.8e, pi_z_J [T1] %2.8e \n',znum,max(abs(gather(zM)-gather(zG)),[],'all'),max(abs(gather(pM)-gather(pG)),[],'all'))
    to=struct(); to.verbose=0;
    [zT,pT]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,ones(1,J),mu1,sd1,znum,J,nS,to);
    fo=struct(); fo.verbose=0; fo.nSigmas=nS;
    [zF,pF]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew+mu1,rho,sd1,znum,J,fo);
    fprintf('znum=%3i: wGM_Tauchen(nmix=1) vs FGPTauchen - z_grid_J [T1] %2.8e, pi_z_J [T1] %2.8e \n',znum,max(abs(gather(zT)-gather(zF)),[],'all'),max(abs(gather(pT)-gather(pF)),[],'all'))
end

%% 2. The frozen identity: constant parameters plus a stationary initial condition
% Both stationary counterparts take the mixture directly, so no folding is needed here.
fprintf('\n--- the frozen identity, stationary initial condition (holds at every age) --- \n')
fmew=calib.frozen.mew; frho=calib.frozen.rho;
fp=calib.frozen.mixprobs_i; fmu=calib.frozen.mu_i; fsd=calib.frozen.sigma_i;
r=frho(1);
% the stationary standard deviation of z for the frozen process
fe_var=sum(fp(:,1).*(fmu(:,1).^2+fsd(:,1).^2))-sum(fp(:,1).*fmu(:,1))^2;
statsigmaz=sqrt(fe_var/(1-r^2));
for c_c=1:length(znums)
    znum=znums(c_c);
    nS=min(sqrt(znum-1),4);
    ko=struct(); ko.verbose=0; ko.method='even'; ko.nMoments=4; ko.nSigmas=nS;
    ko.initialj0sigmaz=statsigmaz;
    [zK,pK]=discretizeLifeCycleAR1wGM_KFTT(fmew,frho,fp,fmu,fsd,znum,J,ko);
    so=struct(); so.verbose=0; so.method='even'; so.nMoments=4; so.nSigmas=nS;
    [zS,pS]=discretizeAR1wGM_FarmerToda(fmew(1),r,fp(:,1),fmu(:,1),fsd(:,1),znum,so);
    zK=gather(zK); pK=gather(pK); zS=gather(zS); pS=gather(pS);
    fprintf('znum=%3i: wGM_KFTT vs AR1wGM_FarmerToda - z_grid at every age [T1] %2.8e, pi_z at every age [T1] %2.8e \n',znum,max(abs(zK-repmat(zS(:),1,J)),[],'all'),max(abs(pK-repmat(pS,1,1,J-1)),[],'all'))
    to=struct(); to.verbose=0; to.initialj0sigmaz=statsigmaz;
    [zT,pT]=discretizeLifeCycleAR1wGM_Tauchen(fmew,frho,fp,fmu,fsd,znum,J,nS,to);
    [zR,pR]=discretizeAR1wGM_Tauchen(fmew(1),r,fp(:,1),fmu(:,1),fsd(:,1),znum,nS,struct('verbose',0));
    zT=gather(zT); pT=gather(pT); zR=gather(zR); pR=gather(pR);
    fprintf('znum=%3i: wGM_Tauchen vs AR1wGM_Tauchen - z_grid at every age [T1] %2.8e, pi_z at every age [T1] %2.8e \n',znum,max(abs(zT-repmat(zR(:),1,J)),[],'all'),max(abs(pT-repmat(pR,1,1,J-1)),[],'all'))
end

%% 3. The transient form of the frozen identity
% Default z(0)=0, constant parameters, growing J. The last age approaches the stationary answer,
% and the check is that it gets tighter as J grows - a single large J proves nothing.
fprintf('\n--- the frozen identity, transient form (last age only, as J grows) --- \n')
znumT=15; nS=min(sqrt(znumT-1),4);
[zR,~]=discretizeAR1wGM_Tauchen(0,r,fp(:,1),fmu(:,1),fsd(:,1),znumT,nS,struct('verbose',0));
zR=gather(zR);
Jlist=[10,25,50,calib.Jbig];
dlast=zeros(1,length(Jlist));
for j_c=1:length(Jlist)
    Jt=Jlist(j_c);
    to=struct(); to.verbose=0;
    [zT,~]=discretizeLifeCycleAR1wGM_Tauchen(zeros(1,Jt),r*ones(1,Jt),repmat(fp(:,1),1,Jt),repmat(fmu(:,1),1,Jt),repmat(fsd(:,1),1,Jt),znumT,Jt,nS,to);
    zT=gather(zT);
    dlast(j_c)=max(abs(zT(:,Jt)-zR(:)));
    fprintf('J=%4i: the last-age grid against the stationary grid, difference %2.3e \n',Jt,dlast(j_c))
end
fprintf('the difference shrinks monotonically as J grows [T2], this should be one: %i \n',all(diff(dlast)<=10^(-14)))

%% 4. rho(j)=0 at every age
fprintf('\n--- rho(j)=0 at every age --- \n')
znum0=11;
p=calib.vary.mixprobs_i; mu=calib.vary.mu_i; sd=calib.vary.sigma_i;
for k_c=1:2
    if k_c==1
        [z0,pz0,j0]=discretizeLifeCycleAR1wGM_KFTT(mew,zeros(1,J),p,mu,sd,znum0,J,struct('verbose',0));
        nm='wGM_KFTT';
    else
        [z0,pz0,j0]=discretizeLifeCycleAR1wGM_Tauchen(mew,zeros(1,J),p,mu,sd,znum0,J,3,struct('verbose',0));
        nm='wGM_Tauchen';
    end
    z0=gather(z0); pz0=gather(pz0); j0=gather(j0);
    maxrowdiff=0;
    for j_c=1:J-1
        s=pz0(:,:,j_c);
        maxrowdiff=max(maxrowdiff,max(abs(s-s(1,:)),[],'all'));
    end
    fprintf('%-12s: every transition matrix has identical rows [T1], this should be zero: %2.8e \n',nm,maxrowdiff)
    [~,v0,~]=MarkovChainMoments_FHorz(z0,pz0,j0);
    fprintf('%-12s: the variance at every age is the innovation variance [T2], this should be small: %2.3e \n',nm,max(abs(v0(:)'-calib.vary.evar)))
end

%% 5. The J-1 shape, otheroutputs, and the two centring conventions
fprintf('\n--- the J-1 shape, otheroutputs, and the centring conventions --- \n')
znumS=15;
[~,pKa,~,ooK]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znumS,J,struct('verbose',0));
[zTa,pTa,~,ooT]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znumS,J,3,struct('verbose',0));
fprintf('wGM_KFTT: size(pi_z_J,3) is J-1 [T0], this should be zero: %i \n',size(pKa,3)~=J-1)
fprintf('wGM_Tauchen: size(pi_z_J,3) is J-1 [T0], this should be zero: %i \n',size(pTa,3)~=J-1)
% sigma_z does not involve the mixture MEAN, so the two commands must agree on it exactly even
% though they disagree about centring. That separates the two questions cleanly.
fprintf('the two commands agree on otheroutputs.sigma_z [T1], this should be zero: %2.8e \n',max(abs(gather(ooK.sigma_z(:))-gather(ooT.sigma_z(:)))))
fprintf('wGM_Tauchen: otheroutputs.sigma_z matches the cumulant recursion [T1], this should be zero: %2.8e \n',max(abs(gather(ooT.sigma_z(:))'-sqrt(calib.vary.varz))))
% And here is the centring difference itself, measured rather than assumed. The Tauchen command
% centres on E(z_j); the KFTT one on a recursion with no E(e) term. On a mean-zero mixture these
% coincide, which is why P3 could not see it; this calibration is not mean zero.
mewz_noE=zeros(1,J); mewz_noE(1)=mew(1);
for j_c=2:J
    mewz_noE(j_c)=mew(j_c)+rho(j_c)*mewz_noE(j_c-1);
end
gmidT=(zTa(end,:)+zTa(1,:))/2;
fprintf('wGM_Tauchen centres on E(z_j) [T1], this should be zero: %2.8e \n',max(abs(gather(gmidT)-calib.vary.mewz)))
fprintf('the two conventions differ by up to %2.4f on this calibration, against a sd(z) of about %2.4f \n',max(abs(calib.vary.mewz-mewz_noE)),max(sqrt(calib.vary.varz)))
fprintf('   (E(e)/(1-rho) is the offset, so persistence makes it large; on a mean-zero mixture it \n')
fprintf('    would be exactly zero and this block would not be able to see the difference at all) \n')
output.centringgap=max(abs(calib.vary.mewz-mewz_noE));

%% 6. Scalar versus age-vector width, and the GMQ forwarding regression (B27)
fprintf('\n--- scalar versus age-vector width (the B28 regression for these commands) --- \n')
o1=struct('verbose',0,'nSigmas',3); o2=struct('verbose',0,'nSigmas',3*ones(J,1));
[zA,pA,jA]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,15,J,o1);
[zB,pB,jB]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,15,J,o2);
fprintf('wGM_KFTT: z_grid_J [T0], this should be zero: %2.8e \n',max(abs(gather(zA)-gather(zB)),[],'all'))
fprintf('wGM_KFTT: pi_z_J [T0], this should be zero: %2.8e \n',max(abs(gather(pA)-gather(pB)),[],'all'))
fprintf('wGM_KFTT: jequaloneDistz [T0], this should be zero: %2.8e \n',max(abs(gather(jA)-gather(jB))))
[zC,pC,jC]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,15,J,3,struct('verbose',0));
[zD,pD,jD]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,15,J,3*ones(1,J),struct('verbose',0));
fprintf('wGM_Tauchen: z_grid_J [T0], this should be zero: %2.8e \n',max(abs(gather(zC)-gather(zD)),[],'all'))
fprintf('wGM_Tauchen: pi_z_J [T0], this should be zero: %2.8e \n',max(abs(gather(pC)-gather(pD)),[],'all'))
fprintf('wGM_Tauchen: jequaloneDistz [T0], this should be zero: %2.8e \n',max(abs(gather(jC)-gather(jD))))

% B27, moved here from P6 where it was sitting in the wrong block. method='GMQ' together with
% initialj0sigmaz used to crash inside discretizeAR1_FarmerToda, which has no GMQ case; the fix
% translates GMQ to 'even' for that one forwarded call, period 0 being a plain normal.
fprintf('\n--- method=GMQ together with initialj0sigmaz (the B27 regression) --- \n')
kg=struct(); kg.verbose=0; kg.method='GMQ'; kg.initialj0sigmaz=0.3;
[zG,pG,jG,ooG]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,11,J,kg);
zG=gather(zG); pG=gather(pG); jG=gather(jG);
fprintf('it runs [T0], this should be zero: %i \n',any(size(zG)~=[11,J])+any(size(pG)~=[11,11,J-1]))
fprintf('rows of pi_z_J sum to one [T1], this should be zero: %2.8e \n',max(abs(sum(pG,2)-1),[],'all'))
fprintf('jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',abs(sum(jG)-1))
ke=kg; ke.method='even';
[~,~,~,ooE]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,11,J,ke);
fprintf('the period 0 distribution matches the method=''even'' call [T0], this should be zero: %2.8e \n',max(abs(gather(ooG.jequalzeroDistz)-gather(ooE.jequalzeroDistz))))

%% 7. The computed grid width, against the tail-mass rule recomputed here
% Both mixture Tauchen commands now set their default width by solving for the w that leaves the
% same tail mass outside the grid as four standard deviations does for a normal, rather than capping
% at a constant 4. The rule is recomputed here from the calibration, independently of the command,
% and compared against the width the command actually used - which is readable off the grid, since
% a Tauchen grid spans exactly E(z_j) +- Tauchen_q*sd(z_j).
%
% Two things are checked and they are different claims. That a one-component mixture returns exactly
% 4, which is what makes the change backward compatible for a normal. And that a fat-tailed mixture
% returns the width the rule predicts.
fprintf('\n--- the computed default grid width, against the tail-mass rule --- \n')
epstail=2*(1-normcdf(4));
znumW=31;
% (i) a one-component mixture must give exactly 4
[zg1,~]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,ones(1,J),calib.vary.mu_i(1,:),calib.vary.sigma_i(1,:),znumW,J,[],struct('verbose',0));
zg1=gather(zg1);
sd1=zeros(1,J); sd1(1)=calib.vary.sigma_i(1,1);
for j_c=2:J
    sd1(j_c)=sqrt(rho(j_c)^2*sd1(j_c-1)^2+calib.vary.sigma_i(1,j_c)^2);
end
w1=(zg1(end,:)-zg1(1,:))/2./sd1;
fprintf('a one-component mixture gives width min(sqrt(znum-1),4)=%2.4f [T1], this should be zero: %2.8e \n',min(sqrt(znumW-1),4),max(abs(w1-min(sqrt(znumW-1),4))))
% (ii) the fat-tailed mixture must give what the rule predicts, recomputed here by bisection
[zgM,~]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znumW,J,[],struct('verbose',0));
zgM=gather(zgM);
wM=(zgM(end,:)-zgM(1,:))/2./sqrt(calib.vary.varz);
wpred=zeros(1,J);
for j_c=1:J
    se=sqrt(calib.vary.evar(j_c)); me=calib.vary.emean(j_c);
    wlo=0.5; whi=40;
    for b_c=1:200
        wmid=(wlo+whi)/2;
        tm=sum(p(:,j_c).*((1-normcdf((me+wmid*se-mu(:,j_c))./sd(:,j_c)))+normcdf((me-wmid*se-mu(:,j_c))./sd(:,j_c))));
        if tm>epstail
            wlo=wmid;
        else
            whi=wmid;
        end
    end
    wpred(j_c)=min(sqrt(znumW-1),(wlo+whi)/2);
end
fprintf('the fat-tailed mixture gives the width the rule predicts [T1], this should be zero: %2.8e \n',max(abs(wM-wpred)))
fprintf('   the rule asks for %2.2f to %2.2f standard deviations across the ages, against the 4 it \n',min(wpred),max(wpred))
fprintf('   would have used before; the innovation excess kurtosis runs %2.2f to %2.2f \n',min(calib.vary.exkurtz),max(calib.vary.exkurtz))

%% 8. Error paths
fprintf('\n--- error paths --- \n')
try
    discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,11,J,struct('verbose',0,'method','not-a-real-method'));
    fprintf('an unsupported method name: ran without error [T0], this should not happen \n')
catch ME
    fprintf('an unsupported method name: errors, and the message names the option [T0], this should be one: %i \n',contains(ME.message,'kfttoptions.method'))
end
try
    badp=p; badp(1,3)=badp(1,3)+0.2; % this column no longer sums to one
    discretizeLifeCycleAR1wGM_Tauchen(mew,rho,badp,mu,sd,11,J,3,struct('verbose',0));
    fprintf('mixture weights not summing to one: ran without error [T0], this should not happen \n')
catch ME
    fprintf('mixture weights not summing to one: errors, as it should. The message is: %s \n',ME.message)
end

output.done=1;

end
