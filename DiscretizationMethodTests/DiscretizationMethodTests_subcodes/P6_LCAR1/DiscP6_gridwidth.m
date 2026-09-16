function output=DiscP6_gridwidth(calib,znums,figure_c)
% P6: grid width against grid size, for the three life-cycle AR(1) methods
%
% The same sweep as DiscP2_gridwidth, and this is the block where the width question is sharpest,
% for two reasons.
%
% First, discretizeLifeCycleAR1_FellaGallipoliPan is a Rouwenhorst construction, so sqrt(znum-1) is
% not a default for it, it is a requirement - and its shipped default caps the width at 4, which
% P6 measured as costing it the exact variance match from znum>17 onwards. Its column of this table
% is therefore not a hyperparameter sweep at all: it is a picture of how fast the exactness decays
% as the width is moved away from the one value that works.
%
% Second, sigma_z here is not one number but an age profile, and it rises with age - by a factor of
% three over the life cycle on this calibration. Every method sets its width as
% nSigmas*sigmaz(j), so the grid widens with age automatically. What the sweep varies is the
% multiplier that applies at every age.
%
% Nothing here is asserted. It is a measurement to choose a default from.

fprintf('\n========== P6: grid width against grid size ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho; sigma=calib.vary.sigma;
sigmazT=calib.vary.sigmaz;
% the driftless profile, for the FGP command which takes no mew
sz0=zeros(1,J); sz0(1)=sigma(1);
for j_c=2:J
    sz0(j_c)=sqrt(rho(j_c)^2*sz0(j_c-1)^2+sigma(j_c)^2);
end
fprintf('calibration: sd(z) rises from %2.4f at age 1 to %2.4f at age %i, peaking at %2.4f \n',sigmazT(1),sigmazT(J),J,max(sigmazT))

zlist=[5,9,15,31,51];
wlist=[1.5,2,2.5,3,4,5,7];
nz=length(zlist); nw=length(wlist);
eK=NaN(nz,nw); eG=NaN(nz,nw); eT=NaN(nz,nw);

for z_c=1:nz
    znum=zlist(z_c);
    for w_c=1:nw
        w=wlist(w_c);
        ko=struct(); ko.method='even'; ko.nSigmas=w;
        [zg,pz,j1]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,ko);
        [~,v,~]=MarkovChainMoments_FHorz(gather(zg),gather(pz),gather(j1));
        eK(z_c,w_c)=max(abs(v(:)'-sigmazT.^2));

        go=struct(); go.nSigmas=w;
        [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,go);
        [~,v,~]=MarkovChainMoments_FHorz(gather(zg),gather(pz),gather(j1));
        eG(z_c,w_c)=max(abs(v(:)'-sz0.^2));

        to=struct(); to.nSigmas=w;
        [zg,pz,j1]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,to);
        [~,v,~]=MarkovChainMoments_FHorz(gather(zg),gather(pz),gather(j1));
        eT(z_c,w_c)=max(abs(v(:)'-sigmazT.^2));
    end
end

for tb=1:3
    switch tb
        case 1
            A=eK; nm='worst-age variance error, discretizeLifeCycleAR1_KFTT';
        case 2
            A=eG; nm='worst-age variance error, discretizeLifeCycleAR1_FellaGallipoliPan';
        case 3
            A=eT; nm='worst-age variance error, discretizeLifeCycleAR1_FellaGallipoliPanTauchen';
    end
    fprintf('\n--- %s (rows znum, columns width in sigma_z) --- \n',nm)
    fprintf('%7s',' ');
    for w_c=1:nw
        fprintf('%11.1f',wlist(w_c));
    end
    % ALL THREE COMMANDS DEFAULT DIFFERENTLY, and the three values are worth having side by side:
    %    KFTT          min(sqrt(znum-1),4)   capped, like the stationary Tauchen commands
    %    FGP           sqrt(znum-1)          uncapped, and REQUIRED - the exact variance match is a
    %                                        property of that specific width, so it is not a choice
    %    FGP-Tauchen   min(sqrt(znum-1),3)   capped at THREE, the only 3 anywhere in the family
    % That last one is not a typo in this subcode: it is what the command ships, on the reasoning in
    % its own comment that Tauchen "would anyway typically just put zeros outside +-3 sigma". The
    % sweep below is the measurement of whether 3 is the right place to stop.
    fprintf('%13s%16s \n','default','sqrt(znum-1)');
    for z_c=1:nz
        fprintf('%7i',zlist(z_c));
        for w_c=1:nw
            fprintf('%11.2e',A(z_c,w_c));
        end
        switch tb
            case 1 % KFTT
                fprintf('%13.2f%16.2f \n',min(sqrt(zlist(z_c)-1),4),sqrt(zlist(z_c)-1));
            case 2 % FGP
                fprintf('%13.2f%16.2f \n',sqrt(zlist(z_c)-1),sqrt(zlist(z_c)-1));
            case 3 % FGP-Tauchen
                fprintf('%13.2f%16.2f \n',min(sqrt(zlist(z_c)-1),3),sqrt(zlist(z_c)-1));
        end
    end
end

fprintf('\n--- where the best width sits, against each command''s own default --- \n')
for z_c=1:nz
    [bk,bki]=min(eK(z_c,:));
    [bg,bgi]=min(eG(z_c,:));
    [bt,bti]=min(eT(z_c,:));
    fprintf('znum=%3i: KFTT best at %4.1f (%2.2e), default min(sqrt(znum-1),4)=%5.2f \n',zlist(z_c),wlist(bki),bk,min(sqrt(zlist(z_c)-1),4))
    fprintf('%10s FGP best at %4.1f (%2.2e), default sqrt(znum-1)=%5.2f (required, not a choice) \n',' ',wlist(bgi),bg,sqrt(zlist(z_c)-1))
    fprintf('%10s FGP-Tauchen best at %4.1f (%2.2e), default min(sqrt(znum-1),3)=%5.2f \n',' ',wlist(bti),bt,min(sqrt(zlist(z_c)-1),3))
end
fprintf('\nThe FGP row is the one to read differently: its best width should sit exactly at \n')
fprintf('sqrt(znum-1) and nowhere else, because that is what its construction requires. If the \n')
fprintf('sweep shows a minimum somewhere else, the swept grid simply does not contain the right \n')
fprintf('value - not that a different width is better. \n')

output.eK=eK; output.eG=eG; output.eT=eT;
output.zlist=zlist; output.wlist=wlist;

%% Figure
figure(figure_c)
subplot(1,3,1)
semilogy(wlist,max(eK',10^(-18)))
xlabel('width, in sigma_z'); ylabel('worst-age variance error'); title('KFTT')
legend(arrayfun(@(x) ['znum=',num2str(x)],zlist,'UniformOutput',false),'Location','best')
subplot(1,3,2)
semilogy(wlist,max(eG',10^(-18)))
xlabel('width, in sigma_z'); title('FGP (width is required, not chosen)')
subplot(1,3,3)
semilogy(wlist,max(eT',10^(-18)))
xlabel('width, in sigma_z'); title('FGP-Tauchen')

end
