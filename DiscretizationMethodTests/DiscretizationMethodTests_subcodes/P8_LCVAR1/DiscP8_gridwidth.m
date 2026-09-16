function output=DiscP8_gridwidth(calib,figure_c)
% P8: grid width against grid size, for the life-cycle VAR(1)
%
% THE CAP HERE IS THREE, not four. discretizeLifeCycleVAR1_Tauchen defaults nSigmas to
% min(sqrt(znum-1),3), which it shares with discretizeLifeCycleAR1_FellaGallipoliPanTauchen and with
% nothing else in the family - every stationary Tauchen command caps at 4, Farmer-Toda and
% Fella-Gallipoli-Pan do not cap at all, and the two mixture Tauchen commands use the tail-mass rule.
% So three is the odd value out, and this subcode is the measurement of whether it belongs.
%
% WHY THE ANSWER NEED NOT BE THE SAME AS P2's. The lever P2 identified is that Tauchen pays for width
% in spacing: at width w and znum points the spacing is 2w/(znum-1) standard deviations, so widening
% at fixed znum resolves the conditional distribution more coarsely. In a VAR that trade-off is
% sharper, because the joint state count is znum^M rather than znum - buying resolution costs
% quadratically here where it costs linearly in P2. That alone is an argument for a tighter cap.
%
% Nothing here is asserted. It is a measurement to choose a default from.

fprintf('\n========== P8: grid width against grid size ========== \n')

output=struct();
J=calib.J; M=calib.M;
Mew_J=calib.vary.Mew_J; Rho_J=calib.vary.Rho_J; SigmaSq_J=calib.vary.SigmaSq_J;
sigmazT=calib.vary.sigmaz; mewzT=calib.vary.mewz; ccT=calib.vary.crosscorr;

zlist=[5,7,9,13,17];
wlist=[1.5,2,2.5,3,4,5,7];
nz=length(zlist); nw=length(wlist);
vE=NaN(nz,nw); mE=NaN(nz,nw); cE=NaN(nz,nw);

for z_c=1:nz
    znum=zlist(z_c);
    for w_c=1:nw
        to=struct(); to.verbose=0; to.nSigmas=wlist(w_c)*ones(M,1);
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum,J,to);
        z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);
        d=jequaloneDistz;
        wm=0; wv=0; wc=0;
        for j_c=1:J
            zvals=CreateGridvals(znum*ones(M,1),z_grid_J(:,j_c),1);
            mj=(d'*zvals)';
            dev=zvals-mj';
            Vj=(dev.*d)'*dev;
            wm=max(wm,max(abs(mj-mewzT(:,j_c))));
            wv=max(wv,max(abs(diag(Vj)-sigmazT(:,j_c).^2)));
            wc=max(wc,abs(Vj(1,2)/sqrt(Vj(1,1)*Vj(2,2))-ccT(j_c)));
            if j_c<J
                d=pi_z_J(:,:,j_c)'*d;
            end
        end
        mE(z_c,w_c)=wm; vE(z_c,w_c)=wv; cE(z_c,w_c)=wc;
    end
end

for tb=1:3
    switch tb
        case 1
            A=vE; nm='worst-age variance error';
        case 2
            A=cE; nm='worst-age cross-correlation error';
        case 3
            A=mE; nm='worst-age mean error';
    end
    fprintf('\n--- %s, discretizeLifeCycleVAR1_Tauchen (rows znum per variable, columns width in sd) --- \n',nm)
    fprintf('%7s',' ');
    for w_c=1:nw
        fprintf('%11.1f',wlist(w_c));
    end
    fprintf('%13s%16s \n','default','sqrt(znum-1)');
    for z_c=1:nz
        fprintf('%7i',zlist(z_c));
        for w_c=1:nw
            fprintf('%11.2e',A(z_c,w_c));
        end
        fprintf('%13.2f%16.2f \n',min(sqrt(zlist(z_c)-1),3),sqrt(zlist(z_c)-1));
    end
end

fprintf('\n--- where the best width sits, against the shipped default of min(sqrt(znum-1),3) --- \n')
for z_c=1:nz
    [bv,bvi]=min(vE(z_c,:));
    [bc,bci]=min(cE(z_c,:));
    wdef=min(sqrt(zlist(z_c)-1),3);
    [~,di]=min(abs(wlist-wdef));
    fprintf('znum=%3i: variance best at width %4.1f (%2.2e), cross-correlation best at %4.1f (%2.2e); default %4.2f gives %2.2e on variance \n',zlist(z_c),wlist(bvi),bv,wlist(bci),bc,wdef,vE(z_c,di))
end
fprintf('\nRead down a column for the truncation floor at a fixed width, and across a row for the \n')
fprintf('resolution cost of widening. If the best width keeps rising with znum the cap is costing \n')
fprintf('accuracy at every grid size; if it settles at or below 3, the cap is in the right place. \n')
fprintf('Note the joint state count is znum^2 here, so a column of this table costs quadratically \n')
fprintf('more than the same column in P2 - which is itself part of the argument for a tighter cap. \n')

output.vE=vE; output.cE=cE; output.mE=mE; output.zlist=zlist; output.wlist=wlist;

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(wlist,max(vE',1e-18))
xlabel('width in standard deviations'); ylabel('worst-age variance error')
title('LifeCycleVAR1\_Tauchen: variance vs width')
legend(arrayfun(@(z) ['znum=',num2str(z)],zlist,'UniformOutput',false),'Location','best')
subplot(1,2,2)
semilogy(wlist,max(cE',1e-18))
xlabel('width in standard deviations'); ylabel('worst-age cross-correlation error')
title('cross-correlation vs width')
sgtitle('P8: grid width')

end
