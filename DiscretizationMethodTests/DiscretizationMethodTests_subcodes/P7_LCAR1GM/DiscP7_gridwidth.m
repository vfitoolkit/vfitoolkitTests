function output=DiscP7_gridwidth(calib,znums,figure_c)
% P7: grid width against grid size, for the life-cycle gaussian-mixture methods
%
% The third and last unmeasured width cap in the toolkit. discretizeLifeCycleAR1wGM_KFTT defaults
% nSigmas to min(sqrt(znum-1),4); discretizeLifeCycleAR1wGM_Tauchen takes Tauchen_q and defaults to
% min(sqrt(znum-1),4) as well, matching discretizeAR1_Tauchen.
%
% What the earlier blocks established, and what this one is checking against: the Farmer-Toda family
% wants width and is limited only by its own solve collapsing (P4 removed a cap of 2 that was giving
% the wrong SIGN on excess kurtosis); the Tauchen family pays for width in grid spacing and its
% optimum settles around 3 to 4 (P2, and P6 for the life-cycle version). A mixture has fatter tails
% than a normal, so if anything wants MORE width than its gaussian counterpart it should be these.
%
% EXCESS KURTOSIS IS THE MOMENT, for the same reason as P3 and P4: KFTT forces the conditional
% moments it targets whatever the grid, so its variance says little about width. Both are reported.
%
% Nothing here is asserted. It is a measurement.

fprintf('\n========== P7: grid width against grid size ========== \n')

output=struct();
J=calib.J;
mew=calib.vary.mew; rho=calib.vary.rho;
p=calib.vary.mixprobs_i; mu=calib.vary.mu_i; sd=calib.vary.sigma_i;
varzT=calib.vary.varz; exkurtT=calib.vary.exkurtz;
fprintf('calibration: excess kurtosis of z runs %2.4f to %2.4f \n',exkurtT(1),exkurtT(J))

zlist=[5,9,15,31,51];
% The sweep runs past 7 deliberately. On the run of 2026-08-28 the wGM_Tauchen kurtosis error was
% still falling at width 7 (0.056 at znum=51 and dropping), so the optimum was not bracketed and all
% that could be said was ">=7". The tail-mass rule that now sets the default predicts about 6.3 for
% this calibration, so whether the true optimum is near 6 or out at 10 is exactly what distinguishes
% the rule being right from it running low.
wlist=[1.5,2,2.5,3,4,5,7,10,14];
nz=length(zlist); nw=length(wlist);
kK=NaN(nz,nw); kT=NaN(nz,nw); vK=NaN(nz,nw); vT=NaN(nz,nw); fbK=NaN(nz,nw);

for z_c=1:nz
    znum=zlist(z_c);
    for w_c=1:nw
        w=wlist(w_c);
        ko=struct(); ko.verbose=0; ko.nSigmas=w;
        [zg,pz,j1,oo]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,ko);
        zg=gather(zg); pz=gather(pz); j1=gather(j1);
        [~,v,~,dd]=MarkovChainMoments_FHorz(zg,pz,j1);
        ek=zeros(1,J);
        for j_c=1:J
            d=gather(dd(:,j_c)); d=d/sum(d); g=zg(:,j_c);
            m1=sum(d.*g); vv=sum(d.*(g-m1).^2);
            ek(j_c)=sum(d.*(g-m1).^4)/vv^2-3;
        end
        vK(z_c,w_c)=max(abs(v(:)'-varzT));
        kK(z_c,w_c)=max(abs(ek-exkurtT));
        fbK(z_c,w_c)=mean(oo.nMoments_grid(:)<4);

        to=struct(); to.verbose=0;
        [zg,pz,j1]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,w,to);
        zg=gather(zg); pz=gather(pz); j1=gather(j1);
        [~,v,~,dd]=MarkovChainMoments_FHorz(zg,pz,j1);
        ek=zeros(1,J);
        for j_c=1:J
            d=gather(dd(:,j_c)); d=d/sum(d); g=zg(:,j_c);
            m1=sum(d.*g); vv=sum(d.*(g-m1).^2);
            ek(j_c)=sum(d.*(g-m1).^4)/vv^2-3;
        end
        vT(z_c,w_c)=max(abs(v(:)'-varzT));
        kT(z_c,w_c)=max(abs(ek-exkurtT));
    end
end

for tb=1:5
    switch tb
        case 1
            A=kK; nm='worst-age excess kurtosis error, discretizeLifeCycleAR1wGM_KFTT';
        case 2
            A=kT; nm='worst-age excess kurtosis error, discretizeLifeCycleAR1wGM_Tauchen';
        case 3
            A=vK; nm='worst-age variance error, discretizeLifeCycleAR1wGM_KFTT';
        case 4
            A=vT; nm='worst-age variance error, discretizeLifeCycleAR1wGM_Tauchen';
        case 5
            A=100*fbK; nm='percent of conditional distributions matching fewer than 4 moments, KFTT';
    end
    fprintf('\n--- %s (rows znum, columns width in sd(z)) --- \n',nm)
    fprintf('%7s',' ');
    for w_c=1:nw
        fprintf('%11.1f',wlist(w_c));
    end
    % The two commands do NOT share a default. KFTT uses min(sqrt(znum-1),4). Tauchen uses the
    % tail-mass rule, which solves for the width at which the mixture leaves the same mass outside
    % the grid that a normal leaves beyond four standard deviations - so it depends on the mixture
    % and on the age, not on znum, and cannot be printed as a column here. What it actually chose
    % for this calibration is reported by DiscP7_crosstests.
    fprintf('%10s%16s \n','default','sqrt(znum-1)');
    for z_c=1:nz
        fprintf('%7i',zlist(z_c));
        for w_c=1:nw
            fprintf('%11.2e',A(z_c,w_c));
        end
        if tb==2 || tb==4 % the Tauchen tables: the default is the tail-mass rule, not a function of znum
            fprintf('%10s%16.2f \n','tail-mass',sqrt(zlist(z_c)-1));
        else              % the KFTT tables
            fprintf('%10.2f%16.2f \n',min(sqrt(zlist(z_c)-1),4),sqrt(zlist(z_c)-1));
        end
    end
end

% Cells where the solve fell back are excluded from the KFTT best-width line, for the reason
% DiscP4_gridwidth records: a row that matched only its conditional mean has an unconstrained
% variance, so its contribution to the kurtosis is an artefact of the failure.
fbmax=10;
fprintf('\n--- where the best width sits for excess kurtosis --- \n')
fprintf('(KFTT cells with more than %g%% fallback are excluded) \n',fbmax)
for z_c=1:nz
    elig=(100*fbK(z_c,:)<=fbmax);
    kKe=kK(z_c,:); kKe(~elig)=NaN;
    [bt,bti]=min(kT(z_c,:));
    if any(elig)
        [bk,bki]=min(kKe);
        fprintf('znum=%3i: KFTT best at %4.1f (%2.2e), Tauchen best at %4.1f (%2.2e); KFTT default min(sqrt(znum-1),4)=%4.2f, Tauchen default from the tail-mass rule \n',zlist(z_c),wlist(bki),bk,wlist(bti),bt,min(sqrt(zlist(z_c)-1),4))
    else
        fprintf('znum=%3i: KFTT has NO eligible width at this fallback threshold; Tauchen best at %4.1f (%2.2e) \n',zlist(z_c),wlist(bti),bt)
    end
end

output.kK=kK; output.kT=kT; output.vK=vK; output.vT=vT; output.fbK=fbK;
output.zlist=zlist; output.wlist=wlist;

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(wlist,max(kK',10^(-18)))
xlabel('width, in sd(z)'); ylabel('worst-age excess kurtosis error'); title('wGM\_KFTT')
legend(arrayfun(@(x) ['znum=',num2str(x)],zlist,'UniformOutput',false),'Location','best')
subplot(1,2,2)
semilogy(wlist,max(kT',10^(-18)))
xlabel('width, in sd(z)'); ylabel('worst-age excess kurtosis error'); title('wGM\_Tauchen')

end
