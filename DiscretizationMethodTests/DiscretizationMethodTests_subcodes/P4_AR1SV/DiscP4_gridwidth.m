function output=DiscP4_gridwidth(calib,znums,figure_c)
% P4: grid width against grid size, for stochastic volatility
%
% WHY THIS EXISTS SEPARATELY FROM THE SWEEP IN DiscP4_crosstests. That one varies nSigmas at a
% single grid size, which was enough to be the B24 regression - it proves the option is live and
% that the grid spans exactly +-nSigmas*sigmaz. It is NOT enough to choose a default, because the
% question a default has to answer is how the best width moves with znum, and a single znum cannot
% see that. This sweeps both.
%
% WHAT THE MEASUREMENT IS FOR. discretizeAR1wSV_FarmerToda defaults nSigmas to
% min(sqrt((znum-1)/2),2) - capped at 2, and the narrowest default anywhere in the toolkit. At that
% width P4 measures the excess kurtosis of z at -0.52 against a truth of +0.83: the wrong sign. The
% question is whether the cap should go, and to what.
%
% EXCESS KURTOSIS IS THE MOMENT, not variance. Farmer-Toda's maximum entropy step forces the
% conditional moments it targets whatever the grid, so the variance is nearly width-insensitive and
% says little; the fat tails that stochastic volatility exists to produce live beyond the grid edge,
% and only kurtosis sees them. Both are reported, and the contrast between the two tables is itself
% the argument for reading the kurtosis one.
%
% Nothing here is asserted. It is a measurement to choose a default from.

fprintf('\n========== P4: grid width against grid size ========== \n')

output=struct();
rho=calib.rho; phi=calib.phi; sigmau=calib.sigmau; sigmae=calib.sigmae;
varzT=calib.z.var; exkurtT=calib.z.exkurt;
fprintf('calibration: rho=%g, phi=%g, and z has variance %2.6f, excess kurtosis %2.4f \n',rho,phi,varzT,exkurtT)

xnum=9;   % held fixed: this is a question about the z grid, and P4's main sweep already varies xnum
zlist=[5,9,15,31,51,101];
wlist=[1.5,2,2.5,3,4,5,7,10];
nz=length(zlist); nw=length(wlist);
kF=NaN(nz,nw); kT=NaN(nz,nw); vF=NaN(nz,nw); vT=NaN(nz,nw); fbF=NaN(nz,nw);

for z_c=1:nz
    znum=zlist(z_c);
    if xnum*znum>calib.Ncap
        fprintf('znum=%i: SKIPPED, xnum*znum=%i exceeds the ceiling of %i \n',znum,xnum*znum,calib.Ncap)
        continue
    end
    for w_c=1:nw
        w=wlist(w_c);
        opts=struct(); opts.nSigmas=w; opts.verbose=0;
        [zg,pz,oo]=discretizeAR1wSV_FarmerToda(rho,phi,sigmau,sigmae,xnum,znum,opts);
        zg=gather(zg); pz=gather(pz);
        sd=ones(xnum*znum,1)/(xnum*znum);
        for i_c=1:10000
            sn=pz'*sd;
            if max(abs(sn-sd))<10^(-14), sd=sn; break; end
            sd=sn;
        end
        sd=sd/sum(sd);
        Pxz=reshape(sd,[xnum,znum]); zb=zg(xnum+1:end); zm=sum(Pxz,1)';
        m=sum(zm.*zb); v=sum(zm.*(zb-m).^2);
        vF(z_c,w_c)=abs(v-varzT);
        kF(z_c,w_c)=abs(sum(zm.*(zb-m).^4)/v^2-3-exkurtT);
        fbF(z_c,w_c)=mean(oo.nMoments_grid(:)<2);

        [zg,pz]=discretizeAR1wSV_Tauchen(rho,phi,sigmau,sigmae,xnum,znum,w,struct());
        zg=gather(zg); pz=gather(pz);
        sd=ones(xnum*znum,1)/(xnum*znum);
        for i_c=1:10000
            sn=pz'*sd;
            if max(abs(sn-sd))<10^(-14), sd=sn; break; end
            sd=sn;
        end
        sd=sd/sum(sd);
        Pxz=reshape(sd,[xnum,znum]); zb=zg(xnum+1:end); zm=sum(Pxz,1)';
        m=sum(zm.*zb); v=sum(zm.*(zb-m).^2);
        vT(z_c,w_c)=abs(v-varzT);
        kT(z_c,w_c)=abs(sum(zm.*(zb-m).^4)/v^2-3-exkurtT);
    end
end

for tb=1:5
    switch tb
        case 1
            A=kF; nm='excess kurtosis error, discretizeAR1wSV_FarmerToda';
        case 2
            A=kT; nm='excess kurtosis error, discretizeAR1wSV_Tauchen';
        case 3
            A=vF; nm='variance error, discretizeAR1wSV_FarmerToda';
        case 4
            A=vT; nm='variance error, discretizeAR1wSV_Tauchen';
        case 5
            A=100*fbF; nm='percent of rows matching fewer than 2 conditional moments, wSV_FarmerToda';
    end
    fprintf('\n--- %s (rows znum, columns width in sigma_z, xnum=%i) --- \n',nm,xnum)
    fprintf('%7s',' ');
    for w_c=1:nw
        fprintf('%11.1f',wlist(w_c));
    end
    % The two commands ship DIFFERENT default widths, so the default column has to follow whichever
    % command this table is for. discretizeAR1wSV_FarmerToda uses sqrt(znum-1) uncapped (its entropy
    % step matches the conditional moments whatever the grid is, so width costs it tails and failed
    % solves rather than accuracy); discretizeAR1wSV_Tauchen caps at 4 (it pays for width in spacing).
    fprintf('%10s%16s \n','default','sqrt(znum-1)');
    for z_c=1:nz
        fprintf('%7i',zlist(z_c));
        for w_c=1:nw
            fprintf('%11.2e',A(z_c,w_c));
        end
        if tb==2 || tb==4 % the discretizeAR1wSV_Tauchen tables
            fprintf('%10.2f%16.2f \n',min(sqrt(zlist(z_c)-1),4),sqrt(zlist(z_c)-1));
        else              % the discretizeAR1wSV_FarmerToda tables
            fprintf('%10.2f%16.2f \n',sqrt(zlist(z_c)-1),sqrt(zlist(z_c)-1));
        end
    end
end

% CELLS WHERE THE SOLVE FAILED ARE EXCLUDED, and that is not a detail. On the run of 2026-08-27 the
% naive minimum for FarmerToda at znum=9 was width 5 with a kurtosis error of 7.7e-02 - produced by
% a solve that fell back to one moment on 61.7%% of rows. A row that matched only the conditional
% mean has an unconstrained conditional variance, so its contribution to the kurtosis is an accident,
% and reporting it as the best width would recommend a setting on the strength of a failure. Only
% cells below the threshold are eligible, and how many were dropped is printed rather than hidden.
fbmax=10; % percent
fprintf('\n--- where the best width sits for excess kurtosis --- \n')
fprintf('(FarmerToda cells with more than %g%% of rows falling back to one moment are excluded: their \n',fbmax)
fprintf(' kurtosis is an artefact of the failed solve, not a property of the width) \n')
for z_c=1:nz
    if all(isnan(kF(z_c,:)))
        continue
    end
    elig=(100*fbF(z_c,:)<=fbmax);
    kFe=kF(z_c,:); kFe(~elig)=NaN;
    [bt,bti]=min(kT(z_c,:));
    if any(elig)
        [bf,bfi]=min(kFe);
        fprintf('znum=%3i: FarmerToda best at %4.1f (%2.2e, from %i of %i eligible widths), Tauchen best at %4.1f (%2.2e); shipped defaults are FarmerToda sqrt(znum-1)=%5.2f and Tauchen min(sqrt(znum-1),4)=%4.2f \n',zlist(z_c),wlist(bfi),bf,sum(elig),nw,wlist(bti),bt,sqrt(zlist(z_c)-1),min(sqrt(zlist(z_c)-1),4))
    else
        fprintf('znum=%3i: FarmerToda has NO eligible width - the solve falls back on more than %g%% of rows at every width tested; Tauchen best at %4.1f (%2.2e) \n',zlist(z_c),fbmax,wlist(bti),bt)
    end
end
fprintf('\nThe shipped default caps at 2 and so never exceeds it at any znum. If the best width keeps \n')
fprintf('rising with znum, the cap is costing tails at every grid size and should go; if it settles, \n')
fprintf('the question is only what it should settle at. Read the fallback table alongside: a width \n')
fprintf('that is too wide for the number of points shows up there before it shows up in the moments. \n')

output.kF=kF; output.kT=kT; output.vF=vF; output.vT=vT; output.fbF=fbF;
output.zlist=zlist; output.wlist=wlist; output.xnum=xnum;

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(wlist,max(kF',10^(-18)))
xlabel('grid half-width, in sigma_z'); ylabel('excess kurtosis error'); title('wSV\_FarmerToda')
legend(arrayfun(@(x) ['znum=',num2str(x)],zlist,'UniformOutput',false),'Location','best')
subplot(1,2,2)
semilogy(wlist,max(kT',10^(-18)))
xlabel('grid half-width, in sigma_z'); ylabel('excess kurtosis error'); title('wSV\_Tauchen')

end
