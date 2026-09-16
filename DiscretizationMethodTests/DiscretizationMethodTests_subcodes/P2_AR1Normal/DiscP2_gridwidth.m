function output=DiscP2_gridwidth(calib,znums,figure_c)
% P2: grid width against grid size, for the three stationary AR(1) methods
%
% WHAT THIS IS FOR. All three methods place their grid at +-width*sigma_z, and they disagree about
% what width should be:
%    Rouwenhorst  sqrt(znum-1)*sigma_z, and it has NO choice - the exact variance match is a
%                 property of that specific width, so it is not a hyperparameter at all
%    Tauchen      Tauchen_q*sigma_z, a free hyperparameter, historically 2 or 3, and until now it
%                 did not grow with znum at all; it now defaults to min(sqrt(znum-1),4)
%    Farmer-Toda  nSigmas*sigma_z, a free hyperparameter, whose default used to be capped at 3 and
%                 is now sqrt(znum-1) with no cap
% The two free hyperparameters do NOT share a default, and the split is the thing to check here.
% Tauchen pays for width in spacing, since a wider grid with the same znum resolves the conditional
% distribution more coarsely, so its default is capped. Farmer-Toda's max-entropy step matches the
% conditional moments whatever the grid is, so width costs it tails and failed solves rather than
% accuracy, and its default is uncapped. sqrt(znum-1) is 7.07 at znum=51 and 10.0 at znum=101, so
% the cap bites hard at the top of the sweep and this subcode measures whether it should.
%
% THE TRADE-OFF, stated before the numbers so the table can be read against it. At width w and znum
% points the spacing is 2*w/(znum-1) standard deviations. Holding w fixed, spacing falls like
% 1/znum. Setting w=sqrt(znum-1), spacing is 2/sqrt(znum-1) and falls only like 1/sqrt(znum) - so a
% width that grows buys tail coverage at the price of resolution, and past some znum the resolution
% loss should dominate. Where that happens, and whether it happens at all before znum=100, is what
% the sweep measures.
%
% Nothing here is asserted. It is a measurement to choose a default from.

fprintf('\n========== P2: grid width against grid size ========== \n')

output=struct();

% TWO CALIBRATIONS, because the answer should depend on persistence and it would be a mistake to
% set a default from one. The truncation term depends on how much mass sits beyond w standard
% deviations, and the spacing term on how finely the grid resolves the conditional distribution
% whose width is sigma rather than sigma_z; the ratio sigma/sigma_z is sqrt(1-rho^2), so a more
% persistent process has a conditional distribution that is narrow relative to its grid. Whether
% the measured optimum moves with rho is exactly the question.
calibnames={'moderate','drift'};
for cal_c=1:2
cn=calibnames{cal_c};
mew=calib.(cn).mew; rho=calib.(cn).rho; sigma=calib.(cn).sigma;
sigmaz=sigma/sqrt(1-rho^2); varz=sigmaz^2; acz=rho;
fprintf('\n=== calibration %s: mew=%g, rho=%g, sigma=%g, so sd(z)=%2.6f, sigma/sigma_z=%2.4f === \n',cn,mew,rho,sigma,sigmaz,sqrt(1-rho^2))

zlist=[5,9,15,31,51,101];
wlist=[1.5,2,2.5,3,4,5,7,10];
nz=length(zlist); nw=length(wlist);
errT=NaN(nz,nw); errF=NaN(nz,nw); acT=NaN(nz,nw); acF=NaN(nz,nw);
errR=NaN(nz,1); acR=NaN(nz,1); wR=NaN(nz,1);

for z_c=1:nz
    znum=zlist(z_c);
    % Rouwenhorst has no width option: its width IS sqrt(znum-1), and that is the whole reason it
    % is exact. It appears in the table as a single reference point per row.
    [zg,pz]=discretizeAR1_Rouwenhorst(mew,rho,sigma,znum,struct());
    zg=gather(zg); pz=gather(pz);
    [~,v,a]=MarkovChainMoments(zg,pz);
    errR(z_c)=abs(v-varz); acR(z_c)=abs(a-acz); wR(z_c)=sqrt(znum-1);
    for w_c=1:nw
        w=wlist(w_c);
        [zg,pz]=discretizeAR1_Tauchen(mew,rho,sigma,znum,w,struct());
        zg=gather(zg); pz=gather(pz);
        [~,v,a]=MarkovChainMoments(zg,pz);
        errT(z_c,w_c)=abs(v-varz); acT(z_c,w_c)=abs(a-acz);
        fo=struct(); fo.method='even'; fo.nSigmas=w; fo.verbose=0;
        [zg,pz]=discretizeAR1_FarmerToda(mew,rho,sigma,znum,fo);
        zg=gather(zg); pz=gather(pz);
        [~,v,a]=MarkovChainMoments(zg,pz);
        errF(z_c,w_c)=abs(v-varz); acF(z_c,w_c)=abs(a-acz);
    end
end

fprintf('\n--- variance error, discretizeAR1_Tauchen (rows znum, columns width in sigma_z) --- \n')
fprintf('%7s',' ');
for w_c=1:nw
    fprintf('%11.1f',wlist(w_c));
end
fprintf('%13s%9s \n','sqrt(znum-1)','spacing');
for z_c=1:nz
    fprintf('%7i',zlist(z_c));
    for w_c=1:nw
        fprintf('%11.2e',errT(z_c,w_c));
    end
    fprintf('%13.2f%9.3f \n',wR(z_c),2*wR(z_c)/(zlist(z_c)-1));
end

fprintf('\n--- variance error, discretizeAR1_FarmerToda, method=even (rows znum, columns width) --- \n')
fprintf('%7s',' ');
for w_c=1:nw
    fprintf('%11.1f',wlist(w_c));
end
fprintf(' \n');
for z_c=1:nz
    fprintf('%7i',zlist(z_c));
    for w_c=1:nw
        fprintf('%11.2e',errF(z_c,w_c));
    end
    fprintf(' \n');
end

fprintf('\n--- autocorrelation error, discretizeAR1_Tauchen (rows znum, columns width) --- \n')
fprintf('%7s',' ');
for w_c=1:nw
    fprintf('%11.1f',wlist(w_c));
end
fprintf(' \n');
for z_c=1:nz
    fprintf('%7i',zlist(z_c));
    for w_c=1:nw
        fprintf('%11.2e',acT(z_c,w_c));
    end
    fprintf(' \n');
end

fprintf('\n--- discretizeAR1_Rouwenhorst, which has no width to choose --- \n')
for z_c=1:nz
    fprintf('znum=%3i: width is sqrt(znum-1)=%5.2f, variance error %2.2e, autocorrelation error %2.2e \n',zlist(z_c),wR(z_c),errR(z_c),acR(z_c))
end

% The reading. For each znum, which width in the sweep was best, and how each method's OWN default
% compares against it. Tauchen's default is capped at 4 and Farmer-Toda's is not, so they must be
% read at different widths - scoring Tauchen at sqrt(znum-1) would be scoring a width it never uses.
fprintf('\n--- where the best width sits, against each method''s own default --- \n')
for z_c=1:nz
    [bt,bti]=min(errT(z_c,:));
    [bf,bfi]=min(errF(z_c,:));
    % what each default would give, interpolated onto the sweep by nearest listed width
    wTdef=min(sqrt(zlist(z_c)-1),4); % discretizeAR1_Tauchen: capped
    wFdef=sqrt(zlist(z_c)-1);        % discretizeAR1_FarmerToda: uncapped
    [~,dT]=min(abs(wlist-wTdef));
    [~,dF]=min(abs(wlist-wFdef));
    fprintf('znum=%3i: Tauchen best at width %4.1f (%2.2e); default min(sqrt(znum-1),4)=%5.2f, nearest swept width %4.1f gives %2.2e \n',zlist(z_c),wlist(bti),bt,wTdef,wlist(dT),errT(z_c,dT))
    fprintf('%10s FarmerToda best at width %4.1f (%2.2e); default sqrt(znum-1)=%5.2f, nearest swept width %4.1f gives %2.2e \n',' ',wlist(bfi),bf,wFdef,wlist(dF),errF(z_c,dF))
end
fprintf('\nRead the Tauchen table down a column to see the truncation floor at a fixed width, and \n')
fprintf('across a row to see the resolution cost of widening. Farmer-Toda matches the conditional \n')
fprintf('moments exactly whatever the grid, so its table should be far flatter - what a bad width \n')
fprintf('costs it is tails and failed solves, not the second moment. \n')

output.(cn).errT=errT; output.(cn).errF=errF; output.(cn).acT=acT; output.(cn).acF=acF;
output.(cn).errR=errR; output.(cn).acR=acR; output.(cn).rho=rho;

% The criterion the default is actually chosen on: the WORSE of the two moments. Reporting the
% variance-best and the autocorrelation-best separately invites picking whichever suits, and the
% two point in opposite directions - the variance wants a narrow grid for spacing, the
% autocorrelation wants a wide one for support.
fprintf('\n--- the width minimising the WORSE of the variance and autocorrelation errors, %s --- \n',cn)
fprintf('    (this block is Tauchen only, so it is scored at the Tauchen default) \n')
for z_c=1:nz
    worst=max(errT(z_c,:),acT(z_c,:));
    [bw,bwi]=min(worst);
    wTdef=min(sqrt(zlist(z_c)-1),4);
    [~,dT]=min(abs(wlist-wTdef));
    fprintf('znum=%3i: best width %4.1f (worse-moment error %2.2e); default min(sqrt(znum-1),4)=%5.2f (uncapped this would be %5.2f); at the nearest swept width %4.1f the worse moment is %2.2e \n',zlist(z_c),wlist(bwi),bw,wTdef,sqrt(zlist(z_c)-1),wlist(dT),worst(dT))
end

%% Figure
figure(figure_c)
subplot(2,2,2*cal_c-1)
semilogy(wlist,max(errT',10^(-18)))
xlabel('half-width, in sigma_z'); ylabel('variance error'); title(['Tauchen, ',cn,' (rho=',num2str(rho),')'])
if cal_c==1
    legend(arrayfun(@(x) ['znum=',num2str(x)],zlist,'UniformOutput',false),'Location','best')
end
subplot(2,2,2*cal_c)
semilogy(wlist,max(acT',10^(-18)))
xlabel('half-width, in sigma_z'); ylabel('autocorrelation error'); title(['Tauchen autocorr, ',cn])

end

output.zlist=zlist; output.wlist=wlist; output.wR=wR; output.calibnames=calibnames;

end
