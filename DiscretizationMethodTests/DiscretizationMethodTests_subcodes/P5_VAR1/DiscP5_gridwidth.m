function output=DiscP5_gridwidth(calib,znums,figure_c)
% P5: grid width against grid size, for a VAR(1)
%
% The same sweep as DiscP2_gridwidth, in the multivariate case, and it asks one extra question that
% the univariate blocks cannot. discretizeVAR1_Tauchen applies Tauchen_q to EACH variable's own
% standard deviation, so a single width parameter has to serve variables that may differ in spread;
% discretizeVAR1_FarmerToda builds its grid in a rotated basis and maps back, so nSigmas is a width
% in the transformed space rather than in z. Whether a common default width means the same thing in
% both is not obvious, and this measures it.
%
% There is no Rouwenhorst here - it has no multivariate form in the toolkit - so the reference point
% that anchors the univariate tables is missing, and the sweep stands on its own.
%
% Nothing here is asserted. It is a measurement to choose a default from.

fprintf('\n========== P5: grid width against grid size ========== \n')

output=struct();
Mew=calib.full.Mew; Rho=calib.full.Rho; SigmaSq=calib.full.SigmaSq;
SigmaSqzT=calib.full.SigmaSqz;
fprintf('calibration: M=2, sd(z) = [%2.4f %2.4f], cross-covariance %2.6f \n',calib.full.sigmaz(1),calib.full.sigmaz(2),SigmaSqzT(1,2))

zlist=[5,9,15,31];   % znum^2 states, so 31 is 961 and the ceiling of %i bites above that
wlist=[1.5,2,2.5,3,4,5,7];
nz=length(zlist); nw=length(wlist);
vT=NaN(nz,nw); cT=NaN(nz,nw); vF=NaN(nz,nw); cF=NaN(nz,nw);

for z_c=1:nz
    znum=zlist(z_c);
    if znum^2>calib.Ncap
        fprintf('znum=%i: SKIPPED, znum^2=%i exceeds the ceiling of %i \n',znum,znum^2,calib.Ncap)
        continue
    end
    for w_c=1:nw
        w=wlist(w_c);
        [zg,pz]=discretizeVAR1_Tauchen(Mew,Rho,SigmaSq,znum,w,struct());
        zg=gather(zg); pz=gather(pz);
        sd=ones(znum^2,1)/(znum^2);
        for i_c=1:10000
            sn=pz'*sd;
            if max(abs(sn-sd))<10^(-14), sd=sn; break; end
            sd=sn;
        end
        sd=sd/sum(sd);
        zv=CreateGridvals([znum;znum],zg,1);
        m=(sd'*zv)'; dz=zv-m'; V=(dz.*sd)'*dz;
        vT(z_c,w_c)=max(abs(diag(V)-diag(SigmaSqzT)));
        cT(z_c,w_c)=abs(V(1,2)-SigmaSqzT(1,2));

        fo=struct(); fo.method='even'; fo.nSigmas=w; fo.parallel=1; fo.verbose=0;
        [zg,pz]=discretizeVAR1_FarmerToda(Mew,Rho,SigmaSq,znum,fo);
        zg=gather(zg); pz=gather(pz);
        sd=ones(znum^2,1)/(znum^2);
        for i_c=1:10000
            sn=pz'*sd;
            if max(abs(sn-sd))<10^(-14), sd=sn; break; end
            sd=sn;
        end
        sd=sd/sum(sd);
        m=(sd'*zg)'; dz=zg-m'; V=(dz.*sd)'*dz;
        vF(z_c,w_c)=max(abs(diag(V)-diag(SigmaSqzT)));
        cF(z_c,w_c)=abs(V(1,2)-SigmaSqzT(1,2));
    end
end

for tb=1:4
    switch tb
        case 1
            A=vT; nm='variance error, discretizeVAR1_Tauchen';
        case 2
            A=cT; nm='cross-covariance error, discretizeVAR1_Tauchen';
        case 3
            A=vF; nm='variance error, discretizeVAR1_FarmerToda';
        case 4
            A=cF; nm='cross-covariance error, discretizeVAR1_FarmerToda';
    end
    fprintf('\n--- %s (rows znum, columns width in sigma_z) --- \n',nm)
    fprintf('%7s',' ');
    for w_c=1:nw
        fprintf('%11.1f',wlist(w_c));
    end
    % THE TWO COMMANDS ARE NOT SYMMETRIC HERE. discretizeVAR1_FarmerToda defaults nSigmas to
    % sqrt(znum(1)-1). discretizeVAR1_Tauchen has NO default at all - Tauchen_q is a required
    % argument, and passing [] makes it an empty vector which then fails its own length check
    % against the number of variables in the VAR. So there is a default width to report for one of
    % these and not the other, and the column says so rather than printing a number that is not used.
    fprintf('%13s%16s \n','default','sqrt(znum-1)');
    for z_c=1:nz
        fprintf('%7i',zlist(z_c));
        for w_c=1:nw
            fprintf('%11.2e',A(z_c,w_c));
        end
        if tb==1 || tb==2 % the discretizeVAR1_Tauchen tables
            fprintf('%13s%16.2f \n','(none)',sqrt(zlist(z_c)-1));
        else              % the discretizeVAR1_FarmerToda tables
            fprintf('%13.2f%16.2f \n',sqrt(zlist(z_c)-1),sqrt(zlist(z_c)-1));
        end
    end
end

fprintf('\n--- where the best width sits, against each command''s own default --- \n')
for z_c=1:nz
    if all(isnan(vT(z_c,:)))
        continue
    end
    [bt,bti]=min(vT(z_c,:));
    [bf,bfi]=min(vF(z_c,:));
    fprintf('znum=%3i: Tauchen best at width %4.1f (%2.2e), and it has no default - Tauchen_q is required; FarmerToda best at width %4.1f (%2.2e), default sqrt(znum-1)=%5.2f \n',zlist(z_c),wlist(bti),bt,wlist(bfi),bf,sqrt(zlist(z_c)-1))
end

output.vT=vT; output.cT=cT; output.vF=vF; output.cF=cF;
output.zlist=zlist; output.wlist=wlist;

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(wlist,max(vT',10^(-18)))
xlabel('grid half-width, in sigma_z'); ylabel('variance error'); title('VAR1\_Tauchen')
legend(arrayfun(@(x) ['znum=',num2str(x)],zlist,'UniformOutput',false),'Location','best')
subplot(1,2,2)
semilogy(wlist,max(vF',10^(-18)))
xlabel('grid half-width, in sigma_z'); ylabel('variance error'); title('VAR1\_FarmerToda')

end
