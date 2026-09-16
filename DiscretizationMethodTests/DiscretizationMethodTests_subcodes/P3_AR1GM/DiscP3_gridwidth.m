function output=DiscP3_gridwidth(calib,znums,figure_c)
% P3: grid width against grid size, for a gaussian-mixture AR(1)
%
% The same sweep as DiscP2_gridwidth, on a process with fat tails. It is here rather than only in
% P2 because the whole question - how wide should the grid be - has a different answer when the
% innovation is a mixture: the tails beyond three standard deviations carry real mass, so a width
% that is fine for a gaussian may not be. This is the same mechanism P4 found for stochastic
% volatility, where excess kurtosis converged to the TRUNCATED distribution's value rather than the
% truth, and the width, not the number of points, was the lever.
%
% Rouwenhorst does not exist for a mixture, so the two methods here are discretizeAR1wGM_FarmerToda
% and discretizeAR1wGM_Tauchen. The moment that matters is excess kurtosis, not variance:
% Farmer-Toda matches the conditional moments it targets whatever the grid, so variance says little
% about the width, where kurtosis says everything.
%
% Nothing here is asserted. It is a measurement to choose a default from.

fprintf('\n========== P3: grid width against grid size ========== \n')

output=struct();
mew=calib.mew; rho=calib.rho;
mixprobs=calib.mixprobs; mu_i=calib.mu; sigma_i=calib.sigma;
varzT=calib.z.var; skewT=calib.z.skew; exkurtT=calib.z.exkurt;
fprintf('calibration: rho=%g, and z has variance %2.6f, skewness %2.4f, excess kurtosis %2.4f \n',rho,varzT,skewT,exkurtT)

zlist=[5,9,15,31,51,101];
wlist=[1.5,2,2.5,3,4,5,7,10];
nz=length(zlist); nw=length(wlist);
kF=NaN(nz,nw); kT=NaN(nz,nw); vF=NaN(nz,nw); vT=NaN(nz,nw);

for z_c=1:nz
    znum=zlist(z_c);
    for w_c=1:nw
        w=wlist(w_c);
        fo=struct(); fo.method='even'; fo.nSigmas=w; fo.nMoments=4; fo.verbose=0;
        [zg,pz]=discretizeAR1wGM_FarmerToda(mew,rho,mixprobs,mu_i,sigma_i,znum,fo);
        zg=gather(zg); pz=gather(pz);
        [~,~,~,sd]=MarkovChainMoments(zg,pz);
        m=sum(sd(:).*zg(:)); v=sum(sd(:).*(zg(:)-m).^2);
        vF(z_c,w_c)=abs(v-varzT);
        kF(z_c,w_c)=abs(sum(sd(:).*(zg(:)-m).^4)/v^2-3-exkurtT);
        [zg,pz]=discretizeAR1wGM_Tauchen(mew,rho,mixprobs,mu_i,sigma_i,znum,w,struct());
        zg=gather(zg); pz=gather(pz);
        [~,~,~,sd]=MarkovChainMoments(zg,pz);
        m=sum(sd(:).*zg(:)); v=sum(sd(:).*(zg(:)-m).^2);
        vT(z_c,w_c)=abs(v-varzT);
        kT(z_c,w_c)=abs(sum(sd(:).*(zg(:)-m).^4)/v^2-3-exkurtT);
    end
end

for tb=1:4
    switch tb
        case 1
            A=kF; nm='excess kurtosis error, discretizeAR1wGM_FarmerToda';
        case 2
            A=kT; nm='excess kurtosis error, discretizeAR1wGM_Tauchen';
        case 3
            A=vF; nm='variance error, discretizeAR1wGM_FarmerToda';
        case 4
            A=vT; nm='variance error, discretizeAR1wGM_Tauchen';
    end
    fprintf('\n--- %s (rows znum, columns width in sigma_z) --- \n',nm)
    fprintf('%7s',' ');
    for w_c=1:nw
        fprintf('%11.1f',wlist(w_c));
    end
    % NEITHER command defaults to sqrt(znum-1), so the default column follows the command this table
    % is for. discretizeAR1wGM_FarmerToda keeps Toda's persistence branch: sqrt(2*(znum-1)) when
    % rho<=1-2/(znum-1) and sqrt(znum-1) otherwise, so it is the one width in the family that can be
    % WIDER than sqrt(znum-1). discretizeAR1wGM_Tauchen uses the tail-mass rule, which solves for the
    % width at which the mixture leaves the same mass outside the grid that a normal leaves beyond
    % four standard deviations - so it depends on the mixture, not on znum, and cannot be a column
    % here. What it chose for this calibration is reported by DiscP3_AR1wGM_Tauchen.
    fprintf('%13s%16s \n','default','sqrt(znum-1)');
    for z_c=1:nz
        fprintf('%7i',zlist(z_c));
        for w_c=1:nw
            fprintf('%11.2e',A(z_c,w_c));
        end
        if tb==2 || tb==4 % the discretizeAR1wGM_Tauchen tables
            fprintf('%13s%16.2f \n','tail-mass',sqrt(zlist(z_c)-1));
        else              % the discretizeAR1wGM_FarmerToda tables
            if rho <= 1-2/(zlist(z_c)-1)
                fprintf('%13.2f%16.2f \n',sqrt(2*(zlist(z_c)-1)),sqrt(zlist(z_c)-1));
            else
                fprintf('%13.2f%16.2f \n',sqrt(zlist(z_c)-1),sqrt(zlist(z_c)-1));
            end
        end
    end
end

fprintf('\n--- where the best width sits for excess kurtosis, against each command''s own default --- \n')
for z_c=1:nz
    [bf,bfi]=min(kF(z_c,:));
    [bt,bti]=min(kT(z_c,:));
    if rho <= 1-2/(zlist(z_c)-1)
        wFdef=sqrt(2*(zlist(z_c)-1)); brname='sqrt(2*(znum-1))';
    else
        wFdef=sqrt(zlist(z_c)-1); brname='sqrt(znum-1)   ';
    end
    fprintf('znum=%3i: FarmerToda best at width %4.1f (%2.2e), default %s=%5.2f; Tauchen best at width %4.1f (%2.2e), default from the tail-mass rule \n',zlist(z_c),wlist(bfi),bf,brname,wFdef,wlist(bti),bt)
end
fprintf('\nIf the best width for kurtosis keeps rising with znum rather than settling, that is the fat \n')
fprintf('tails talking, and it is an argument for a width rule that grows - which sqrt(znum-1) does. \n')
fprintf('If it settles, a constant beyond some point would be the better rule. \n')

output.kF=kF; output.kT=kT; output.vF=vF; output.vT=vT;
output.zlist=zlist; output.wlist=wlist;

%% Figure
figure(figure_c)
subplot(1,2,1)
semilogy(wlist,max(kF',10^(-18)))
xlabel('grid half-width, in sigma_z'); ylabel('excess kurtosis error'); title('wGM\_FarmerToda')
legend(arrayfun(@(x) ['znum=',num2str(x)],zlist,'UniformOutput',false),'Location','best')
subplot(1,2,2)
semilogy(wlist,max(kT',10^(-18)))
xlabel('grid half-width, in sigma_z'); ylabel('excess kurtosis error'); title('wGM\_Tauchen')

end
