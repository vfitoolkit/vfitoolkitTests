function output=DiscP4_crosstests(calib,znums)
% P4: the identities that hold within the stochastic-volatility block
%
% No figures, per the convention in the other banks.

fprintf('\n========== P4: cross-tests ========== \n')

output=struct();
rho=calib.rho; phi=calib.phi; sigmau=calib.sigmau;
xnum=5; znum=15;

%% 1. The volatility block must be an ordinary gaussian AR(1)
% x is a gaussian AR(1) by construction, so the x marginal of the joint chain has exact truth:
% mean xBar, variance sigmae^2/(1-phi^2), autocorrelation phi. This also pins the ORDERING - if x
% did not vary fastest in the joint index, the reshape below would pick out the wrong marginal and
% these moments would be wrong rather than merely imprecise.
for sigmae=[0.3,0.1]
    optsF=struct(); optsF.nSigmas=3;
    [zgF,pzF]=discretizeAR1wSV_FarmerToda(rho,phi,sigmau,sigmae,xnum,znum,optsF);
    sX=(sigmae^2)/(1-phi^2); xB=2*log(sigmau)-sX/2;
    statdist=ones(xnum*znum,1)/(xnum*znum);
    for i_c=1:10000
        sn=pzF'*statdist;
        if max(abs(sn-statdist))<10^(-14), statdist=sn; break; end
        statdist=sn;
    end
    statdist=statdist/sum(statdist);
    Pxz=reshape(statdist,[xnum,znum]); % x fastest
    xg=zgF(1:xnum); xmarg=sum(Pxz,2);
    mx=sum(xmarg.*xg); vx=sum(xmarg.*(xg-mx).^2);
    fprintf('sigmae=%g: the x marginal has mean xBar [T2], this should be small: %2.3e \n',sigmae,abs(mx-xB))
    fprintf('sigmae=%g: the x marginal has variance sigmae^2/(1-phi^2) [T2], this should be small: %2.3e \n',sigmae,abs(vx-sX))
    fprintf('sigmae=%g: the x marginal sums to one [T1], this should be zero: %2.8e \n',sigmae,abs(sum(xmarg)-1))
end

%% 2. As sigmae -> 0 the volatility stops moving and z becomes an ordinary AR(1)
% NOTE ON WHAT IS COMPARED. An earlier version of this check compared the z GRID against the plain
% AR(1) grid and reported 0.000e+00 at every sigmae, which looked like clean convergence and was
% in fact vacuous: the command sets sigmaz=sqrt(exp(xBar+sigmaX/2)/(1-rho^2)) with
% xBar=2*log(sigmau)-sigmaX/2, so exp(xBar+sigmaX/2)=sigmau^2 identically and THE z GRID DOES NOT
% DEPEND ON sigmae AT ALL. The two grids were equal by construction, not by convergence. The
% reduction has to be tested on the transition probabilities, which is what is done here.
[zgP,pzP]=discretizeAR1_Tauchen(0,rho,sigmau,znum,3,struct());
prev=NaN;
for sigmae=[0.3,0.1,0.03,0.01]
    [zgT,pzT]=discretizeAR1wSV_Tauchen(rho,phi,sigmau,sigmae,xnum,znum,3,struct());
    zg=zgT(xnum+1:end);
    fprintf('sigmae=%g: the z grid equals the plain AR(1) grid [T1], this should be zero: %2.8e \n',sigmae,max(abs(zg-zgP)))
    fprintf('   (that is exact at every sigmae by construction, not convergence - see the note above) \n')
    % Collapse the joint transition onto z by averaging over the x marginal, and compare THAT
    statdist=ones(xnum*znum,1)/(xnum*znum);
    for i_c=1:10000
        sn=pzT'*statdist;
        if max(abs(sn-statdist))<10^(-14), statdist=sn; break; end
        statdist=sn;
    end
    statdist=statdist/sum(statdist);
    Pxz=reshape(statdist,[xnum,znum]);
    pz_marginal=zeros(znum,znum);
    for zi_c=1:znum
        w=Pxz(:,zi_c); w=w/sum(w); % distribution over x given z_i
        for xi_c=1:xnum
            ii=xi_c+(zi_c-1)*xnum;
            rowz=sum(reshape(pzT(ii,:),[xnum,znum]),1); % marginalize the destination over x
            pz_marginal(zi_c,:)=pz_marginal(zi_c,:)+w(xi_c)*rowz;
        end
    end
    d=max(abs(pz_marginal-pzP),[],'all');
    fprintf('sigmae=%g: the z-marginal transition approaches the plain AR(1) one, difference %2.3e \n',sigmae,d)
    if ~isnan(prev)
        fprintf('sigmae=%g: and the difference is smaller than at the previous sigmae [T2], this should be one: %i \n',sigmae,d<=prev)
    end
    prev=d;
end

%% 3. rho=0: z becomes iid, and its excess kurtosis has a one-line closed form
% With rho=0, z_t = u_t, a scale mixture of normals, whose excess kurtosis is 3*exp(sigmaX)-3.
% This is the cleanest possible check that the volatility machinery produces the right tails,
% because it removes the AR(1) averaging that otherwise thins them.
sigmae=calib.sigmae;
sX=(sigmae^2)/(1-phi^2);
ek_true=3*exp(sX)-3;
for cmd_c=1:2
    if cmd_c==1
        optsF=struct(); optsF.nSigmas=4;
        [zg0,pz0]=discretizeAR1wSV_FarmerToda(0,phi,sigmau,sigmae,9,31,optsF);
        nm='wSV_FarmerToda';
    else
        [zg0,pz0]=discretizeAR1wSV_Tauchen(0,phi,sigmau,sigmae,9,31,4,struct());
        nm='wSV_Tauchen';
    end
    statdist=ones(9*31,1)/(9*31);
    for i_c=1:10000
        sn=pz0'*statdist;
        if max(abs(sn-statdist))<10^(-14), statdist=sn; break; end
        statdist=sn;
    end
    statdist=statdist/sum(statdist);
    Pxz=reshape(statdist,[9,31]);
    zg=zg0(10:end); zmarg=sum(Pxz,1)';
    mz=sum(zmarg.*zg); vz=sum(zmarg.*(zg-mz).^2);
    ek=sum(zmarg.*(zg-mz).^4)/vz^2-3;
    fprintf('rho=0, %s: excess kurtosis %2.4f vs the closed form 3*exp(sigmaX)-3 = %2.4f, error [T2] %2.3e \n',nm,ek,ek_true,abs(ek-ek_true))
end

%% 4. The x block against a standalone discretization of the same AR(1)
% discretizeAR1wSV_FarmerToda discretizes its x process with discretizeVAR1_FarmerToda at M=1.
% P2 established that this does NOT agree with discretizeAR1_FarmerToda on the same problem - same
% grid, different transition matrix, because both apply the same absolute prior floor kappa=1e-8 to
% densities built on different scales. This is the first place that difference propagates into a
% user-facing result, so it is measured rather than left implicit.
sigmae=calib.sigmae;
sX=(sigmae^2)/(1-phi^2); xB=2*log(sigmau)-sX/2;
optsFT=struct(); optsFT.method='even'; optsFT.nSigmas=3; optsFT.verbose=0;
[xgA,pxA,~]=discretizeAR1_FarmerToda(xB*(1-phi),phi,sigmae,xnum,optsFT);
optsV=struct(); optsV.method='even'; optsV.nSigmas=3;
[xgB,pxB]=discretizeVAR1_FarmerToda(xB*(1-phi),phi,sigmae^2,xnum,optsV);
fprintf('the x process, AR1_FarmerToda vs VAR1_FarmerToda(M=1): x_grid [T1], this should be zero: %2.8e \n',max(abs(xgA(:)-xgB(:))))
fprintf('the x process, AR1_FarmerToda vs VAR1_FarmerToda(M=1): pi_x differs by %2.3e (reported; the \n',max(abs(pxA-pxB),[],'all'))
fprintf('   kappa=1e-8 prior floor is applied on different scales in the two, see P2) \n')

%% 5. The grid width must actually respond to nSigmas (regression test for B24)
% B24, found by this block on 2026-08-25: discretizeAR1wSV_FarmerToda set farmertodaoptions.nSigmas
% to 2 in order to discretize its x block - which is correct, Farmer-Toda use 2 sigmas for the
% volatility process - but wrote that into the shared options struct and never restored it, so the
% z grid was built at nSigmas=2 whatever the user asked for. The option was computed, validated
% against a warning threshold, and then discarded before its only use.
%
% The consequence was not a small loss of accuracy. The excess kurtosis of z came out at about
% -0.5 at every grid size against a truth of +0.83 - the WRONG SIGN - and got worse, not better,
% as znum rose, because the grid width never changed. So the first check here is the mechanical
% one: the z grid must span exactly +-nSigmas*sigmaz. That is an identity, not an approximation.
fprintf('\n--- the z grid width against nSigmas (regression test for B24) --- \n')
sigmaz=sqrt(calib.z.var); % = sigmau/sqrt(1-rho^2), since E[exp(x)]=exp(xBar+sigmaX/2)=sigmau^2
xnumB=9; znumB=31;
nSlist=[1.5,2,2.5,3,4];
ekB=zeros(1,length(nSlist));
for s_c=1:length(nSlist)
    optsB=struct(); optsB.nSigmas=nSlist(s_c);
    [zgB,pzB]=discretizeAR1wSV_FarmerToda(rho,phi,sigmau,calib.sigmae,xnumB,znumB,optsB);
    zblock=zgB(xnumB+1:end);
    fprintf('nSigmas=%3.1f: the z grid spans exactly +-nSigmas*sigmaz [T1], this should be zero: %2.8e \n',nSlist(s_c),max(abs([zblock(1);zblock(end)]-[-1;1]*nSlist(s_c)*sigmaz)))
    statdist=ones(xnumB*znumB,1)/(xnumB*znumB);
    for i_c=1:10000
        sn=pzB'*statdist;
        if max(abs(sn-statdist))<10^(-14), statdist=sn; break; end
        statdist=sn;
    end
    statdist=statdist/sum(statdist);
    Pxz=reshape(statdist,[xnumB,znumB]);
    zmarg=sum(Pxz,1)';
    mz=sum(zmarg.*zblock); vz=sum(zmarg.*(zblock-mz).^2);
    ekB(s_c)=sum(zmarg.*(zblock-mz).^4)/vz^2-3;
    fprintf('nSigmas=%3.1f: excess kurtosis %2.4f vs truth %2.4f, error %2.3e \n',nSlist(s_c),ekB(s_c),calib.z.exkurt,abs(ekB(s_c)-calib.z.exkurt))
end
% The substantive one. Truth is positive, and a discretization that cannot even get the sign of
% the moment the process exists to produce is not usable for it. Before B24 was fixed this was
% negative at EVERY entry of the list above, including nSigmas=4, because the option did nothing.
fprintf('at the widest grid the excess kurtosis is at least positive [T2], this should be one: %i \n',ekB(end)>0)
[~,bestat]=min(abs(ekB-calib.z.exkurt));
fprintf('the most accurate spacing in this sweep is nSigmas=%3.1f (reported, not asserted: wider is \n',nSlist(bestat))
fprintf('   not uniformly better at fixed znum, since widening the grid also coarsens the spacing) \n')

output.ek_nSigmas=ekB; output.nSlist=nSlist;
output.done=1;

end
