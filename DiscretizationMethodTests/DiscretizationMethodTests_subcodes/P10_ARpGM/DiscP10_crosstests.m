function DiscP10_crosstests(calib,znum)
% P10: cross-tests for the two AR(p) commands
%
% BOTH COMMANDS ARE NEW, so there is no prior version to regress against. What there is instead is a
% set of special cases where an AR(p) command must reproduce something already in the toolkit and
% separately tested. Those are the oracles, and they are what this block is for.
%
% The reduction identities all turn on the same thing: an AR(p) with p=1 IS an AR(1), and an AR(p)
% whose last coefficient is zero IS an AR(p-1). A command that got the lag bookkeeping wrong would
% still return a valid stochastic matrix with plausible moments; only a reduction to something
% independently implemented catches it.

fprintf('\n========== P10: cross-tests ========== \n')
mew=calib.mew; sigma=calib.sigma;
mixprobs_i=calib.mixprobs_i; mu_i=calib.mu_i; sigma_i=calib.sigma_i;

%% 1. p=1 against discretizeAR1_FarmerToda
% THE METHOD AND WIDTH ARE PINNED EXPLICITLY rather than left to the defaults, and that is not
% tidiness. The two commands branch their defaults on different persistence measures - rho for the
% p=1 command, the largest modulus companion eigenvalue for the AR(p) one - and those agree only for
% rho>=0 (B34). Pinning both options makes this an identity about the ALGORITHM rather than about
% whether two default rules happen to coincide, which is what it is supposed to test.
rho1=0.7;
fprintf('\n--- 1. p=1 against discretizeAR1_FarmerToda, rho=%g --- \n',rho1)
for mm=[2,3,4]
    o1=struct(); o1.verbose=0; o1.method='even'; o1.nSigmas=sqrt(znum-1); o1.nMoments=mm; o1.parallel=1;
    [zgA,pzA]=discretizeAR1_FarmerToda(mew,rho1,sigma,znum,o1);
    [zgB,pzB]=discretizeARp_FarmerToda(mew,rho1,sigma,znum,o1);
    zgA=gather(zgA); pzA=gather(pzA); zgB=gather(zgB); pzB=gather(pzB);
    fprintf('nMoments=%i: z_grid agrees [T0], this should be zero: %2.8e \n',mm,max(abs(zgA-zgB)))
    fprintf('nMoments=%i: pi_z agrees [T0], this should be zero: %2.8e \n',mm,max(abs(pzA-pzB),[],'all'))
end
% ...and a NEGATIVE rho, where the two commands' defaults are known to differ. With both options
% pinned the identity must still hold exactly; if it does, B34 is purely about the defaults and not
% about the algorithm.
% rho=-0.9, NOT -0.7. The divergence B34 describes needs |rho|>0.8 with rho<0: the p=1 command
% reads rho<=0.8 and takes the low-persistence branch, while this one reads maxabseig=|rho|>0.8 and
% takes the high-persistence one. At rho=-0.7 the modulus is ALSO below 0.8, both commands pick the
% same branch, and the probe measures nothing - which is what the run of 2026-09-17 reported, a
% difference of 0.000e+00 underneath a line of text claiming the two diverge.
rhoneg=-0.9;
o1=struct(); o1.verbose=0; o1.method='even'; o1.nSigmas=sqrt(znum-1); o1.nMoments=2; o1.parallel=1;
[zgA,pzA]=discretizeAR1_FarmerToda(mew,rhoneg,sigma,znum,o1);
[zgB,pzB]=discretizeARp_FarmerToda(mew,rhoneg,sigma,znum,o1);
fprintf('rho=%g with both options pinned: z_grid agrees [T0], this should be zero: %2.8e \n',rhoneg,max(abs(gather(zgA)-gather(zgB))))
fprintf('rho=%g with both options pinned: pi_z agrees [T0], this should be zero: %2.8e \n',rhoneg,max(abs(gather(pzA)-gather(pzB)),[],'all'))
% The defaults themselves: a record of the divergence B34 describes, printed rather than asserted,
% because which one is right is the open question.
od=struct(); od.verbose=0; od.parallel=1;
[~,pzC]=discretizeAR1_FarmerToda(mew,rhoneg,sigma,znum,od);
[~,pzD]=discretizeARp_FarmerToda(mew,rhoneg,sigma,znum,od);
defaultgap=max(abs(gather(pzC)-gather(pzD)),[],'all');
fprintf('at rho=%g on DEFAULTS the two differ by %2.3e \n',rhoneg,defaultgap)
fprintf('   This is B34. The p=1 command branches its method default on SIGNED rho and this one on \n')
fprintf('   |rho|, so at rho=-0.9 the first picks gauss-hermite and the second even. Which is right \n')
fprintf('   is the open question; that they differ is not, so it is asserted rather than printed. \n')
fprintf('the defaults DO diverge at negative persistent rho, which is what B34 is about [T0], this should be one: %i \n',defaultgap>1e-10)

%% 2. p=1 against discretizeAR1wGM_FarmerToda
fprintf('\n--- 2. p=1 against discretizeAR1wGM_FarmerToda --- \n')
og=struct(); og.verbose=0; og.method='even'; og.nSigmas=sqrt(znum-1); og.nMoments=4; og.parallel=1;
[zgE,pzE]=discretizeAR1wGM_FarmerToda(mew,rho1,mixprobs_i,mu_i,sigma_i,znum,og);
[zgF,pzF]=discretizeARpwGM_FarmerToda(mew,rho1,mixprobs_i,mu_i,sigma_i,znum,og);
fprintf('z_grid agrees [T0], this should be zero: %2.8e \n',max(abs(gather(zgE)-gather(zgF))))
fprintf('pi_z agrees [T0], this should be zero: %2.8e \n',max(abs(gather(pzE)-gather(pzF)),[],'all'))

%% 3. a one-component mixture is a gaussian
% The GM command with nmix=1 must reproduce the gaussian command exactly. This is the check that the
% mixture machinery does not quietly change the answer in the case where it should not.
fprintf('\n--- 3. nmix=1 against the gaussian AR(p) command --- \n')
Rho2=calib.AR2.Rho;
o3=struct(); o3.verbose=0; o3.method='even'; o3.nSigmas=sqrt(znum-1); o3.nMoments=2; o3.parallel=1;
[zgG,pzG]=discretizeARp_FarmerToda(mew,Rho2,sigma,znum,o3);
[zgH,pzH]=discretizeARpwGM_FarmerToda(mew,Rho2,1,0,sigma,znum,o3);
fprintf('p=2, nmix=1 with mu=0 and the same sigma: z_grid [T0], this should be zero: %2.8e \n',max(abs(gather(zgG)-gather(zgH))))
fprintf('p=2, nmix=1 with mu=0 and the same sigma: pi_z [T1], this should be zero: %2.8e \n',max(abs(gather(pzG)-gather(pzH)),[],'all'))
fprintf('   (T1 rather than T0: the two reach the same conditional moments by different arithmetic - \n')
fprintf('   the gaussian command writes TBar down in closed form, the mixture one sums over components.) \n')

%% 4. the last coefficient zero collapses to one fewer lag
% An AR(3) with Rho(3)=0 is an AR(2) in the first two dimensions, with the third dimension along for
% the ride. The chains are not the same SIZE, so what is compared is the implied dynamics: the
% marginal over the first two dimensions of the AR(3) chain must match the AR(2) chain.
fprintf('\n--- 4. Rho(p)=0 collapses to AR(p-1) --- \n')
o4=struct(); o4.verbose=0; o4.method='even'; o4.nSigmas=sqrt(znum-1); o4.nMoments=2; o4.parallel=1;
[zg2,pz2]=discretizeARp_FarmerToda(mew,[Rho2(1),Rho2(2)],sigma,znum,o4);
[zg3,pz3]=discretizeARp_FarmerToda(mew,[Rho2(1),Rho2(2),0],sigma,znum,o4);
zg2=gather(zg2); pz2=gather(pz2); zg3=gather(zg3); pz3=gather(pz3);
fprintf('the shared one-dimensional grid is the same [T1], this should be zero: %2.8e \n',max(abs(zg3(1:znum)-zg2(1:znum))))
% The AR(3) transition from (a1,a2,a3) must not depend on a3 when Rho(3)=0, and the probabilities
% over the new z must equal the AR(2) ones from (a1,a2).
gap=0; dep=0;
for i1=1:znum
    for i2=1:znum
        ii2=i1+znum*(i2-1); % AR(2) joint index, dimension 1 fastest
        p2row=pz2(ii2,(i1-1)*znum+(1:znum)); % destinations from (i1,i2) are znum*(i1-1)+(1:znum)
        for i3=1:znum
            ii3=i1+znum*(i2-1)+znum^2*(i3-1);
            base3=znum*mod(ii3-1,znum^2);
            p3row=pz3(ii3,base3+(1:znum));
            gap=max(gap,max(abs(p3row-p2row)));
            if i3>1
                ii3a=i1+znum*(i2-1);
                base3a=znum*mod(ii3a-1,znum^2);
                dep=max(dep,max(abs(p3row-pz3(ii3a,base3a+(1:znum)))));
            end
        end
    end
end
fprintf('the AR(3) conditional over the new z equals the AR(2) one [T1], this should be zero: %2.8e \n',gap)
fprintf('and it does not depend on the third lag at all [T1], this should be zero: %2.8e \n',dep)

%% 5. the joint index ordering
% P5 found the two stationary VAR commands order their joint index OPPOSITELY, so this is measured
% rather than assumed. With Rho=[rho,0] the two dimensions are independent draws one period apart,
% and the transition factors in a way that says which dimension moves fastest.
fprintf('\n--- 5. the joint index ordering --- \n')
o5=struct(); o5.verbose=0; o5.method='even'; o5.nSigmas=sqrt(znum-1); o5.nMoments=2; o5.parallel=1;
[zgI,pzI]=discretizeARp_FarmerToda(mew,[0.5,0],sigma,znum,o5);
zgI=gather(zgI); pzI=gather(pzI);
zvI=CreateGridvals(znum*ones(2,1),zgI,1);
% Under dimension-1-fastest, state ii has dimension 1 equal to zgI(mod(ii-1,znum)+1). Check that
% CreateGridvals and the command agree by reading the conditional mean back off the transition: from
% each state the mean of the new z must be mew+0.5*(dimension 1 of that state).
worst=0;
for ii=1:znum^2
    base=znum*mod(ii-1,znum);
    pm=pzI(ii,base+(1:znum))*zgI(1:znum);
    worst=max(worst,abs(pm-(mew+0.5*zvI(ii,1))));
end
fprintf('the conditional mean read off pi_z equals mew+rho*(dimension 1) at every state [T1], this should be zero: %2.8e \n',worst)
fprintf('   which is CreateGridvals'' ordering, dimension 1 varying FASTEST - the same convention \n')
fprintf('   discretizeVAR1_Tauchen uses and discretizeVAR1_FarmerToda does not (see P5). \n')

%% 6. error paths
fprintf('\n--- 6. error paths --- \n')
oe=struct(); oe.verbose=0; oe.parallel=1;
try
    discretizeARp_FarmerToda(mew,[0.9,0.5],sigma,znum,oe); % sum>1, explosive
    fprintf('a non-stationary AR(p) was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('a non-stationary AR(p) errors, and the message says stationary, this should be one: %i \n',contains(ME.message,'stationar'))
end
try
    discretizeARp_FarmerToda(mew,Rho2,sigma,[znum,znum+2],oe);
    fprintf('unequal znum per dimension was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('unequal znum per dimension errors, and the message says so, this should be one: %i \n',contains(ME.message,'same number'))
end
try
    discretizeARp_FarmerToda(mew,[0.2,0.1,0.1,0.1,0.1,0.1],sigma,4,oe);
    fprintf('p=6 was accepted, but CreateGridvals handles at most five [T0], this should be one: 0 \n')
catch ME
    fprintf('p=6 errors, and the message names p rather than n_x, this should be one: %i \n',contains(ME.message,'p>5')||contains(ME.message,'p=6'))
end
try
    discretizeARpwGM_FarmerToda(mew,Rho2,[0.6;0.5],mu_i,sigma_i,znum,oe);
    fprintf('mixture probabilities not summing to one were accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('mixture probabilities not summing to one error, this should be one: %i \n',contains(ME.message,'sum to one'))
end

end
