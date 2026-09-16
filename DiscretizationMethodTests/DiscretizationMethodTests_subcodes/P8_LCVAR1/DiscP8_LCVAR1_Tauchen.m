function output=DiscP8_LCVAR1_Tauchen(calib,znums,figure_c)
% P8: discretizeLifeCycleVAR1_Tauchen
%
% THE LAST UNCOVERED COMMAND. Before this subcode, discretizeLifeCycleVAR1_Tauchen was the only
% user-facing command in DiscretizationMethods that the bank never called - the string did not
% appear anywhere in the test tree. It is also the command that carried B7, B11 and B14, all three
% found by reading rather than by running, so nothing in it has ever been executed under a check.
%
% NO SIBLING METHOD. Every other block races two or three commands against each other. Here there is
% one, so the accuracy claim rests entirely on the analytic age profiles from DiscSetup_LCVAR1 and
% on the cross-test identities in DiscP8_crosstests, where discretizeVAR1_Tauchen and
% discretizeLifeCycleAR1_FellaGallipoliPanTauchen stand in as independent implementations of special
% cases this command has to reproduce.
%
% WHAT TAUCHEN CAN AND CANNOT DO HERE. The grid is +-nSigmas*sigmaz(j) per variable with nSigmas
% defaulting to min(sqrt(znum-1),3), and the transition probabilities come from the joint normal cdf
% over the bin rectangles. So the variance is biased DOWN by truncation at every age, and the bias
% is a property of the width rather than of znum - refining the grid at a fixed width drives the
% error to a floor rather than to zero, exactly as P1 and P2 found for the univariate Tauchen. The
% bars below are set accordingly: they are regression bars on a biased estimator, not accuracy claims.

fprintf('\n========== P8: discretizeLifeCycleVAR1_Tauchen ========== \n')

output=struct();
J=calib.J; M=calib.M;
Mew_J=calib.vary.Mew_J; Rho_J=calib.vary.Rho_J; SigmaSq_J=calib.vary.SigmaSq_J;
mewzT=calib.vary.mewz; SigmaSqzT=calib.vary.SigmaSqz; sigmazT=calib.vary.sigmaz;
acT=calib.vary.autocorr; ccT=calib.vary.crosscorr;
nz=length(znums);
err_mean=zeros(1,nz); err_var=zeros(1,nz); err_cov=zeros(1,nz); err_ac=zeros(1,nz);
runtime=zeros(1,nz);

for c_c=1:nz
    znum=znums(c_c);
    to=struct(); to.verbose=0;

    tic;
    [z_grid_J,pi_z_J,jequaloneDistz,oo]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum,J,to);
    twarm=toc;
    if twarm>calib.timethreshold
        nreps=1;
    else
        nreps=calib.nreps;
    end
    treps=zeros(1,nreps);
    for r_c=1:nreps
        tic;
        [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum,J,to);
        treps(r_c)=toc;
    end
    runtime(c_c)=median(treps);
    z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);
    Nz=znum^M;

    % --- invariants
    % z_grid_J is STACKED, sum(znum)-by-J, which is the other toolkit convention from the joint
    % (prod(znum))-by-M form discretizeVAR1_FarmerToda returns. B8 recorded that the header used to
    % document the joint shape; this asserts the shape that actually ships.
    fprintf('znum=%3i: size of z_grid_J is sum(znum)-by-J [T0], this should be zero: %i \n',znum,any(size(z_grid_J)~=[M*znum,J]))
    fprintf('znum=%3i: size of pi_z_J is prod(znum)-by-prod(znum)-by-(J-1) [T0], this should be zero: %i \n',znum,any(size(pi_z_J)~=[Nz,Nz,J-1]))
    fprintf('znum=%3i: size of jequaloneDistz is prod(znum)-by-1 [T0], this should be zero: %i \n',znum,any(size(jequaloneDistz)~=[Nz,1]))
    asc=1;
    for j_c=1:J
        for m_c=1:M
            asc=asc*issorted(z_grid_J((m_c-1)*znum+1:m_c*znum,j_c),'strictascend');
        end
    end
    fprintf('znum=%3i: every variable''s grid is strictly ascending at every age [T0], this should be one: %i \n',znum,asc)
    fprintf('znum=%3i: pi_z_J is in [0,1] [T0], this should be zero: %i \n',znum,any(pi_z_J(:)<0)+any(pi_z_J(:)>1))
    % B14: the third dimension is J-1, not J, so there is no padding slice to exempt - every slice
    % asserted, no exclusions.
    fprintf('znum=%3i: rows of pi_z_J sum to one at every one of the J-1 slices, this should be below 1e-07: %2.8e \n',znum,max(abs(sum(pi_z_J,2)-1),[],'all'))
    fprintf('znum=%3i: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',znum,abs(sum(jequaloneDistz)-1))
    fprintf('znum=%3i: jequaloneDistz is non-negative [T0], this should be zero: %i \n',znum,any(jequaloneDistz<0))
    fprintf('znum=%3i: no NaN or Inf [T0], this should be zero: %i \n',znum,any(~isfinite(z_grid_J(:)))+any(~isfinite(pi_z_J(:)))+any(~isfinite(jequaloneDistz)))

    % --- the command's own recursion against the one recomputed in the setup
    % otheroutputs is where the command reports the moments it USED to place the grid. If those are
    % wrong the grid is centred in the wrong place, which is exactly what B11 was: initialj0mewz was
    % read and then discarded, so period 1's mean was wrong whenever z0 was non-zero.
    fprintf('znum=%3i: otheroutputs.mewz matches the recursion [T1], this should be zero: %2.8e \n',znum,max(abs(gather(oo.mewz)-mewzT),[],'all'))
    fprintf('znum=%3i: otheroutputs.SigmaSqz matches the recursion [T1], this should be zero: %2.8e \n',znum,max(abs(gather(oo.SigmaSqz)-SigmaSqzT),[],'all'))
    fprintf('znum=%3i: otheroutputs.sigma_z is the sqrt of the diagonal of SigmaSqz [T1], this should be zero: %2.8e \n',znum,max(abs(gather(oo.sigma_z)-sigmazT),[],'all'))

    % --- the age profiles the chain actually delivers
    % Propagate jequaloneDistz forward. CreateGridvals with the third argument 1 unstacks the grid
    % into (znum^M)-by-M with VARIABLE 1 VARYING FASTEST, which is the toolkit's joint-index
    % convention; DiscP8_crosstests is where that ordering is tested rather than assumed.
    mz=zeros(M,J); vz=zeros(M,J); cz=zeros(1,J); acj=nan(M,J);
    d=jequaloneDistz;
    zprev=[]; dprev=[];
    for j_c=1:J
        zvals=CreateGridvals(znum*ones(M,1),z_grid_J(:,j_c),1);
        mz(:,j_c)=(d'*zvals)';
        dev=zvals-mz(:,j_c)';
        Vz=(dev.*d)'*dev;
        vz(:,j_c)=diag(Vz);
        cz(j_c)=Vz(1,2)/sqrt(Vz(1,1)*Vz(2,2));
        if j_c>1
            % lag-one cross moment: E[z_j z_{j-1}'] under the joint of consecutive ages
            Ezz=zeros(M,M);
            Pj=pi_z_J(:,:,j_c-1);
            for m1=1:M
                for m2=1:M
                    Ezz(m1,m2)=sum((dprev.*zprev(:,m2)).*(Pj*zvals(:,m1)));
                end
            end
            Clag=Ezz-mz(:,j_c)*mz(:,j_c-1)';
            acj(:,j_c)=diag(Clag)./(sqrt(vz(:,j_c)).*sqrt(vz(:,j_c-1)));
        end
        if j_c<J
            zprev=zvals; dprev=d;
            d=pi_z_J(:,:,j_c)'*d; % slice j_c is the transition from age j_c to age j_c+1
        end
    end
    [err_mean(c_c),imw]=max(abs(mz(:)-mewzT(:)));
    [~,jmw]=ind2sub([M,J],imw);
    err_var(c_c)=max(abs(vz-sigmazT.^2),[],'all');
    err_cov(c_c)=max(abs(cz-ccT));
    err_ac(c_c)=max(abs(acj(:,2:end)-acT(:,2:end)),[],'all');
    fprintf('znum=%3i: worst-age mean error [T2] %2.3e at age %i, variance error [T2] %2.3e \n',znum,err_mean(c_c),jmw,err_var(c_c))
    % THESE THREE LINES FOUND B33, and they stay because they are what localises this class of bug.
    % The 2026-09-16 run had the worst-age mean error flat at ~3.1e-02 across the whole znum sweep
    % and converging to the same 3.10e-02 across the width sweep - a constant approached from both
    % directions, which truncation cannot produce. Age 1 was exactly zero and the error grew with
    % age, so it was accumulating through the transitions rather than sitting in jequaloneDistz.
    % That is an off-by-one signature, and it was: discretizeLifeCycleVAR1_Tauchen indexed the
    % transition's Mew_J, Rho_J and SigmaSq_J at jj rather than jj+1, so the chain's mean followed
    % m(j+1)=Mew(:,j)+Rho(:,:,j)*m(j). Iterating that wrong recursion by hand reproduces 3.10e-02 at
    % age 15, 3.54e-03 at age 2 and 2.18e-02 at age 25, against the 3.071e-02, 3.535e-03 and
    % 2.147e-02 measured - the residual being the discretization's own error on top of the bias.
    %
    % The bug is invisible without age-varying parameters, which is why nothing else caught it: the
    % frozen cross-test holds Mew and Rho constant, so indexing them at jj or jj+1 is the same thing.
    fprintf('znum=%3i: mean error at age 1 alone [T2] %2.3e, at age 2 %2.3e, at the last age %2.3e \n',znum,max(abs(mz(:,1)-mewzT(:,1))),max(abs(mz(:,2)-mewzT(:,2))),max(abs(mz(:,J)-mewzT(:,J))))
    % jequaloneDistz is built by putting N(mewz(:,1),SigmaSqz(:,:,1)) on the period-1 grid, so its
    % own mean is a separate thing from the grid's centre and can be checked directly.
    zv1=CreateGridvals(znum*ones(M,1),z_grid_J(:,1),1);
    fprintf('znum=%3i: jequaloneDistz has mean mewz(:,1) [T2], this should be below %g: %2.3e \n',znum,0.01,max(abs((jequaloneDistz'*zv1)'-mewzT(:,1))))
    fprintf('znum=%3i: worst-age cross-correlation error [T2] %2.3e, autocorrelation error [T2] %2.3e \n',znum,err_cov(c_c),err_ac(c_c))
    fprintf('znum=%3i: runtime %2.6f s (nreps=%i) \n',znum,runtime(c_c),nreps)

    output.sweep(c_c).mz=mz; output.sweep(c_c).vz=vz; output.sweep(c_c).cz=cz; output.sweep(c_c).acj=acj;
end
output.err_mean=err_mean; output.err_var=err_var; output.err_cov=err_cov; output.err_ac=err_ac;
output.runtime=runtime; output.znums=znums;

%% Convergence, and the floor
% Refining the grid at a fixed width cannot remove the truncation bias, so the thing to assert is
% that the error FALLS and then flattens, not that it goes to zero.
%
% THE MEAN IS NOT ASSERTED TO FALL - it is asserted to be EXACT. With B33 fixed the mean error came
% out at 8e-17 to 1.5e-16 across the whole sweep, which is machine precision, and that is not an
% accident: Tauchen's grid is symmetric about mewz(:,j) and the transition now uses the right age's
% parameters, so the discretized mean is exact by symmetry however coarse the grid is. Asserting it
% FALLS would compare 1.5e-16 against 8.3e-17 - two numbers that are both zero to within rounding,
% where which is smaller is noise. That assertion duly failed on 2026-09-16 with nothing wrong, the
% same vacuous-ordering mistake this bank had just fixed in P9. So the claim is the one with content.
fprintf('\nworst-age mean error is at machine precision at every grid size [T2], this should be below %g: %2.3e \n',1e-12,max(err_mean))
fprintf('   (before B33 was fixed this was 2.7e-02 to 3.1e-02 and flat in znum; a mean error that does \n')
fprintf('   not shrink with the grid is the signature that found it.) \n')
fprintf('worst-age variance error falls over the sweep [T2], this should be one: %i \n',err_var(end)<err_var(1))
fprintf('   mean          :');
for c_c=1:nz
    fprintf(' %2.1e',err_mean(c_c));
end
fprintf(' \n   variance      :');
for c_c=1:nz
    fprintf(' %2.1e',err_var(c_c));
end
fprintf(' \n   cross-corr    :');
for c_c=1:nz
    fprintf(' %2.1e',err_cov(c_c));
end
fprintf(' \n   autocorr      :');
for c_c=1:nz
    fprintf(' %2.1e',err_ac(c_c));
end
fprintf(' \n')

%% How znum is passed
% THE SHAPE OF znum IS NOT COSMETIC HERE, and this block exists because reading the command
% suggested it matters. The default width is set by
%     tauchenoptions.nSigmas = min(sqrt(znum-1)',3);
% at a point BEFORE the scalar-to-vector normalisation further down, and it is later used as
%     q_sigmaz = tauchenoptions.nSigmas.*sigmaz;
% with sigmaz of size M-by-J. A scalar znum leaves nSigmas scalar, which broadcasts against anything.
% A ROW vector znum makes nSigmas M-by-1, which broadcasts correctly against M-by-J. A COLUMN vector
% znum - which is the shape the command's own header documents, "(M x 1)" - makes nSigmas 1-by-M,
% and 1-by-M against M-by-J has no valid expansion unless J happens to equal M.
%
% So the prediction under test is: scalar works, row works, column errors. If all three agree
% instead, the reading was wrong and this block says so. Either way it is written as a comparison
% against the scalar call rather than as an assumption about which one is right.
znum_s=znums(min(2,nz));
to=struct(); to.verbose=0;
[zg_s,pz_s]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,J,to);
zg_s=gather(zg_s); pz_s=gather(pz_s);
fprintf('\nznum passed as a scalar: it runs [T0], this should be one: %i \n',all(isfinite(zg_s(:))))
try
    [zg_r,pz_r]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s*ones(1,M),J,to);
    zg_r=gather(zg_r); pz_r=gather(pz_r);
    fprintf('znum passed as a ROW vector: it runs, and equals the scalar call: z_grid_J [T0], this should be zero: %2.8e \n',max(abs(zg_r-zg_s),[],'all'))
    fprintf('znum passed as a ROW vector: pi_z_J equals the scalar call [T0], this should be zero: %2.8e \n',max(abs(pz_r-pz_s),[],'all'))
catch ME
    fprintf('znum passed as a ROW vector: ERRORED, which was not expected. The message is: %s \n',ME.message)
end
try
    [zg_c,pz_c]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s*ones(M,1),J,to);
    zg_c=gather(zg_c); pz_c=gather(pz_c);
    fprintf('znum passed as a COLUMN vector, the shape the header documents: it runs, and equals the scalar call: z_grid_J [T0], this should be zero: %2.8e \n',max(abs(zg_c-zg_s),[],'all'))
    fprintf('znum passed as a COLUMN vector: pi_z_J equals the scalar call [T0], this should be zero: %2.8e \n',max(abs(pz_c-pz_s),[],'all'))
catch ME
    fprintf('znum passed as a COLUMN vector, the shape the header documents: ERRORED. The message is: %s \n',ME.message)
    fprintf('   If that message is about incompatible array sizes, it is the nSigmas transpose described above, \n')
    fprintf('   and the command cannot be called the way its own header says to call it. \n')
end

%% nSigmas
% The default is min(sqrt(znum-1),3), a cap of THREE - the same value discretizeLifeCycleAR1_-
% FellaGallipoliPanTauchen uses and one less than the four every stationary Tauchen command caps at.
% DiscP8_gridwidth measures whether three is the right place to stop; here the only claim is that the
% option is read, which is the B5/B6/B10/B11 failure mode.
for nS=[2,3,4]
    ton=struct(); ton.verbose=0; ton.nSigmas=nS*ones(M,1);
    [zgn,~]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,J,ton);
    zgn=gather(zgn);
    w1=(zgn(znum_s,1)-zgn(1,1))/2; % half-width of variable 1 at age 1
    fprintf('nSigmas=%i: variable 1''s half-width at age 1 is nSigmas*sigmaz [T1], this should be zero: %2.8e \n',nS,abs(w1-nS*sigmazT(1,1)))
end

%% initialj0mewz, which is B11
% B11 was that this option was read into z0 and then never used, so period 1's mean was Mew(:,1)
% rather than Mew(:,1)+Rho(:,:,1)*z0. It is silent at the default z0=0, which is why it survived, so
% the check has to use a NON-ZERO z0 or it tests nothing.
z0=[0.5;-0.3];
t0=struct(); t0.verbose=0; t0.initialj0mewz=z0;
[zg0,~,~,oo0]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,J,t0);
mewz0=zeros(M,J); mewz0(:,1)=Mew_J(:,1)+Rho_J(:,:,1)*z0;
for j_c=2:J
    mewz0(:,j_c)=Mew_J(:,j_c)+Rho_J(:,:,j_c)*mewz0(:,j_c-1);
end
fprintf('\ninitialj0mewz: the mean profile is Mew(:,1)+Rho(:,:,1)*z0 recursed forward [T1], this should be zero: %2.8e \n',max(abs(gather(oo0.mewz)-mewz0),[],'all'))
fprintf('initialj0mewz: and it DIFFERS from the z0=0 profile, so the option is not being ignored [T0], this should be one: %i \n',max(abs(gather(oo0.mewz)-mewzT),[],'all')>1e-8)
fprintf('initialj0mewz: the grid centre at age 1 moved with it [T1], this should be zero: %2.8e \n',abs((gather(zg0(1,1))+gather(zg0(znum_s,1)))/2-mewz0(1,1)))
fprintf('initialj0mewz: the covariance profile is UNCHANGED, since z0 is a point mass [T1], this should be zero: %2.8e \n',max(abs(gather(oo0.SigmaSqz)-SigmaSqzT),[],'all'))

%% initialj1mewz and initialj1SigmaSqz
% These overwrite period 1 outright rather than deriving it, and they can be set independently. The
% two one-sided forms are the interesting ones: setting only the mean zeroes the covariance, which is
% the degenerate period 1 that B7 had to handle, and setting only the covariance zeroes the mean.
m1=[0.2;-0.1]; S1=[0.09,0.02;0.02,0.05];
t1=struct(); t1.verbose=0; t1.initialj1mewz=m1; t1.initialj1SigmaSqz=S1;
[~,~,~,oo1]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,J,t1);
fprintf('\ninitialj1mewz and initialj1SigmaSqz together: period 1 mean is taken as given [T1], this should be zero: %2.8e \n',max(abs(gather(oo1.mewz(:,1))-m1)))
fprintf('initialj1mewz and initialj1SigmaSqz together: period 1 covariance is taken as given [T1], this should be zero: %2.8e \n',max(abs(gather(oo1.SigmaSqz(:,:,1))-S1),[],'all'))
mewz1=zeros(M,J); mewz1(:,1)=m1; SigmaSqz1=zeros(M,M,J); SigmaSqz1(:,:,1)=S1;
for j_c=2:J
    mewz1(:,j_c)=Mew_J(:,j_c)+Rho_J(:,:,j_c)*mewz1(:,j_c-1);
    SigmaSqz1(:,:,j_c)=Rho_J(:,:,j_c)*SigmaSqz1(:,:,j_c-1)*Rho_J(:,:,j_c)'+SigmaSq_J(:,:,j_c);
end
fprintf('initialj1mewz and initialj1SigmaSqz together: ages 2 onward recurse from it [T1], this should be zero: %2.8e \n',max(abs(gather(oo1.mewz)-mewz1),[],'all'))
fprintf('initialj1mewz and initialj1SigmaSqz together: and the covariance likewise [T1], this should be zero: %2.8e \n',max(abs(gather(oo1.SigmaSqz)-SigmaSqz1),[],'all'))

t2=struct(); t2.verbose=0; t2.initialj1mewz=m1;
[~,~,jd2,oo2]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,J,t2);
fprintf('initialj1mewz alone: period 1 covariance is set to zero, a point mass [T1], this should be zero: %2.8e \n',max(abs(gather(oo2.SigmaSqz(:,:,1))),[],'all'))
% B7: a singular period-1 covariance is where mvncdf cannot be used, and the command falls back to a
% point mass on the median grid point. A degenerate jequaloneDistz must still be a distribution.
jd2=gather(jd2);
fprintf('initialj1mewz alone: jequaloneDistz still sums to one on the degenerate period 1 [T1], this should be zero: %2.8e \n',abs(sum(jd2)-1))
fprintf('initialj1mewz alone: and it is a point mass, so exactly one state carries all the mass [T0], this should be one: %i \n',sum(jd2>1e-12)==1)

t3=struct(); t3.verbose=0; t3.initialj1SigmaSqz=S1;
[~,~,~,oo3]=discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,J,t3);
fprintf('initialj1SigmaSqz alone: period 1 mean is set to zero [T1], this should be zero: %2.8e \n',max(abs(gather(oo3.mewz(:,1)))))
fprintf('initialj1SigmaSqz alone: period 1 covariance is taken as given [T1], this should be zero: %2.8e \n',max(abs(gather(oo3.SigmaSqz(:,:,1))-S1),[],'all'))

%% Error paths
te=struct(); te.verbose=0;
try
    discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,2,J,te);
    fprintf('\nznum=2 was accepted, but the command documents a minimum of 3 [T0], this should be one: 0 \n')
catch ME
    fprintf('\nznum below 3 errors as documented, this should be one: %i \n',contains(ME.message,'znum'))
end
try
    discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSq_J,znum_s,1,te);
    fprintf('J=1 was accepted, but the command documents a minimum of 2 [T0], this should be one: 0 \n')
catch ME
    fprintf('J below 2 errors as documented, this should be one: %i \n',contains(ME.message,'horizon')||contains(ME.message,'J'))
end
try
    discretizeLifeCycleVAR1_Tauchen(Mew_J(:,1:end-1),Rho_J,SigmaSq_J,znum_s,J,te);
    fprintf('a Mew_J with the wrong number of columns was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('a Mew_J with the wrong number of columns errors, and the message names Mew, this should be one: %i \n',contains(ME.message,'Mew'))
end
try
    SigmaSqbad=SigmaSq_J; SigmaSqbad(1,1,3)=-1;
    discretizeLifeCycleVAR1_Tauchen(Mew_J,Rho_J,SigmaSqbad,znum_s,J,te);
    fprintf('a negative variance was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('a negative variance errors, and the message names the age, this should be one: %i \n',contains(ME.message,'SigmaSq')||contains(ME.message,'variance'))
end

%% REGRESSION BARS ON THE ACCURACY NUMBERS ABOVE
% Regression bars, not accuracy claims, for the reason set out at the top: Tauchen's variance is
% biased down by truncation at nSigmas=3 and refining znum drives the error to that floor rather
% than to zero. Each family gets the pair DiscP7 established - the worst over the sweep against a
% loose bar, which catches a gross regression anywhere, and the finest grid against a tighter one.
fprintf('\nworst mean error over the whole sweep [T2], this should be below %g: %2.3e \n',0.5,max(err_mean))
fprintf('mean error at the finest grid [T2], this should be below %g: %2.3e \n',0.05,err_mean(end))
fprintf('worst variance error over the whole sweep [T2], this should be below %g: %2.3e \n',0.5,max(err_var))
fprintf('variance error at the finest grid [T2], this should be below %g: %2.3e \n',0.1,err_var(end))
fprintf('worst cross-correlation error over the whole sweep [T2], this should be below %g: %2.3e \n',0.5,max(err_cov))
fprintf('worst autocorrelation error over the whole sweep [T2], this should be below %g: %2.3e \n',0.5,max(err_ac))

%% Figure
figure(figure_c)
subplot(1,3,1)
semilogy(znums,max(err_mean,1e-18),'o-',znums,max(err_var,1e-18),'s-',znums,max(err_cov,1e-18),'d-',znums,max(err_ac,1e-18),'^-')
set(gca,'XScale','log')
xlabel('znum (per variable)'); ylabel('worst-age |error|')
title('LifeCycleVAR1\_Tauchen: accuracy'); legend('mean','variance','cross-corr','autocorr','Location','best')
subplot(1,3,2)
plot(1:J,sigmazT(1,:),'k-','LineWidth',1.5); hold on
for c_c=1:nz
    plot(1:J,sqrt(output.sweep(c_c).vz(1,:)),'--')
end
hold off
xlabel('age j'); ylabel('sd of variable 1'); title('sd profile, variable 1 (black is truth)')
subplot(1,3,3)
loglog(znums,runtime,'o-')
xlabel('znum (per variable)'); ylabel('seconds'); title('runtime')
sgtitle('P8: discretizeLifeCycleVAR1\_Tauchen')

end
