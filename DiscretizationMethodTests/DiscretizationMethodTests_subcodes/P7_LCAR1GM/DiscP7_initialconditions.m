function output=DiscP7_initialconditions(calib)
% P7: the initial-condition options, including the gaussian-mixture period 1
%
% This command has the largest initial-condition surface in the toolkit: the four scalar options
% that P6 covered, plus a period 1 that can itself be a gaussian mixture, signalled by giving
% initialj1sigmaz as a VECTOR (that was B6). discretizeLifeCycleAR1wGM_Tauchen was written to the
% same convention, so both are exercised here.
%
% jequaloneDistz is the object these options determine, and three separate things have to hold: it
% sums to one, it lives on z_grid_J(:,1), and its moments are those of the period 1 distribution.
% The middle one is not automatic - B16 was the gaussian KFTT building jequaloneDistz on a
% different grid from z_grid_J(:,1), which puts mass on the wrong points and which nothing
% downstream would notice.

fprintf('\n========== P7: initial conditions and jequaloneDistz ========== \n')

output=struct();
J=calib.J; znum=15;
mew=calib.vary.mew; rho=calib.vary.rho;
p=calib.vary.mixprobs_i; mu=calib.vary.mu_i; sd=calib.vary.sigma_i;
r1=rho(1); mw1=mew(1);
e1mean=calib.vary.emean(1); e1var=calib.vary.evar(1);

m0=0.3; s0=0.25; m1=0.4; s1=0.35;
% For each configuration: the label, the option fields, and what period 1's mean and variance must
% be. Note the KFTT and Tauchen commands differ on whether the mixture mean enters the CENTRING,
% but not on what the distribution of z(1) actually is, so these targets serve both.
cfg={};
cfg{1}={'default',           struct(),                                         mw1+e1mean,        e1var};
cfg{2}={'initialj0mewz',     struct('initialj0mewz',m0),                        mw1+r1*m0+e1mean,  e1var};
cfg{3}={'initialj0sigmaz',   struct('initialj0sigmaz',s0),                      mw1+e1mean,        r1^2*s0^2+e1var};
cfg{4}={'both initialj0',    struct('initialj0mewz',m0,'initialj0sigmaz',s0),   mw1+r1*m0+e1mean,  r1^2*s0^2+e1var};
cfg{5}={'initialj1mewz',     struct('initialj1mewz',m1),                        m1,                0};
cfg{6}={'initialj1sigmaz',   struct('initialj1sigmaz',s1),                      0,                 s1^2};
cfg{7}={'both initialj1',    struct('initialj1mewz',m1,'initialj1sigmaz',s1),   m1,                s1^2};

cmdname={'LifeCycleAR1wGM_KFTT','LifeCycleAR1wGM_Tauchen'};
for k_c=1:2
    fprintf('\n--- %s --- \n',cmdname{k_c})
    for c_c=1:7
        opts=cfg{c_c}{2}; lbl=cfg{c_c}{1}; mtrue=cfg{c_c}{3}; vtrue=cfg{c_c}{4};
        opts.verbose=0;
        if k_c==1
            [z_grid_J,~,jequaloneDistz]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,opts);
        else
            [z_grid_J,~,jequaloneDistz]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,3,opts);
        end
        z_grid_J=gather(z_grid_J); jequaloneDistz=gather(jequaloneDistz);
        g1=z_grid_J(:,1);
        m=sum(jequaloneDistz(:).*g1);
        v=sum(jequaloneDistz(:).*(g1-m).^2);
        fprintf('%-18s: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',lbl,abs(sum(jequaloneDistz)-1))
        fprintf('%-18s: its mean is %+2.4f against %+2.4f [T2], error %2.3e \n',lbl,m,mtrue,abs(m-mtrue))
        fprintf('%-18s: its variance is %2.5f against %2.5f [T2], error %2.3e \n',lbl,v,vtrue,abs(v-vtrue))
        if vtrue==0
            fprintf('%-18s: variance zero, so it is a point mass [T0], this should be one: %i \n',lbl,sum(jequaloneDistz>10^(-12))==1)
        end
        output.(sprintf('cmd%i',k_c)).cfg(c_c).label=lbl;
        output.(sprintf('cmd%i',k_c)).cfg(c_c).moments=[m,v];
    end
end

%% Period 1 as a gaussian mixture (B6), which is unique to this pair of commands
% A VECTOR initialj1sigmaz means period 1 is itself a mixture, with initialj1mixprobs and
% initialj1mu giving the weights and means. The distribution of z(1) is then that mixture directly.
fprintf('\n--- period 1 as a gaussian mixture (the B6 convention) --- \n')
j1p=[0.6;0.4]; j1m=[0.2;-0.3]; j1s=[0.15;0.35];
j1mean=sum(j1p.*j1m);
j1var=sum(j1p.*(j1m.^2+j1s.^2))-j1mean^2;
fprintf('the requested period 1 mixture has mean %+2.4f and variance %2.5f \n',j1mean,j1var)
for k_c=1:2
    o=struct(); o.verbose=0;
    o.initialj1sigmaz=j1s; o.initialj1mixprobs=j1p; o.initialj1mu=j1m;
    if k_c==1
        [zg,~,j1d]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,o);
    else
        [zg,~,j1d]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,3,o);
    end
    zg=gather(zg); j1d=gather(j1d);
    g1=zg(:,1); m=sum(j1d(:).*g1); v=sum(j1d(:).*(g1-m).^2);
    fprintf('%-24s: sums to one [T1], this should be zero: %2.8e \n',cmdname{k_c},abs(sum(j1d)-1))
    fprintf('%-24s: mean %+2.4f against %+2.4f [T2], error %2.3e \n',cmdname{k_c},m,j1mean,abs(m-j1mean))
    fprintf('%-24s: variance %2.5f against %2.5f [T2], error %2.3e \n',cmdname{k_c},v,j1var,abs(v-j1var))
end

%% The documented error case
% A vector initialj1mewz with a scalar or absent initialj1sigmaz would be "a normal distribution
% with a vector mean", which is not a thing. Erroring beats guessing.
fprintf('\n--- a vector initialj1mewz without a vector initialj1sigmaz must error --- \n')
for k_c=1:2
    o=struct(); o.verbose=0; o.initialj1mewz=[0.2;-0.3];
    try
        if k_c==1
            discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,o);
        else
            discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,3,o);
        end
        fprintf('%-24s: ran without error [T0], this should not happen \n',cmdname{k_c})
    catch ME
        fprintf('%-24s: errors, as it should. The message is: %s \n',cmdname{k_c},ME.message)
    end
end

%% The initial condition must propagate through the whole age profile, not just age 1
% Setting initialj0mewz moves mewz(1), and the recursion carries that forward as
% rho(2)*rho(3)*...*rho(j) times the shift. A command that applied the initial condition to
% jequaloneDistz but not to the grids would pass every age-1 check above and fail this one.
fprintf('\n--- the initial condition must propagate through the whole age profile --- \n')
shift_true=zeros(1,J);
shift_true(1)=r1*m0;
for j_c=2:J
    shift_true(j_c)=rho(j_c)*shift_true(j_c-1);
end
for k_c=1:2
    oA=struct(); oA.verbose=0;
    oB=struct(); oB.verbose=0; oB.initialj0mewz=m0;
    if k_c==1
        [zA,pA,jA]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,oA);
        [zB,pB,jB]=discretizeLifeCycleAR1wGM_KFTT(mew,rho,p,mu,sd,znum,J,oB);
    else
        [zA,pA,jA]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,3,oA);
        [zB,pB,jB]=discretizeLifeCycleAR1wGM_Tauchen(mew,rho,p,mu,sd,znum,J,3,oB);
    end
    [mA,~,~]=MarkovChainMoments_FHorz(gather(zA),gather(pA),gather(jA));
    [mB,~,~]=MarkovChainMoments_FHorz(gather(zB),gather(pB),gather(jB));
    d=abs((mB(:)'-mA(:)')-shift_true);
    [worstd,worstj]=max(d);
    fprintf('%-24s: the mean profile shifts by exactly prod(rho)*initialj0mewz [T2], worst age j=%i, error %2.3e \n',cmdname{k_c},worstj,worstd)
end
fprintf('   (the shift decays from %2.4f at age 1 to %2.4f at age %i, so this is a shape test) \n',shift_true(1),shift_true(J),J)

output.shift_true=shift_true;

end
