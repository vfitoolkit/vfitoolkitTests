function output=DiscP6_initialconditions(calib)
% P6: the initial-condition options, and jequaloneDistz
%
% Its own subcode because this is the largest option surface in the bank - four options across two
% commands, each alone and in the documented pairs - and because it is where B5, B6, B7, B10 and
% B11 all landed. The default is z(0)=0 as a point, so z(1) is N(mew(1),sigma(1)^2). The options
% replace that, either by giving period 0 a mean and/or a standard deviation, or by setting
% period 1 directly.
%
% jequaloneDistz is tested here rather than in the method subcodes, because it is the object the
% initial conditions actually determine. Three things have to hold and they are different claims:
% it sums to one; it lives on z_grid_J(:,1) rather than on some other grid; and its moments are
% mewz(1) and sigmaz(1). discretizeLifeCycleAR1_KFTT warns internally when its own grid check
% fails, so a silent run is itself part of the result.

fprintf('\n========== P6: initial conditions and jequaloneDistz ========== \n')

output=struct();
J=calib.J; znum=15;
mew=calib.vary.mew; rho=calib.vary.rho; sigma=calib.vary.sigma;

% The four options, and the pairs. Each row is: a label, then the option struct fields, then the
% period-1 mean and standard deviation the recursion says should result.
%   default:            z0=0 point           -> mewz(1)=mew(1),          sigmaz(1)=sigma(1)
%   initialj0mewz=m:    z0=m point           -> mewz(1)=mew(1)+rho(1)*m, sigmaz(1)=sigma(1)
%   initialj0sigmaz=s:  z0~N(0,s^2)          -> mewz(1)=mew(1),          sigmaz(1)=sqrt(rho(1)^2*s^2+sigma(1)^2)
%   both j0:            z0~N(m,s^2)          -> mewz(1)=mew(1)+rho(1)*m, sigmaz(1)=sqrt(rho(1)^2*s^2+sigma(1)^2)
%   initialj1mewz=m:    period 1 set direct  -> mewz(1)=m,               sigmaz(1)=0
%   initialj1sigmaz=s:  period 1 set direct  -> mewz(1)=0,               sigmaz(1)=s
%   both j1:            period 1 set direct  -> mewz(1)=m,               sigmaz(1)=s
m0=0.3; s0=0.25; m1=0.4; s1=0.35;
r1=rho(1); mw1=mew(1); sg1=sigma(1);
cfg={};
cfg{1}={'default',            struct(),                                                     mw1,              sg1};
cfg{2}={'initialj0mewz',      struct('initialj0mewz',m0),                                   mw1+r1*m0,        sg1};
cfg{3}={'initialj0sigmaz',    struct('initialj0sigmaz',s0),                                 mw1,              sqrt(r1^2*s0^2+sg1^2)};
cfg{4}={'both initialj0',     struct('initialj0mewz',m0,'initialj0sigmaz',s0),              mw1+r1*m0,        sqrt(r1^2*s0^2+sg1^2)};
cfg{5}={'initialj1mewz',      struct('initialj1mewz',m1),                                   m1,               0};
cfg{6}={'initialj1sigmaz',    struct('initialj1sigmaz',s1),                                 0,                s1};
cfg{7}={'both initialj1',     struct('initialj1mewz',m1,'initialj1sigmaz',s1),              m1,               s1};

cmdname={'LifeCycleAR1_KFTT','LifeCycleAR1_FellaGallipoliPanTauchen'};
for k_c=1:2
    fprintf('\n--- %s --- \n',cmdname{k_c})
    for c_c=1:7
        opts=cfg{c_c}{2}; lbl=cfg{c_c}{1}; mtrue=cfg{c_c}{3}; strue=cfg{c_c}{4};
        if k_c==1
            [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,opts);
        else
            [z_grid_J,pi_z_J,jequaloneDistz]=discretizeLifeCycleAR1_FellaGallipoliPanTauchen(mew,rho,sigma,znum,J,opts);
        end
        z_grid_J=gather(z_grid_J); pi_z_J=gather(pi_z_J); jequaloneDistz=gather(jequaloneDistz);
        g1=z_grid_J(:,1);
        m=sum(jequaloneDistz(:).*g1);
        s=sqrt(sum(jequaloneDistz(:).*(g1-m).^2));
        fprintf('%-18s: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',lbl,abs(sum(jequaloneDistz)-1))
        fprintf('%-18s: its mean is %+2.4f against %+2.4f [T2], error %2.3e \n',lbl,m,mtrue,abs(m-mtrue))
        fprintf('%-18s: its std dev is %2.4f against %2.4f [T2], error %2.3e \n',lbl,s,strue,abs(s-strue))
        % A degenerate period 1 (sigmaz(1)=0) is a point mass, and the only honest test of that is
        % that exactly one grid point carries all the mass.
        if strue==0
            fprintf('%-18s: sigmaz(1)=0, so it is a point mass [T0], this should be one: %i \n',lbl,sum(jequaloneDistz>10^(-12))==1)
        end
        output.(sprintf('cmd%i',k_c)).cfg(c_c).label=lbl;
        output.(sprintf('cmd%i',k_c)).cfg(c_c).moments=[m,s];
        output.(sprintf('cmd%i',k_c)).cfg(c_c).truth=[mtrue,strue];
    end
end

%% FellaGallipoliPan supports only initialj0sigmaz, and takes no mew at all
% It is the odd one out of the three, so its one supported option is tested on its own terms: the
% driftless recursion, where the period 1 standard deviation is sqrt(rho(1)^2*s0^2+sigma(1)^2).
fprintf('\n--- LifeCycleAR1_FellaGallipoliPan (initialj0sigmaz is its only initial-condition option) --- \n')
for c_c=[1,3]
    opts=cfg{c_c}{2}; lbl=cfg{c_c}{1};
    if c_c==1
        strue=sg1;
    else
        strue=sqrt(r1^2*s0^2+sg1^2);
    end
    [z_grid_J,~,jequaloneDistz]=discretizeLifeCycleAR1_FellaGallipoliPan(rho,sigma,znum,J,opts);
    z_grid_J=gather(z_grid_J); jequaloneDistz=gather(jequaloneDistz);
    g1=z_grid_J(:,1);
    m=sum(jequaloneDistz(:).*g1);
    s=sqrt(sum(jequaloneDistz(:).*(g1-m).^2));
    fprintf('%-18s: jequaloneDistz sums to one [T1], this should be zero: %2.8e \n',lbl,abs(sum(jequaloneDistz)-1))
    fprintf('%-18s: its mean is %+2.4f against 0 (no drift) [T2], error %2.3e \n',lbl,m,abs(m))
    fprintf('%-18s: its std dev is %2.4f against %2.4f [T2], error %2.3e \n',lbl,s,strue,abs(s-strue))
end

%% The whole age profile must shift with the initial condition, not just age 1
% Setting initialj0mewz moves mewz(1), and the recursion then propagates that through every later
% age as rho(2)*rho(3)*...*rho(j) times the shift. If a command applied the initial condition to
% jequaloneDistz but not to the grids, age 1 would look right and every later age would be wrong -
% which is exactly the kind of error the age-1-only checks above cannot see.
fprintf('\n--- the initial condition must propagate through the whole age profile --- \n')
optsA=struct(); optsB=struct(); optsB.initialj0mewz=m0;
[zA,pA,jA]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,optsA);
[zB,pB,jB]=discretizeLifeCycleAR1_KFTT(mew,rho,sigma,znum,J,optsB);
[mA,~,~]=MarkovChainMoments_FHorz(gather(zA),gather(pA),gather(jA));
[mB,~,~]=MarkovChainMoments_FHorz(gather(zB),gather(pB),gather(jB));
shift_true=zeros(1,J);
shift_true(1)=r1*m0;
for j_c=2:J
    shift_true(j_c)=rho(j_c)*shift_true(j_c-1);
end
d=abs((mB(:)'-mA(:)')-shift_true);
[worstd,worstj]=max(d);
fprintf('the mean profile shifts by exactly prod(rho)*initialj0mewz [T2], worst age j=%i, error %2.3e \n',worstj,worstd)
fprintf('   (the shift decays from %2.4f at age 1 to %2.4f at age %i, so this is not a constant offset) \n',shift_true(1),shift_true(J),J)
output.shift_true=shift_true;
output.shift_err=d;

end
