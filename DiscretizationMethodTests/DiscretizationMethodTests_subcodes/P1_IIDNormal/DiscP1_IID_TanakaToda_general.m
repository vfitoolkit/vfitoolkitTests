function output=DiscP1_IID_TanakaToda_general(calib,figure_c)
% P1: the generalised options of discretizeIID_TanakaToda
%
% discretizeIID_TanakaToda is the one command in the family that is meant to grow beyond the normal -
% its normal-only twin discretizeIIDNormal_TanakaToda exists precisely so it can. It has just gained
% distribution, distparams, truncate, e_grid, targetmoments, prior and masspoints. This subcode is
% the first exercise of all of them.
%
% THE REGRESSION CHECK IS FREE, AND IT IS THE MOST VALUABLE THING HERE. Every new option defaults
% off, so the default path must be bit-for-bit what it was. Normally that needs baseline numbers
% saved from before the change; here it does not, because discretizeIIDNormal_TanakaToda is a frozen
% copy of exactly that default path and is NOT being generalised. So the copy identity IS the
% regression test, it needs no maintenance, and it keeps working forever. It is run below across all
% four methods and all four nMoments rather than at defaults only, which is what DiscP1_crosstests
% already does.
%
% TRUTH COMES FROM CLOSED FORMS, deliberately. The command computes the moments of every non-normal
% family by quadrature on the density - one code path, so it composes with truncation. If this
% subcode also used quadrature it would be checking the routine against itself, so every target here
% is written down in closed form instead.

fprintf('\n========== P1: discretizeIID_TanakaToda, the generalised options ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma;

%% 1. THE REGRESSION: the default path has not moved
% All four methods, all four nMoments, three grid sizes. If any of the new options leaked into the
% default path, this is where it shows.
fprintf('\n--- 1. the default path against the frozen copy discretizeIIDNormal_TanakaToda --- \n')
methods={'even','gauss-legendre','clenshaw-curtis','gauss-hermite'};
worstg=0; worstp=0; ncheck=0;
for m_c=1:4
    for mm=1:4
        for ee=[5,9,15]
            oa=struct(); oa.verbose=0; oa.parallel=1; oa.method=methods{m_c}; oa.nMoments=mm;
            [ega,pea]=discretizeIID_TanakaToda(mew,sigma,ee,oa);
            [egb,peb]=discretizeIIDNormal_TanakaToda(mew,sigma,ee,oa);
            worstg=max(worstg,max(abs(gather(ega)-gather(egb))));
            worstp=max(worstp,max(abs(gather(pea)-gather(peb))));
            ncheck=ncheck+1;
        end
    end
end
fprintf('across %i combinations of method, nMoments and enum: e_grid [T0], this should be zero: %2.8e \n',ncheck,worstg)
fprintf('across %i combinations of method, nMoments and enum: pi_e [T0], this should be zero: %2.8e \n',ncheck,worstp)
fprintf('   (discretizeIIDNormal_TanakaToda is the normal-only twin and is not being generalised, so \n')
fprintf('   it is a frozen reference. Any drift in the default path breaks this and nothing else.) \n')

%% 2. the new otheroutputs fields
fprintf('\n--- 2. otheroutputs --- \n')
o2=struct(); o2.verbose=0; o2.parallel=1; o2.nMoments=4;
[~,pe2,oo2]=discretizeIID_TanakaToda(mew,sigma,15,o2);
fprintf('nMoments is reported and is at most what was requested [T0], this should be one: %i \n',oo2.nMoments<=4)
fprintf('momentError is finite at the accepted solution [T0], this should be one: %i \n',isfinite(oo2.momentError))
fprintf('TBar is the normal closed form [0,sigma^2,0,3sigma^4] [T1], this should be zero: %2.8e \n',max(abs(oo2.TBar(:)-[0;sigma^2;0;3*sigma^4])))
fprintf('mew_eff is mew on the default path [T1], this should be zero: %2.8e \n',abs(oo2.mew_eff-mew))
fprintf('q has one entry per grid point [T0], this should be zero: %i \n',length(oo2.q)~=15)
% THE FALLBACK MUST BE VISIBLE. A grid far too wide for its node count fails silently and falls back
% on the prior; the note that prompted this work records a case where that put the standard deviation
% out by a factor of four while every probability still summed to one. nMoments is how a caller sees it.
% A GRID TOO NARROW, not too wide. The first version of this used nSigmas=25 at enum=5 and did not
% demonstrate anything: a wide grid can still carry the variance by putting a little weight far out -
% at +-12.5 sigma it needs p=1/(2*12.5^2)=0.0032, which is perfectly feasible - so the solve matched
% three moments and returned the right standard deviation, underneath a line of text claiming it had
% been poisoned. A narrow grid is the real hazard: on [-0.5 sigma, +0.5 sigma] the largest variance
% ANY distribution can have is 0.25 sigma^2, so the target is infeasible, the solve must fall back,
% and the standard deviation cannot exceed half the truth however many points are used.
%
% The grid is passed through e_grid rather than nSigmas so the command's own nSigmas<1.2 warning
% does not fire - here the narrowness is the point of the test, not a mistake to be warned about.
o2b=struct(); o2b.verbose=0; o2b.parallel=1; o2b.nMoments=4;
o2b.e_grid=(mew+linspace(-0.5*sigma,0.5*sigma,5))';
[eg2b,pe2b,oo2b]=discretizeIID_TanakaToda(mew,sigma,5,o2b);
eg2b=gather(eg2b); pe2b=gather(pe2b);
sdb=sqrt(sum(pe2b.*(eg2b-sum(pe2b.*eg2b)).^2));
fprintf('on a grid too narrow to carry the variance the solve falls back [T0], this should be one: %i \n',oo2b.nMoments<2)
fprintf('   it reports nMoments=%i, and the standard deviation comes out %2.4f against the true %2.4f - \n',oo2b.nMoments,sdb,sigma)
fprintf('   it cannot exceed half the truth on this grid, whatever the solver does. \n')
fprintf('the standard deviation really is badly wrong, which is the hazard [T0], this should be one: %i \n',sdb<0.6*sigma)
fprintf('and yet every probability still sums to one, so nMoments is the ONLY way a caller sees it [T1], this should be zero: %2.8e \n',abs(sum(pe2b)-1))
fprintf('and pi_e still sums to one even then [T1], this should be zero: %2.8e \n',abs(sum(gather(pe2b))-1))

%% 3. e_grid
fprintf('\n--- 3. e_grid --- \n')
eg_even=(mew+linspace(-3*sigma,3*sigma,15))';
o3=struct(); o3.verbose=0; o3.parallel=1; o3.method='even'; o3.nSigmas=3; o3.nMoments=2;
[egA,peA]=discretizeIID_TanakaToda(mew,sigma,15,o3);
o3b=struct(); o3b.verbose=0; o3b.parallel=1; o3b.e_grid=eg_even; o3b.nMoments=2;
[egB,peB]=discretizeIID_TanakaToda(mew,sigma,15,o3b);
fprintf('the grid method=even builds IS the grid being passed in [T1], this should be zero: %2.8e \n',max(abs(gather(egA)-eg_even)))
fprintf('passing it explicitly reproduces the constructed call: e_grid [T0], this should be zero: %2.8e \n',max(abs(gather(egB)-gather(egA))))
fprintf('passing it explicitly reproduces the constructed call: pi_e [T0], this should be zero: %2.8e \n',max(abs(gather(peB)-gather(peA))))
% a non-uniform grid, where a constant bin width would not do
uu=linspace(-1,1,15)';
eg_nu=mew+3*sigma*sinh(2*uu)/sinh(2);
o3c=struct(); o3c.verbose=0; o3c.parallel=1; o3c.e_grid=eg_nu; o3c.nMoments=2;
[egC,peC]=discretizeIID_TanakaToda(mew,sigma,15,o3c);
egC=gather(egC); peC=gather(peC);
mC=sum(peC.*egC); vC=sum(peC.*(egC-mC).^2);
fprintf('on a NON-UNIFORM grid it still matches the mean [T2], this should be below 1e-08: %2.3e \n',abs(mC-mew))
fprintf('on a NON-UNIFORM grid it still matches the variance [T2], this should be below 1e-08: %2.3e \n',abs(vC-sigma^2))
fprintf('and it is returned as a column matching what was given [T0], this should be zero: %2.8e \n',max(abs(egC-eg_nu)))

%% 4. targetmoments
fprintf('\n--- 4. targetmoments --- \n')
% Handing it the normal's own moments must reproduce the normal. THE FIRST ENTRY IS THE MEAN, not
% zero - the moment functions are centred on it - and mew is deliberately non-zero here so getting
% that convention backwards is visible rather than silent.
o4=struct(); o4.verbose=0; o4.parallel=1; o4.method='even'; o4.nSigmas=3; o4.nMoments=4;
[egD,peD]=discretizeIID_TanakaToda(mew,sigma,15,o4);
o4b=o4; o4b.targetmoments=[mew,sigma^2,0,3*sigma^4];
[egE,peE]=discretizeIID_TanakaToda(mew,sigma,15,o4b);
fprintf('targetmoments set to the normal''s own moments reproduces it: e_grid [T0], this should be zero: %2.8e \n',max(abs(gather(egE)-gather(egD))))
fprintf('targetmoments set to the normal''s own moments reproduces it: pi_e [T1], this should be below %g: %2.8e \n',1e-06,max(abs(gather(peE)-gather(peD))))
fprintf('   (T1 with a bar rather than a zero, and the bar is the SOLVER''S, not machine precision. \n')
fprintf('   The two calls reach the same targets by different routes - the closed form and the given \n')
fprintf('   vector - so each runs its own fminunc, and the entropy solve accepts once its moment \n')
fprintf('   error is below 1e-05. Two independent solves of the same problem therefore agree to about \n')
fprintf('   that, not to 1e-16. The run of 2026-09-17 measured 1.33e-09 against a 1e-10 bar.) \n')
% and the moments it was asked for are the moments it delivers
o4c=o4; o4c.targetmoments=[mew,sigma^2,0.4*sigma^3,3.5*sigma^4]; % a skewed, fat-tailed target
[egF,peF,ooF]=discretizeIID_TanakaToda(mew,sigma,21,o4c);
egF=gather(egF); peF=gather(peF);
mF=sum(peF.*egF); vF=sum(peF.*(egF-mF).^2);
s3F=sum(peF.*(egF-mF).^3); s4F=sum(peF.*(egF-mF).^4);
% THE BARS ARE THE SOLVER'S TOLERANCE, not something tighter that sounds impressive. The entropy
% solve stops once norm(momentError) is below 1e-05 on its SCALED moment conditions, so a delivered
% moment cannot be expected to beat that - asking for 1e-07 asks for more precision than the routine
% promises, and the run of 2026-09-17 duly failed on 9.4e-07 and 1.1e-06 with nothing wrong.
mombar=1e-05;
fprintf('a skewed target: mean delivered [T2], this should be below %g: %2.3e \n',mombar,abs(mF-mew))
fprintf('a skewed target: variance delivered [T2], this should be below %g: %2.3e \n',mombar,abs(vF-sigma^2))
if ooF.nMoments>=3
    fprintf('a skewed target: 3rd central moment delivered [T2], this should be below %g: %2.3e \n',mombar,abs(s3F-0.4*sigma^3))
end
if ooF.nMoments>=4
    fprintf('a skewed target: 4th central moment delivered [T2], this should be below %g: %2.3e \n',mombar,abs(s4F-3.5*sigma^4))
end
fprintf('a skewed target: it reports matching %i moments \n',ooF.nMoments)
fprintf('mew_eff is targetmoments(1), the MEAN and not zero [T1], this should be zero: %2.8e \n',abs(ooF.mew_eff-mew))

%% 5. the distribution families, against closed-form truth
fprintf('\n--- 5. distribution, against closed forms --- \n')
enum5=21;
for d_c=1:3
    od=struct(); od.verbose=0; od.parallel=1; od.nMoments=4; od.method='even';
    if d_c==1
        od.distribution='uniform'; od.distparams.lb=-1; od.distparams.ub=2; nm='uniform[-1,2]';
        Tm=0.5; Tv=(2-(-1))^2/12; Tsk=0; Tek=-6/5;
    elseif d_c==2
        od.distribution='exponential'; od.distparams.lambda=1.5; nm='exponential(1.5)';
        Tm=1/1.5; Tv=1/1.5^2; Tsk=2; Tek=6;
    else
        od.distribution='lognormal'; od.distparams.mu=0.1; od.distparams.sigma=0.4; nm='lognormal(0.1,0.4)';
        s2=0.4^2;
        Tm=exp(0.1+s2/2); Tv=(exp(s2)-1)*exp(2*0.1+s2);
        Tsk=(exp(s2)+2)*sqrt(exp(s2)-1); Tek=exp(4*s2)+2*exp(3*s2)+3*exp(2*s2)-6;
    end
    [egd,ped,ood]=discretizeIID_TanakaToda([],[],enum5,od);
    egd=gather(egd); ped=gather(ped);
    md=sum(ped.*egd); vd=sum(ped.*(egd-md).^2);
    skd=sum(ped.*(egd-md).^3)/vd^1.5; ekd=sum(ped.*(egd-md).^4)/vd^2-3;
    fprintf('%s: pi_e sums to one [T1], this should be zero: %2.8e \n',nm,abs(sum(ped)-1))
    fprintf('%s: the grid is strictly ascending [T0], this should be one: %i \n',nm,issorted(egd,'strictascend'))
    fprintf('%s: mean %2.6f against truth %2.6f, error [T2] %2.3e \n',nm,md,Tm,abs(md-Tm))
    fprintf('%s: variance %2.6f against truth %2.6f, error [T2] %2.3e \n',nm,vd,Tv,abs(vd-Tv))
    fprintf('%s: skewness %+2.4f against truth %+2.4f, excess kurtosis %+2.4f against truth %+2.4f (matched %i moments) \n',nm,skd,Tsk,ekd,Tek,ood.nMoments)
    % the TARGETS the command computed by quadrature, against the closed form. This separates a bad
    % quadrature from a bad entropy solve - if the targets are right and the moments are not, the
    % solve fell back; if the targets are wrong, the quadrature is.
    fprintf('%s: the target variance it computed by quadrature matches the closed form [T2], this should be below 1e-08: %2.3e \n',nm,abs(ood.TBar(2)-Tv))
    fprintf('%s: and the target mean [T2], this should be below 1e-08: %2.3e \n',nm,abs(ood.mew_eff-Tm))
end
% mew and sigma really are unused for a non-normal family
ou1=struct(); ou1.verbose=0; ou1.parallel=1; ou1.distribution='uniform'; ou1.distparams.lb=0; ou1.distparams.ub=1;
[egu1,peu1]=discretizeIID_TanakaToda([],[],15,ou1);
[egu2,peu2]=discretizeIID_TanakaToda(99,77,15,ou1);
fprintf('for a non-normal family mew and sigma are ignored entirely [T0], this should be zero: %2.8e \n',max(abs(gather(egu1)-gather(egu2)))+max(abs(gather(peu1)-gather(peu2))))

%% 6. truncate
fprintf('\n--- 6. truncate --- \n')
% TRUNCATING A UNIFORM GIVES ANOTHER UNIFORM, so this is an identity rather than an accuracy check -
% much the stronger statement, and it needs no closed form for a truncated density.
ot1=struct(); ot1.verbose=0; ot1.parallel=1; ot1.nMoments=4; ot1.method='even';
ot1.distribution='uniform'; ot1.distparams.lb=0; ot1.distparams.ub=1; ot1.truncate=[0.25,0.75];
[egt1,pet1]=discretizeIID_TanakaToda([],[],15,ot1);
ot2=struct(); ot2.verbose=0; ot2.parallel=1; ot2.nMoments=4; ot2.method='even';
ot2.distribution='uniform'; ot2.distparams.lb=0.25; ot2.distparams.ub=0.75;
[egt2,pet2]=discretizeIID_TanakaToda([],[],15,ot2);
fprintf('U[0,1] truncated to [0.25,0.75] equals U[0.25,0.75]: e_grid [T1], this should be zero: %2.8e \n',max(abs(gather(egt1)-gather(egt2))))
fprintf('U[0,1] truncated to [0.25,0.75] equals U[0.25,0.75]: pi_e [T1], this should be below %g: %2.8e \n',1e-06,max(abs(gather(pet1)-gather(pet2))))
fprintf('   (a bar, not a zero, for the same reason as section 4: two routes to the same targets means \n')
fprintf('   two entropy solves, and they agree to the solver''s tolerance. The grid above IS exact.) \n')
% truncating to the natural support changes nothing
ot3=struct(); ot3.verbose=0; ot3.parallel=1; ot3.nMoments=4; ot3.method='even';
ot3.distribution='exponential'; ot3.distparams.lambda=1.5;
[egt3,pet3]=discretizeIID_TanakaToda([],[],15,ot3);
ot4=ot3; ot4.truncate=[0,Inf];
[egt4,pet4]=discretizeIID_TanakaToda([],[],15,ot4);
fprintf('truncating an exponential to [0,Inf), its own support, changes nothing [T1], this should be zero: %2.8e \n',max(abs(gather(egt3)-gather(egt4)))+max(abs(gather(pet3)-gather(pet4))))
% ONE-SIDED truncation, which is the one place truncate and nSigmas interact: the infinite side
% falls back to nSigmas standard deviations, the finite side is taken from the bound.
ot5=struct(); ot5.verbose=0; ot5.parallel=1; ot5.nMoments=2; ot5.method='even'; ot5.nSigmas=3;
ot5.truncate=[mew-0.5*sigma,Inf];
[egt5,pet5,oot5]=discretizeIID_TanakaToda(mew,sigma,15,ot5);
egt5=gather(egt5); pet5=gather(pet5);
fprintf('one-sided truncate [a,Inf): the grid starts exactly at a [T1], this should be zero: %2.8e \n',abs(egt5(1)-(mew-0.5*sigma)))
fprintf('one-sided truncate [a,Inf): the upper edge is mean+nSigmas*sd, not infinite [T0], this should be one: %i \n',isfinite(egt5(end)))
% and the truncated normal's mean, in closed form, which is the check with content
al=(mew-0.5*sigma-mew)/sigma;
phial=exp(-0.5*al^2)/sqrt(2*pi);
Phial=0.5*erfc(-al/sqrt(2));
Ttm=mew+sigma*phial/(1-Phial);
fprintf('one-sided truncate: the target mean matches the truncated-normal closed form [T2], this should be below 1e-08: %2.3e \n',abs(oot5.mew_eff-Ttm))

%% 7. masspoints
fprintf('\n--- 7. masspoints --- \n')
om=struct(); om.verbose=0; om.parallel=1; om.nMoments=2; om.method='even'; om.nSigmas=3;
[egm0,pem0]=discretizeIID_TanakaToda(mew,sigma,15,om);
egm0=gather(egm0); pem0=gather(pem0);
% an atom OFF the grid: a new node appears and the grid stays sorted
atomv=mew+5*sigma; atomp=0.2;
om1=om; om1.masspoints=[atomv,atomp];
[egm1,pem1]=discretizeIID_TanakaToda(mew,sigma,15,om1);
egm1=gather(egm1); pem1=gather(pem1);
fprintf('an atom off the grid adds one node [T0], this should be one: %i \n',length(egm1)==length(egm0)+1)
fprintf('the grid is still strictly ascending [T0], this should be one: %i \n',issorted(egm1,'strictascend'))
fprintf('pi_e still sums to one [T1], this should be zero: %2.8e \n',abs(sum(pem1)-1))
fprintf('the new node carries exactly the atom''s probability [T1], this should be zero: %2.8e \n',abs(pem1(egm1==atomv)-atomp))
fprintf('the continuous part is scaled by 1-sum(atoms) [T1], this should be zero: %2.8e \n',max(abs(pem1(egm1~=atomv)-(1-atomp)*pem0)))
% an atom ON an existing node: it merges, and the node carries MORE than the atom's probability
om2=om; om2.masspoints=[egm0(8),atomp];
[egm2,pem2]=discretizeIID_TanakaToda(mew,sigma,15,om2);
egm2=gather(egm2); pem2=gather(pem2);
fprintf('an atom on an existing node adds no node [T0], this should be one: %i \n',length(egm2)==length(egm0))
fprintf('pi_e still sums to one [T1], this should be zero: %2.8e \n',abs(sum(pem2)-1))
fprintf('that node carries MORE than the atom''s probability [T0], this should be one: %i \n',pem2(8)>atomp)
fprintf('   it carries %2.4f against the atom''s %2.4f. THIS IS NOT A BUG: the continuous part also \n',pem2(8),atomp)
fprintf('   puts mass on that node and the grid cannot tell the two apart. The note that prompted \n')
fprintf('   this work records the same thing for the GMR zeta shock, atom 0.560 and node 0.714. \n')
fprintf('and it is exactly the atom plus the scaled continuous mass [T1], this should be zero: %2.8e \n',abs(pem2(8)-(atomp+(1-atomp)*pem0(8))))
% atoms carrying ALL the mass: the continuous part is scaled to zero and only the atoms remain
om3=om; om3.masspoints=[-1,0.3;0,0.3;1,0.4];
[egm3,pem3]=discretizeIID_TanakaToda(mew,sigma,15,om3);
egm3=gather(egm3); pem3=gather(pem3);
fprintf('atoms summing to one leave all the mass on the atoms [T1], this should be zero: %2.8e \n',abs(sum(pem3(ismember(egm3,[-1;0;1])))-1))
fprintf('and every other node has zero probability [T1], this should be zero: %2.8e \n',max(abs(pem3(~ismember(egm3,[-1;0;1])))))

%% 8. precedence
fprintf('\n--- 8. precedence --- \n')
% GRID: e_grid beats truncate beats nSigmas.
op=struct(); op.verbose=0; op.parallel=1; op.nMoments=2; op.method='even'; op.nSigmas=3;
op.distribution='uniform'; op.distparams.lb=0; op.distparams.ub=1; op.truncate=[0.2,0.8];
op.e_grid=linspace(0.3,0.7,11)';
[egp1,~]=discretizeIID_TanakaToda([],[],11,op);
fprintf('e_grid beats truncate for the grid [T1], this should be zero: %2.8e \n',max(abs(gather(egp1)-linspace(0.3,0.7,11)')))
% TARGETS: targetmoments beats distribution.
op2=struct(); op2.verbose=0; op2.parallel=1; op2.nMoments=2; op2.method='even';
op2.distribution='uniform'; op2.distparams.lb=0; op2.distparams.ub=1;
op2.targetmoments=[0.5,0.01,0,3*0.01^2];
[~,~,oop2]=discretizeIID_TanakaToda([],[],15,op2);
fprintf('targetmoments beats distribution for the targets [T1], this should be zero: %2.8e \n',abs(oop2.TBar(2)-0.01))
fprintf('   (the uniform on [0,1] has variance 1/12=%2.6f, so if distribution had won it would show.) \n',1/12)
% PRIOR: an explicit prior beats the density.
op3=struct(); op3.verbose=0; op3.parallel=1; op3.nMoments=2; op3.method='even'; op3.nSigmas=3;
[~,~,oop3a]=discretizeIID_TanakaToda(mew,sigma,15,op3);
op3b=op3; op3b.prior=ones(1,15);
[~,~,oop3b]=discretizeIID_TanakaToda(mew,sigma,15,op3b);
fprintf('an explicit prior beats the density [T0], this should be one: %i \n',max(abs(oop3b.q(:)-oop3a.q(:)))>1e-10)
fprintf('and the prior returned is the one that was given [T1], this should be zero: %2.8e \n',max(abs(oop3b.q(:)-ones(15,1))))

%% 9. error paths
fprintf('\n--- 9. error paths --- \n')
oe=struct(); oe.verbose=0; oe.parallel=1;
try
    oe1=oe; oe1.e_grid=linspace(-1,1,10)';
    discretizeIID_TanakaToda(mew,sigma,15,oe1);
    fprintf('an e_grid whose length disagrees with enum was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('an e_grid whose length disagrees with enum errors, this should be one: %i \n',contains(ME.message,'e_grid'))
end
try
    oe2=oe; oe2.distribution='weibull';
    discretizeIID_TanakaToda(mew,sigma,15,oe2);
    fprintf('an unsupported distribution was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('an unsupported distribution errors, and the message lists the supported ones, this should be one: %i \n',contains(ME.message,'lognormal'))
end
try
    oe3=oe; oe3.distribution='uniform'; oe3.distparams.lb=0;
    discretizeIID_TanakaToda([],[],15,oe3);
    fprintf('a uniform with no ub was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('a uniform missing distparams.ub errors, and the message names it, this should be one: %i \n',contains(ME.message,'distparams'))
end
try
    oe4=oe; oe4.distribution='uniform'; oe4.distparams.lb=0; oe4.distparams.ub=1; oe4.method='gauss-hermite';
    discretizeIID_TanakaToda([],[],15,oe4);
    fprintf('gauss-hermite on a non-normal was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('gauss-hermite on a non-normal errors, and the message says why, this should be one: %i \n',contains(ME.message,'gaussian'))
end
try
    oe5=oe; oe5.masspoints=[0,0.6;1,0.7];
    discretizeIID_TanakaToda(mew,sigma,15,oe5);
    fprintf('masspoints summing above one were accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('masspoints summing above one error, this should be one: %i \n',contains(ME.message,'masspoints'))
end
try
    oe6=oe; oe6.targetmoments=[0,1,0];
    discretizeIID_TanakaToda(mew,sigma,15,oe6);
    fprintf('a targetmoments of the wrong length was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('a targetmoments of the wrong length errors, this should be one: %i \n',contains(ME.message,'targetmoments'))
end

%% 10. the four defects found by review on 2026-09-17
% NONE OF SECTIONS 1 TO 9 WOULD HAVE CAUGHT ANY OF THESE, and the reason is worth recording because
% it is a property of how the sections are written rather than an oversight in any one of them.
% Every section above pins method='even' so that the option under test is the only thing varying -
% which is right, and is exactly what routed around the gauss-hermite branch where three of the four
% defects lived. The fourth was a value never checked: section 2 asserted momentError was finite and
% nothing more. A test suite that isolates variables will not find bugs in the paths it isolates
% away from, so those paths need their own section. This is it.
fprintf('\n--- 10. the paths that sections 1 to 9 isolate away from --- \n')

% (1) gauss-hermite cannot respect a finite truncation, and now says so instead of ignoring it.
% Before the fix, truncate=[-1,1] at enum=9 gave a grid spanning [-4.51,4.51] with eight of nine
% nodes OUTSIDE the interval, carrying 27% of the probability at points where the truncated density
% is zero - and the moments still came out right, because the entropy fit compensates.
try
    oq=struct(); oq.verbose=0; oq.parallel=1; oq.method='gauss-hermite'; oq.truncate=[-1,1];
    discretizeIID_TanakaToda(mew,sigma,9,oq);
    fprintf('gauss-hermite with a finite truncate was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('gauss-hermite with a finite truncate errors, and the message says why, this should be one: %i \n',contains(ME.message,'unbounded'))
end
% ...and the invariant that would have caught it in the first place: with a finite truncation every
% grid point must lie inside it. Checked on each method that IS allowed with one.
for m_c=1:3
    mn={'even','gauss-legendre','clenshaw-curtis'};
    ot=struct(); ot.verbose=0; ot.parallel=1; ot.method=mn{m_c}; ot.truncate=[mew-sigma,mew+sigma]; ot.nMoments=2;
    [egt,~]=discretizeIID_TanakaToda(mew,sigma,9,ot);
    egt=gather(egt);
    fprintf('%-16s with truncate=[a,b]: every grid point lies inside [a,b] [T0], this should be one: %i \n',mn{m_c},all(egt>=mew-sigma-1e-12 & egt<=mew+sigma+1e-12))
end

% (2) targetmoments with NO method set. The default was gauss-hermite, which needs mew and sigma -
% and targetmoments is precisely the case where they may be empty. It now defaults to 'even'.
ot2=struct(); ot2.verbose=0; ot2.parallel=1; ot2.nMoments=2;
ot2.targetmoments=[0.2,0.09,0,3*0.09^2];
[egv,pev,oov]=discretizeIID_TanakaToda([],[],11,ot2);
egv=gather(egv); pev=gather(pev);
mv=sum(pev.*egv); vv=sum(pev.*(egv-mv).^2);
fprintf('\ntargetmoments with no method and empty mew/sigma: it runs [T0], this should be one: %i \n',all(isfinite(pev)))
fprintf('targetmoments with no method: the mean is delivered [T2], this should be below 1e-05: %2.3e \n',abs(mv-0.2))
fprintf('targetmoments with no method: the variance is delivered [T2], this should be below 1e-05: %2.3e \n',abs(vv-0.09))
fprintf('targetmoments with no method: mew_eff is targetmoments(1) [T1], this should be zero: %2.8e \n',abs(oov.mew_eff-0.2))

% (3) targetmoments WITH distribution='normal' named explicitly. The old prior rule made this the
% breaking case and the unnamed one the working case - the opposite way round from expectation.
ot3=ot2; ot3.distribution='normal';
[egw,pew]=discretizeIID_TanakaToda([],[],11,ot3);
fprintf('\ntargetmoments with distribution=''normal'' named: it runs [T0], this should be one: %i \n',all(isfinite(gather(pew))))
fprintf('and naming the distribution changes nothing, since targetmoments overrides it [T0], this should be zero: %2.8e \n',max(abs(gather(pew)-pev))+max(abs(gather(egw)-egv)))

% (4) momentError must be the error of the fit that was KEPT. The ladder accepts only when the error
% is below 1e-05, so momentError<=1e-05 whenever nMoments>=2 is exactly the discriminator: under the
% old code a fall-back-to-two reported the REJECTED three- or four-moment error, which is above
% 1e-05 by definition - that is why it was rejected. The sweep is wide enough to produce fallbacks,
% and the count is asserted so this cannot pass by never falling back at all.
fprintf('\n');
nfall=0; nge2=0; worstacc=0;
for ee=[5,7,9,11,15,21]
    for nS=[2,3,4,6,8]
        om=struct(); om.verbose=0; om.parallel=1; om.nMoments=4; om.method='even'; om.nSigmas=nS;
        [~,~,oom]=discretizeIID_TanakaToda(mew,sigma,ee,om);
        if oom.nMoments<4
            nfall=nfall+1;
        end
        if oom.nMoments>=2
            nge2=nge2+1;
            worstacc=max(worstacc,oom.momentError);
        end
    end
end
fprintf('over 30 configurations, %i fell back below the 4 moments requested and %i kept 2 or more \n',nfall,nge2)
fprintf('at least one configuration fell back, so the next check is not vacuous [T0], this should be one: %i \n',nfall>0)
fprintf('momentError is the error of the ACCEPTED fit everywhere nMoments>=2 [T0], this should be below %g: %2.3e \n',1e-05,worstacc)
fprintf('   (the ladder accepts only below 1e-05, so a reported error above it could only be a \n')
fprintf('   REJECTED solve''s - which is what the fall-back branches used to report.) \n')

% and the validation that makes the empty-mew case say what is wrong
try
    discretizeIID_TanakaToda([],[],9);
    fprintf('\nan empty mew and sigma with no targetmoments was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('\nan empty mew/sigma with no targetmoments errors, and the message says so, this should be one: %i \n',contains(ME.message,'cannot be empty'))
end

output.regression_grid=worstg; output.regression_pi=worstp;
output.nfallbacks=nfall; output.worstacceptederror=worstacc;

%% Figure
figure(figure_c)
subplot(1,2,1)
plot(egm1,pem1,'o-'); hold on
plot(egm0,(1-atomp)*pem0,'s--'); hold off
xlabel('e'); ylabel('probability'); title('an atom off the grid')
legend('with atom','continuous part, scaled','Location','best')
subplot(1,2,2)
od=struct(); od.verbose=0; od.parallel=1; od.nMoments=4; od.method='even';
od.distribution='exponential'; od.distparams.lambda=1.5;
[ege,pee]=discretizeIID_TanakaToda([],[],21,od);
plot(gather(ege),gather(pee),'o-')
xlabel('e'); ylabel('probability'); title('exponential(1.5) on 21 points')
sgtitle('P1: discretizeIID\_TanakaToda, generalised')

end
