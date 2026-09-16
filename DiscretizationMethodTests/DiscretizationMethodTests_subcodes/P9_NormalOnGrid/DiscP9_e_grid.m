function output=DiscP9_e_grid(calib,figure_c)
% P9: the e_grid option, and the three-command comparison it makes possible
%
% THE OPTION IS NEW, added for this block. discretizeIIDNormal_Tauchen and
% discretizeIIDNormal_TanakaToda both now accept an e_grid in their options, so three separate
% commands can be handed the same grid and asked for probabilities on it. Nothing here existed
% before, so the first job is the option itself and only then the comparison.
%
% THE THREE ARE NOT TRYING TO DO THE SAME THING, and the comparison is only worth reading once that
% is stated. Tauchen and MVNormal both integrate the normal density over the bin around each grid
% point - same rule, different code - so they must agree to machine precision, and that is a [T1]
% identity. Tanaka-Toda does something else entirely: it solves a maximum entropy problem that
% MATCHES MOMENTS on the given grid. On a grid wide enough to carry them it should therefore be
% better on the moments and worse as an approximation to the bin probabilities. Neither is "right";
% which one you want depends on whether the grid is there to represent a density or to carry moments.

fprintf('\n========== P9: the e_grid option ========== \n')

output=struct();
mew=calib.mew; sigma=calib.sigma; enum=calib.enum;
gnames={'even','nonuniform','offcentre'};

%% 1. The option is read at all
% On the EVEN grid, Tauchen's own construction coincides with what is being passed in, so handing it
% that grid must reproduce the no-e_grid call exactly. This is the one grid where an exact identity
% is available, and it is the check that the option is wired to the same arithmetic rather than to a
% parallel code path.
eg_even=calib.grid.even;
t0=struct(); t0.parallel=1;
[egA,peA]=discretizeIIDNormal_Tauchen(mew,sigma,enum,3,t0);
te=struct(); te.parallel=1; te.e_grid=eg_even;
[egB,peB]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],te);
egA=gather(egA); peA=gather(peA); egB=gather(egB); peB=gather(peB);
fprintf('\n--- 1. handing Tauchen the grid it would have built anyway --- \n')
fprintf('the even grid IS what Tauchen_q=3 constructs [T1], this should be zero: %2.8e \n',max(abs(egA-eg_even)))
fprintf('e_grid reproduces the no-e_grid call: e_grid [T0], this should be zero: %2.8e \n',max(abs(egB-egA)))
fprintf('e_grid reproduces the no-e_grid call: pi_e [T0], this should be zero: %2.8e \n',max(abs(peB-peA)))
fprintf('and Tauchen_q was ignored, since it was passed as empty and the call still worked [T0], this should be one: %i \n',all(isfinite(peB)))

%% 2. Tauchen's e_grid against MVNormal at M=1, on all three grids
% Same rule, separate implementations, so this is an identity rather than a comparison. The
% non-uniform and off-centre grids are where it has content: an evenly spaced grid would let a
% command get away with a constant bin width, and these two will not.
fprintf('\n--- 2. discretizeIIDNormal_Tauchen(e_grid) against MVNormal_ProbabilitiesOnGrid(M=1) --- \n')
mvo=struct(); mvo.parallel=1; mvo.verbose=0;
gapTM=zeros(1,3);
for g_c=1:3
    eg=calib.grid.(gnames{g_c});
    tg=struct(); tg.parallel=1; tg.e_grid=eg;
    [~,pT]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],tg);
    pM=MVNormal_ProbabilitiesOnGrid(eg,mew,sigma^2,enum,mvo);
    pT=gather(pT); pM=gather(pM);
    gapTM(g_c)=max(abs(pT(:)-pM(:)));
    fprintf('%11s grid: the two agree [T1], this should be zero: %2.8e \n',gnames{g_c},gapTM(g_c))
    fprintf('%11s grid: Tauchen''s probabilities sum to one [T1], this should be zero: %2.8e \n',gnames{g_c},abs(sum(pT)-1))
    fprintf('%11s grid: MVNormal''s sum to one [T1], this should be zero: %2.8e \n',gnames{g_c},abs(sum(pM(:))-1))
end

%% 3. All three commands, on all three grids
% The moment comparison. Tanaka-Toda is matching moments by construction so it should win on them;
% the other two are integrating the density so they should win on being the density. Both numbers are
% reported for all three, and neither is asserted to be smaller than the other - the ORDERING is what
% is asserted, and only where theory says it must hold.
fprintf('\n--- 3. the three commands on the same grid --- \n')
fprintf('%11s %14s %12s %12s %12s \n','grid','command','mean err','var err','sum-1');
errm=zeros(3,3); errv=zeros(3,3);
for g_c=1:3
    eg=calib.grid.(gnames{g_c});
    for m_c=1:3
        if m_c==1
            tg=struct(); tg.parallel=1; tg.e_grid=eg;
            [~,pp]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],tg); nm='Tauchen';
        elseif m_c==2
            sg=struct(); sg.parallel=1; sg.verbose=0; sg.e_grid=eg; sg.nMoments=2;
            [~,pp]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,sg); nm='TanakaToda';
        else
            pp=MVNormal_ProbabilitiesOnGrid(eg,mew,sigma^2,enum,mvo); nm='MVNormal';
        end
        pp=gather(pp); pp=pp(:);
        mh=sum(pp.*eg); vh=sum(pp.*(eg-mh).^2);
        errm(g_c,m_c)=abs(mh-mew); errv(g_c,m_c)=abs(vh-sigma^2);
        fprintf('%11s %14s %12.2e %12.2e %12.2e \n',gnames{g_c},nm,errm(g_c,m_c),errv(g_c,m_c),abs(sum(pp)-1));
    end
end
% Tanaka-Toda with nMoments=2 SOLVES for the mean and variance on this grid, so unless the solve
% fell back it must beat a rule that never looked at them. This is the one ordering theory pins.
% ...but only where there is something to order. On a grid symmetric about the mean, Tauchen's mean
% is exact BY SYMMETRY - the 2026-09-16 run had it at 5.55e-17 on the even and non-uniform grids -
% so comparing the two means compares two numbers that are both zero to within rounding, and which
% comes out smaller is noise. The 2026-09-16 run duly failed that comparison on the non-uniform grid
% with 1.94e-15 against 5.55e-17: a fail with no content. So the ordering is asserted against
% max(Tauchen's error, machine precision), which lets a tie at the floor pass and still catches
% Tanaka-Toda genuinely losing. The off-centre grid is the case that carries the comparison, and
% there it is not close: 2.87e-12 against 3.22e-03.
mprec=1e-14;
for g_c=1:3
    fprintf('%11s grid: TanakaToda is at least as good as Tauchen on the mean, or both are at the machine-precision floor [T2], this should be one: %i \n',gnames{g_c},errm(g_c,2)<=max(errm(g_c,1),mprec))
    fprintf('%11s grid: TanakaToda beats Tauchen on the variance [T2], this should be one: %i \n',gnames{g_c},errv(g_c,2)<=errv(g_c,1))
end
fprintf('the off-centre grid is the one where the mean comparison has content: Tauchen %2.2e against TanakaToda %2.2e [T2], this should be one: %i \n',errm(3,1),errm(3,2),errm(3,2)<errm(3,1)/100)

%% 4. TanakaToda's e_grid against its own 'even' construction
% The same identity as check 1, for the other command: on the grid its method='even' would have
% built, passing that grid explicitly must reproduce the call exactly. The translation this relies on
% is that a user grid takes the 'even' prior - W=ones - which is the only one of the four built-in
% methods whose weights are not tied to a specific set of nodes.
nS=3;
sa=struct(); sa.parallel=1; sa.verbose=0; sa.method='even'; sa.nSigmas=nS; sa.nMoments=2;
[egE,peE]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,sa);
egE=gather(egE); peE=gather(peE);
sb=struct(); sb.parallel=1; sb.verbose=0; sb.e_grid=egE; sb.nMoments=2;
[egF,peF]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,sb);
egF=gather(egF); peF=gather(peF);
fprintf('\n--- 4. TanakaToda: e_grid against its own method=''even'' on the same grid --- \n')
fprintf('e_grid returns the grid it was given [T0], this should be zero: %2.8e \n',max(abs(egF-egE)))
fprintf('e_grid reproduces the method=''even'' call: pi_e [T0], this should be zero: %2.8e \n',max(abs(peF-peE)))
% ...and it must NOT reproduce gauss-hermite, whose grid is different, or the option is being ignored
% and the check above is vacuous.
sc=struct(); sc.parallel=1; sc.verbose=0; sc.method='gauss-hermite'; sc.nMoments=2;
[egG,~]=discretizeIIDNormal_TanakaToda(mew,sigma,enum,sc);
fprintf('and the gauss-hermite grid is genuinely different from the even one [T0], this should be one: %i \n',max(abs(gather(egG)-egE))>1e-06)

%% 5. Error paths on the new option
fprintf('\n--- 5. error paths on e_grid --- \n')
try
    tb=struct(); tb.parallel=1; tb.e_grid=calib.grid.even(1:end-1);
    discretizeIIDNormal_Tauchen(mew,sigma,enum,[],tb);
    fprintf('Tauchen: an e_grid of the wrong length was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('Tauchen: an e_grid of the wrong length errors, and the message names e_grid, this should be one: %i \n',contains(ME.message,'e_grid'))
end
try
    tb=struct(); tb.parallel=1; tb.e_grid=flipud(calib.grid.even);
    discretizeIIDNormal_Tauchen(mew,sigma,enum,[],tb);
    fprintf('Tauchen: a descending e_grid was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('Tauchen: a descending e_grid errors, and the message says ascending, this should be one: %i \n',contains(ME.message,'ascend'))
end
try
    sb2=struct(); sb2.parallel=1; sb2.verbose=0; sb2.e_grid=calib.grid.even(1:end-1);
    discretizeIIDNormal_TanakaToda(mew,sigma,enum,sb2);
    fprintf('TanakaToda: an e_grid of the wrong length was accepted [T0], this should be one: 0 \n')
catch ME
    fprintf('TanakaToda: an e_grid of the wrong length errors, and the message names e_grid, this should be one: %i \n',contains(ME.message,'e_grid'))
end
% A row vector must be accepted and treated as a column, which is what the other commands do.
tr=struct(); tr.parallel=1; tr.e_grid=calib.grid.even';
[egR,peR]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],tr);
fprintf('Tauchen: an e_grid passed as a ROW vector is accepted and returned as a column [T0], this should be zero: %2.8e \n',max(abs(gather(egR)-calib.grid.even)))
fprintf('Tauchen: and gives the same probabilities as the column form [T0], this should be zero: %2.8e \n',max(abs(gather(peR)-peB)))

%% 6. parallel=0 vs 1 vs 2 on a user grid
% The e_grid branch is duplicated across the cpu and gpu paths in discretizeIIDNormal_Tauchen, so it
% is exactly the kind of place the two can drift apart. P1 checks this for the constructed grid; this
% checks it for a user grid, which is separate code.
eg=calib.grid.nonuniform;
o0=struct(); o0.parallel=0; o0.e_grid=eg;
o1=struct(); o1.parallel=1; o1.e_grid=eg;
[g0,p0]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],o0);
[g1,p1]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],o1);
fprintf('\n--- 6. the e_grid branch on cpu against gpu --- \n')
fprintf('parallel=0 vs 1 on a user grid: pi_e [T0], this should be zero: %2.8e \n',max(abs(gather(p0)-gather(p1))))
if gpuDeviceCount>0
    o2=struct(); o2.parallel=2; o2.e_grid=eg;
    [g2,p2]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],o2);
    fprintf('parallel=2 returns gpuArrays [T0], this should be one: %i \n',isa(g2,'gpuArray')&&isa(p2,'gpuArray'))
    fprintf('parallel=1 vs 2 on a user grid: e_grid [T0], this should be zero: %2.8e \n',max(abs(gather(g1)-gather(g2))))
    fprintf('parallel=1 vs 2 on a user grid: pi_e [T0], this should be zero: %2.8e \n',max(abs(gather(p1)-gather(p2))))
    fprintf('   (both branches run the same erfc expression, so any residual is the cpu and gpu erfc \n')
    fprintf('   libraries disagreeing in the last bit, not a different formula.) \n')
else
    fprintf('parallel=2 skipped: no gpu on this machine \n')
end

output.errm=errm; output.errv=errv; output.gapTM=gapTM;

%% Figure
figure(figure_c)
subplot(1,2,1)
for g_c=1:3
    eg=calib.grid.(gnames{g_c});
    tg=struct(); tg.parallel=1; tg.e_grid=eg;
    [~,pT]=discretizeIIDNormal_Tauchen(mew,sigma,enum,[],tg);
    plot(eg,gather(pT),'o-'); hold on
end
hold off
xlabel('e'); ylabel('probability'); title('Tauchen on three given grids')
legend(gnames,'Location','best')
subplot(1,2,2)
bar(errm)
set(gca,'XTickLabel',gnames); set(gca,'YScale','log')
ylabel('|mean error|'); title('mean error by command')
legend({'Tauchen','TanakaToda','MVNormal'},'Location','best')
sgtitle('P9: a normal distribution on a given grid')

end
