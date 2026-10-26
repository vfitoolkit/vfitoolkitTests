function output=EZFHorz_CrossTests_feasibility(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% DELIBERATE TRIGGER for two feasibility conventions of the Epstein-Zin raws (both horizons; this is
% the FHorz half, EZInfHorz_CrossTests_feasibility is the InfHorz half). Written 2026-10-27 to run RED
% on the raws as they stood, then green once they are fixed. The bank's other models have every
% state feasible and strictly signed returns, so nothing else in the bank can reach either convention.
%
% A: A STATE WITH NO FEASIBLE CHOICE (V=-Inf). a_grid is extended below the natural borrowing limit
%    with integer-spaced nodes from -12 to -1, so the lowest states cannot afford any aprime at the
%    terminal age, and the dead region grows backwards through the ages. The vNM solve of the same
%    budget (exoticpreferences not set) gives the reference dead set. Every Epstein-Zin raw carries
%    a -Inf continuation through the transform (ezc4*V')^ezc5 UNTOUCHED (the transform is applied to
%    finite entries only) and then through ezc3*beta*(.). For consumption units and positive utility
%    (ezc3=ezc4=1) -Inf stays -Inf and the dead state is correctly avoided. For NEGATIVE utility
%    (ezc3=ezc4=-1) the untouched -Inf comes out of ezc3*beta*(.) as +Inf, so a choice that leads to
%    a dead state is valued +Inf and is CHOSEN, and the +Inf then propagates as -Inf one step back,
%    and so on: V holds +Inf entries, and the dead set no longer matches the budget. The fix is one
%    line per site: with ezc4=-1 a -Inf continuation must become +Inf BEFORE the expectation
%    ((-V')^(1+risk) of a -Inf is +Inf), so that ezc3*beta*(.) returns it to -Inf.
%    Checks: no +Inf and no NaN in V; the -Inf set equals the vNM dead set; the dead set is the same
%    in all three cases (same budget); ValueFnFromPolicy agrees on the finite set and reproduces
%    the dead set; and the same census on the lowmemory, divide-and-conquer and grid-interpolation
%    legs (the fix is a sweep over every raw, and these legs are different files).
%
% B: A RETURN OF EXACTLY ZERO AT A FEASIBLE CHOICE. The floor ReturnFns give consumption below a
%    subsistence level xbar (feasible, c>0) a return of floorval; for negative utility the return is
%    satiated at floorval from above. Solved twice, floorval=0 and floorval=+-1e-10. Whatever the
%    right semantics of a zero return, V must be CONTINUOUS in it: the two solves must have the same
%    -Inf set and agree on V to the scale of 1e-10. Every Epstein-Zin raw currently rewrites a return
%    of exactly zero to -Inf (temp2(ReturnMatrix==0)=-Inf), so the floorval=0 solve treats every
%    sub-subsistence choice as infeasible and, where the floor is forced (a=0, low income), the state
%    becomes dead: the -Inf sets differ, which is the red line. The companion convention
%    temp(V==0)=0 (a continuation of exactly zero contributes nothing to the expectation) is also
%    exercised here: it is the solver's marker for an uninitialised V0, and in consumption units
%    with EIS<1 a genuine V=0 (zero composite good) is a catastrophe the certainty equivalent must
%    collapse to. The 1e-10 solve shows what the limit looks like; the zero solve must match it.
%    Anti-vacuity: the share of states whose chosen consumption is below the floor must be well
%    above zero, or the floor was never hit and the test says nothing.
%
% Each Epstein-Zin solve is wrapped in try/catch, so a crash prints as a failed check and the rest
% of the bank runs. The reference vNM solves are not wrapped: an error there must be seen.

N_z=prod(n_z);
DF=DiscountFactorParamNames;

ReturnFn_cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) EZReturnFn_cons_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_posU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_positiveUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_negU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZReturnFn_negativeUtils_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_cons_floor=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension,xbar,floorval) EZReturnFn_cons_nod_z_floor(aprime,a,z,r,w,kappa_j,agej,Jr,pension,xbar,floorval);
ReturnFn_posU_floor=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,xbar,floorval) EZReturnFn_positiveUtils_nod_z_floor(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,xbar,floorval);
ReturnFn_negU_floor=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,kappa,floorval) EZReturnFn_negativeUtils_nod_z_floor(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension,kappa,floorval);

fprintf('\n===== Epstein-Zin feasibility conventions: -Inf states (A) and zero returns (B) ===== \n')

%% A: states with no feasible choice
n_neg=12;
a_gridNBL=[(-n_neg:1:-1)'; a_grid];
n_aNBL=length(a_gridNBL);
fprintf('A: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid (n_a=%i) \n',n_neg,-n_neg,n_aNBL)
legnames={'lowmemory=1','divide-and-conquer','grid interpolation','DC+GI'};
Vdead=cell(3,1);
for ezcase=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    if ezcase==1
        casestr='cons-units';
        vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        Paramsref=Params; Paramsref.ezsigma=Params.ezgamma; ReturnFnref=ReturnFn_negU;
    elseif ezcase==2
        casestr='positive utils';
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=1; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        Paramsref=Params; ReturnFnref=ReturnFn_posU;
    else
        casestr='negative utils';
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=0; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        Paramsref=Params; ReturnFnref=ReturnFn_negU;
    end
    % the vNM reference dead set (same budget)
    vfN=struct();
    [Vn,~]=ValueFnIter_Case1_FHorz(0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFnref,Paramsref,DF,[],vfN);
    deadN=gather(Vn==-Inf);
    fprintf('A [EZ %s] vNM reference census: finite %i, -Inf %i, +Inf %i, NaN %i of %i \n',casestr,gather(sum(isfinite(Vn(:)))),sum(deadN(:)),gather(sum(Vn(:)==Inf)),gather(sum(isnan(Vn(:)))),numel(Vn))
    fprintf('A [EZ %s] vNM reference dead set is non-empty, count should be WELL ABOVE zero: %i \n',casestr,sum(deadN(:)))
    % the Epstein-Zin solve, plain path
    ok=0; msg='';
    try
        [V1,P1]=ValueFnIter_Case1_FHorz(0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,[],vfEZ);
        ok=1;
    catch ME
        msg=ME.message;
    end
    fprintf('A [EZ %s] plain solve ran, this should be one: %d \n',casestr,ok)
    if ok~=1
        fprintf('   error: %s \n',msg)
        Vdead{ezcase}=deadN; % so the cross-case line below still runs
        continue
    end
    V1=gather(V1);
    fin=isfinite(V1);
    fprintf('A [EZ %s] census: finite %i, -Inf %i, +Inf %i, NaN %i of %i \n',casestr,sum(fin(:)),sum(V1(:)==-Inf),sum(V1(:)==Inf),sum(isnan(V1(:))),numel(V1))
    fprintf('A [EZ %s] finite set is non-empty (else the max-over-finite lines below are vacuous), this should be one: %d \n',casestr,sum(fin(:))>0)
    fprintf('A [EZ %s] V has no +Inf, this should be one: %d \n',casestr,sum(V1(:)==Inf)==0)
    fprintf('A [EZ %s] V has no NaN, this should be one: %d \n',casestr,sum(isnan(V1(:)))==0)
    fprintf('A [EZ %s] -Inf set equals the vNM dead set, this should be one: %d \n',casestr,isequal(V1==-Inf,deadN))
    fprintf('A [EZ %s] states dead in exactly one of the two, DIAGNOSTIC: %i \n',casestr,sum(xor(V1(:)==-Inf,deadN(:))))
    Vdead{ezcase}=(V1==-Inf);
    Vfp=gather(ValueFnFromPolicy_FHorz(P1,0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,vfEZ));
    fprintf('A [EZ %s] ValueFnFromPolicy on the finite set, this should be zero: %.3e \n',casestr,max(abs(Vfp(fin)-V1(fin))))
    fprintf('A [EZ %s] ValueFnFromPolicy -Inf set equals the solver''s, this should be one: %d \n',casestr,isequal(Vfp==-Inf,V1==-Inf))
    % the other legs: same census, same dead set, and (off the GI legs) the same V
    VGI=[];
    for lg=1:4
        vfL=vfEZ;
        if lg==1
            vfL.lowmemory=1;
        elseif lg==2
            vfL.divideandconquer=1;
        elseif lg==3
            vfL.gridinterplayer=1; vfL.ngridinterp=5;
        else
            vfL.divideandconquer=1; vfL.gridinterplayer=1; vfL.ngridinterp=5;
        end
        ok=0; msg='';
        try
            [VL,PL]=ValueFnIter_Case1_FHorz(0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,[],vfL);
            ok=1;
        catch ME
            msg=ME.message;
        end
        fprintf('A [EZ %s] %s leg ran, this should be one: %d \n',casestr,legnames{lg},ok)
        if ok~=1
            fprintf('   error: %s \n',msg)
            continue
        end
        VL=gather(VL);
        fprintf('A [EZ %s] %s leg: V has no +Inf and no NaN, this should be one: %d \n',casestr,legnames{lg},(sum(VL(:)==Inf)+sum(isnan(VL(:))))==0)
        fprintf('A [EZ %s] %s leg: -Inf set equals the plain leg''s, this should be one: %d \n',casestr,legnames{lg},isequal(VL==-Inf,V1==-Inf))
        if lg<=2
            fprintf('A [EZ %s] %s leg: V on the finite set vs the plain leg, this should be zero: %.3e \n',casestr,legnames{lg},max(abs(VL(fin)-V1(fin))))
        else
            VfpL=gather(ValueFnFromPolicy_FHorz(PL,0,n_aNBL,n_z,N_j,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,vfL));
            finL=isfinite(VL);
            fprintf('A [EZ %s] %s leg: ValueFnFromPolicy on the finite set, this should be zero: %.3e \n',casestr,legnames{lg},max(abs(VfpL(finL)-VL(finL))))
            if lg==3
                VGI=VL;
            else
                fprintf('A [EZ %s] DC+GI vs GI, V on the finite set, this should be zero: %.3e \n',casestr,max(abs(VL(finL)-VGI(finL))))
            end
        end
    end
end
fprintf('A: dead set identical across the three Epstein-Zin cases (same budget), this should be one: %d \n',isequal(Vdead{1},Vdead{2}) && isequal(Vdead{1},Vdead{3}))

%% B: a return of exactly zero at a feasible choice
Pf=Params;
Pf.xbar=0.5;  % subsistence consumption: at a=0 with kappa_j=0.5 the whole wage is below it, so the floor is FORCED there
Pf.kappa=0.5; % negative utils: u(c)+kappa reaches zero at c=1/kappa=2, so the satiation floor is hit at high assets
epsval=10^(-10);
fprintf('B: xbar=%g (cons-units and positive utils), kappa=%g (negative utils), floorval 0 vs %g \n',Pf.xbar,Pf.kappa,epsval)
for ezcase=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    if ezcase==1
        casestr='cons-units';
        vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
        ReturnFn=ReturnFn_cons_floor;
        epssigned=epsval;
    elseif ezcase==2
        casestr='positive utils';
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=1; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU_floor;
        epssigned=epsval;
    else
        casestr='negative utils';
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=0; vfEZ.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU_floor;
        epssigned=-epsval;
    end
    P0=Pf; P0.floorval=0;
    Pe=Pf; Pe.floorval=epssigned;
    oke=0; msg='';
    try
        [Ve,Pe_pol]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn,Pe,DF,[],vfEZ);
        oke=1;
    catch ME
        msg=ME.message;
    end
    fprintf('B [EZ %s] floorval=%g solve ran, this should be one: %d \n',casestr,epssigned,oke)
    if oke~=1
        fprintf('   error: %s \n',msg)
    end
    ok0=0; msg='';
    try
        [V0,P0_pol]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn,P0,DF,[],vfEZ);
        ok0=1;
    catch ME
        msg=ME.message;
    end
    fprintf('B [EZ %s] floorval=0 solve ran, this should be one: %d \n',casestr,ok0)
    if ok0~=1
        fprintf('   error: %s \n',msg)
    end
    if oke~=1 || ok0~=1
        continue
    end
    Ve=gather(Ve); V0=gather(V0);
    % anti-vacuity: the floor must actually be chosen somewhere in the 1e-10 solve
    PolicyValse=gather(PolicyInd2Val_FHorz(Pe_pol,0,n_a,n_z,N_j,[],a_grid,vfEZ));
    aprimevals=reshape(PolicyValse,[n_a,N_z,N_j]);
    cvals=zeros(n_a,N_z,N_j);
    for jj=1:N_j
        if Params.agej(jj)<Params.Jr
            cvals(:,:,jj)=(1+Params.r)*a_grid+Params.w*Params.kappa_j(jj)*z_grid'-aprimevals(:,:,jj);
        else
            cvals(:,:,jj)=(1+Params.r)*a_grid+Params.pension-aprimevals(:,:,jj);
        end
    end
    if ezcase<=2
        floorhit=(cvals>0)&(cvals<Pf.xbar);
    else
        floorhit=(cvals>0)&((cvals.^(1-Params.ezsigma))/(1-Params.ezsigma)+Pf.kappa>=0);
    end
    fprintf('B [EZ %s] share of (a,z,j) whose chosen consumption sits on the floor (floorval=%g solve), should be WELL ABOVE zero: %2.4f \n',casestr,epssigned,mean(floorhit(:)))
    fprintf('B [EZ %s] floorval=0 census: finite %i, -Inf %i, +Inf %i, NaN %i, exactly 0: %i of %i \n',casestr,sum(isfinite(V0(:))),sum(V0(:)==-Inf),sum(V0(:)==Inf),sum(isnan(V0(:))),sum(V0(:)==0),numel(V0))
    fprintf('B [EZ %s] floorval=0 solve: V has no +Inf and no NaN, this should be one: %d \n',casestr,(sum(V0(:)==Inf)+sum(isnan(V0(:))))==0)
    fprintf('B [EZ %s] -Inf set of the floorval=0 solve equals that of the floorval=%g solve, this should be one: %d \n',casestr,epssigned,isequal(V0==-Inf,Ve==-Inf))
    fprintf('B [EZ %s] states dead in exactly one of the two, DIAGNOSTIC: %i \n',casestr,sum(xor(V0(:)==-Inf,Ve(:)==-Inf)))
    fin=isfinite(V0)&isfinite(Ve);
    fprintf('B [EZ %s] common finite set is non-empty, this should be one: %d \n',casestr,sum(fin(:))>0)
    fprintf('B [EZ %s] V(floorval=0) vs V(floorval=%g) on the common finite set, this should be zero: %.3e \n',casestr,epssigned,max(abs(V0(fin)-Ve(fin))))
    fprintf('B [EZ %s] Policy(floorval=0) vs Policy(floorval=%g), raw index difference, DIAGNOSTIC: %.3e at %i of %i \n',casestr,epssigned,max(abs(double(P0_pol(:))-double(Pe_pol(:)))),gather(sum(P0_pol(:)~=Pe_pol(:))),numel(P0_pol))
    Vfp0=gather(ValueFnFromPolicy_FHorz(P0_pol,0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn,P0,DF,vfEZ));
    fin0=isfinite(V0);
    fprintf('B [EZ %s] ValueFnFromPolicy on the floorval=0 solve (finite set), this should be zero: %.3e \n',casestr,max(abs(Vfp0(fin0)-V0(fin0))))
end

%%
output=struct();

end
