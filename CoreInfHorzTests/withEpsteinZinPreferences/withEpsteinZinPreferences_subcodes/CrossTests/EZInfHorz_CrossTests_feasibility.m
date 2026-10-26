function output=EZInfHorz_CrossTests_feasibility(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% DELIBERATE TRIGGER for two feasibility conventions of the Epstein-Zin raws: the InfHorz half of
% EZFHorz_CrossTests_feasibility (read that header for the full account). Written 2026-10-27 to run
% RED on the raws as they stood, then green once they are fixed.
%
% A: A STATE WITH NO FEASIBLE CHOICE (V=-Inf). a_grid extended below the natural borrowing limit
%    (integer nodes -12..-1): the lowest states at low z cannot afford any aprime, and the dead set
%    is the closure of those under the dynamics (the InfHorz vNM raws, fixed for -Inf in V on
%    2026-10-22..27, give the reference). Negative utility (ezc3=ezc4=-1) turns a -Inf continuation
%    into +Inf through ezc3*beta*(.), so V holds +Inf entries and the iteration alternates sign
%    rather than converging: the solve is bounded by vfoptions.maxiter=2000 and wrapped in
%    try/catch so that the bank finishes either way.
%    Checks: no +Inf and no NaN in V; the -Inf set equals the vNM dead set; the dead set is the
%    same in all three cases; ValueFnFromPolicy agrees on the finite set and reproduces the dead
%    set; the lowmemory leg gives the same census and the same V.
%
% B: A RETURN OF EXACTLY ZERO AT A FEASIBLE CHOICE (subsistence floor, floorval=0 vs +-1e-10):
%    continuity of V in the return. The zero-return rewrite temp2(ReturnMatrix==0)=-Inf makes the
%    forced-floor states dead in the floorval=0 solve, so the -Inf sets differ: the red line.

N_a=prod(n_a);
N_z=prod(n_z);
DF=DiscountFactorParamNames;

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_cons_floor=@(aprime,a,z,r,w,xbar,floorval) EZReturnFn_cons_nod_z_floor(aprime,a,z,r,w,xbar,floorval);
ReturnFn_posU_floor=@(aprime,a,z,r,w,ezsigma,xbar,floorval) EZReturnFn_positiveUtils_nod_z_floor(aprime,a,z,r,w,ezsigma,xbar,floorval);
ReturnFn_negU_floor=@(aprime,a,z,r,w,ezsigma,kappa,floorval) EZReturnFn_negativeUtils_nod_z_floor(aprime,a,z,r,w,ezsigma,kappa,floorval);

vftol=10^(-10);

fprintf('\n===== Epstein-Zin feasibility conventions: -Inf states (A) and zero returns (B) ===== \n')

%% A: states with no feasible choice
n_neg=12;
a_gridNBL=[(-n_neg:1:-1)'; a_grid];
n_aNBL=length(a_gridNBL);
fprintf('A: a_grid has %i integer-spaced nodes from %i to -1 below the usual cubic grid (n_a=%i) \n',n_neg,-n_neg,n_aNBL)
Vdead=cell(3,1);
for ezcase=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    vfEZ.tolerance=vftol;
    vfEZ.maxiter=2000; % bounds the red run: on the defective negative-utility raws the iteration alternates sign and never converges
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
    vfN=struct(); vfN.tolerance=vftol;
    [Vn,~]=ValueFnIter_InfHorz(0,n_aNBL,n_z,[],a_gridNBL,z_grid,pi_z,ReturnFnref,Paramsref,DF,[],vfN);
    deadN=gather(Vn==-Inf);
    fprintf('A [EZ %s] vNM reference census: finite %i, -Inf %i, +Inf %i, NaN %i of %i \n',casestr,gather(sum(isfinite(Vn(:)))),sum(deadN(:)),gather(sum(Vn(:)==Inf)),gather(sum(isnan(Vn(:)))),numel(Vn))
    fprintf('A [EZ %s] vNM reference dead set is non-empty, count should be WELL ABOVE zero: %i \n',casestr,sum(deadN(:)))
    ok=0; msg='';
    try
        [V1,P1]=ValueFnIter_InfHorz(0,n_aNBL,n_z,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,[],vfEZ);
        ok=1;
    catch ME
        msg=ME.message;
    end
    fprintf('A [EZ %s] plain solve ran, this should be one: %d \n',casestr,ok)
    if ok~=1
        fprintf('   error: %s \n',msg)
        Vdead{ezcase}=deadN;
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
    okfp=0; msg='';
    try
        Vfp=gather(ValueFnFromPolicy_InfHorz(P1,0,n_aNBL,n_z,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,vfEZ));
        okfp=1;
    catch ME
        msg=ME.message;
    end
    fprintf('A [EZ %s] ValueFnFromPolicy ran, this should be one: %d \n',casestr,okfp)
    if okfp==1
        fprintf('A [EZ %s] ValueFnFromPolicy on the finite set, this should be zero to convergence tolerance: %.3e \n',casestr,max(abs(Vfp(fin)-V1(fin))))
        fprintf('A [EZ %s] ValueFnFromPolicy -Inf set equals the solver''s, this should be one: %d \n',casestr,isequal(Vfp==-Inf,V1==-Inf))
    else
        fprintf('   error: %s \n',msg)
    end
    % the lowmemory leg (a different raw)
    vfL=vfEZ; vfL.lowmemory=1;
    ok=0; msg='';
    try
        [VL,~]=ValueFnIter_InfHorz(0,n_aNBL,n_z,[],a_gridNBL,z_grid,pi_z,ReturnFn,Params,DF,[],vfL);
        ok=1;
    catch ME
        msg=ME.message;
    end
    fprintf('A [EZ %s] lowmemory=1 leg ran, this should be one: %d \n',casestr,ok)
    if ok==1
        VL=gather(VL);
        fprintf('A [EZ %s] lowmemory=1 leg: V has no +Inf and no NaN, this should be one: %d \n',casestr,(sum(VL(:)==Inf)+sum(isnan(VL(:))))==0)
        fprintf('A [EZ %s] lowmemory=1 leg: -Inf set equals the plain leg''s, this should be one: %d \n',casestr,isequal(VL==-Inf,V1==-Inf))
        fprintf('A [EZ %s] lowmemory=1 leg: V on the finite set vs the plain leg, this should be zero: %.3e \n',casestr,max(abs(VL(fin)-V1(fin))))
    else
        fprintf('   error: %s \n',msg)
    end
end
fprintf('A: dead set identical across the three Epstein-Zin cases (same budget), this should be one: %d \n',isequal(Vdead{1},Vdead{2}) && isequal(Vdead{1},Vdead{3}))

%% B: a return of exactly zero at a feasible choice
Pf=Params;
Pf.xbar=0.5;  % subsistence consumption: at a=0 the lowest wage w*z_min is below it, so the floor is FORCED there
Pf.kappa=0.5; % negative utils: u(c)+kappa reaches zero at c=1/kappa=2, so the satiation floor is hit at high assets
epsval=10^(-10);
fprintf('B: xbar=%g (cons-units and positive utils), kappa=%g (negative utils), floorval 0 vs %g \n',Pf.xbar,Pf.kappa,epsval)
for ezcase=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    vfEZ.tolerance=vftol;
    vfEZ.maxiter=2000;
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
        [Ve,Pe_pol]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Pe,DF,[],vfEZ);
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
        [V0,P0_pol]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,P0,DF,[],vfEZ);
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
    aprimevals=reshape(gather(PolicyInd2Val_InfHorz(Pe_pol,0,n_a,n_z,[],a_grid,vfEZ)),[N_a,N_z]);
    cvals=(1+Params.r)*a_grid+Params.w*z_grid'-aprimevals;
    if ezcase<=2
        floorhit=(cvals>0)&(cvals<Pf.xbar);
    else
        floorhit=(cvals>0)&((cvals.^(1-Params.ezsigma))/(1-Params.ezsigma)+Pf.kappa>=0);
    end
    fprintf('B [EZ %s] share of (a,z) whose chosen consumption sits on the floor (floorval=%g solve), should be WELL ABOVE zero: %2.4f \n',casestr,epssigned,mean(floorhit(:)))
    fprintf('B [EZ %s] floorval=0 census: finite %i, -Inf %i, +Inf %i, NaN %i, exactly 0: %i of %i \n',casestr,sum(isfinite(V0(:))),sum(V0(:)==-Inf),sum(V0(:)==Inf),sum(isnan(V0(:))),sum(V0(:)==0),numel(V0))
    fprintf('B [EZ %s] floorval=0 solve: V has no +Inf and no NaN, this should be one: %d \n',casestr,(sum(V0(:)==Inf)+sum(isnan(V0(:))))==0)
    fprintf('B [EZ %s] -Inf set of the floorval=0 solve equals that of the floorval=%g solve, this should be one: %d \n',casestr,epssigned,isequal(V0==-Inf,Ve==-Inf))
    fprintf('B [EZ %s] states dead in exactly one of the two, DIAGNOSTIC: %i \n',casestr,sum(xor(V0(:)==-Inf,Ve(:)==-Inf)))
    fin=isfinite(V0)&isfinite(Ve);
    fprintf('B [EZ %s] common finite set is non-empty, this should be one: %d \n',casestr,sum(fin(:))>0)
    fprintf('B [EZ %s] V(floorval=0) vs V(floorval=%g) on the common finite set, this should be zero to convergence tolerance: %.3e \n',casestr,epssigned,max(abs(V0(fin)-Ve(fin))))
    fprintf('B [EZ %s] Policy(floorval=0) vs Policy(floorval=%g), raw index difference, DIAGNOSTIC: %.3e at %i of %i \n',casestr,epssigned,max(abs(double(P0_pol(:))-double(Pe_pol(:)))),gather(sum(P0_pol(:)~=Pe_pol(:))),numel(P0_pol))
    okfp=0; msg='';
    try
        Vfp0=gather(ValueFnFromPolicy_InfHorz(P0_pol,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,P0,DF,vfEZ));
        okfp=1;
    catch ME
        msg=ME.message;
    end
    if okfp==1
        fin0=isfinite(V0);
        fprintf('B [EZ %s] ValueFnFromPolicy on the floorval=0 solve (finite set), this should be zero to convergence tolerance: %.3e \n',casestr,max(abs(Vfp0(fin0)-V0(fin0))))
    else
        fprintf('B [EZ %s] ValueFnFromPolicy on the floorval=0 solve ran, this should be one: 0 \n',casestr)
        fprintf('   error: %s \n',msg)
    end
end

%%
output=struct();

end
