function output=EZInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Cross-tests for the Epstein-Zin INTERFACE guards, the model shapes that InfHorz Epstein-Zin
% does not support, and the 2A probe.
%
% The checks that expect an error are wrapped in try/catch, so that a guard which stops firing
% shows up as a failed check rather than as a silently-accepted bad input. Each also asserts a
% PHRASE of the message, so that a refusal that is replaced by a different failure (a crash
% inside a raw, say) shows up too: see the audit-refusal-checks-as-a-set rule. The phrases, and
% where each lives (all on 2026-10-22):
%   1  'EZriskaversion'                  ValueFnIter_InfHorz_EpsteinZin (and the ValueFnFromPolicy mirror)
%   2  'EZeis'                           ValueFnIter_InfHorz_EpsteinZin
%   3  'less than one' / 'less than zero' / 'equal to one'   ValueFnIter_InfHorz_EpsteinZin (the ceis==1 one is new, copied from FHorz)
%   4  'EZutils must be 0'               ValueFnIter_InfHorz_EpsteinZin (new)
%   5  'EZoneminusbeta must be 0 or 1'   ValueFnIter_InfHorz_EpsteinZin (new)
%   6  'EZmortalityriskaversion'         ValueFnIter_InfHorz_EpsteinZin (new)
%   7  'survivalprobability'             ValueFnIter_InfHorz_EpsteinZin
%   8  'n_e' / 'n_semiz' / 'endogenousvfi'   ValueFnIter_InfHorz (the main dispatcher)
%   9  'exogenous state z'               ValueFnIter_InfHorz_EpsteinZin (new; was V=0 silently)
%   10 'grid interpolation layer'        ValueFnIter_InfHorz_EpsteinZin (new; was silently ignored)
%   10b 'lowmemory'                      ValueFnIter_InfHorz_EpsteinZin (new; lowmemory=2 fell through to an undefined VKron)
%   12 'not a recognized option'         ValueFnIter_InfHorz (the exoticpreferences whitelist)
%   13 'EZriskaversion' / 'grid interpolation layer'   ValueFnFromPolicy_InfHorz_EpsteinZin
%
% Test 11 is the 2A PROBE: two endogenous states have no guard on the Epstein-Zin path, and are
% expected to WORK, because the Epstein-Zin raws only ever use prod(n_a) (unlike the
% quasi-hyperbolic raws, which fail inside repelem on a vector n_a and so are guarded). It is
% written as a real test: an a2-ignored model must reproduce the one-endogenous-state solve.
% It is wrapped in try/catch like everything else here, so that if it does NOT work the bank
% still finishes; a flag then reads 0 and the decision is to turn this test into a guard.
%
% A ONE-OFF block confirming that vfoptions.maxiter and vfoptions.verbose reach the Epstein-Zin
% raws ran once (2026-10-27) and was then deleted; the numbers are in the driver header.

N_a=prod(n_a);
N_z=prod(n_z);
N_d=prod(n_d);

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_cons_d=@(d,aprime,a,z,r,w,varphi) EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi);

vftol=10^(-10);

% The bank's standard cons-units spec, the base for most of the bad inputs below
vfbase=struct();
vfbase.exoticpreferences='EpsteinZin';
vfbase.EZutils=0;
vfbase.EZriskaversion='ezgamma';
vfbase.EZeis='ezphi';
vfbase.tolerance=vftol;

fprintf('\n===== Epstein-Zin interface guards ===== \n')

%% 0. Positive control: the full cons-units spec solves
% Without this line, every refusal below passing would be consistent with Epstein-Zin being
% rejected outright.
ok0=0; msg='';
try
    [V0,Policy0]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vfbase);
    ok0=gather(all(isfinite(V0(:))));
catch ME
    msg=ME.message;
end
fprintf('0. the full cons-units Epstein-Zin spec solves and gives a finite V, this should be one: %d \n',ok0)
if ok0~=1
    fprintf('   unexpected error: %s \n',msg)
end

%% 1. EZriskaversion not declared at all must error
ok=0; msg='';
try
    vf=vfbase; vf=rmfield(vf,'EZriskaversion');
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'EZriskaversion')); msg=ME.message;
end
fprintf('1. missing EZriskaversion errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 2. cons-units (EZutils=0) without EZeis must error
ok=0; msg='';
try
    vf=vfbase; vf=rmfield(vf,'EZeis');
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'EZeis')); msg=ME.message;
end
fprintf('2. EZutils=0 without EZeis errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 3. cons-units parameter restrictions: gamma<1, phi<=0, phi==1
ok=0; msg='';
try
    P=Params; P.ezgamma=0.5;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vfbase); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'less than one')); msg=ME.message;
end
fprintf('3a. cons-units with ezgamma<1 errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
ok=0; msg='';
try
    P=Params; P.ezphi=-0.5;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vfbase); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'less than zero')); msg=ME.message;
end
fprintf('3b. cons-units with ezphi<=0 errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
ok=0; msg='';
try
    P=Params; P.ezphi=1;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vfbase); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'equal to one')); msg=ME.message;
end
fprintf('3c. cons-units with ezphi=1 errors (new 2026-10-22, copied from the FHorz dispatcher), this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 4. EZutils=2 must error (it was never implemented; it used to die on an undefined variable inside a raw)
ok=0; msg='';
try
    vf=vfbase; vf.EZutils=2;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'EZutils must be 0')); msg=ME.message;
end
fprintf('4. EZutils=2 errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 5. EZoneminusbeta=2 must error (the (1-sj*beta) form needs ages; it used to be ignored)
ok=0; msg='';
try
    vf=vfbase; vf.EZoneminusbeta=2;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'EZoneminusbeta must be 0 or 1')); msg=ME.message;
end
fprintf('5. EZoneminusbeta=2 errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 6. EZmortalityriskaversion must error (no survival probabilities in infinite horizon; it used to be ignored)
ok=0; msg='';
try
    vf=vfbase; vf.EZmortalityriskaversion='ezgamma';
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'EZmortalityriskaversion')); msg=ME.message;
end
fprintf('6. EZmortalityriskaversion errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 7. survivalprobability and WarmGlowBequestsFn must error
ok=0; msg='';
try
    vf=vfbase; vf.survivalprobability='sj'; P=Params; P.sj=0.98;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,P,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'survivalprobability')); msg=ME.message;
end
fprintf('7a. survivalprobability errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
ok=0; msg='';
try
    vf=vfbase; vf.WarmGlowBequestsFn=@(aprime) aprime;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'survivalprobability')); msg=ME.message;
end
fprintf('7b. WarmGlowBequestsFn errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 8. iid e, semiz and endogenousvfi are refused by ValueFnIter_InfHorz itself
ok=0; msg='';
try
    vf=vfbase; vf.n_e=3; vf.e_grid=[0.9;1;1.1]; vf.pi_e=[1/3;1/3;1/3];
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'n_e')); msg=ME.message;
end
fprintf('8a. iid e with Epstein-Zin errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
ok=0; msg='';
try
    vf=vfbase; vf.n_semiz=2; vf.semiz_grid=[0;1]; vf.SemiExoStateFn=@(n,nprime,dsemiz) 0.5;
    Vt=ValueFnIter_InfHorz(2,n_a,n_z,[0;1],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'n_semiz')); msg=ME.message;
end
fprintf('8b. semiz with Epstein-Zin errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
ok=0; msg='';
try
    vf=vfbase; vf.endogenousvfi=1;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'endogenousvfi')); msg=ME.message;
end
fprintf('8c. endogenousvfi=1 with Epstein-Zin errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 9. No exogenous state (n_z=0) must error
% Until 2026-10-22 this returned V=0 and Policy=0 without a word: every Epstein-Zin InfHorz raw
% loops over z, and with N_z=0 that loop never ran, currdist was 0 on the first pass, and the
% zeros came back. A silent zero would fail this check (no error is caught), which is the point.
ReturnFn_noz=@(aprime,a,r,w) EZReturnFn_cons_nod_z(aprime,a,1,r,w);
ok=0; msg='';
try
    Vt=ValueFnIter_InfHorz(0,n_a,0,[],a_grid,[],[],ReturnFn_noz,Params,DiscountFactorParamNames,[],vfbase); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'exogenous state z')); msg=ME.message;
end
fprintf('9. n_z=0 with Epstein-Zin errors with the interface message (not a silent V=0), this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 10. The grid interpolation layer must error (it was silently ignored until 2026-10-22)
% When Tier 2 (an Epstein-Zin InfHorz GI implementation) lands, this test is retired and figs
% 3-4 take its place.
ok=0; msg='';
try
    vf=vfbase; vf.gridinterplayer=1; vf.ngridinterp=5;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'grid interpolation layer')); msg=ME.message;
end
fprintf('10. gridinterplayer=1 with Epstein-Zin errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
ok=0; msg='';
try
    vf=vfbase; vf.lowmemory=2;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'lowmemory')); msg=ME.message;
end
fprintf('10b. lowmemory=2 with Epstein-Zin errors with the interface message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 11. The 2A probe: two endogenous states, the second ignored, must reproduce the one-state solve
% n_a=[n_a,3], a2 takes values 1,2,3 that the return fn never reads. V must be identical across
% the a2 slices and equal to the 1A solve; the a1prime policy likewise. The a2prime channel is a
% pure tie (every a2prime gives the same value), so it resolves to the first index and is printed
% as a diagnostic, not graded. ValueFnFromPolicy_InfHorz is run on the 2A Policy too (its
% Epstein-Zin subfn takes a vector n_a the same way the plain one does).
n_a2=3;
a2_grid=[1;2;3];
n_a_2A=[n_a,n_a2];
a_grid_2A=[a_grid; a2_grid];
ReturnFn_cons_2A=@(aprime,a2prime,a,a2,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU_2A=@(aprime,a2prime,a,a2,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU_2A=@(aprime,a2prime,a,a2,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_cons_d_2A=@(d,aprime,a2prime,a,a2,z,r,w,varphi) EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi);
for ezcase=1:3
    vf11=struct();
    vf11.exoticpreferences='EpsteinZin';
    vf11.tolerance=vftol;
    if ezcase==1
        casestr='cons-units';
        vf11.EZutils=0; vf11.EZriskaversion='ezgamma'; vf11.EZeis='ezphi';
        ReturnFn1A=ReturnFn_cons; ReturnFn2A=ReturnFn_cons_2A;
    elseif ezcase==2
        casestr='positive utils';
        vf11.EZutils=1; vf11.EZpositiveutility=1; vf11.EZriskaversion='ezrisk';
        ReturnFn1A=ReturnFn_posU; ReturnFn2A=ReturnFn_posU_2A;
    else
        casestr='negative utils';
        vf11.EZutils=1; vf11.EZpositiveutility=0; vf11.EZriskaversion='ezrisk';
        ReturnFn1A=ReturnFn_negU; ReturnFn2A=ReturnFn_negU_2A;
    end
    [V1A,P1A]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn1A,Params,DiscountFactorParamNames,[],vf11);
    ok=0; msg='';
    try
        [V2A,P2A]=ValueFnIter_InfHorz(0,n_a_2A,n_z,[],a_grid_2A,z_grid,pi_z,ReturnFn2A,Params,DiscountFactorParamNames,[],vf11);
        ok=1;
        fprintf('11. [EZ %s] 2A (a2 ignored) solve ran, this should be one: %d \n',casestr,ok)
        fprintf('11. [EZ %s] size(V) should be [%d %d %d], is [%s] \n',casestr,N_a,n_a2,N_z,num2str(size(V2A)))
        fprintf('11. [EZ %s] size(Policy) should be [2 %d %d %d], is [%s] \n',casestr,N_a,n_a2,N_z,num2str(size(P2A)))
        Vdev=V2A-reshape(V1A,[N_a,1,N_z]); % implicit expansion over the a2 slices
        fprintf('11. [EZ %s] 2A V vs 1A V (every a2 slice), this should be zero: %2.8e \n',casestr,max(abs(Vdev(:))))
        Pdev=reshape(P2A(1,:,:,:),[N_a,n_a2,N_z])-reshape(P1A,[N_a,1,N_z]);
        fprintf('11. [EZ %s] 2A a1prime policy vs 1A policy (every a2 slice), this should be zero: %.3e \n',casestr,max(abs(Pdev(:))))
        fprintf('11. [EZ %s] a2prime channel resolves the tie to the first index, DIAGNOSTIC: share equal to 1 is %2.4f \n',casestr,gather(mean(P2A(2,:)==1)))
        V2Afp=ValueFnFromPolicy_InfHorz(P2A,0,n_a_2A,n_z,[],a_grid_2A,z_grid,pi_z,ReturnFn2A,Params,DiscountFactorParamNames,vf11);
        fprintf('11. [EZ %s] 2A ValueFnFromPolicy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V2Afp(:)-V2A(:))))
    catch ME
        msg=ME.message;
        fprintf('11. [EZ %s] 2A (a2 ignored) solve ran, this should be one: %d \n',casestr,ok)
        fprintf('   error: %s \n',msg)
        fprintf('   (if this is a shape failure inside a raw, the decision is to GUARD 2A on the Epstein-Zin path, as QH does, not to extend the raws) \n')
    end
end
% and with d, cons-units (the with-d raw: Policy rows are d, a1prime, a2prime)
[V1A,P1A]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons_d,Params,DiscountFactorParamNames,[],vfbase);
ok=0; msg='';
try
    [V2A,P2A]=ValueFnIter_InfHorz(n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,pi_z,ReturnFn_cons_d_2A,Params,DiscountFactorParamNames,[],vfbase);
    ok=1;
    fprintf('11. [EZ cons-units, with d] 2A (a2 ignored) solve ran, this should be one: %d \n',ok)
    fprintf('11. [EZ cons-units, with d] size(Policy) should be [3 %d %d %d], is [%s] \n',N_a,n_a2,N_z,num2str(size(P2A)))
    Vdev=V2A-reshape(V1A,[N_a,1,N_z]);
    fprintf('11. [EZ cons-units, with d] 2A V vs 1A V (every a2 slice), this should be zero: %2.8e \n',max(abs(Vdev(:))))
    Pdev=reshape(P2A(1:2,:,:,:),[2,N_a,n_a2,N_z])-reshape(P1A,[2,N_a,1,N_z]);
    fprintf('11. [EZ cons-units, with d] 2A (d,a1prime) policy vs 1A policy (every a2 slice), this should be zero: %.3e \n',max(abs(Pdev(:))))
    fprintf('11. [EZ cons-units, with d] d index in range, this should be one: %d \n',gather(all(all(P2A(1,:)>=1)) && all(all(P2A(1,:)<=N_d))))
catch ME
    msg=ME.message;
    fprintf('11. [EZ cons-units, with d] 2A (a2 ignored) solve ran, this should be one: %d \n',ok)
    fprintf('   error: %s \n',msg)
end

%% 12. A misspelt exoticpreferences must hit the whitelist, not fall through to a standard solve
ok=0; msg='';
try
    vf=vfbase; vf.exoticpreferences='EpsteinZinn';
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=double(contains(ME.message,'not a recognized option')); msg=ME.message;
end
fprintf('12. misspelt exoticpreferences errors with the whitelist message, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 13. ValueFnFromPolicy_InfHorz under Epstein-Zin: its own guards
if ok0==1
    ok=0; msg='';
    try
        vf=vfbase; vf=rmfield(vf,'EZriskaversion');
        Vt=ValueFnFromPolicy_InfHorz(Policy0,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,vf); %#ok<NASGU>
    catch ME
        ok=double(contains(ME.message,'EZriskaversion')); msg=ME.message;
    end
    fprintf('13a. ValueFnFromPolicy with missing EZriskaversion errors, this should be one: %d \n',ok)
    fprintf('   message: %s \n',msg)
    ok=0; msg='';
    try
        vf=vfbase; vf.gridinterplayer=1; vf.ngridinterp=5;
        Vt=ValueFnFromPolicy_InfHorz(Policy0,0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_cons,Params,DiscountFactorParamNames,vf); %#ok<NASGU>
    catch ME
        ok=double(contains(ME.message,'grid interpolation layer')); msg=ME.message;
    end
    fprintf('13b. ValueFnFromPolicy with gridinterplayer=1 errors, this should be one: %d \n',ok)
    fprintf('   message: %s \n',msg)
else
    fprintf('13. skipped: the positive control did not solve, so there is no Policy to hand to ValueFnFromPolicy, this should be one: 0 \n')
end

%% (A ONE-OFF block confirming that maxiter and verbose reach the Epstein-Zin raws stood here until
% the 2026-10-27 run; its numbers are recorded in the driver header and it was deleted, per the
% driver. Do not re-add it.)

%%
output=struct();

end
