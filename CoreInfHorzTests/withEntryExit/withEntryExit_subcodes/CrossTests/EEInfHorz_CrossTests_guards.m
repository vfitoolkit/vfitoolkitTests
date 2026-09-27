function output=EEInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline)
% Interface guards for entry and exit, and probes of the paths this bank does not cover.
%
% Everything is wrapped in try/catch for two reasons: a guard that stops firing must show up as
% a failed check rather than as a silently-accepted bad input, and several probes below are
% expected to fail today, so running them bare would abort the bank.
%
% Probes 6 and 7 are the CPU raws. They are deliberately out of this bank's scope (the bank is
% GPU only), but 6 of the 11 entry-exit raws live there and it is worth having the diary record
% what happens when they are asked for, so that their status is not simply unknown.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,p,alpha,tau,cf,empcap) ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap);
ReturnFn_d=@(d,aprime,a,z,p,alpha,tau,cf,empcap,psi) ReturnFn_EE_d(d,aprime,a,z,p,alpha,tau,cf,empcap,psi);

fprintf('\n===== entry-exit interface guards ===== \n')

%% 0. Positive control: the supported shape solves
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf);
    Vt=gather(Vt);
    ok=any(isfinite(Vt(:))) && ~any(isnan(Vt(:)));
catch ME
    msg=ME.message;
end
fprintf('0. endogenousexit=1 solves with no NaN, this should be one: %d \n',ok)
if ok~=1
    fprintf('   unexpected error or NaN: %s \n',msg)
end

%% 1. endogenousexit=1 without a ReturnToExitFn must error
ok=0; msg='';
try
    vf=struct(); vf.parallel=2; vf.endogenousexit=1; % but no vf.ReturnToExitFn
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('1. endogenousexit=1 without ReturnToExitFn errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 2. endogenousexit=2 without exitprobabilities must error
ok=0; msg='';
try
    vf=vfoptionsbaseline; vf.endogenousexit=2;
    vf.ReturnToExitFn=@(d,aprime,a,z,tau,minexit) ReturnToExitFn_EE_daprime(d,aprime,a,z,tau,minexit);
    % but no vf.exitprobabilities
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('2. endogenousexit=2 without exitprobabilities errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 3. lowmemory=1 must error
% ValueFnIter_InfHorz_EndogExit says so explicitly: 'endogenousexit does allow for
% vfoptions.lowmemory=1' (the message reads as a typo for 'does not allow').
ok=0; msg='';
try
    vf=vfoptionsbaseline; vf.lowmemory=1;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
catch ME
    ok=1; msg=ME.message;
end
fprintf('3. lowmemory=1 errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 4. StationaryDist with agententryandexit=1 but no EntryExitParamNames must error
ok=0; msg='';
try
    vf=vfoptionsbaseline;
    [~,Policyt,ExitPolicyt]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf);
    Paramst=Params; Paramst.zeta=1-gather(ExitPolicyt);
    StationaryDistt=StationaryDist_InfHorz(Policyt,0,n_a,n_z,pi_z,simoptionsbaseline,Paramst); % no EntryExitParamNames
catch ME
    ok=1; msg=ME.message;
end
fprintf('4. entry-exit StationaryDist without EntryExitParamNames errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 5. PROBE: a semi-endogenous shock together with d
% The dispatcher only looks for vfoptions.SemiEndogShock inside its n_d(1)==0 branch, so with d
% present the field is silently IGNORED and the plain path runs. There is no 'with d'
% semi-endogenous raw. This probe records that: if the two answers agree, the field was ignored.
ok=0; msg=''; samesame=NaN;
try
    vf=vfoptionsbaseline;
    vf.SemiEndogShock=(1-0.5)*repmat(shiftdim(gather(pi_z),-1),[N_a,1,1])+0.5*repmat(shiftdim(eye(N_z),-1),[N_a,1,1]);
    Vs=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vf);
    vfp=vfoptionsbaseline;
    Vpl=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfp);
    samesame=max(abs(gather(Vs(:))-gather(Vpl(:))));
    ok=1;
catch ME
    msg=ME.message;
end
% This started as a probe and became a guard. The first run showed it ran happily and returned a
% difference of exactly zero from the plain path, i.e. SemiEndogShock was being silently ignored
% whenever there was a d, because the dispatcher only reads it inside its n_d(1)==0 branch. It
% now errors instead, so this is a 'should be one' like the other guards.
fprintf('5. semi-endogenous shock WITH d errors (it used to be silently ignored), this should be one: %i \n',ok==0)
if ok==1
    fprintf('   RAN INSTEAD OF ERRORING; difference from the plain path: %.3e (zero means it was ignored again) \n',samesame)
else
    fprintf('   message: %s \n',msg)
end

%% 6. PROBE: endogenousexit=1 on the CPU (parallel=0)
% Out of this bank's scope, but 3 of the 11 raws are here and this records their status.
ok=0; msg='';
try
    vf=vfoptionsbaseline; vf.parallel=0;
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
    ok=1;
catch ME
    msg=ME.message;
end
fprintf('6. PROBE endogenousexit=1 with parallel=0 (CPU) ran: %d \n',ok)
if ok~=1
    fprintf('   message: %s \n',msg)
end

%% 7. PROBE: endogenousexit=2 on the CPU (parallel=0)
% Expected to misbehave: the parallel==0 and parallel==1 branches of
% ValueFnIter_InfHorz_EndogExit2 are commented out, so VKron is never assigned and the failure
% is an 'undefined variable' rather than a clean 'not implemented'.
ok=0; msg='';
try
    vf=vfoptionsbaseline; vf.parallel=0; vf.endogenousexit=2;
    vf.ReturnToExitFn=@(d,aprime,a,z,tau,minexit) ReturnToExitFn_EE_daprime(d,aprime,a,z,tau,minexit);
    vf.exitprobabilities={'exitprob_endog','exitprob_exog'};
    vf.endogenousexitcontinuationcost={'continuationcost'};
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf); %#ok<NASGU>
    ok=1;
catch ME
    msg=ME.message;
end
fprintf('7. PROBE endogenousexit=2 with parallel=0 (CPU) ran: %d \n',ok)
if ok~=1
    fprintf('   message: %s \n',msg)
end

%%
output=struct();

end
