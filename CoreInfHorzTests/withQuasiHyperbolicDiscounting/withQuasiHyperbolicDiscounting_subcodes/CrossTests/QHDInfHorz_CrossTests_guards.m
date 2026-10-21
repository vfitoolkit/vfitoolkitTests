function output=QHDInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for the quasi-hyperbolic INTERFACE guards, and probes of the model shapes that
% InfHorz quasi-hyperbolic does not support.
%
% Everything here is wrapped in try/catch, for two reasons. A guard that stops firing must show
% up as a failed check rather than as a silently-accepted bad input; and the two KNOWN ISSUE
% probes at the end are expected to error today, so running them bare would abort the bank.
%
% Guards 1-5 must all pass. Probes 6 and 7 are diagnostics: they record the current behaviour of
% shapes the dispatcher does not support, so that if support is ever added the diary shows it.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_nod=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

fprintf('\n===== quasi-hyperbolic interface guards ===== \n')

%% 0. Positive control: the char form of QHadditionaldiscount works
% The counterpart to guard 2. Without this line, guard 2 passing would be consistent with the
% field being rejected outright.
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0'; % char, the supported form
    [Vt,~,~,~]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf);
    ok=all(isfinite(gather(Vt(:))));
catch ME
    msg=ME.message;
end
fprintf('0. char QHadditionaldiscount solves and gives a finite V, this should be one: %d \n',ok)
if ok~=1
    fprintf('   unexpected error: %s \n',msg)
end

%% 1. QHadditionaldiscount not declared at all must error
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('1. missing QHadditionaldiscount errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 2. The cell form {'beta0'} must error (the interface is char only)
% Before the char-only change this was the ONLY accepted form, and the FHorz quasi-hyperbolic
% banks still declare it. This check is what stops the cell form creeping back in here.
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount={'beta0'};
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('2. cell-form QHadditionaldiscount errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 3. A non-scalar beta0 must error
% Age-dependent beta0 was deliberately removed. Without this guard a 1-by-N vector would flow
% into beta0beta=beta0*beta and broadcast into a plausible-looking but wrong answer, which is
% the one outcome worth spending an error message to prevent.
ok=0; msg='';
try
    Params_vec=Params;
    Params_vec.beta0=[0.9,0.8];
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0';
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params_vec,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('3. non-scalar beta0 errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 4. A misspelled quasi_hyperbolic must error
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='naive'; % lowercase
    vf.QHadditionaldiscount='beta0';
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('4. misspelled quasi_hyperbolic errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
% Note on what this check is really worth: the ROUTER decides Naive vs Sophisticated with
% ~isfield(vfoptions,'quasi_hyperbolic') || strcmp(...,'Naive'), so a misspelling is routed as
% Sophisticated (three outputs) and only then rejected by the dispatcher. The error is correct;
% the number of outputs the router had already committed to is not. Asking for four outputs
% with a misspelling would fail on the output count, not on the spelling.

%% 5. Grid interpolation layer with two endogenous states must error
% The quasi-hyperbolic grid interpolation path is implemented for scalar n_a only, and says so.
n_a2=3;
a2_grid=[1;2;3]; % values are irrelevant, a2 is ignored by the return function
n_a_2A=[n_a,n_a2];
a_grid_2A=[a_grid; a2_grid];
ReturnFn_2A=@(aprime,a2prime,a,a2,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz_2Aignored(aprime,a2prime,a,a2,z,r,w,sigma);

ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0';
    vf.gridinterplayer=1; vf.ngridinterp=5;
    Vt=ValueFnIter_InfHorz(0,n_a_2A,n_z,[],a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('5. two endogenous states with the grid interpolation layer errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)

%% 5b. Sophisticated + grid interpolation with no maxaprimediff must error
% Added 2026-10-21. maxaprimediff has no default on this path because no single value is right
% across calibrations: the width the search window needs is the quasi-hyperbolic-vs-exponential
% policy gap, which widens as beta0 falls (measured 53 -> 71 for a 0.10 fall at n_a=251), and a
% window that is too narrow pins the choice at the window edge and returns a wrong answer with no
% error. Silent-and-wrong is the outcome worth an error, so the user must choose the value.
% The counterpart positive control is the whole rest of the bank: every Sophisticated + GI solve
% in it sets maxaprimediff explicitly and runs, so this refusal cannot be passing vacuously.
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Sophisticated';
    vf.QHadditionaldiscount='beta0';
    vf.gridinterplayer=1; vf.ngridinterp=5; % deliberately NO maxaprimediff
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('5b. Sophisticated + grid interpolation with no maxaprimediff errors, this should be one: %d \n',ok)
fprintf('   message: %s \n',msg)
% And the Naive counterpart must NOT error: Naive keeps a default, because its quasi-hyperbolic
% step is one extra argmax after V_std has converged rather than part of the iteration.
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0';
    vf.gridinterplayer=1; vf.ngridinterp=5; % deliberately NO maxaprimediff
    Vt=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf);
    ok=all(isfinite(gather(Vt(:))));
catch ME
    msg=ME.message;
end
fprintf('5c. Naive + grid interpolation with no maxaprimediff still solves, this should be one: %d \n',ok)
if ok~=1
    fprintf('   unexpected error: %s \n',msg)
end

%% 6. PROBE (KNOWN ISSUE): a two-dimensional exogenous state
% ValueFnIter_InfHorz.m converts z_grid to z_gridvals at line ~289, and every other solver
% branch is handed z_gridvals. The two exotic-preference branches (QuasiHyperbolic at :356/:359
% and EpsteinZin at :364) are handed the raw stacked z_grid instead. For a single z variable the
% two are the same column vector, so nothing shows; with two z variables z_gridvals is
% prod(n_z)-by-2 while z_grid is sum(n_z)-by-1, and CreateReturnFnMatrix_Disc reads
% z_gridvals(:,2). Expected to error today. The fix is one word at each of the two call sites.
n_z_2d=[3,2];
z1_grid=z_grid(1:3);
z2_grid=[0.9;1.1];
z_grid_2d=[z1_grid;z2_grid];
pi_z1=pi_z(1:3,1:3)./sum(pi_z(1:3,1:3),2);
pi_z2=[0.7,0.3;0.3,0.7];
pi_z_2d=kron(pi_z2,pi_z1); % ordering irrelevant here, the probe never gets as far as using it
ReturnFn_2z=@(aprime,a,z1,z2,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z1*z2,r,w,sigma);

ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0';
    Vt=ValueFnIter_InfHorz(0,n_a,n_z_2d,[],a_grid,z_grid_2d,pi_z_2d,ReturnFn_2z,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('6. PROBE: two-dimensional n_z. Errors today (=1) because the dispatcher is passed z_grid, not z_gridvals: %d \n',ok)
fprintf('   message: %s \n',msg)
fprintf('   If this ever prints 0, multi-dimensional z now works and this probe should become a real test. \n')

%% 7. PROBE (KNOWN ISSUE): two endogenous states WITHOUT the grid interpolation layer
% There is no guard on this path, unlike guard 5. The plain-discretization branch of the
% dispatcher passes n_a straight into ValueFnIter_InfHorz_nod_raw, whose second input is N_a,
% so a vector n_a reaches code that expects a scalar. Recorded here rather than asserted,
% because what it does today is not a specification.
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0';
    Vt=ValueFnIter_InfHorz(0,n_a_2A,n_z,[],a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DiscountFactorParamNames,[],vf);
    msg='(no error: it returned something)';
catch ME
    ok=1; msg=ME.message;
end
fprintf('7. PROBE: two endogenous states without the grid interpolation layer. Errored: %d \n',ok)
fprintf('   message: %s \n',msg)
fprintf('   Two endogenous states are not supported on any InfHorz quasi-hyperbolic path. Guard 5 rejects \n')
fprintf('   the grid-interpolation case explicitly; this path has no such guard. \n')

%%
output=struct();

end
