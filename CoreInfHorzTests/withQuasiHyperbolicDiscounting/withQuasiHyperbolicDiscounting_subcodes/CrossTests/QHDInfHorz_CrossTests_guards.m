function output=QHDInfHorz_CrossTests_guards(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for the quasi-hyperbolic INTERFACE guards, and probes of the model shapes that
% InfHorz quasi-hyperbolic does not support.
%
% The checks that expect an error are wrapped in try/catch, so that a guard which stops firing
% shows up as a failed check rather than as a silently-accepted bad input. Test 6 expects a
% SOLVE, so it is deliberately not wrapped: an error there must abort and be seen.
%
% Guards 0-5c and tests 6-7 must all pass. 6 and 7 were KNOWN ISSUE probes recording unsupported
% behaviour until 2026-10-21, when both were fixed: multi-dimensional z now solves (the exotic-
% preference branches were being handed z_grid where every other branch gets z_gridvals), and two
% endogenous states now raise the interface guard on every path rather than only under the grid
% interpolation layer.

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

%% 6. A two-dimensional exogenous state must solve, and agree with the equivalent single z
% This was PROBE 6 until 2026-10-21, and it errored with 'Index exceeds matrix dimension'.
% ValueFnIter_InfHorz.m converts z_grid to z_gridvals and hands z_gridvals to every solver branch,
% but the two exotic-preference branches (QuasiHyperbolic, EpsteinZin) were handed the raw stacked
% z_grid. With one z variable the two forms are the same column, so nothing showed; with two,
% z_gridvals is prod(n_z)-by-2 while z_grid is sum(n_z)-by-1, and the return-matrix construction
% reads column 2. Now fixed at the call sites, so this is a real test.
%
% The check: two independent z variables entering the return function only through their product
% z1*z2 is the SAME model as one z whose grid holds those products in the joint order, with the
% same joint transition matrix. So V and Policy must agree exactly. This is not a comparison of
% identical-by-construction paths: the two-dimensional solve goes through CreateGridvals with
% n_z=[3,2] and a two-z return function, the single-z solve through a one-z return function, and
% if the joint ordering convention were wrong the 2-d V would come back permuted.
n_z_2d=[3,2];
z1_grid=z_grid(1:3);
z2_grid=[0.9;1.1];
z_grid_2d=[z1_grid;z2_grid];
pi_z1=pi_z(1:3,1:3)./sum(pi_z(1:3,1:3),2);
pi_z2=[0.7,0.3;0.3,0.7];
% z1 is the FASTEST-moving index in the joint state, so pi_z2 is the left kron factor. The ordering
% matters now that the solve runs: get it wrong and the two models below are genuinely different.
pi_z_2d=kron(pi_z2,pi_z1);
ReturnFn_2z=@(aprime,a,z1,z2,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z1*z2,r,w,sigma);

vf6=struct(); vf6.exoticpreferences='QuasiHyperbolic'; vf6.quasi_hyperbolic='Naive';
vf6.QHadditionaldiscount='beta0';
[V6,Policy6]=ValueFnIter_InfHorz(0,n_a,n_z_2d,[],a_grid,z_grid_2d,pi_z_2d,ReturnFn_2z,Params,DiscountFactorParamNames,[],vf6);

% The equivalent single z: the products, in the same joint order (z1 fastest).
zprod_grid=kron(z2_grid,ones(3,1)).*kron(ones(2,1),z1_grid);
[V6one,Policy6one]=ValueFnIter_InfHorz(0,n_a,prod(n_z_2d),[],a_grid,zprod_grid,pi_z_2d,ReturnFn_nod,Params,DiscountFactorParamNames,[],vf6);

fprintf('6. two-dimensional n_z vs the equivalent single z, V, this should be zero: %2.8e \n',max(abs(reshape(V6,[N_a*prod(n_z_2d),1])-reshape(V6one,[N_a*prod(n_z_2d),1]))))
fprintf('6. two-dimensional n_z vs the equivalent single z, Policy, this should be zero: %.3e \n',max(abs(double(reshape(Policy6,[],1))-double(reshape(Policy6one,[],1)))))
% Discriminator: the two-dimensional solve must actually depend on z2. If z2 were being dropped
% (which is what reading the wrong grid column could do silently once the shape error is gone),
% V would collapse onto the z2-free model and the line above could still pass for the wrong reason.
[V6flat,~]=ValueFnIter_InfHorz(0,n_a,n_z_2d,[],a_grid,[z1_grid;[1;1]],pi_z_2d,ReturnFn_2z,Params,DiscountFactorParamNames,[],vf6);
fprintf('6. z2 actually reaches the return fn: max|V(z2 spread) - V(z2 degenerate)| should NOT be zero: %2.8f \n',max(abs(reshape(V6,[N_a*prod(n_z_2d),1])-reshape(V6flat,[N_a*prod(n_z_2d),1]))))

%% 7. Two endogenous states must error with the INTERFACE message, on every path
% This was PROBE 7 until 2026-10-21. Two endogenous states were refused cleanly with the grid
% interpolation layer (guard 5) but fell through on the plain and refined paths to fail inside a
% raw with 'In repelem(A,N1,...,Nn), Nj must be a scalar or a vector of length size(A,j)' -- the
% same outcome, an unrecognisable message. The guard now sits with the other interface checks, so
% both entry points raise it, and guard 5 above is the grid-interpolation half of this same check.
% Matching on a phrase, not the whole string: see the audit-refusal-checks-as-a-set rule, and note
% that a refusal lifting must show up as a FAILURE here rather than as silence.
ok=0; msg='';
try
    vf=struct(); vf.exoticpreferences='QuasiHyperbolic'; vf.quasi_hyperbolic='Naive';
    vf.QHadditionaldiscount='beta0';
    Vt=ValueFnIter_InfHorz(0,n_a_2A,n_z,[],a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DiscountFactorParamNames,[],vf);
catch ME
    ok=1; msg=ME.message;
end
fprintf('7. two endogenous states without the grid interpolation layer errors, this should be one: %d \n',ok)
fprintf('7. and it is the interface guard, not a failure inside a raw, this should be one: %d \n',double(contains(msg,'single endogenous state')))
fprintf('   message: %s \n',msg)

%%
output=struct();

end
