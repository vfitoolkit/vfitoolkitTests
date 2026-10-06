function output=CoreInfHorzTPath_localsearch(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Checks vfoptions.localsearch=1, the local-search value function step for InfHorz transition paths.
%
% WHAT THE FEATURE IS. With vfoptions.localsearch=1 the backward value function step maximises over
% aprime only within a window of 2*vfoptions.nlocalsearch+1 grid points around an
% aprimeReferencePolicy, instead of over all N_a points. The window SLIDES at the grid ends (same
% device as divide-and-conquer's loweredge=min(maxindex1,n_a-maxgap)) so it is always exactly 2n+1
% distinct points inside [1,N_a]. It is an approximation unless the window happens to contain the
% unrestricted optimum, and is meant as a cheap warm start for the first few passes of a path solve.
%
% WHAT THIS PART DOES NOT DO. No general eqm. The feature is tested here purely as a value function
% solver, through ValueFnOnTransPath_InfHorz (which evaluates V and Policy along a GIVEN price path)
% and through one direct call to the single-step dispatcher. Working it into the general eqm
% transition path is separate and comes later.
%
% THE METHOD. Four kinds of check, because an approximation cannot be held to an identity in general:
%   EXACT where the restriction is not binding. Two cases, and they are the backbone of this part:
%     - nlocalsearch=(N_a-1)/2 makes the window the whole grid at every state (loweredge collapses to
%       1), so the restricted solve must reproduce the unrestricted path BIT FOR BIT. The return
%       matrix is built by a different command (CreateReturnFnMatrix_Disc_DC1_nod rather than
%       CreateReturnFnMatrix_Disc) and the discounted EV is gathered by linear index rather than
%       broadcast, but both are the same two doubles added in the same order, so exact is the right
%       bar and anything else is a bug.
%     - a reference pinned to the UNRESTRICTED optimum. The window always contains its own centre, so
%       the restricted max is then the unrestricted max at every n. This is the only check that can
%       tell the difference between an aprimeReferencePolicy that is used and one that is ignored:
%       every other check here would pass just as well if the input were silently discarded in favour
%       of the default.
%   MONOTONE. V is non-decreasing in nlocalsearch at every state and every period. Pointwise, not on
%       a norm: the windows are nested in n (shown in the raw's header), the Bellman operator is
%       monotone in the continuation value, and the whole backward sweep starts from the same V_final,
%       so it follows by induction over t. A violation anywhere is a defect.
%   THE WINDOW ITSELF. |aprime-reference| <= 2n everywhere and <= n at interior states. The two bars
%       differ because of the sliding: at a reference within n of either end the window is no longer
%       centred on it. Clamping instead would hold <=n everywhere but would repeat indices and so
%       search fewer than 2n+1 distinct points near the ends.
%   THE ERRORS. Five documented refusals (nlocalsearch unset, nlocalsearch too large for the grid,
%       combination with the grid interpolation layer, combination with divide-and-conquer, and any
%       model other than scalar-a/no-d/z). Each is checked by its message and not merely by the fact
%       that something was thrown, so a refusal cannot pass for the wrong reason.
%
% The approximation error itself is reported but not checked: max|V_n-V_full| and the share of Policy
% entries that match, per n, together with the largest number of grid points the UNRESTRICTED policy
% moves from the current a index (the window the model actually needs). Those numbers are the whole
% case for or against the method, and a bar on them would be invented rather than derived.
%
% SPEED, AND WHY IT NEEDS THE BIG GRID. The whole point of the method is that the step costs
% N_a/(2n+1) less, and at n_a=101 that cannot be seen: the step is kernel-launch bound there and every
% sweep, restricted or not, comes in at the timer's resolution. Section (9) therefore repeats the
% exercise on the bank's n_a_big=1001 grid, where the predicted factor at n=1 is 334 and the
% unrestricted step is arithmetic bound. The measured ratio there is a LOWER bound on the saving: the
% restricted side is still launch bound at 1001 as well, so what is being measured is partly the floor
% and not the arithmetic. The times are measurements and carry no checks - a bar on a runtime is a bar
% on the machine - but the two exactness identities are repeated at 1001 so the index arithmetic is
% confirmed at a size where it is doing real work.

vfoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

% The bank's d grid is kept aside, used only by the with-d refusal check at the end
n_d_withd=n_d; d_grid_withd=d_grid;

n_d=0; d_grid=[]; % local search is currently implemented for no d

N_a=prod(n_a);
N_z=prod(n_z);

if mod(N_a,2)==0
    error('This part needs an odd N_a: the window is 2n+1 points wide, so only an odd grid can be covered exactly by one n')
end
nwide=(N_a-1)/2; % the window is then 1,...,N_a at every state

% V_final and Policy_final: the stationary model at Params, which IS the terminal point of both paths
% (PricePath.r ends at Params.r and ParamPath.w ends at Params.w), so the path is a genuine transition
% rather than a path that has to jump at T.
[V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

%% The unrestricted path, which everything below is measured against
% GPU work is queued asynchronously, so every timing here waits on the device before stopping the
% clock. Without that a restricted sweep 'takes' whatever it costs to ENQUEUE it, which is the one
% number that is guaranteed not to depend on the size of the problem.
gpudev=gpuDevice;
tic;
[VPath0,PolicyPath0]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
wait(gpudev)
time0=toc;
VPath0=gather(VPath0);
Pol0=gather(reshape(PolicyPath0(1,:,:,1:T-1),[N_a,N_z,T-1])); % periods 1,...,T-1 are the ones computed; T is just Policy_final copied in
V0=VPath0(:,:,1:T-1);

aindex=(1:1:N_a)'; % the default reference: one aprime index per state, the current a index, so the window sits around staying put
maxmove=max(abs(Pol0-aindex),[],'all'); % aindex is [N_a,1] and Pol0 is [N_a,N_z,T-1], so this is per state and period
fprintf('localsearch: the unrestricted policy moves at most %i grid points from the current a index on this path, so that is the window this model would need, and the unrestricted sweep took %5.1f seconds \n',maxmove,time0)

%% The ladder in nlocalsearch
nladder=unique([0,1,2,5,10,25,nwide]); % unique also sorts, which the monotonicity check needs; the trim keeps it honest if the bank's n_a ever shrinks
nladder=nladder(nladder<=nwide);
Vgap=zeros(1,length(nladder));
Polsame=zeros(1,length(nladder));
timen=zeros(1,length(nladder));
monotone=1; windowall=1; windowint=1;
Vprev=[];
Vn0=[]; Poln0=[]; Vn1=[]; Vnwide=[]; Polnwide=[];
for n_c=1:length(nladder)
    vfoptionsLS=vfoptions;
    vfoptionsLS.localsearch=1;
    vfoptionsLS.nlocalsearch=nladder(n_c);
    tic;
    [VPathn,PolicyPathn]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLS);
    wait(gpudev)
    timen(n_c)=toc;
    VPathn=gather(VPathn);
    Vn=VPathn(:,:,1:T-1);
    Poln=gather(reshape(PolicyPathn(1,:,:,1:T-1),[N_a,N_z,T-1]));

    Vgap(n_c)=max(abs(Vn-V0),[],'all');
    Polsame(n_c)=mean(Poln(:)==Pol0(:));

    % Monotone in n, pointwise
    if ~isempty(Vprev)
        monotone=monotone*(min(Vn-Vprev,[],'all')>=0);
    end
    Vprev=Vn;

    % The window, against the default reference
    dev=abs(Poln-aindex);
    windowall=windowall*(max(dev,[],'all')<=2*nladder(n_c));
    interior=(aindex-nladder(n_c)>=1 & aindex+nladder(n_c)<=N_a);
    windowint=windowint*(max(dev(interior,:,:),[],'all')<=nladder(n_c)); % interior is non-empty for every n in the ladder, since the ladder stops at (N_a-1)/2, where it is the single middle state

    fprintf('localsearch: nlocalsearch=%3i: max|V-Vfull|=%.3e, share of Policy entries equal to unrestricted %6.2f, sweep took %5.1f seconds \n',nladder(n_c),Vgap(n_c),Polsame(n_c),timen(n_c))

    if nladder(n_c)==0
        Vn0=Vn; Poln0=Poln;
    elseif nladder(n_c)==1
        Vn1=Vn;
    end
    if nladder(n_c)==nwide
        Vnwide=Vn; Polnwide=Poln;
    end
end

%% (1) The whole-grid window must reproduce the unrestricted path exactly
fprintf('localsearch: (1) nlocalsearch=%i covers the grid, V vs the unrestricted path, this should be zero: %.3e \n',nwide,max(abs(Vnwide-V0),[],'all'))
fprintf('localsearch: (1) nlocalsearch=%i covers the grid, Policy vs the unrestricted path, this should be zero: %.3e \n',nwide,max(abs(Polnwide-Pol0),[],'all'))

%% (2) For scale: a narrow window must NOT give the unrestricted answer, or (1) is empty
fprintf('localsearch: (2) for scale, nlocalsearch=1 vs the unrestricted path, V, should NOT be zero: %.3e \n',max(abs(Vn1-V0),[],'all'))
fprintf('localsearch: (2) for scale, nlocalsearch=0 vs the unrestricted path, V, should NOT be zero: %.3e \n',max(abs(Vn0-V0),[],'all'))

%% (3) nlocalsearch=0 with the default reference leaves aprime at the current a index
% The window is then the single point the reference names, so this is what the default reference IS.
fprintf('localsearch: (3) nlocalsearch=0 with the default reference, Policy vs the current a index, this should be zero: %.3e \n',max(abs(Poln0-aindex),[],'all'))

%% (4) V non-decreasing in nlocalsearch, (5)-(6) the window
fprintf('localsearch: (4) V is non-decreasing in nlocalsearch at every state and period, should be one: %i \n',monotone)
fprintf('localsearch: (5) |aprime-reference|<=2*nlocalsearch at every state, should be one: %i \n',windowall)
fprintf('localsearch: (6) |aprime-reference|<=nlocalsearch at every interior state, should be one: %i \n',windowint)

%% (7) The reference is honoured: pin it to the unrestricted optimum and the restriction stops binding
% One single step, called straight on the dispatcher, because ValueFnOnTransPath_InfHorz takes no
% aprimeReferencePolicy input (it builds the default internally). The dispatcher expects grids on the
% GPU and vfoptions already fully filled in, both of which its callers do for it.
ParamsS=Params;
ParamsS.r=PricePath.r(1);
ParamsS.w=ParamPath.w(1);
ReturnFnParamNames=ReturnFnParamNamesFn(ReturnFn,n_d,n_a,n_z,0,vfoptions,ParamsS);
a_gridG=gpuArray(a_grid);
pi_zG=gpuArray(pi_z);
z_gridvals=CreateGridvals(n_z,gpuArray(z_grid),1);
VnextS=gpuArray(reshape(V_final,[N_a,N_z]));

vfoptionsS.exoticpreferences='None';
vfoptionsS.lowmemory=0;
vfoptionsS.gridinterplayer=0;
vfoptionsS.divideandconquer=0;
vfoptionsS.localsearch=0;
[V1,P1]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsS);

% P1 is [1,N_a,N_z], which is exactly the shape aprimeReferencePolicy is documented to take
vfoptionsS.localsearch=1;
for nls=[0,1,5]
    vfoptionsS.nlocalsearch=nls;
    [VLS,PLS]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, P1, vfoptionsS);
    fprintf('localsearch: (7) reference pinned to the unrestricted optimum, nlocalsearch=%i, V, this should be zero: %.3e \n',nls,max(abs(gather(VLS)-gather(V1)),[],'all'))
    fprintf('localsearch: (7) reference pinned to the unrestricted optimum, nlocalsearch=%i, Policy, this should be zero: %.3e \n',nls,max(abs(gather(PLS)-gather(P1)),[],'all'))
end

% And the same single step with the DEFAULT reference, to show (7) is not just two ways of writing the
% unrestricted solve: the restriction does bind when the reference is the current a index.
vfoptionsS.nlocalsearch=1;
[VLSd,~]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsS);
fprintf('localsearch: (7) for scale, the same step with the default reference instead, V, should NOT be zero: %.3e \n',max(abs(gather(VLSd)-gather(V1)),[],'all'))

%% (8) The five documented refusals, each checked by its own message
% Each is tested at the level that raises it: nlocalsearch unset and nlocalsearch too large for the
% grid are the command's own validation, the three combination refusals are the dispatcher's. Going
% through the command for the last three would not work -- a gridinterplayer Policy_final has three
% channels and a with-d one has two, so the command would throw on the size of its own Policy path
% before the dispatcher was ever reached, and the check would pass on the wrong error.
err=zeros(1,5);

vfoptionsE=vfoptions; vfoptionsE.localsearch=1; % no nlocalsearch
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(1)=contains(ME.message,'must set vfoptions.nlocalsearch');
end

vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=N_a; % too large for the grid
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(2)=contains(ME.message,'too large for the grid');
end

vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.gridinterplayer=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
catch ME
    err(3)=contains(ME.message,'combined with vfoptions.gridinterplayer');
end

vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.divideandconquer=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
catch ME
    err(4)=contains(ME.message,'combined with vfoptions.divideandconquer');
end

% With d. Nothing after the refusal is reached, so the d grid is all that has to be real here.
vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d_withd,n_a,n_z,gpuArray(d_grid_withd), a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
catch ME
    err(5)=contains(ME.message,'only implemented for one endogenous state');
end

fprintf('localsearch: (8) nlocalsearch unset is refused by name, should be one: %i \n',err(1))
fprintf('localsearch: (8) nlocalsearch too large for the grid is refused by name, should be one: %i \n',err(2))
fprintf('localsearch: (8) combination with the grid interpolation layer is refused by name, should be one: %i \n',err(3))
fprintf('localsearch: (8) combination with divide-and-conquer is refused by name, should be one: %i \n',err(4))
fprintf('localsearch: (8) a model with d is refused by name, should be one: %i \n',err(5))

%% (9) The same thing on the big grid, where speed can actually be measured
% n_a_big=1001, so the predicted saving at nlocalsearch=1 is 1001/3=334. Everything here is the same
% model on a finer asset grid; the price and parameter paths, z and T are unchanged.
% A note on memory, since the last rung is the widest window the grid allows: the raw MATERIALISES its
% index arrays (aprimeindexes, a_grid(aprimeindexes), aprimez), so at that rung it holds about three
% times what the unrestricted step holds to compute the same answer. That is irrelevant at the window
% sizes the method is for, where the arrays are (2n+1)/N_a of a row, and it is why the wide rung is
% the first thing to drop if this part is ever run on a much larger grid.
N_a_big=prod(n_a_big);
if mod(N_a_big,2)==0
    error('This check needs an odd N_a_big, for the same reason as N_a above')
end
nwide_big=(N_a_big-1)/2;

[V_finalB,Policy_finalB]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

% One untimed sweep first, so the timed ones are not paying for compilation and the first allocation
% of a 1001-by-1001-by-5 return matrix.
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
wait(gpudev)

tic;
[VPathB0,PolicyPathB0]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
wait(gpudev)
time0B=toc;
VPathB0=gather(VPathB0);
VB0=VPathB0(:,:,1:T-1);
PolB0=gather(reshape(PolicyPathB0(1,:,:,1:T-1),[N_a_big,N_z,T-1]));
aindexB=(1:1:N_a_big)';
fprintf('localsearch: big grid (n_a=%i), the unrestricted policy moves at most %i grid points from the current a index, and the unrestricted sweep took %6.2f seconds \n',N_a_big,max(abs(PolB0-aindexB),[],'all'),time0B)

nladderB=[1,2,5,10,nwide_big];
for n_c=1:length(nladderB)
    vfoptionsLS=vfoptions;
    vfoptionsLS.localsearch=1;
    vfoptionsLS.nlocalsearch=nladderB(n_c);
    tic;
    [VPathBn,PolicyPathBn]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLS);
    wait(gpudev)
    timeBn=toc;
    VPathBn=gather(VPathBn);
    VBn=VPathBn(:,:,1:T-1);
    PolBn=gather(reshape(PolicyPathBn(1,:,:,1:T-1),[N_a_big,N_z,T-1]));
    fprintf('localsearch: big grid, nlocalsearch=%4i: %6.2f seconds, measured speed-up %7.1fx against a predicted %7.1fx, max|V-Vfull|=%.3e, share of Policy entries equal to unrestricted %6.2f \n',nladderB(n_c),timeBn,time0B/timeBn,N_a_big/(2*nladderB(n_c)+1),max(abs(VBn-VB0),[],'all'),mean(PolBn(:)==PolB0(:)))

    if nladderB(n_c)==nwide_big
        fprintf('localsearch: (9) big grid, nlocalsearch=%i covers the grid, V vs the unrestricted path, this should be zero: %.3e \n',nwide_big,max(abs(VBn-VB0),[],'all'))
        fprintf('localsearch: (9) big grid, nlocalsearch=%i covers the grid, Policy vs the unrestricted path, this should be zero: %.3e \n',nwide_big,max(abs(PolBn-PolB0),[],'all'))
    end
end

output=struct();

end
