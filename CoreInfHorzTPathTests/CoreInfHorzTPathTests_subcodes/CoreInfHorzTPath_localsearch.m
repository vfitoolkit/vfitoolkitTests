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
% WITH THE GRID INTERPOLATION LAYER, sections (10) to (14). There the window is still 2*nlocalsearch+1
% COARSE points around the same coarse reference, with ngridinterp points put between each consecutive
% pair, and the search is one max over all of those fine points. There is no coarse layer at all,
% which is where the saving comes from: the GI raw pays a full N_a coarse max before its second layer,
% and this pays nothing.
% That makes the two solvers DIFFERENT, not two routes to one answer, and the checks are written to
% that. At nlocalsearch=1 centred on the midpoint the GI raw itself would use, the two search exactly
% the same fine points and the answer is an identity -- check (10), which is the backbone. With the
% window covering the whole grid the search is the entire fine grid, which is the exact maximum the
% interpolation layer approximates, so it can only be weakly BETTER than the GI raw -- check (11),
% which also measures how often the GI two-layer search mislocates.
% Policy is compared as the DECODED fine grid point, (1+ngridinterp)*(P1-1)+L2, never as the raw
% channels: a point sitting on a coarse node has two legal spellings and neither raw is canonical
% about which it gives, so a channel comparison would report a difference that is not one.
%
% nlocalsearch=0 is refused everywhere, which is why there is no point-search rung here.
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
%
% AND THE BENCHMARK THAT ACTUALLY MATTERS, section (15). Every speed number here is against the
% UNRESTRICTED sweep, which no one would use on a large model. Divide-and-conquer attacks the same
% coarse max, is already in the toolkit, and is exact under monotonicity of the policy in a without
% needing a reference at all. (15) times it, at several level1n since the toolkit says the default is
% not its best setting, and checks that it gives the same V so the times are comparable. If DC
% already banks the saving then local search has nothing left to win and the reference work should be
% dropped; (15) states the prediction either way before the run decides it.

vfoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

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
nladder=unique([1,2,5,10,25,nwide]); % unique also sorts, which the monotonicity check needs; the trim keeps it honest if the bank's n_a ever shrinks
nladder=nladder(nladder<=nwide);
Vgap=zeros(1,length(nladder));
Polsame=zeros(1,length(nladder));
timen=zeros(1,length(nladder));
monotone=1; windowall=1; windowint=1;
Vprev=[];
Vn1=[]; Vnwide=[]; Polnwide=[];
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

    if nladder(n_c)==1
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

%% The single-step dispatcher, used by (3), (7), (10) and (11)
% ValueFnOnTransPath_InfHorz takes no aprimeReferencePolicy input, so every check that is about the
% reference itself has to go one level down. The dispatcher expects grids on the GPU and vfoptions
% already fully filled in, both of which its callers do for it.
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

%% (3) The default reference IS the current a index
% Passing [] and passing the current a index explicitly must give the same thing, bit for bit. This
% is what makes the default a documented object rather than whatever the dispatcher happens to do,
% and nlocalsearch=0 is no longer available to show it the short way.
vfoptionsD=vfoptionsS;
vfoptionsD.localsearch=1;
vfoptionsD.nlocalsearch=1;
[VD1,PD1]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsD);
[VD2,PD2]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, repmat((1:1:N_a),[1,1,N_z]), vfoptionsD);
fprintf('localsearch: (3) the default reference vs the current a index passed explicitly, V, this should be zero: %.3e \n',max(abs(gather(VD1)-gather(VD2)),[],'all'))
fprintf('localsearch: (3) the default reference vs the current a index passed explicitly, Policy, this should be zero: %.3e \n',max(abs(gather(PD1)-gather(PD2)),[],'all'))

%% (4) V non-decreasing in nlocalsearch, (5)-(6) the window
fprintf('localsearch: (4) V is non-decreasing in nlocalsearch at every state and period, should be one: %i \n',monotone)
fprintf('localsearch: (5) |aprime-reference|<=2*nlocalsearch at every state, should be one: %i \n',windowall)
fprintf('localsearch: (6) |aprime-reference|<=nlocalsearch at every interior state, should be one: %i \n',windowint)

%% (7) The reference is honoured: pin it to the unrestricted optimum and the restriction stops binding
% P1 (above) is [1,N_a,N_z], which is exactly the shape aprimeReferencePolicy is documented to take.
% The window always contains its own centre, so pinning the reference to the unrestricted optimum
% makes the restricted max the unrestricted max at every n.
vfoptionsS.localsearch=1;
for nls=[1,5]
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
% Each is tested at the level that raises it: the three nlocalsearch refusals are the command's own
% validation, the divide-and-conquer one is the dispatcher's.
% There is no with-d refusal to check any more. It was here until d was implemented, and it is now
% doPart(5) that covers that case -- the same thing that happened to the gridinterplayer refusal when
% the interpolation layer was implemented. A refusal check has to be deleted the moment the thing it
% refuses starts working, or it fails for the best possible reason.
err=zeros(1,4);

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

vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=0; % a single grid point, so no search at all
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(3)=contains(ME.message,'must be at least 1');
end

vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.divideandconquer=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
catch ME
    err(4)=contains(ME.message,'combined with vfoptions.divideandconquer');
end

fprintf('localsearch: (8) nlocalsearch unset is refused by name, should be one: %i \n',err(1))
fprintf('localsearch: (8) nlocalsearch too large for the grid is refused by name, should be one: %i \n',err(2))
fprintf('localsearch: (8) nlocalsearch=0 is refused by name, should be one: %i \n',err(3))
fprintf('localsearch: (8) combination with divide-and-conquer is refused by name, should be one: %i \n',err(4))

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

%% (10) Grid interpolation layer: reproducing the GI raw exactly at nlocalsearch=1
% The GI raw takes a full max over all N_a coarse points, turns the winner into a midpoint, and
% searches the 2*ngridinterp+3 fine points one coarse step either side. Its coarse max is the SAME
% max the plain raw does -- same return matrix command, same entireRHS -- so P1 above IS that
% maxindex, and the GI midpoint is P1 clamped to [2,N_a-1]. A window of nlocalsearch=1 around that
% has loweredge=midpoint-1 and therefore spans exactly the fine points the GI raw looks at.
% Same points, same EVinterp, same ReturnFn: an identity, not a comparison. It is the one check that
% validates the fine indexing, the Policy decode and the L2flag derivation all at once.
ngi=5; % the value the rest of this bank uses for ngridinterp
vfoptionsG=vfoptionsS;
vfoptionsG.localsearch=0;
vfoptionsG.gridinterplayer=1;
vfoptionsG.ngridinterp=ngi;
[VG,PG]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsG);

refGI=max(min(P1,N_a-1),2); % the midpoint the GI raw uses
vfoptionsLG=vfoptionsG;
vfoptionsLG.localsearch=1;
vfoptionsLG.nlocalsearch=1;
[VLG,PLG]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refGI, vfoptionsLG);

% The decoded fine grid point is what the two encodings agree on even where the channels do not: the
% LS raw always names a point sitting on a coarse node as (that node, L2=1), while the GI raw may
% name it (previous node, L2=ngridinterp+2). Same aprime, two spellings.
fineG =(1+ngi)*(gather(PG(1,:,:)) -1)+gather(PG(2,:,:));
fineLG=(1+ngi)*(gather(PLG(1,:,:))-1)+gather(PLG(2,:,:));
fprintf('localsearch+GI: (10) nlocalsearch=1 centred where the GI raw centres, V, this should be zero: %.3e \n',max(abs(gather(VLG)-gather(VG)),[],'all'))
fprintf('localsearch+GI: (10) same, the decoded aprime fine grid point, this should be zero: %.3e \n',max(abs(fineLG-fineG),[],'all'))
% HOW TO READ A FAILURE OF THE FLAG LINE. The LS raw tests BOTH endpoints of the cell the optimum
% landed in; the GI raw only ever tests the far one, because its near endpoint is the midpoint and
% the midpoint is its coarse argmax and so feasible. The one case where they can part is when the
% clamp to [2,N_a-1] bit, so the midpoint is NOT the argmax, AND the ReturnFn is -Inf there. That
% cannot happen in this model, and if it ever does the LS raw is the one that is right.
fprintf('localsearch+GI: (10) same, the L2flag channel, this should be zero: %.3e \n',max(abs(gather(PLG(3,:,:))-gather(PG(3,:,:))),[],'all'))
fprintf('localsearch+GI: for information, the raw Policy channels differ at %i of %i states, which is the two spellings of one fine point and not a different aprime \n',sum(gather(PLG(1,:,:))~=gather(PG(1,:,:)),'all'),N_a*N_z)

%% (11) The whole-grid window is the exact fine-grid optimum, so it cannot fall below the GI raw
% With the window covering every coarse point, the search is every point of the fine grid, which is
% the exact maximum the interpolation layer is an approximation to. The GI raw only ever looks one
% coarse step either side of its coarse argmax, so it must be weakly worse. The comparison is exact
% rather than tolerant: where both find the same point the arithmetic is bit-identical.
% How often it is STRICTLY better is a measurement of how often the two-layer search mislocates.
vfoptionsLGw=vfoptionsG;
vfoptionsLGw.localsearch=1;
vfoptionsLGw.nlocalsearch=nwide;
[VLGw,~]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsLGw);
fprintf('localsearch+GI: (11) the whole-grid window is the exact max over the entire fine grid, so it cannot fall below the GI raw, should be one: %i \n',min(gather(VLGw)-gather(VG),[],'all')>=0)
fprintf('localsearch+GI: (11) for information, it is strictly above the GI raw at %i of %i states, which is how often the GI two-layer search mislocates the fine optimum \n',sum(gather(VLGw)>gather(VG),'all'),N_a*N_z)

%% (12) The ladder again, along a path, with the grid interpolation layer on
% Measured against the whole-grid window of the same family rather than against the GI raw, because
% by (11) the GI raw is not the exact object here and the wide window is.
vfoptionsGI=vfoptions;
vfoptionsGI.gridinterplayer=1;
vfoptionsGI.ngridinterp=ngi;
[V_finalGI,Policy_finalGI]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsGI);

vfoptionsGIw=vfoptionsGI;
vfoptionsGIw.localsearch=1;
vfoptionsGIw.nlocalsearch=nwide;
[VPathGIw,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalGI, Policy_finalGI, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsGIw);
VPathGIw=gather(VPathGIw);
VGIw=VPathGIw(:,:,1:T-1);

nladderGI=unique([1,2,5,10,25,nwide]);
nladderGI=nladderGI(nladderGI<=nwide);
monotoneGI=1; windowallGI=1; windowintGI=1;
VprevGI=[];
for n_c=1:length(nladderGI)
    vfoptionsGIn=vfoptionsGI;
    vfoptionsGIn.localsearch=1;
    vfoptionsGIn.nlocalsearch=nladderGI(n_c);
    [VPathGIn,PolicyPathGIn]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalGI, Policy_finalGI, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsGIn);
    VPathGIn=gather(VPathGIn);
    VGIn=VPathGIn(:,:,1:T-1);
    % Channel 1 is the lower coarse grid point, which is where the window constrains the answer
    PolGIn=gather(reshape(PolicyPathGIn(1,:,:,1:T-1),[N_a,N_z,T-1]));

    if ~isempty(VprevGI)
        monotoneGI=monotoneGI*(min(VGIn-VprevGI,[],'all')>=0);
    end
    VprevGI=VGIn;

    devGI=abs(PolGIn-aindex);
    windowallGI=windowallGI*(max(devGI,[],'all')<=2*nladderGI(n_c));
    interiorGI=(aindex-nladderGI(n_c)>=1 & aindex+nladderGI(n_c)<=N_a);
    windowintGI=windowintGI*(max(devGI(interiorGI,:,:),[],'all')<=nladderGI(n_c));

    fprintf('localsearch+GI: nlocalsearch=%3i: %5i fine points searched, max|V-Vwide|=%.3e \n',nladderGI(n_c),2*nladderGI(n_c)*(1+ngi)+1,max(abs(VGIn-VGIw),[],'all'))
end
fprintf('localsearch+GI: (12) V is non-decreasing in nlocalsearch at every state and period, should be one: %i \n',monotoneGI)
fprintf('localsearch+GI: (12) the lower coarse point of aprime is within 2*nlocalsearch of the reference at every state, should be one: %i \n',windowallGI)
fprintf('localsearch+GI: (12) and within nlocalsearch at every interior state, should be one: %i \n',windowintGI)

%% (13) lowmemory is refused with local search, with and without the interpolation layer
errLS=zeros(1,2);
vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.lowmemory=1;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    errLS(1)=contains(ME.message,'lowmemory is not used with vfoptions.localsearch');
end
vfoptionsE=vfoptionsGI; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.lowmemory=1;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalGI, Policy_finalGI, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    errLS(2)=contains(ME.message,'lowmemory is not used with vfoptions.localsearch');
end
fprintf('localsearch: (13) lowmemory with local search is refused by name, should be one: %i \n',errLS(1))
fprintf('localsearch+GI: (13) lowmemory with local search and the interpolation layer is refused by name, should be one: %i \n',errLS(2))

%% (14) Speed on the big grid, with the interpolation layer on
% The unrestricted GI step costs its full coarse layer plus its second layer, about N_a+2*ngi+3
% points per state. The local search costs 2*n*(1+ngi)+1 and no coarse layer at all, so at n=1 and
% ngi=5 the predicted factor is 1014/13, about 78. Expect the measured number to be well below that
% and in fact below the 4.6x the no-GI section reports, for the reason section (9) gives: the
% restricted side is launch bound, and the GI step has strictly more kernel launches than the plain
% one. If it comes back above 4.6x then that reading of the floor is wrong, which is worth knowing.
[V_finalBGI,Policy_finalBGI]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsGI);

ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsGI);
wait(gpudev)
tic;
[VPathB0GI,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsGI);
wait(gpudev)
time0BGI=toc;
VPathB0GI=gather(VPathB0GI);
VB0GI=VPathB0GI(:,:,1:T-1); % (15) uses this, so the sweep is not paid for twice
fprintf('localsearch+GI: big grid (n_a=%i, ngridinterp=%i), the unrestricted GI sweep took %6.2f seconds \n',N_a_big,ngi,time0BGI)

% THE ONE-PASS VERSUS TWO-LAYER COMPARISON, which is what decides which encoding to keep.
% One pass searches 2*n*(1+ngi)+1 fine points, growing at 2*(1+ngi) per unit of nlocalsearch.
% Two layers search (2n+1) coarse plus a FIXED 2*ngi+3 fine, growing at 2. So the point counts cross
% over at n=1+3/(2*ngi), below 2 for any sensible ngi, and the gap widens from there: at n=10 and
% ngi=5 it is 121 points against 34. Against that, two layers pay an extra kernel launch set and an
% exposed reduction latency per step, since layer 2 cannot start until layer 1 retires.
% The ratio of TIMES against the ratio of POINTS is the number to read. If time tracks points then
% the kernel is bandwidth bound and two layers should be kept; if the time ratio stays near 1 while
% the point ratio grows, the step is latency bound at this size and the counts do not matter.
for nls=[1,2,5,10]
    vfoptionsBGI=vfoptionsGI;
    vfoptionsBGI.localsearch=1;
    vfoptionsBGI.nlocalsearch=nls;
    vfoptionsBGI.localsearchtwolayer=0;
    tic;
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsBGI);
    wait(gpudev)
    timeBGI=toc;
    vfoptionsBGI2=vfoptionsBGI;
    vfoptionsBGI2.localsearchtwolayer=1;
    tic;
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsBGI2);
    wait(gpudev)
    timeBGI2=toc;
    nfine1=2*nls*(1+ngi)+1;
    nfine2=(2*nls+1)+(2*ngi+3);
    fprintf('localsearch+GI: big grid, nlocalsearch=%4i: ONE PASS %6.2fs over %5i points (%6.1fx vs unrestricted), TWO LAYERS %6.2fs over %5i points (%6.1fx); two-layer takes %5.2f of the one-pass time for %5.2f of the points \n',nls,timeBGI,nfine1,time0BGI/timeBGI,timeBGI2,nfine2,time0BGI/timeBGI2,timeBGI2/timeBGI,nfine2/nfine1)
end

%% (15) Divide-and-conquer: the competitor that is already in the toolkit
% Every speed number above compares local search against the UNRESTRICTED sweep, which is not the
% benchmark that matters. Divide-and-conquer attacks the same N_a-by-N_a-by-N_z coarse max, is already
% implemented and tested, and is exact under its own assumption -- monotonicity of the policy in a --
% without needing a reference at all. If it already banks the saving, then a local search with a real
% reference policy would be competing for something that is not there to win, and is not worth
% building. That is what this section is for.
%
% level1n is swept rather than left at its default of floor(sqrt(n_a)), because the toolkit itself
% says in ValueFnOnTransPath_InfHorz that the default is not the best setting on a large model, and
% DC should be met at its best.
%
% A PREDICTION, so the run adjudicates it rather than confirming whatever is convenient. The DC raw
% LOOPS over its level1n-1 intervals, serially, inside every one of the T-1 steps. So its kernel
% launch count scales with level1n: at level1n=31 that is about 30*99 launches against a handful for
% the local search. Section (9) showed the unrestricted step sits only about 4x above the launch
% floor on this grid, so DC has very little room, and I expect its best to come at a SMALL level1n
% and to land no better than, and possibly worse than, the local search at nlocalsearch=1. If DC
% instead comes back comfortably faster, the reference work should be dropped.
%
% The two identities keep the comparison honest: DC must give the same V, or it is not solving the
% same problem and the times mean nothing. V and not Policy, because a coarse-level tie can send DC
% and the full search to different optima of equal value, which is a known and accepted artefact of
% this bank and not something to re-litigate inside a timing section.
level1nlist=[15,21,31,51];

%% No grid interpolation layer
vfoptionsLSb=vfoptions;
vfoptionsLSb.localsearch=1;
vfoptionsLSb.nlocalsearch=1;
tic;
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLSb);
wait(gpudev)
timeLSb=toc;
fprintf('DC comparison: big grid, no GI: unrestricted %6.2f s, local search at nlocalsearch=1 %6.2f s \n',time0B,timeLSb)

dcVgap=0; bestDC=Inf; bestlevel1n=0;
for l_c=1:length(level1nlist)
    vfoptionsDC=vfoptions;
    vfoptionsDC.divideandconquer=1;
    vfoptionsDC.level1n=level1nlist(l_c);
    tic;
    [VPathDC,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsDC);
    wait(gpudev)
    timeDC=toc;
    VPathDC=gather(VPathDC);
    dcVgap=max(dcVgap,max(abs(VPathDC(:,:,1:T-1)-VB0),[],'all'));
    if timeDC<bestDC
        bestDC=timeDC; bestlevel1n=level1nlist(l_c);
    end
    fprintf('DC comparison: big grid, no GI: divide-and-conquer at level1n=%3i took %6.2f s, a %5.2fx speed-up on the unrestricted sweep \n',level1nlist(l_c),timeDC,time0B/timeDC)
end
fprintf('DC comparison: big grid, no GI: the best divide-and-conquer is level1n=%i at %6.2f s, and the local search is %5.2fx ITS time (below one means the local search wins) \n',bestlevel1n,bestDC,timeLSb/bestDC)
fprintf('DC comparison: (15) divide-and-conquer V vs the unrestricted path, worst over the level1n swept, this should be zero: %.3e \n',dcVgap)

%% With the grid interpolation layer
vfoptionsLSbGI=vfoptionsGI;
vfoptionsLSbGI.localsearch=1;
vfoptionsLSbGI.nlocalsearch=1;
tic;
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLSbGI);
wait(gpudev)
timeLSbGI=toc;
fprintf('DC comparison: big grid, with GI: unrestricted %6.2f s, local search at nlocalsearch=1 %6.2f s \n',time0BGI,timeLSbGI)

dcVgapGI=0; bestDCGI=Inf; bestlevel1nGI=0;
for l_c=1:length(level1nlist)
    vfoptionsDCGI=vfoptionsGI;
    vfoptionsDCGI.divideandconquer=1;
    vfoptionsDCGI.level1n=level1nlist(l_c);
    tic;
    [VPathDCGI,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsDCGI);
    wait(gpudev)
    timeDCGI=toc;
    VPathDCGI=gather(VPathDCGI);
    dcVgapGI=max(dcVgapGI,max(abs(VPathDCGI(:,:,1:T-1)-VB0GI),[],'all'));
    if timeDCGI<bestDCGI
        bestDCGI=timeDCGI; bestlevel1nGI=level1nlist(l_c);
    end
    fprintf('DC comparison: big grid, with GI: divide-and-conquer at level1n=%3i took %6.2f s, a %5.2fx speed-up on the unrestricted GI sweep \n',level1nlist(l_c),timeDCGI,time0BGI/timeDCGI)
end
fprintf('DC comparison: big grid, with GI: the best divide-and-conquer is level1n=%i at %6.2f s, and the local search is %5.2fx ITS time (below one means the local search wins) \n',bestlevel1nGI,bestDCGI,timeLSbGI/bestDCGI)
fprintf('DC comparison: (15) divide-and-conquer with GI, V vs the unrestricted GI path, worst over the level1n swept, this should be zero: %.3e \n',dcVgapGI)

%% (16) The third output
% LS1_nod: the new reference is just the chosen aprime, so it must equal Policy.
% LS1_GI1_nod: the answer is a fine point and the reference is coarse, so it must be the nearest
% coarse point to it, recomputed here from Policy channels 1 and 2 rather than taken on trust.
vfoptions16=vfoptionsS;
vfoptions16.gridinterplayer=0;
vfoptions16.localsearch=1;
vfoptions16.nlocalsearch=2;
[V16,P16,ref16]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions16);
fprintf('localsearch: (16) the new reference vs Policy, no d and no GI, this should be zero: %.3e \n',max(abs(gather(ref16)-gather(P16)),[],'all'))

vfoptions16GI=vfoptions16;
vfoptions16GI.gridinterplayer=1;
vfoptions16GI.ngridinterp=ngi;
[~,P16GI,ref16GI]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions16GI);
nearestcoarse=gather(P16GI(1,:,:))+round((gather(P16GI(2,:,:))-1)/(1+ngi));
fprintf('localsearch: (16) the new reference vs the nearest coarse point to Policy, no d with GI, this should be zero: %.3e \n',max(abs(gather(ref16GI)-nearestcoarse),[],'all'))

%% (17) Feeding the new reference back in cannot lower V
% The window re-centred on the previous answer still contains it, so V can only improve. This is the
% step of the reference iteration a warm start would do, and the direction is all that is claimed:
% it is NOT an identity, because the re-centred window may contain a better point than the one it
% was centred on.
[V17,~,~]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref16, vfoptions16);
fprintf('localsearch: (17) the new reference fed back in, V cannot fall below the V it came from, should be one: %i \n',min(gather(V17)-gather(V16),[],'all')>=0)
fprintf('localsearch: (17) for information, feeding it back raised V at %i of %i states \n',sum(gather(V17)>gather(V16),'all'),N_a*N_z)

%% (18) The two GI encodings: same window, one exact inside it, the other two-pass
% vfoptions.localsearchtwolayer=1 maxes over the window's COARSE points, clamps that argmax to the
% window's interior, and runs the interpolation layer one coarse cell either side. The clamp makes
% layer 2 reach exactly from coarse loweredge to coarse loweredge+2n, so BOTH ENCODINGS RESTRICT
% aprime TO IDENTICALLY THE SAME FINE SET. They differ only in how they search it, which gives an
% ordering that holds with no assumption at all:
%   one pass is the EXACT max over that set, so it can never fall below two layers;
%   two layers can MISLOCATE -- the fine optimum need not lie in a cell touching the coarse argmax --
%   and how often it does is the accuracy price of the cheaper search.
% That is the same relationship section (11) measured between the whole-grid window and the GI raw,
% and it is measured the same way: a weak-dominance flag plus a count of the strict cases.
%
% AT nlocalsearch=1 THEY MUST AGREE EXACTLY. The clamp forces the midpoint to the window centre, so
% layer 2's 2*ngi+3 points are the one-pass raw's 2*(1+ngi)+1 points IN THE SAME ORDER, making the
% max a bit-identical reduction. Held as an identity, which also pins down that the clamp is right:
% get it wrong by one and the two search sets stop coinciding.
nladder18=unique([1,2,5,10]);
nladder18=nladder18(nladder18<=nwide);
dom18=1; mislocated18=0; total18=0; flagdiff18=0;
V18one_at1=[]; V18two_at1=[]; P18one_at1=[]; P18two_at1=[];
for n_c=1:length(nladder18)
    vfoptions18=vfoptionsGI;
    vfoptions18.localsearch=1;
    vfoptions18.nlocalsearch=nladder18(n_c);
    vfoptions18.localsearchtwolayer=0;
    [VP18one,PP18one]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalGI, Policy_finalGI, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions18);
    vfoptions18.localsearchtwolayer=1;
    [VP18two,PP18two]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalGI, Policy_finalGI, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions18);
    V18one=gather(VP18one); V18one=V18one(:,:,1:T-1);
    V18two=gather(VP18two); V18two=V18two(:,:,1:T-1);
    dom18=dom18*(min(V18one-V18two,[],'all')>=0);
    mislocated18=mislocated18+sum(V18one>V18two,'all');
    total18=total18+numel(V18one);
    % The DECODED fine point, because the two encodings spell it differently: one pass keeps L2 in
    % [1,1+ngi] while two layers inherit the standard GI spelling where L2 reaches ngi+2, so the same
    % fine point has two names. Channels 1 and 2 are the lower coarse point and L2.
    P18one=gather(reshape((1+ngi)*(PP18one(1,:,:,1:T-1)-1)+PP18one(2,:,:,1:T-1),[N_a,N_z,T-1]));
    P18two=gather(reshape((1+ngi)*(PP18two(1,:,:,1:T-1)-1)+PP18two(2,:,:,1:T-1),[N_a,N_z,T-1]));
    % Channel 3 is the L2 flag. The two raws locate the cell endpoints differently (one pass tests
    % both endpoints of the cell the optimum landed in, two layers test the two ends of the
    % interpolation window, exactly as the standard GI raw does), so a difference here is a known
    % structural difference and is reported rather than asserted.
    F18one=gather(reshape(PP18one(3,:,:,1:T-1),[N_a,N_z,T-1]));
    F18two=gather(reshape(PP18two(3,:,:,1:T-1),[N_a,N_z,T-1]));
    flagdiff18=flagdiff18+sum(F18one~=F18two,'all');
    if nladder18(n_c)==1
        V18one_at1=V18one; V18two_at1=V18two; P18one_at1=P18one; P18two_at1=P18two;
    end
end
fprintf('localsearch+GI: (18) one pass is the exact max over the window, so it never falls below two layers, should be one: %i \n',dom18)
fprintf('localsearch+GI: (18) for information, two layers mislocated the fine optimum at %i of %i state-periods, over nlocalsearch=%s \n',mislocated18,total18,mat2str(nladder18))
fprintf('localsearch+GI: (18) for information, the L2 flag differs at %i of %i state-periods (the two raws locate the cell endpoints differently) \n',flagdiff18,total18)
fprintf('localsearch+GI: (18) at nlocalsearch=1 the two encodings search the same points in the same order, V, this should be zero: %.3e \n',max(abs(V18one_at1-V18two_at1),[],'all'))
fprintf('localsearch+GI: (18) at nlocalsearch=1, the decoded fine policy, this should be zero: %.3e \n',max(abs(P18one_at1-P18two_at1),[],'all'))

%% (19) localsearchtwolayer without the interpolation layer is refused by name
% There is only one encoding without the layer, so the option would be silently ignored.
err19=0;
vfoptions19=vfoptions;
vfoptions19.gridinterplayer=0;
vfoptions19.localsearch=1;
vfoptions19.nlocalsearch=2;
vfoptions19.localsearchtwolayer=1;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions19);
catch ME
    err19=contains(ME.message,'only means anything with vfoptions.gridinterplayer=1');
end
fprintf('localsearch: (19) localsearchtwolayer without the interpolation layer is refused by name, should be one: %i \n',err19)

%% (20) THE GI POLICY SPELLING INVARIANT, stressed at the top of the grid
% The agent distribution code always forms the upper grid point as P1+1, so the GI convention is
% P1<=N_a-1 and the top grid point is spelled (N_a-1, L2=ngridinterp+2), never (N_a, L2=1). Nothing
% in the value fn checks can see a violation: they all compare the DECODED fine point, which is the
% same number under either spelling. The failure surfaces only once a Policy reaches the
% distribution, as an out-of-range index inside sparse() in
% AgentDist_InfHorz_TPath_SingleStep_nProbs_raw -- which is where it surfaced, in doPart(6), on the
% first GI solve ever run inside a general eqm path.
% So it is checked directly here, and stressed: the reference is PINNED TO THE TOP OF THE GRID, so
% loweredge sits at N_a-2n for every state and the window's top coarse point IS N_a. Any state whose
% optimum is the top grid point then produces the illegal spelling unless the raw re-spells it.
% Both encodings are checked. The two-layer one is immune by construction -- its coarse argmax is
% clamped to the window interior, so its midpoint is at most N_a-1 -- but checking it is what makes
% that claim a measurement rather than an argument.
refTop=N_a*ones(1,N_a,N_z);
P1max20=0; P1min20=Inf; inwindow20=1;
nladder20=[1,2,5];
nladder20=nladder20(nladder20<=nwide);
for nls=nladder20
    for tl=[0,1]
        vfoptions20=vfoptionsG;
        vfoptions20.localsearch=1;
        vfoptions20.nlocalsearch=nls;
        vfoptions20.localsearchtwolayer=tl;
        [~,P20]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,[], a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refTop, vfoptions20);
        P1=gather(P20(1,:,:)); L2=gather(P20(2,:,:));
        P1max20=max(P1max20,max(P1,[],'all'));
        P1min20=min(P1min20,min(P1,[],'all'));
        % The re-spelling must not move the answer out of the window. loweredge is N_a-2*nls here.
        fine20=(1+ngi)*(P1-1)+L2;
        lo20=(1+ngi)*((N_a-2*nls)-1)+1;
        hi20=(1+ngi)*((N_a-2*nls)+2*nls-1)+1;
        inwindow20=inwindow20*all(fine20>=lo20 & fine20<=hi20,'all');
    end
end
fprintf('localsearch+GI: (20) the coarse channel never reaches N_a, so the upper point P1+1 stays in range, should be one: %i \n',P1max20<=N_a-1)
fprintf('localsearch+GI: (20) the coarse channel never falls below 1, should be one: %i \n',P1min20>=1)
fprintf('localsearch+GI: (20) the decoded fine point still lies inside the window after any re-spelling, should be one: %i \n',inwindow20)
fprintf('localsearch+GI: (20) for information, the largest coarse channel seen was %i against N_a=%i \n',P1max20,N_a)

output=struct();

end
