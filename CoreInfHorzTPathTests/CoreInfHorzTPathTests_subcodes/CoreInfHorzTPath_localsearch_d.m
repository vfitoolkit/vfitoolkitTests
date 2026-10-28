function output=CoreInfHorzTPath_localsearch_d(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Checks vfoptions.localsearch=1 on a model WITH a decision variable. The counterpart of
% CoreInfHorzTPath_localsearch, which covers the no-d case; read that header first.
%
% WHAT IS DIFFERENT WITH d. The reference is conditional on d: aprimeReferencePolicy is
% [N_d,N_a,N_z], one aprime index per (d,a,z). It has to be, because d enters the budget and so the
% optimal aprime moves with d -- a reference that did not condition on d would need a window
% spanning the whole range of aprime*(d), the same failure as a reference that did not condition on
% z. Divide-and-conquer with d already does this, its loweredge being n_d-by-1-by-n_z.
% d itself is NOT windowed: it is searched in full, because d is often unordered (an occupation, a
% sector) and nothing can tell which kind of d a user has.
% So the reference is NOT shaped like Policy here. Policy holds (d*,aprime*) and gives aprime* only
% AT the optimal d, which is not what a search keeping d free needs. The no-d case being
% Policy-shaped is a coincidence of l_d=0.
%
% THE THIRD OUTPUT. The raws return aprimeReferencePolicyNew, the aprime chosen for EACH d -- an
% object no other raw in the toolkit produces, since the unrestricted raw maxes over a folded
% (d,aprime) dimension and keeps only the joint optimum. Two things follow, and they are the
% backbone of this part:
%   - at the WIDE window it is the unrestricted per-d aprime optimum. Fed back at nlocalsearch=1 it
%     must reproduce the wide-window answer EXACTLY, because the optimum is in the re-centred window
%     and is the max of it. That is the property a warm start would rely on. Section (6).
%   - it is also exactly the midpoint the GI raw derives, since the GI raw gets its midpoint from a
%     full coarse max over aprime FOR EACH d. So it gives the GI identity for free. Section (8).
%
% THE LINE MOST LIKELY TO BE WRONG is the gather of aprime at the chosen d, so section (5) drives a
% reference that VARIES with d. Under the default reference the reference does not vary with d at
% all and that check could not fail however the gather was written.

vfoptions=struct();

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

if mod(N_a,2)==0
    error('This part needs an odd N_a: the window is 2n+1 points wide, so only an odd grid can be covered exactly by one n')
end
nwide=(N_a-1)/2; % the window is then 1,...,N_a at every (d,a,z)
ngi=5; % the value the rest of this bank uses for ngridinterp

[V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

gpudev=gpuDevice;

%% The unrestricted path, which the path-level checks are measured against
[VPath0,PolicyPath0]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
VPath0=gather(VPath0);
V0=VPath0(:,:,1:T-1);
Pol0d=gather(reshape(PolicyPath0(1,:,:,1:T-1),[N_a,N_z,T-1])); % d
Pol0a=gather(reshape(PolicyPath0(2,:,:,1:T-1),[N_a,N_z,T-1])); % aprime

aindex=(1:1:N_a)';
fprintf('localsearch d: the unrestricted policy moves at most %i grid points from the current a index on this path \n',max(abs(Pol0a-aindex),[],'all'))

%% The ladder in nlocalsearch, along the path
nladder=unique([1,2,5,10,25,nwide]);
nladder=nladder(nladder<=nwide);
monotone=1;
Vprev=[]; Vn1=[]; Vnwide=[]; Polnwided=[]; Polnwidea=[];
for n_c=1:length(nladder)
    vfoptionsLS=vfoptions;
    vfoptionsLS.localsearch=2;
    vfoptionsLS.nlocalsearch=nladder(n_c);
    [VPathn,PolicyPathn]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLS);
    VPathn=gather(VPathn);
    Vn=VPathn(:,:,1:T-1);

    if ~isempty(Vprev)
        monotone=monotone*(min(Vn-Vprev,[],'all')>=0);
    end
    Vprev=Vn;

    fprintf('localsearch d: nlocalsearch=%3i: %5i candidates per state instead of %6i, max|V-Vfull|=%.3e \n',nladder(n_c),N_d*(2*nladder(n_c)+1),N_d*N_a,max(abs(Vn-V0),[],'all'))

    if nladder(n_c)==1
        Vn1=Vn;
    end
    if nladder(n_c)==nwide
        Vnwide=Vn;
        Polnwided=gather(reshape(PolicyPathn(1,:,:,1:T-1),[N_a,N_z,T-1]));
        Polnwidea=gather(reshape(PolicyPathn(2,:,:,1:T-1),[N_a,N_z,T-1]));
    end
end

% NOTE ON vfoptions.localsearch=2 BELOW. ValueFnOnTransPath_InfHorz IGNORES localsearch=1 and warns:
% it produces the value and policy path the welfare numbers and the agent distribution are built
% from, and it has no verification sweep to certify a restricted answer against, so a restricted
% answer there would simply be trusted. localsearch=2 is how a caller forces it anyway, and it is
% what this part has to use -- otherwise every check below would compare an unrestricted sweep
% against an unrestricted sweep and pass while testing nothing.
% Calls that go to the single-step dispatcher, or to TransitionPath_InfHorz, keep localsearch=1;
% those commands REFUSE the value 2 by name, so the two cannot be confused.
%% (1) The whole-grid window must reproduce the unrestricted path
% The channels are compared separately so that a d-only or an aprime-only defect is located rather
% than merely reported. The two maxima this raw takes break exact ties differently from the single
% folded max of the unrestricted raw -- lowest d among the per-d maxima, rather than first in the
% folded order -- so Policy here can in principle differ on a tie while V cannot. Same class as the
% accepted DC tie artefact of this bank.
fprintf('localsearch d: (1) nlocalsearch=%i covers the grid, V vs the unrestricted path, this should be zero: %.3e \n',nwide,max(abs(Vnwide-V0),[],'all'))
fprintf('localsearch d: (1) nlocalsearch=%i covers the grid, Policy d vs the unrestricted path, this should be zero: %.3e \n',nwide,max(abs(Polnwided-Pol0d),[],'all'))
fprintf('localsearch d: (1) nlocalsearch=%i covers the grid, Policy aprime vs the unrestricted path, this should be zero: %.3e \n',nwide,max(abs(Polnwidea-Pol0a),[],'all'))

%% (2) For scale: a narrow window must NOT give the unrestricted answer, or (1) is empty
fprintf('localsearch d: (2) for scale, nlocalsearch=1 vs the unrestricted path, V, should NOT be zero: %.3e \n',max(abs(Vn1-V0),[],'all'))

%% (3) V non-decreasing in nlocalsearch
fprintf('localsearch d: (3) V is non-decreasing in nlocalsearch at every state and period, should be one: %i \n',monotone)

%% The single-step dispatcher, used by (4), (5), (6), (7) and (8)
ParamsS=Params;
ParamsS.r=PricePath.r(1);
ParamsS.w=ParamPath.w(1);
ReturnFnParamNames=ReturnFnParamNamesFn(ReturnFn,n_d,n_a,n_z,0,vfoptions,ParamsS);
a_gridG=gpuArray(a_grid);
pi_zG=gpuArray(pi_z);
z_gridvals=CreateGridvals(n_z,gpuArray(z_grid),1);
if isscalar(n_d)
    d_gridvals=gpuArray(d_grid);
else
    d_gridvals=CreateGridvals(n_d,gpuArray(d_grid),1);
end
VnextS=gpuArray(reshape(V_final,[N_a,N_z]));

vfoptionsS.exoticpreferences='None';
vfoptionsS.lowmemory=0;
vfoptionsS.gridinterplayer=0;
vfoptionsS.divideandconquer=0;
vfoptionsS.localsearch=0;

%% (4) The default reference IS the current a index, for every d
vfoptionsD=vfoptionsS;
vfoptionsD.localsearch=1;
vfoptionsD.nlocalsearch=1;
[VD1,PD1]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsD);
[VD2,PD2]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, repmat((1:1:N_a),[N_d,1,N_z]), vfoptionsD);
fprintf('localsearch d: (4) the default reference vs the current a index passed explicitly for every d, V, this should be zero: %.3e \n',max(abs(gather(VD1)-gather(VD2)),[],'all'))
fprintf('localsearch d: (4) the default reference vs the current a index passed explicitly for every d, Policy, this should be zero: %.3e \n',max(abs(gather(PD1)-gather(PD2)),[],'all'))

%% (5) The window, against a reference that VARIES with d
% ref(d,a,z)=clamp(a+d-ceil(N_d/2),1,N_a). Under the default reference the reference is the same for
% every d, so reading it at the wrong d would give the right answer by accident and this check could
% not fail. With it varying in d, the invariant can only hold if the raw gathers aprime at the d it
% actually chose.
nls5=2;
refd=min(max(repmat((1:1:N_a),[N_d,1,N_z])+((1:1:N_d)'-ceil(N_d/2)),1),N_a); % [N_d,N_a,N_z]
vfoptions5=vfoptionsS;
vfoptions5.localsearch=1;
vfoptions5.nlocalsearch=nls5;
[~,P5]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refd, vfoptions5);
dstar5=reshape(gather(P5(1,:,:)),[N_a,N_z]);
apstar5=reshape(gather(P5(2,:,:)),[N_a,N_z]);
refdg=gather(refd);
refatd=refdg(dstar5+N_d*(0:1:N_a-1)'+N_d*N_a*(0:1:N_z-1)); % the reference at the CHOSEN d
devd=abs(apstar5-refatd);
interiord=(refatd-nls5>=1 & refatd+nls5<=N_a);
fprintf('localsearch d: (5) with a d-varying reference, |aprime-reference at the chosen d|<=2*nlocalsearch everywhere, should be one: %i \n',max(devd,[],'all')<=2*nls5)
fprintf('localsearch d: (5) and <=nlocalsearch wherever that reference is interior, should be one: %i \n',max(devd(interiord),[],'all')<=nls5)

%% (6) The third output at the wide window, fed back at nlocalsearch=1, is exact
% At the wide window the third output is the unrestricted aprime optimum for each d. A window of 1
% around it still contains it, and it is the max of that window, so the answer cannot move. This is
% the property a warm start from a previous pass would rely on, and it is the reason the third
% output exists.
vfoptionsW=vfoptionsS;
vfoptionsW.localsearch=1;
vfoptionsW.nlocalsearch=nwide;
[VW,PW,refW]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsW);
vfoptions6=vfoptionsS;
vfoptions6.localsearch=1;
vfoptions6.nlocalsearch=1;
[V6,P6]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refW, vfoptions6);
fprintf('localsearch d: (6) the wide window new reference fed back at nlocalsearch=1, V, this should be zero: %.3e \n',max(abs(gather(V6)-gather(VW)),[],'all'))
fprintf('localsearch d: (6) the wide window new reference fed back at nlocalsearch=1, Policy, this should be zero: %.3e \n',max(abs(gather(P6)-gather(PW)),[],'all'))

%% (7) The two outputs agree: Policy aprime is the new reference read at the chosen d
dW=reshape(gather(PW(1,:,:)),[N_a,N_z]);
apW=reshape(gather(PW(2,:,:)),[N_a,N_z]);
refWg=gather(refW);
fprintf('localsearch d: (7) Policy aprime vs the new reference read at the chosen d, this should be zero: %.3e \n',max(abs(apW-refWg(dW+N_d*(0:1:N_a-1)'+N_d*N_a*(0:1:N_z-1))),[],'all'))

%% (8) With the grid interpolation layer: reproducing the GI raw exactly at nlocalsearch=1
% The GI raw derives its midpoint from a full coarse max over aprime FOR EACH d, which is exactly
% what refW is. Clamped to [2,N_a-1] it IS the GI raw midpoint.
% At nlocalsearch=1 the localsearch raw's interior clamp, [loweredge+1, loweredge+2n-1], collapses to
% the single point loweredge+1 -- the reference itself -- so handing it refGI PINS its first layer to
% the GI raw's midpoint, per d, and leaves its second layer to do the GI raw's second layer over the
% same 2*ngridinterp+3 fine points. Same points, same EVinterp, same ReturnFn: an identity.
% The decoded fine point is kept as the primary comparison because a point on a coarse node has two
% legal spellings. They should nonetheless agree channel for channel here, since the localsearch raw
% re-spells with the standard GI raw's own lines, and the extra check below holds that outright.
refGI=max(min(refW,N_a-1),2);
vfoptionsG=vfoptionsS;
vfoptionsG.gridinterplayer=1;
vfoptionsG.ngridinterp=ngi;
[VG,PG]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsG);
vfoptionsLG=vfoptionsG;
vfoptionsLG.localsearch=1;
vfoptionsLG.nlocalsearch=1;
[VLG,PLG]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refGI, vfoptionsLG);
fineG =(1+ngi)*(gather(PG(2,:,:)) -1)+gather(PG(3,:,:));
fineLG=(1+ngi)*(gather(PLG(2,:,:))-1)+gather(PLG(3,:,:));
fprintf('localsearch d + GI: (8) nlocalsearch=1 centred where the GI raw centres, V, this should be zero: %.3e \n',max(abs(gather(VLG)-gather(VG)),[],'all'))
fprintf('localsearch d + GI: (8) same, the d channel, this should be zero: %.3e \n',max(abs(gather(PLG(1,:,:))-gather(PG(1,:,:))),[],'all'))
fprintf('localsearch d + GI: (8) same, the decoded aprime fine grid point, this should be zero: %.3e \n',max(abs(fineLG-fineG),[],'all'))
fprintf('localsearch d + GI: (8) same, the L2flag channel, this should be zero: %.3e \n',max(abs(gather(PLG(4,:,:))-gather(PG(4,:,:))),[],'all'))
fprintf('localsearch d + GI: (8) same, the raw coarse channel itself and not just the decoded point, this should be zero: %.3e \n',max(abs(gather(PLG(2,:,:))-gather(PG(2,:,:))),[],'all'))

%% (9) At the WHOLE-GRID window the localsearch GI raw IS the standard GI raw, per d
% nwide=(N_a-1)/2 puts loweredge at 1, so layer 1 maxes over the whole coarse grid for each d --
% which is what the standard GI raw's first layer does -- and the interior clamp becomes [2, N_a-1],
% exactly that raw's own midpoint clamp. So the two must agree bit for bit on every channel, and this
% is the check that pins down the clamp against an independent implementation at the one window where
% it can be.
% Until the one-pass encoding was dropped this section instead compared an EXACT search of the whole
% fine grid against the GI raw, and so measured how often the GI two-layer search mislocates the fine
% optimum. With no exact fine-grid search left in the toolkit that figure is not recoverable.
vfoptionsLGw=vfoptionsG;
vfoptionsLGw.localsearch=1;
vfoptionsLGw.nlocalsearch=nwide;
[VLGw,PLGw]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsLGw);
fineLGw=(1+ngi)*(gather(PLGw(2,:,:))-1)+gather(PLGw(3,:,:));
fprintf('localsearch d + GI: (9) the whole-grid window reproduces the standard GI raw, V, this should be zero: %.3e \n',max(abs(gather(VLGw)-gather(VG)),[],'all'))
fprintf('localsearch d + GI: (9) same, the d channel, this should be zero: %.3e \n',max(abs(gather(PLGw(1,:,:))-gather(PG(1,:,:))),[],'all'))
fprintf('localsearch d + GI: (9) same, the decoded aprime fine grid point, this should be zero: %.3e \n',max(abs(fineLGw-fineG),[],'all'))
fprintf('localsearch d + GI: (9) same, the L2flag channel, this should be zero: %.3e \n',max(abs(gather(PLGw(4,:,:))-gather(PG(4,:,:))),[],'all'))

%% (10) The documented refusals, each checked by its own message
err=zeros(1,6);
vfoptionsE=vfoptions; vfoptionsE.localsearch=2; % no nlocalsearch
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(1)=contains(ME.message,'must set vfoptions.nlocalsearch');
end
vfoptionsE=vfoptions; vfoptionsE.localsearch=2; vfoptionsE.nlocalsearch=N_a;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(2)=contains(ME.message,'too large for the grid');
end
vfoptionsE=vfoptions; vfoptionsE.localsearch=2; vfoptionsE.nlocalsearch=0;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(3)=contains(ME.message,'must be at least 1');
end
% DIVIDE-AND-CONQUER USED TO BE REFUSED HERE, and this used to be the check on that refusal. It is
% now ACCEPTED AND IGNORED on a restricted sweep, silently, exactly as lowmemory is just below -- so
% the check is inverted: the call must SUCCEED. It is kept as a try/catch rather than deleted so that
% a reinstated refusal is REPORTED as a failed check instead of aborting the part, and note that it
% passes no level1n, which only an ignored option permits.
% This is the SIXTH refusal check in this bank to outlive the thing it refused, and the first one
% caught before the run rather than by it: it surfaced while writing a test for the combination,
% which is a different and better instrument than waiting for a refusal to stop firing.
% WHY IT HAD TO GO. The seed sweeps HONOUR divideandconquer -- DC is a memory option here, and with d
% it is the difference between a 2 GB return matrix and 63 MB -- so a user who set it for that reason
% hit the refusal on the first RESTRICTED sweep and could not solve at all, which left the entire DC
% seed route unreachable from the general eqm path it was built for. Section (15) is the identity.
vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.divideandconquer=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
    err(4)=1;
catch
    err(4)=0;
end
% A reference of the no-d shape, [1,N_a,N_z], handed to a model that has a d
vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, repmat((1:1:N_a),[1,1,N_z]), vfoptionsE);
catch ME
    err(6)=contains(ME.message,'one aprime index per state the search conditions on');
end
fprintf('localsearch d: (10) nlocalsearch unset is refused by name, should be one: %i \n',err(1))
fprintf('localsearch d: (10) nlocalsearch too large for the grid is refused by name, should be one: %i \n',err(2))
fprintf('localsearch d: (10) nlocalsearch=0 is refused by name, should be one: %i \n',err(3))
fprintf('localsearch d: (10) divide-and-conquer on a restricted sweep is accepted and ignored, not refused, should be one: %i \n',err(4))
% (10) USED TO CHECK THAT lowmemory WAS REFUSED HERE. It is not any more: the dispatcher zeroes it on
% a restricted sweep, silently and by design, because a window fits in memory whatever the full sweep
% needs, while the full sweeps honour it. The refusal check is replaced by the identity below rather
% than deleted, which is the stronger statement anyway.
% This is the FIFTH refusal check in this bank to outlive the thing it refused, after the with-d one,
% the GI one, the Anderson one and the lowmemory one in doPart(6) (CoreInfHorzTPath_localsearch section (13)). The lesson is not
% "remember to delete it" -- it is that the refusal strings in a bank have to be audited against the
% toolkit's actual error messages as a SET, because fixing the instance that happens to fail leaves
% the others. CoreInfHorzTPath_localsearch section (13) was updated in the same commit that left this one behind.
vfoptions10=vfoptionsS; vfoptions10.localsearch=1; vfoptions10.nlocalsearch=2;
vfoptions10.lowmemory=0;
[V10a,P10a]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions10);
vfoptions10.lowmemory=1;
[V10b,P10b]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions10);
fprintf('localsearch d: (10) lowmemory=1 on a restricted sweep with a d, V vs lowmemory=0, this should be zero: %.3e \n',max(abs(gather(V10a)-gather(V10b)),[],'all'))
fprintf('localsearch d: (10) same, Policy, this should be zero: %.3e \n',max(abs(gather(P10a)-gather(P10b)),[],'all'))
% And the SEED sweep with a d, which is the one raw that had to grow a lowmemory branch of its own:
% LSseed_raw holds the level 1 matrix [N_d,N_a,N_a,N_z], the largest array in the feature, so the
% z-loop is what makes lowmemory mean anything in this corner. Same arithmetic, less resident.
vfoptions10seed=vfoptionsS; vfoptions10seed.localsearch=0; vfoptions10seed.localsearchseed=1;
vfoptions10seed.lowmemory=0;
[V10c,P10c,ref10c]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions10seed);
vfoptions10seed.lowmemory=1;
[V10d,P10d,ref10d]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions10seed);
fprintf('localsearch d: (10) the seed sweep under lowmemory=1, V vs lowmemory=0, this should be zero: %.3e \n',max(abs(gather(V10c)-gather(V10d)),[],'all'))
fprintf('localsearch d: (10) same, Policy, this should be zero: %.3e \n',max(abs(gather(P10c)-gather(P10d)),[],'all'))
fprintf('localsearch d: (10) same, the per-d reference the seed returns, this should be zero: %.3e \n',max(abs(gather(ref10c)-gather(ref10d)),[],'all'))
fprintf('localsearch d: (10) a reference of the no-d shape on a model with d is refused by name, should be one: %i \n',err(6))

%% (11) Divide-and-conquer gives the same V, on the coarse grid where there is no memory risk
% The identity that makes the timings below comparable. It is done at n_a rather than n_a_big because
% with d the unrestricted return matrix is N_d times larger, so the big grid needs about 360MB for it
% and another for entireRHS. Keeping the identity on the coarse grid means it does not ride on that.
dcVgap=0;
for l_c=[7,10,15]
    vfoptionsDC=vfoptions;
    vfoptionsDC.divideandconquer=1;
    vfoptionsDC.level1n=l_c;
    [VPathDC,~]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsDC);
    VPathDC=gather(VPathDC);
    dcVgap=max(dcVgap,max(abs(VPathDC(:,:,1:T-1)-V0),[],'all'));
end
fprintf('localsearch d: (11) divide-and-conquer V vs the unrestricted path, worst over level1n 7, 10 and 15, this should be zero: %.3e \n',dcVgap)

%% (12) Speed on the big grid, with d and no grid interpolation layer. Diagnostics only.
% With d the unrestricted return matrix is [N_d*N_a,N_a,N_z], which at n_a_big=1001 and N_d=9 is
% about 360MB, and entireRHS is another. It fits, measured below; noted because it is the largest
% allocation anywhere in this part.
% THE PREDICTION, recorded before the run decides it. The no-d part measured about 4.8x and argued
% that the ceiling was the kernel launch floor, not the arithmetic: at n_a=1001 and N_z=5 the
% unrestricted step sat only ~4x above that floor. With N_d=9 it sits nine times further above it
% while the restricted step stays near it, so the measured speed-up should be in the TENS. If it
% comes back near 5x again then the floor is not what is limiting it and I have the wrong model.
[V_finalB,Policy_finalB]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
wait(gpudev)
tic;
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
wait(gpudev)
time0B=toc;
fprintf('localsearch d: big grid (n_a=%i, N_d=%i), the unrestricted sweep took %6.2f seconds \n',prod(n_a_big),N_d,time0B)
for nls=[1,2,5,10]
    vfoptionsB=vfoptions;
    vfoptionsB.localsearch=2;
    vfoptionsB.nlocalsearch=nls;
    tic;
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsB);
    wait(gpudev)
    timeB=toc;
    fprintf('localsearch d: big grid, nlocalsearch=%4i: %6.2f seconds, measured speed-up %7.1fx against a predicted %7.1fx \n',nls,timeB,time0B/timeB,prod(n_a_big)/(2*nls+1))
end
for l_c=[7,10,15]
    vfoptionsDCB=vfoptions;
    vfoptionsDCB.divideandconquer=1;
    vfoptionsDCB.level1n=l_c;
    tic;
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalB, Policy_finalB, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsDCB);
    wait(gpudev)
    timeDCB=toc;
    fprintf('localsearch d: big grid, divide-and-conquer at level1n=%3i took %6.2f seconds, a %5.2fx speed-up on the unrestricted sweep \n',l_c,timeDCB,time0B/timeDCB)
end

%% (13) Speed on the big grid, with d AND the grid interpolation layer. Diagnostics only.
% The fourth corner of the table, and the one the earlier runs left blank.
% THE PREDICTION, recorded before the run decided it: the unrestricted GI step pays its full coarse
% layer, [N_d,N_a,N_a,N_z], plus a second layer negligible beside it, so about what the unrestricted
% no-GI step costs -- 4.5 to 5.0 seconds. The restricted step searches 2*n*(1+ngridinterp)+1 fine
% points per d, 13 at n=1 and ngi=5 against 3 without GI, so about 4.3x the arithmetic of the no-GI
% local search but still near the launch floor. On the no-d grid that floor cost GI a factor of 1.43
% (0.10s against 0.07s), so expect roughly 0.20s and a measured speed-up around 23x, against an
% arithmetic ratio of (N_a+2*ngi+3)/(2*(1+ngi)+1)=78x.
% MEASURED 2026-10-07: the unrestricted side came in at 4.94s as predicted, but the restricted side
% took 0.31s not 0.20s, so 16.0x and not 23x. The prediction over-applied the launch floor: on the
% no-d grid the GI restricted step was 65k elements, far enough below the floor that GI cost only
% launches, a factor of 1.43. With d it is 585k and starting to pay arithmetic, so the factor is 2.2.
% The decay in n is the sharper result: 16.0, 11.3, 5.0, 2.6 at n=1,2,5,10 against 31.2, 27.7, 18.2,
% 12.0 without GI. The GI search set grows (1+ngridinterp) times faster in n, so by n=10 the GI local
% search (1.89s) is a tie with divide-and-conquer (2.05s). GI pays only at very small windows.
vfoptionsGI=vfoptions;
vfoptionsGI.gridinterplayer=1;
vfoptionsGI.ngridinterp=ngi;
[V_finalBGI,Policy_finalBGI]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsGI);
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsGI);
wait(gpudev)
tic;
ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsGI);
wait(gpudev)
time0BGI=toc;
fprintf('localsearch d + GI: big grid (n_a=%i, N_d=%i, ngridinterp=%i), the unrestricted GI sweep took %6.2f seconds \n',prod(n_a_big),N_d,ngi,time0BGI)
% The restricted step costs (2n+1) coarse points plus a FIXED 2*ngridinterp+3 fine ones PER d, so it
% grows at 2 points per unit of nlocalsearch rather than at 2*(1+ngridinterp). That is why the
% encoding with a coarse layer was kept: the one-pass version measured here decayed 16.1x, 11.0x,
% 5.0x, 2.6x over nlocalsearch=1,2,5,10 -- by n=10 a tie with divide-and-conquer -- against 15.0x,
% 16.4x, 12.7x, 8.9x for this one. See the header of the LS1_GI1 raw.
for nls=[1,2,5,10]
    vfoptionsBGI=vfoptionsGI;
    vfoptionsBGI.localsearch=2;
    vfoptionsBGI.nlocalsearch=nls;
    tic;
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsBGI);
    wait(gpudev)
    timeBGI=toc;
    npts=(2*nls+1)+(2*ngi+3);
    fprintf('localsearch d + GI: big grid, nlocalsearch=%4i: %6.2f seconds, measured speed-up %7.1fx, over %5i points per d (%4i coarse + %4i fine) instead of %6i \n',nls,timeBGI,time0BGI/timeBGI,npts,2*nls+1,2*ngi+3,prod(n_a_big)+2*ngi+3)
end
for l_c=[7,10,15]
    vfoptionsDCBGI=vfoptionsGI;
    vfoptionsDCBGI.divideandconquer=1;
    vfoptionsDCBGI.level1n=l_c;
    tic;
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalBGI, Policy_finalBGI, Params, n_d,n_a_big,n_z, d_grid,a_grid_big,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsDCBGI);
    wait(gpudev)
    timeDCBGI=toc;
    fprintf('localsearch d + GI: big grid, divide-and-conquer at level1n=%3i took %6.2f seconds, a %5.2fx speed-up on the unrestricted GI sweep \n',l_c,timeDCBGI,time0BGI/timeDCBGI)
end

%% (14) THE GI POLICY SPELLING INVARIANT, stressed at the top of the grid
% The agent distribution code always forms the upper grid point as P1+1, so the GI convention is
% P1<=N_a-1 and the top grid point is spelled (N_a-1, L2=ngridinterp+2), never (N_a, L2=1). Nothing
% in the value fn checks can see a violation: they all compare the DECODED fine point, which is the
% same number under either spelling. The failure surfaces only once a Policy reaches the
% distribution, as an out-of-range index inside sparse() in
% AgentDist_InfHorz_TPath_SingleStep_nProbs_raw -- which is where it surfaced, in doPart(7), on the
% first GI solve ever run inside a general eqm path.
% So it is checked directly here, and stressed: the reference is PINNED TO THE TOP OF THE GRID, so
% loweredge sits at N_a-2n for every state and the window's top coarse point IS N_a. Any state whose
% optimum is the top grid point then produces the illegal spelling unless the raw re-spells it.
% The surviving raw is immune by construction -- its coarse argmax is clamped to the window interior,
% so its midpoint is at most N_a-1 -- but checking it is what makes that claim a measurement rather
% than an argument, and it is the dropped one-pass raw that shipped the violation.
refTop=N_a*ones(N_d,N_a,N_z);
P1max15=0; P1min15=Inf; inwindow15=1;
nladder15=[1,2,5];
nladder15=nladder15(nladder15<=nwide);
for nls=nladder15
    vfoptions15=vfoptionsG;
    vfoptions15.localsearch=1;
    vfoptions15.nlocalsearch=nls;
    [~,P15]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refTop, vfoptions15);
    % Channels are d, lower coarse point, L2, L2flag
    P1=gather(P15(2,:,:)); L2=gather(P15(3,:,:));
    P1max15=max(P1max15,max(P1,[],'all'));
    P1min15=min(P1min15,min(P1,[],'all'));
    fine15=(1+ngi)*(P1-1)+L2;
    lo15=(1+ngi)*((N_a-2*nls)-1)+1;
    hi15=(1+ngi)*((N_a-2*nls)+2*nls-1)+1;
    inwindow15=inwindow15*all(fine15>=lo15 & fine15<=hi15,'all');
end
fprintf('localsearch d + GI: (14) the coarse channel never reaches N_a, so the upper point P1+1 stays in range, should be one: %i \n',P1max15<=N_a-1)
fprintf('localsearch d + GI: (14) the coarse channel never falls below 1, should be one: %i \n',P1min15>=1)
fprintf('localsearch d + GI: (14) the decoded fine point still lies inside the window after any re-spelling, should be one: %i \n',inwindow15)
fprintf('localsearch d + GI: (14) for information, the largest coarse channel seen was %i against N_a=%i \n',P1max15,N_a)

%% (15) DIVIDE-AND-CONQUER WITH LOCAL SEARCH, WITH A d
% The with-d corners of what CoreInfHorzTPath_localsearch section (19) does without one, and the corner that needed real
% arithmetic rather than a returned variable. Two halves, doing opposite things:
%   - A RESTRICTED sweep IGNORES divideandconquer, so the answer is identical bit for bit.
%   - A SEED sweep HONOURS it, and with d and no interpolation layer that meant splitting DC1_raw's
%     level-2 max -- one max over aprime, giving the per-d reference, then one over d, giving Policy
%     -- because a single max over a folded (d,aprime) dimension cannot produce a reference for the
%     N_d-1 values of d it did not choose. Under GI the DC raws already built a per-d midpoint over
%     every a for their own interpolation layer, so there it was one returned variable.
% The seed identity is in V only. DC and the full search break exact ties differently, which moves
% Policy and hence the reference without moving V; the disagreement count is reported instead, and is
% also the only evidence available that DC ran at all, since an exact method leaves no other trace.
% WHAT PINS THE REFERENCE is therefore not an agreement check but an INTERNAL one: read the DC seed's
% own reference at the d its own Policy chose, and it must be that Policy's aprime. That is the check
% on the split, and the one thing that would catch the two maxes being combined in the wrong order.
vfoptions15=vfoptionsS;
vfoptions15.localsearch=1;
vfoptions15.nlocalsearch=2;
[V15,P15d,ref15]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15);
dcVgap15=0; dcPgap15=0; dcRgap15=0;
for l_c=[7,10,15]
    vfoptions15DC=vfoptions15;
    vfoptions15DC.divideandconquer=1;
    vfoptions15DC.level1n=l_c;
    [V15DC,P15DC,ref15DC]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15DC);
    dcVgap15=max(dcVgap15,max(abs(gather(V15DC)-gather(V15)),[],'all'));
    dcPgap15=max(dcPgap15,max(abs(gather(P15DC)-gather(P15d)),[],'all'));
    dcRgap15=max(dcRgap15,max(abs(gather(ref15DC)-gather(ref15)),[],'all'));
end
vfoptions15DC=vfoptions15;
vfoptions15DC.divideandconquer=1; % and with level1n never set at all
[V15DC,P15DC,ref15DC]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15DC);
dcVgap15=max(dcVgap15,max(abs(gather(V15DC)-gather(V15)),[],'all'));
dcPgap15=max(dcPgap15,max(abs(gather(P15DC)-gather(P15d)),[],'all'));
dcRgap15=max(dcRgap15,max(abs(gather(ref15DC)-gather(ref15)),[],'all'));
fprintf('localsearch d: (15) a restricted sweep with d under divideandconquer=1, V vs divideandconquer=0, worst over level1n 7, 10, 15 and unset, this should be zero: %.3e \n',dcVgap15)
fprintf('localsearch d: (15) same, Policy, this should be zero: %.3e \n',dcPgap15)
fprintf('localsearch d: (15) same, the per-d reference, this should be zero: %.3e \n',dcRgap15)

% The same with the interpolation layer. Policy is compared channel for channel, not through the
% decoded fine point: both sides are the SAME raw, so the two legal spellings of a coarse node cannot
% be what differs, and the looser comparison would hide what the check is for.
vfoptions15G=vfoptionsG;
vfoptions15G.localsearch=1;
vfoptions15G.nlocalsearch=2;
[V15G,P15G,ref15G]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15G);
dcVgap15G=0; dcPgap15G=0; dcRgap15G=0;
for l_c=[7,10,15]
    vfoptions15GDC=vfoptions15G;
    vfoptions15GDC.divideandconquer=1;
    vfoptions15GDC.level1n=l_c;
    [V15GDC,P15GDC,ref15GDC]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15GDC);
    dcVgap15G=max(dcVgap15G,max(abs(gather(V15GDC)-gather(V15G)),[],'all'));
    dcPgap15G=max(dcPgap15G,max(abs(gather(P15GDC)-gather(P15G)),[],'all'));
    dcRgap15G=max(dcRgap15G,max(abs(gather(ref15GDC)-gather(ref15G)),[],'all'));
end
fprintf('localsearch d + GI: (15) a restricted GI sweep with d under divideandconquer=1, V vs divideandconquer=0, worst over level1n 7, 10 and 15, this should be zero: %.3e \n',dcVgap15G)
fprintf('localsearch d + GI: (15) same, Policy in every channel, this should be zero: %.3e \n',dcPgap15G)
fprintf('localsearch d + GI: (15) same, the per-d reference, this should be zero: %.3e \n',dcRgap15G)

% THE SEED SWEEP WITH d AND NO GI: DC1_raw's split level-2 max against LSseed_raw, which does the
% same split without DC. This is the pair that section (11) of doPart(7) does for LSseed_raw against
% the folded standard raw; here it is DC1_raw's turn, and the internal reference check is what makes
% it a test of the split rather than of DC.
vfoptions15S=vfoptionsS;
vfoptions15S.localsearch=0;
vfoptions15S.localsearchseed=1;
[V15S,P15S,ref15S]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15S);
seedVgap15=0; seedties15=0; seedref15=0;
for l_c=[7,10,15]
    vfoptions15SDC=vfoptions15S;
    vfoptions15SDC.divideandconquer=1;
    vfoptions15SDC.level1n=l_c;
    [V15SDC,P15SDC,ref15SDC]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15SDC);
    seedVgap15=max(seedVgap15,max(abs(gather(V15SDC)-gather(V15S)),[],'all'));
    seedties15=max(seedties15,sum(gather(P15SDC(1,:,:))~=gather(P15S(1,:,:)),'all'));
    r15=gather(ref15SDC);
    dstar15=gather(P15SDC(1,:,:));
    lin15=dstar15+N_d*(0:1:N_a-1)+N_d*N_a*shiftdim((0:1:N_z-1),-1);
    seedref15=max(seedref15,max(abs(r15(lin15)-gather(P15SDC(2,:,:))),[],'all'));
end
fprintf('localsearch d: (15) the seed sweep with d under divideandconquer=1, V vs divideandconquer=0, worst over level1n 7, 10 and 15, this should be zero: %.3e \n',seedVgap15)
fprintf('localsearch d: (15) the DC seed reference read at the chosen d vs that Policy aprime, worst over level1n 7, 10 and 15, this should be zero: %.3e \n',seedref15)
fprintf('localsearch d: (15) for information, their Policy d channel disagrees at %i of %i states, which is exact ties broken differently and is also the only evidence that DC actually ran \n',seedties15,N_a*N_z)

% AND THE SEED SWEEP WITH d AND GI, where the reference is DC1_GI1_raw's own midpoint per d. The
% clamp relation is what can be checked without re-deriving it: the interpolation layer searches one
% coarse cell either side of the midpoint, so the lower coarse point of the answer is the midpoint or
% one below it at the chosen d, and nothing else. Off by one in the returned midpoint is exactly this.
vfoptions15SG=vfoptionsG;
vfoptions15SG.localsearch=0;
vfoptions15SG.localsearchseed=1;
[V15SG,P15SG,ref15SG]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15SG);
seedVgap15G=0; clamp15G=1; range15G=1;
for l_c=[7,10,15]
    vfoptions15SGDC=vfoptions15SG;
    vfoptions15SGDC.divideandconquer=1;
    vfoptions15SGDC.level1n=l_c;
    [V15SGDC,P15SGDC,ref15SGDC]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15SGDC);
    seedVgap15G=max(seedVgap15G,max(abs(gather(V15SGDC)-gather(V15SG)),[],'all'));
    rG15=gather(ref15SGDC);
    dstarG15=gather(P15SGDC(1,:,:));
    linG15=dstarG15+N_d*(0:1:N_a-1)+N_d*N_a*shiftdim((0:1:N_z-1),-1);
    midG15=rG15(linG15);
    cG15=gather(P15SGDC(2,:,:));
    clamp15G=clamp15G*all(cG15==midG15 | cG15==midG15-1,'all');
    range15G=range15G*all(rG15>=2 & rG15<=N_a-1,'all');
end
fprintf('localsearch d + GI: (15) the GI seed sweep with d under divideandconquer=1, V vs divideandconquer=0, worst over level1n 7, 10 and 15, this should be zero: %.3e \n',seedVgap15G)
fprintf('localsearch d + GI: (15) the DC+GI seed reference at the chosen d is the midpoint of the cell Policy sits in, should be one: %i \n',clamp15G)
fprintf('localsearch d + GI: (15) the DC+GI seed reference lies in [2,N_a-1] for every d, as a midpoint must, should be one: %i \n',range15G)

% THE DC SEED RAWS WITH d UNDER lowmemory=1. The seed sweeps honour lowmemory too -- they are the full
% sweeps, and with d they are the largest arrays anywhere in this feature -- so the per-d reference
% had to be built inside each DC raw's z-loop branch as well, and that branch only runs here. For
% DC1_raw that means the split level-2 max exists twice, once per memory branch, which is exactly the
% shape of the copy-paste drift this bank has found before: identical arithmetic in two places, only
% one of which any check reached. Section (10) covers LSseed_raw's own z-loop branch the same way.
vfoptions15L=vfoptions15S;
vfoptions15L.divideandconquer=1;
vfoptions15L.level1n=10;
vfoptions15L.lowmemory=0;
[V15L0,P15L0,ref15L0]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15L);
vfoptions15L.lowmemory=1;
[V15L1,P15L1,ref15L1]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15L);
fprintf('localsearch d: (15) the DC seed sweep with d under lowmemory=1, V vs lowmemory=0, this should be zero: %.3e \n',max(abs(gather(V15L1)-gather(V15L0)),[],'all'))
fprintf('localsearch d: (15) same, Policy, this should be zero: %.3e \n',max(abs(gather(P15L1)-gather(P15L0)),[],'all'))
fprintf('localsearch d: (15) same, the per-d reference, which is the split max in the z-loop branch, this should be zero: %.3e \n',max(abs(gather(ref15L1)-gather(ref15L0)),[],'all'))

vfoptions15LG=vfoptions15SG;
vfoptions15LG.divideandconquer=1;
vfoptions15LG.level1n=10;
vfoptions15LG.lowmemory=0;
[V15LG0,P15LG0,ref15LG0]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15LG);
vfoptions15LG.lowmemory=1;
[V15LG1,P15LG1,ref15LG1]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions15LG);
fprintf('localsearch d + GI: (15) the DC+GI seed sweep with d under lowmemory=1, V vs lowmemory=0, this should be zero: %.3e \n',max(abs(gather(V15LG1)-gather(V15LG0)),[],'all'))
fprintf('localsearch d + GI: (15) same, Policy in every channel, this should be zero: %.3e \n',max(abs(gather(P15LG1)-gather(P15LG0)),[],'all'))
fprintf('localsearch d + GI: (15) same, the per-d reference, this should be zero: %.3e \n',max(abs(gather(ref15LG1)-gather(ref15LG0)),[],'all'))

% THE DISCRIMINATOR, which the first run of this section showed was missing. Every identity above is
% satisfied whether the DC raw ran or the plain one did -- DC is EXACT, so on a model with no exact
% ties it leaves no trace in V, Policy or the reference, and the tie count came back 0 of 505. So the
% claim "the seed sweep honours DC" rested on reading the dispatcher's fork, not on a measurement.
% There is exactly one observable that separates them: all four DC raws read vfoptions.level1n on
% their first line, and nothing defaults it at this level. So with level1n UNSET a seed sweep under
% divideandconquer=1 must FAIL, and the restricted sweep in the refusal block above must SUCCEED on
% the same options -- which is the pair that pins down which raw each sweep used.
errDC=zeros(1,2);
vfoptionsNL=vfoptions15S;
vfoptionsNL.divideandconquer=1; % and no level1n, deliberately
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsNL);
catch ME
    errDC(1)=contains(ME.message,'level1n');
end
vfoptionsNLG=vfoptions15SG;
vfoptionsNLG.divideandconquer=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsNLG);
catch ME
    errDC(2)=contains(ME.message,'level1n');
end
fprintf('localsearch d: (15) the seed sweep with divideandconquer=1 and no level1n reaches the DC raw and fails on it, which is what proves the seed route honours DC, should be one: %i \n',errDC(1))
fprintf('localsearch d + GI: (15) the same with the interpolation layer, should be one: %i \n',errDC(2))


output=struct();

end
