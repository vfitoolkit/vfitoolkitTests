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
    vfoptionsLS.localsearch=1;
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
% what refW is. Clamped to [2,N_a-1] it IS the GI raw midpoint, and a window of 1 around it spans
% exactly the 2*ngridinterp+3 fine points the GI raw searches, per d. Same points, same EVinterp,
% same ReturnFn: an identity.
% Policy is compared as the decoded fine point, because a point sitting on a coarse node has two
% legal spellings and neither raw is canonical about which it gives.
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

%% (9) The whole-grid window with GI is the exact fine-grid optimum, so it cannot fall below the GI raw
vfoptionsLGw=vfoptionsG;
vfoptionsLGw.localsearch=1;
vfoptionsLGw.nlocalsearch=nwide;
[VLGw,~]=ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsLGw);
fprintf('localsearch d + GI: (9) the whole-grid window is the exact max over the entire fine grid, so it cannot fall below the GI raw, should be one: %i \n',min(gather(VLGw)-gather(VG),[],'all')>=0)
fprintf('localsearch d + GI: (9) for information, it is strictly above the GI raw at %i of %i states \n',sum(gather(VLGw)>gather(VG),'all'),N_a*N_z)

%% (10) The documented refusals, each checked by its own message
err=zeros(1,6);
vfoptionsE=vfoptions; vfoptionsE.localsearch=1; % no nlocalsearch
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(1)=contains(ME.message,'must set vfoptions.nlocalsearch');
end
vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=N_a;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(2)=contains(ME.message,'too large for the grid');
end
vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=0;
try
    ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(3)=contains(ME.message,'must be at least 1');
end
vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.divideandconquer=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
catch ME
    err(4)=contains(ME.message,'combined with vfoptions.divideandconquer');
end
vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1; vfoptionsE.lowmemory=1;
try
    ValueFnIter_InfHorz_TPath_SingleStep(VnextS,n_d,n_a,n_z,d_gridvals, a_gridG, z_gridvals, pi_zG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
catch ME
    err(5)=contains(ME.message,'lowmemory is not used with vfoptions.localsearch');
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
fprintf('localsearch d: (10) combination with divide-and-conquer is refused by name, should be one: %i \n',err(4))
fprintf('localsearch d: (10) lowmemory with local search is refused by name, should be one: %i \n',err(5))
fprintf('localsearch d: (10) a reference of the no-d shape on a model with d is refused by name, should be one: %i \n',err(6))

%% (11) Divide-and-conquer gives the same V, on the coarse grid where there is no memory risk
% The identity that makes the timings below comparable. It is done at n_a rather than n_a_big because
% with d the unrestricted return matrix is N_d times larger and the big grid is near the edge of what
% fits; the timings there are diagnostics and are allowed to fail, a check is not.
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

%% (12) Speed on the big grid, with d. Diagnostics only.
% With d the unrestricted return matrix is [N_d*N_a,N_a,N_z], which at n_a_big=1001 and N_d=9 is
% about 360MB, and entireRHS is another. That is near the edge, so the whole block is allowed to
% fail without taking the part down with it -- the checks above do not depend on it.
% THE PREDICTION, recorded before the run decides it. The no-d part measured about 4.8x and argued
% that the ceiling was the kernel launch floor, not the arithmetic: at n_a=1001 and N_z=5 the
% unrestricted step sat only ~4x above that floor. With N_d=9 it sits nine times further above it
% while the restricted step stays near it, so the measured speed-up should be in the TENS. If it
% comes back near 5x again then the floor is not what is limiting it and I have the wrong model.
try
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
        vfoptionsB.localsearch=1;
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
catch ME
    fprintf('localsearch d: the big grid timing block did not run (%s). No checks depend on it. \n',ME.message)
end

output=struct();

end
