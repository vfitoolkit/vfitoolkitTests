function output=CoreInfHorzTPath_localsearch_e(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE)
% Checks vfoptions.localsearch=1 with z AND an iid e, without and with a decision variable (stage L1 of
% InfHorzTPath_localsearch_noz_e_proposal.md). Written 2026-10-30, before the first run.
% Read CoreInfHorzTPath_localsearch and CoreInfHorzTPath_localsearch_d first: this is those two parts with an e
% dimension, and the oracles are theirs, which are exact wherever the restriction cannot bind.
%
% WHAT IS DIFFERENT WITH e. The reference is per CURRENT e as well: [N_dref,N_a,N_z,N_e], N_dref=max(N_d,1). The
% expectation is over e' and so is the same for every current e, but the return moves with e, and so does the optimal
% aprime; a window shared across e would have to span that whole range. So the line most likely to be wrong is the
% e offset, in the window, in the gather at the chosen d and in Step1's edge count; section (5) drives a reference
% that VARIES with e (and with d) so that a misplaced e offset cannot pass by accident.
%
% THE SEED RAWS. Every e raw has lowmemory 0 (all at once), 1 (loop over e) and 2 (loop over z and e), and every seed
% raw builds the reference inside each branch, so section (7) runs every seed raw under all three. It also carries
% the discriminator the z part found was needed: DC is exact, so only a missing level1n proves the DC raw ran.
%
% Sections:
%   (1) the whole-grid window reproduces the unrestricted path; (2) for scale, a narrow one does not;
%   (3) the reference pinned to the unrestricted optimum gives the unrestricted answer at every n;
%   (4) V non-decreasing in nlocalsearch; (5) the window bound, with a reference varying in e (and d);
%   (6) GI: centred on the GI raw's own midpoint at n=1 it IS the GI raw, and the whole-grid window is too;
%   (7) the seed sweeps (and their DC, lowmemory 1 and 2 versions) against the raws they replace;
%   (8) inside a general eqm path: shooting and Anderson with local search land where the standard solve lands,
%       and an exact re-evaluation of the local search answer is within tolerance;
%   (9) the refusals that remain on this family.

gpudev=gpuDevice;

ngi=5; % the value the rest of this bank uses for ngridinterp

N_a=prod(n_a);
N_z=prod(n_z);
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
N_e=prod(n_e);

if mod(N_a,2)==0
    error('This part needs an odd N_a: the window is 2n+1 points wide, so only an odd grid can be covered exactly by one n')
end
nwide=(N_a-1)/2; % the window is then 1,...,N_a at every state

for dcase=1:2
    % dcase=1 is without d, dcase=2 with d. Everything below is written once and reads the channel layout from this.
    vfoptions=struct();
    vfoptions.n_e=n_e;
    vfoptions.e_grid=e_grid;
    vfoptions.pi_e=pi_e;
    if dcase==1
        lab='localsearch e';
        n_dc=0; d_gridc=[];
        ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
        d_gridvals=[];
        ach=1;      % the aprime channel of Policy without GI
        ach_gi=1;   % the lower coarse point channel of Policy with GI
    else
        lab='localsearch d e';
        n_dc=n_d; d_gridc=d_grid;
        ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
        d_gridvals=gpuArray(d_grid);
        ach=2;
        ach_gi=2;
    end
    N_d=prod(n_dc);
    N_dref=max(N_d,1);

    [V_final,Policy_final]=ValueFnIter_InfHorz(n_dc,n_a,n_z,d_gridc,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

    %% The unrestricted path, which the path-level checks are measured against
    [VPath0,PolicyPath0]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_dc,n_a,n_z, d_gridc,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
    VPath0=reshape(gather(VPath0),[N_a,N_z*N_e,T]);
    V0=VPath0(:,:,1:T-1);
    PolicyPath0=reshape(gather(PolicyPath0),[size(PolicyPath0,1),N_a,N_z*N_e,T]);
    Pol0=PolicyPath0(:,:,:,1:T-1);
    fprintf('%s: the unrestricted policy moves at most %i grid points from the current a index on this path \n',lab,max(abs(reshape(Pol0(ach,:,:,:),[N_a,N_z*N_e,T-1])-(1:1:N_a)'),[],'all'))

    %% The ladder in nlocalsearch, along the path (localsearch=2 forces it in ValueFnOnTransPath_InfHorz, see CoreInfHorzTPath_localsearch_d)
    nladder=unique([1,2,5,10,25,nwide]);
    nladder=nladder(nladder<=nwide);
    monotone=1;
    Vprev=[]; Vn1=[]; Vnwide=[]; Polnwide=[];
    for n_c=1:length(nladder)
        vfoptionsLS=vfoptions;
        vfoptionsLS.localsearch=2;
        vfoptionsLS.nlocalsearch=nladder(n_c);
        [VPathn,PolicyPathn]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_dc,n_a,n_z, d_gridc,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLS);
        VPathn=reshape(gather(VPathn),[N_a,N_z*N_e,T]);
        Vn=VPathn(:,:,1:T-1);
        if ~isempty(Vprev)
            monotone=monotone*(min(Vn-Vprev,[],'all')>=0);
        end
        Vprev=Vn;
        fprintf('%s: nlocalsearch=%3i: %5i candidates per state instead of %6i, max|V-Vfull|=%.3e \n',lab,nladder(n_c),N_dref*(2*nladder(n_c)+1),N_dref*N_a,max(abs(Vn-V0),[],'all'))
        if nladder(n_c)==1
            Vn1=Vn;
        end
        if nladder(n_c)==nwide
            Vnwide=Vn;
            PolicyPathn=reshape(gather(PolicyPathn),[size(PolicyPathn,1),N_a,N_z*N_e,T]);
            Polnwide=PolicyPathn(:,:,:,1:T-1);
        end
    end

    %% (1) The whole-grid window must reproduce the unrestricted path
    % With d the two maxima break exact ties differently from the folded max, so Policy can in principle differ on a
    % tie while V cannot (see CoreInfHorzTPath_localsearch_d section (1)); no tie has been seen on this model.
    fprintf('%s: (1) nlocalsearch=%i covers the grid, V vs the unrestricted path, this should be zero: %.3e \n',lab,nwide,max(abs(Vnwide-V0),[],'all'))
    fprintf('%s: (1) nlocalsearch=%i covers the grid, Policy vs the unrestricted path, this should be zero: %.3e \n',lab,nwide,max(abs(Polnwide-Pol0),[],'all'))

    %% (2) For scale: a narrow window must NOT give the unrestricted answer, or (1) is empty
    fprintf('%s: (2) for scale, nlocalsearch=1 vs the unrestricted path, V, should NOT be zero: %.3e \n',lab,max(abs(Vn1-V0),[],'all'))

    %% (4) V non-decreasing in nlocalsearch
    fprintf('%s: (4) V is non-decreasing in nlocalsearch at every state and period, should be one: %i \n',lab,monotone)

    %% The single-step dispatcher, used by (3), (5), (6), (7) and (9)
    ParamsS=Params;
    ParamsS.r=PricePath.r(1);
    ParamsS.w=ParamPath.w(1);
    ReturnFnParamNames=ReturnFnParamNamesFn(ReturnFn,n_dc,n_a,n_z,0,vfoptions,ParamsS);
    a_gridG=gpuArray(a_grid);
    pi_zG=gpuArray(pi_z);
    pi_eG=gpuArray(pi_e);
    z_gridvals=CreateGridvals(n_z,gpuArray(z_grid),1);
    e_gridvals=CreateGridvals(n_e,gpuArray(e_grid),1);
    VnextS=gpuArray(reshape(V_final,[N_a,N_z,N_e]));

    vfoptionsS.exoticpreferences='None';
    vfoptionsS.lowmemory=0;
    vfoptionsS.gridinterplayer=0;
    vfoptionsS.divideandconquer=0;
    vfoptionsS.localsearch=0;
    vfoptionsS.localsearchseed=0;

    % The linear index of the chosen d into an [N_dref,N_a,N_z,N_e] array, for a Policy whose first channel is d
    % (with no d it is just the state index). Built inline wherever it is needed, from dstar of shape [1,N_a,N_z,N_e].
    stateoffset=N_dref*(0:1:N_a-1)+N_dref*N_a*shiftdim((0:1:N_z-1),-1)+N_dref*N_a*N_z*shiftdim((0:1:N_e-1),-2); % [1,N_a,N_z,N_e]

    % The standard sweep, and the wide window with its reference (the unrestricted per-d coarse argmax)
    [VS,PS]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsS);
    vfoptionsW=vfoptionsS;
    vfoptionsW.localsearch=1;
    vfoptionsW.nlocalsearch=nwide;
    [VW,PW,refW]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsW);
    refWg=gather(refW);
    fprintf('%s: (3) the reference has one entry per (d,a,z,e), should be one: %i \n',lab,numel(refWg)==N_dref*N_a*N_z*N_e && size(refWg,1)==N_dref)
    fprintf('%s: (3) the whole-grid window in one step vs the standard step, V, this should be zero: %.3e \n',lab,max(abs(gather(VW)-gather(VS)),[],'all'))

    %% (3) A reference pinned to the unrestricted optimum gives the unrestricted answer at every n
    % This is the check that the reference is USED: the window around it contains the optimum and the optimum is its
    % max, so nothing can move. Under a reference that was ignored (or read at the wrong e) the n=1 answer would move.
    pin3V=0; pin3P=0;
    for nls=[1,2,5]
        vfoptions3=vfoptionsS;
        vfoptions3.localsearch=1;
        vfoptions3.nlocalsearch=nls;
        [V3,P3]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refW, vfoptions3);
        pin3V=max(pin3V,max(abs(gather(V3)-gather(VW)),[],'all'));
        pin3P=max(pin3P,max(abs(gather(P3)-gather(PW)),[],'all'));
    end
    fprintf('%s: (3) the reference pinned to the unrestricted optimum, worst over nlocalsearch 1, 2 and 5, V, this should be zero: %.3e \n',lab,pin3V)
    fprintf('%s: (3) same, Policy, this should be zero: %.3e \n',lab,pin3P)
    % And Policy's aprime is the new reference read at the chosen d
    PWg=gather(PW);
    if dcase==1
        dstarW=ones(1,N_a,N_z,N_e);
    else
        dstarW=reshape(PWg(1,:,:,:),[1,N_a,N_z,N_e]);
    end
    fprintf('%s: (3) Policy aprime vs the new reference read at the chosen d, this should be zero: %.3e \n',lab,max(abs(reshape(PWg(ach,:,:,:),[1,N_a,N_z,N_e])-refWg(dstarW+stateoffset)),[],'all'))

    %% (5) The window, against a reference that VARIES with e (and with d)
    % ref(d,a,z,e)=clamp(a+(d-ceil(N_d/2))+3*(e-ceil(N_e/2)),1,N_a). With a reference the same for every e, reading it at
    % the wrong e would give the right answer by accident; varying it in e is what makes this check able to fail.
    nls5=2;
    ref5=repmat((1:1:N_a),[N_dref,1,N_z,N_e])+((1:1:N_dref)'-ceil(N_dref/2))+3*(shiftdim((1:1:N_e),-2)-ceil(N_e/2));
    ref5=min(max(ref5,1),N_a); % [N_dref,N_a,N_z,N_e]
    vfoptions5=vfoptionsS;
    vfoptions5.localsearch=1;
    vfoptions5.nlocalsearch=nls5;
    [~,P5,ref5New]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref5, vfoptions5);
    P5g=gather(P5);
    if dcase==1
        dstar5=ones(1,N_a,N_z,N_e);
    else
        dstar5=reshape(P5g(1,:,:,:),[1,N_a,N_z,N_e]);
    end
    refatd=ref5(dstar5+stateoffset);
    dev5=abs(reshape(P5g(ach,:,:,:),[1,N_a,N_z,N_e])-refatd);
    interior5=(refatd-nls5>=1 & refatd+nls5<=N_a);
    devall5=abs(gather(ref5New)-ref5); % every d, not just the chosen one
    interiorall5=(ref5-nls5>=1 & ref5+nls5<=N_a);
    fprintf('%s: (5) with an e-varying reference, |aprime-reference at the chosen d|<=2*nlocalsearch everywhere, should be one: %i \n',lab,max(dev5,[],'all')<=2*nls5)
    fprintf('%s: (5) and <=nlocalsearch wherever that reference is interior, should be one: %i \n',lab,max(dev5(interior5),[],'all')<=nls5)
    fprintf('%s: (5) the new reference for every d, <=nlocalsearch from the old wherever that is interior, should be one: %i \n',lab,max(devall5(interiorall5),[],'all')<=nls5)
    fprintf('%s: (5) the new reference keeps the reference shape, should be one: %i \n',lab,isequal(size(gather(ref5New)),size(ref5)))
    % The window is per e: shift the reference at ONE e only and nothing at the other e may change
    ref5b=ref5;
    ref5b(:,:,:,1)=min(ref5b(:,:,:,1)+7,N_a);
    [V5b,P5b]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref5b, vfoptions5);
    [V5,~]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref5, vfoptions5);
    V5g=gather(V5); V5bg=gather(V5b); P5bg=gather(P5b);
    fprintf('%s: (5) moving the reference at the first e only leaves V at every other e unchanged, this should be zero: %.3e \n',lab,max(abs(V5bg(:,:,2:end)-V5g(:,:,2:end)),[],'all'))
    fprintf('%s: (5) same, Policy, this should be zero: %.3e \n',lab,max(abs(P5bg(:,:,:,2:end)-P5g(:,:,:,2:end)),[],'all'))
    fprintf('%s: (5) and it does change the answer at the first e, should NOT be zero: %.3e \n',lab,max(abs(P5bg(:,:,:,1)-P5g(:,:,:,1)),[],'all'))

    %% (6) With the grid interpolation layer
    vfoptionsG=vfoptionsS;
    vfoptionsG.gridinterplayer=1;
    vfoptionsG.ngridinterp=ngi;
    [VG,PG]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsG);
    PGg=gather(PG);
    % Centred where the GI raw centres: the per-d coarse argmax clamped to [2,N_a-1] is the GI raw's midpoint, and at
    % nlocalsearch=1 the interior clamp collapses to the reference itself, so the first layer is pinned to it
    refGI=max(min(refW,N_a-1),2);
    vfoptionsLG=vfoptionsG;
    vfoptionsLG.localsearch=1;
    vfoptionsLG.nlocalsearch=1;
    [VLG,PLG]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, refGI, vfoptionsLG);
    fprintf('%s + GI: (6) nlocalsearch=1 centred where the GI raw centres, V, this should be zero: %.3e \n',lab,max(abs(gather(VLG)-gather(VG)),[],'all'))
    fprintf('%s + GI: (6) same, Policy in every channel, this should be zero: %.3e \n',lab,max(abs(gather(PLG)-PGg),[],'all'))
    % The whole-grid window: loweredge is 1, so layer 1 is the GI raw's first layer and the clamp is its clamp
    vfoptionsLGw=vfoptionsG;
    vfoptionsLGw.localsearch=1;
    vfoptionsLGw.nlocalsearch=nwide;
    [VLGw,PLGw,refLGw]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsLGw);
    fprintf('%s + GI: (6) the whole-grid window reproduces the standard GI raw, V, this should be zero: %.3e \n',lab,max(abs(gather(VLGw)-gather(VG)),[],'all'))
    fprintf('%s + GI: (6) same, Policy in every channel, this should be zero: %.3e \n',lab,max(abs(gather(PLGw)-PGg),[],'all'))
    fprintf('%s + GI: (6) its new reference is the unclamped coarse argmax, the no-GI wide reference, this should be zero: %.3e \n',lab,max(abs(gather(refLGw)-refWg),[],'all'))
    % GI weakly better than the coarse grid (the fine grid contains the coarse one)
    fprintf('%s + GI: (6) the GI step is weakly better than the coarse step at every state, should be one: %i \n',lab,min(gather(VG)-gather(VS),[],'all')>=-1e-12)
    % The spelling invariant at the top of the grid: the lower coarse point never reaches N_a
    vfoptions6t=vfoptionsG;
    vfoptions6t.localsearch=1;
    vfoptions6t.nlocalsearch=2;
    [~,P6t]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, N_a*ones(N_dref,N_a,N_z,N_e), vfoptions6t);
    P6tg=gather(P6t);
    fprintf('%s + GI: (6) with the reference pinned to the top of the grid the coarse channel stays in [1,N_a-1], should be one: %i \n',lab,all(P6tg(ach_gi,:,:,:)>=1 & P6tg(ach_gi,:,:,:)<=N_a-1,'all'))

    %% (7) The seed sweeps, against the raws they replace
    % Without GI and without d the seed is the standard raw, and the reference is Policy's aprime channel. With d it is
    % LSseed_e_raw (or DC1_e_raw's split max): V identical, Policy identical up to exact ties. With GI it is the
    % standard GI raw with a third output, so V and Policy are identical outright.
    for gi=0:1
        if gi==0
            vfoptions7=vfoptionsS;
            Vstd=gather(VS); Pstd=gather(PS);
            labgi='';
        else
            vfoptions7=vfoptionsG;
            Vstd=gather(VG); Pstd=PGg;
            labgi=' + GI';
        end
        vfoptions7.localsearchseed=1;
        [V7,P7,ref7]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions7);
        V7=gather(V7); P7=gather(P7); ref7=gather(ref7);
        fprintf('%s%s: (7) the seed sweep vs the standard sweep, V, this should be zero: %.3e \n',lab,labgi,max(abs(V7-Vstd),[],'all'))
        if dcase==2 && gi==0
            fprintf('%s%s: (7) for information, the seed and standard Policy differ at %i of %i states, which is exact ties broken differently \n',lab,labgi,sum(any(P7~=Pstd,1),'all'),N_a*N_z*N_e)
        else
            fprintf('%s%s: (7) the seed sweep vs the standard sweep, Policy, this should be zero: %.3e \n',lab,labgi,max(abs(P7-Pstd),[],'all'))
        end
        if gi==0
            fprintf('%s%s: (7) the seed reference is the unrestricted per-d coarse argmax (the whole-grid window reference), this should be zero: %.3e \n',lab,labgi,max(abs(ref7-refWg),[],'all'))
        else
            fprintf('%s%s: (7) the seed reference is the GI raw midpoint, the per-d coarse argmax clamped to [2,N_a-1], this should be zero: %.3e \n',lab,labgi,max(abs(ref7-max(min(refWg,N_a-1),2)),[],'all'))
        end
        for dc=0:1
            vfoptions7dc=vfoptions7;
            vfoptions7dc.divideandconquer=dc;
            vfoptions7dc.level1n=10;
            vfoptions7dc.lowmemory=0;
            [V70,P70,ref70]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions7dc);
            V70=gather(V70); P70=gather(P70); ref70=gather(ref70);
            if dc==1
                labdc=' + DC';
                fprintf('%s%s%s: (7) the DC seed sweep vs the seed sweep, V, this should be zero: %.3e \n',lab,labgi,labdc,max(abs(V70-V7),[],'all'))
                % The internal check on the reference, which does not depend on how DC breaks ties
                if dcase==1
                    dstar7=ones(1,N_a,N_z,N_e);
                else
                    dstar7=reshape(P70(1,:,:,:),[1,N_a,N_z,N_e]);
                end
                ref7atd=ref70(dstar7+stateoffset);
                c7=reshape(P70(ach_gi,:,:,:),[1,N_a,N_z,N_e]);
                if gi==0
                    fprintf('%s%s%s: (7) the DC seed reference read at the chosen d vs that Policy aprime, this should be zero: %.3e \n',lab,labgi,labdc,max(abs(ref7atd-c7),[],'all'))
                else
                    fprintf('%s%s%s: (7) the DC seed reference at the chosen d is the midpoint of the cell Policy sits in, should be one: %i \n',lab,labgi,labdc,all(c7==ref7atd | c7==ref7atd-1,'all'))
                    fprintf('%s%s%s: (7) the DC seed reference lies in [2,N_a-1] for every d, as a midpoint must, should be one: %i \n',lab,labgi,labdc,all(ref70>=2 & ref70<=N_a-1,'all'))
                end
            else
                labdc='';
            end
            for lm=1:2
                vfoptions7lm=vfoptions7dc;
                vfoptions7lm.lowmemory=lm;
                [V7lm,P7lm,ref7lm]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptions7lm);
                fprintf('%s%s%s: (7) the seed sweep under lowmemory=%i vs lowmemory=0, V, Policy and reference, this should be zero: %.3e \n',lab,labgi,labdc,lm,max([max(abs(gather(V7lm)-V70),[],'all'),max(abs(gather(P7lm)-P70),[],'all'),max(abs(gather(ref7lm)-ref70),[],'all')]))
            end
        end
        % THE DISCRIMINATOR: DC is exact, so only a missing level1n proves the seed route reached the DC raw
        vfoptionsNL=vfoptions7;
        vfoptionsNL.divideandconquer=1; % and no level1n, deliberately
        errNL=0;
        try
            ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsNL);
        catch ME
            errNL=contains(ME.message,'level1n');
        end
        fprintf('%s%s: (7) the seed sweep with divideandconquer=1 and no level1n reaches the DC raw and fails on it, should be one: %i \n',lab,labgi,errNL)
    end
    % A RESTRICTED sweep ignores divideandconquer and lowmemory: identical answer
    vfoptions7r=vfoptionsS;
    vfoptions7r.localsearch=1;
    vfoptions7r.nlocalsearch=2;
    [V7r0,P7r0,ref7r0]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref5, vfoptions7r);
    vfoptions7r.divideandconquer=1; % no level1n: only an ignored option permits that
    vfoptions7r.lowmemory=2;
    [V7r1,P7r1,ref7r1]=ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref5, vfoptions7r);
    fprintf('%s: (7) a restricted sweep under divideandconquer=1 and lowmemory=2 vs neither, V, Policy and reference, this should be zero: %.3e \n',lab,max([max(abs(gather(V7r1)-gather(V7r0)),[],'all'),max(abs(gather(P7r1)-gather(P7r0)),[],'all'),max(abs(gather(ref7r1)-gather(ref7r0)),[],'all')]))

    %% (9) The refusals on this family, each checked by its message
    err9=zeros(1,2);
    vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
    try
        % a reference of the no-e shape, [N_dref,N_a,N_z]
        ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, repmat((1:1:N_a),[N_dref,1,N_z]), vfoptionsE);
    catch ME
        err9(1)=contains(ME.message,'one aprime index per state the search conditions on');
    end
    vfoptionsE=vfoptionsS; vfoptionsE.localsearch=2; vfoptionsE.nlocalsearch=1;
    try
        ValueFnIter_InfHorz_TPath_SingleStep_e(VnextS,n_dc,n_a,n_z,n_e,d_gridvals, a_gridG, z_gridvals, e_gridvals, pi_zG, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, [], vfoptionsE);
    catch ME
        err9(2)=contains(ME.message,'must be 0 or 1 here');
    end
    fprintf('%s: (9) a reference without the e dimension is refused by name, should be one: %i \n',lab,err9(1))
    fprintf('%s: (9) localsearch=2 reaching the single step is refused by name, should be one: %i \n',lab,err9(2))
end

%% (8) Inside a general eqm transition path, without and with d
% The model of doPart(7) with an e: K is assets, N is efficiency labour (z*e, or d*z*e with d). Shooting and Anderson,
% each standard and with local search; then the local search answer is re-evaluated with one exact sweep (zero
% damping, maxiter=1), which is the certificate.
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
for dcase=1:2
    vfoptions=struct();
    vfoptions.n_e=n_e;
    vfoptions.e_grid=e_grid;
    vfoptions.pi_e=pi_e;
    simoptions=struct();
    simoptions.n_e=n_e;
    simoptions.e_grid=e_grid;
    simoptions.pi_e=pi_e;
    FnsToEvaluateGE=struct();
    if dcase==1
        lab='localsearch e GE';
        n_dc=0; d_gridc=[];
        ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
        FnsToEvaluateGE.K=@(aprime,a,z,e) a;
        FnsToEvaluateGE.N=@(aprime,a,z,e) z*e;
    else
        lab='localsearch d e GE';
        n_dc=n_d_GE; d_gridc=d_grid_GE;
        ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
        FnsToEvaluateGE.K=@(d,aprime,a,z,e) a;
        FnsToEvaluateGE.N=@(d,aprime,a,z,e) d*z*e;
    end
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

    Params.r=0.04; Params.w=1; % initial guess
    [p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_dc, n_a_GE, n_z, 0, pi_z, d_gridc, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptions, vfoptions);
    fprintf('%s: stationary GE r=%2.8f, w=%2.8f, conditions CapitalMarket=%.3e LabourMarket=%.3e \n',lab,p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)
    ParamsGE=Params;
    ParamsGE.r=p_eqm.r;
    ParamsGE.w=p_eqm.w;
    [V_finalGE,Policy_finalGE]=ValueFnIter_InfHorz(n_dc,n_a_GE,n_z,d_gridc,a_grid_GE,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initialGE=StationaryDist_InfHorz(Policy_finalGE,n_dc,n_a_GE,n_z,pi_z,simoptions,ParamsGE,[]);

    PricePathGE.r=p_eqm.r*ones(1,T);
    PricePathGE.w=p_eqm.w*ones(1,T);
    ParamPathGE=struct();
    ParamPathGE.sigma=Params.sigma*ones(1,T); % constant: a null reform, so the answer is the stationary eqm
    Tbump=floor(2*T/3);
    PricePathBumped=PricePathGE;
    PricePathBumped.r(1:Tbump)=PricePathGE.r(1:Tbump).*(1+0.01+0.04*linspace(0,1,Tbump));
    PricePathBumped.w(1:Tbump)=PricePathGE.w(1:Tbump).*(1+0.05-0.04*linspace(0,1,Tbump));

    GEfloor=max(abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket));
    transpathoptionsS=transpathoptionsbaseline;
    transpathoptionsS.GEnewprice=3;
    transpathoptionsS.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
    transpathoptionsS.maxiter=1000;
    transpathoptionsS.verbose=0;
    transpathoptionsS.toleranceGEprices=10^(-4);
    transpathoptionsS.toleranceGEcondns=5*GEfloor;
    transpathoptionsAnd=transpathoptionsbaseline;
    transpathoptionsAnd.GEnewprice=2;
    transpathoptionsAnd.GEnewprice2.howtoupdate=transpathoptionsS.GEnewprice3.howtoupdate;
    transpathoptionsAnd.maxiter=200;
    transpathoptionsAnd.verbose=0;
    transpathoptionsAnd.toleranceGEprices=transpathoptionsS.toleranceGEprices;
    transpathoptionsAnd.toleranceGEcondns=transpathoptionsS.toleranceGEcondns;

    vfoptionsLS=vfoptions;
    vfoptionsLS.localsearch=1;
    vfoptionsLS.nlocalsearch=1;

    tic;
    PricePathA=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptions, []);
    wait(gpudev); tA=toc;
    tic;
    PricePathB=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsS, simoptions, vfoptionsLS, []);
    wait(gpudev); tB=toc;
    PricePathAnd=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAnd, simoptions, vfoptions, []);
    PricePathAndLS=TransitionPath_InfHorz(PricePathBumped, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsAnd, simoptions, vfoptionsLS, []);
    fprintf('%s: for information, standard shooting %6.2fs, shooting with local search from nlocalsearch=1 %6.2fs \n',lab,tA,tB)
    fprintf('%s: (8) shooting with local search vs the standard shooting solve, r, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathB.r)-gather(PricePathA.r))))
    fprintf('%s: (8) same, w, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathB.w)-gather(PricePathA.w))))
    fprintf('%s: (8) Anderson with local search vs the standard Anderson solve, r, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathAndLS.r)-gather(PricePathAnd.r))))
    fprintf('%s: (8) same, w, this should be close to zero: %.3e \n',lab,max(abs(gather(PricePathAndLS.w)-gather(PricePathAnd.w))))
    fprintf('%s: (8) for information, the standard shooting solve deviates from the stationary eqm by at most r %.3e, w %.3e \n',lab,max(abs(gather(PricePathA.r)-p_eqm.r)),max(abs(gather(PricePathA.w)-p_eqm.w)))
    % The certificate: re-evaluate each local search answer with one exact sweep and zero damping
    transpathoptionsV=transpathoptionsbaseline;
    transpathoptionsV.GEnewprice=3;
    transpathoptionsV.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0; 'LabourMarket','w',0,0};
    transpathoptionsV.maxiter=1;
    transpathoptionsV.verbose=0;
    [~,GEcondnPathV]=TransitionPath_InfHorz(PricePathB, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsV, simoptions, vfoptions, []);
    GEdistV=max(max(abs(gather(GEcondnPathV.CapitalMarket))),max(abs(gather(GEcondnPathV.LabourMarket))));
    [~,GEcondnPathVA]=TransitionPath_InfHorz(PricePathAndLS, ParamPathGE, T, V_finalGE, AgentDist_initialGE, n_dc, n_a_GE, n_z, d_gridc,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsGE, DiscountFactorParamNames, transpathoptionsV, simoptions, vfoptions, []);
    GEdistVA=max(max(abs(gather(GEcondnPathVA.CapitalMarket))),max(abs(gather(GEcondnPathVA.LabourMarket))));
    fprintf('%s: (8) the shooting local search answer re-evaluated with an exact sweep is within toleranceGEcondns, should be one: %i \n',lab,GEdistV<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (8) the Anderson local search answer re-evaluated with an exact sweep is within toleranceGEcondns, should be one: %i \n',lab,GEdistVA<=transpathoptionsS.toleranceGEcondns)
    fprintf('%s: (8) for information, those distances are %.3e and %.3e against toleranceGEcondns %.3e \n',lab,GEdistV,GEdistVA,transpathoptionsS.toleranceGEcondns)
end

output=struct();

end
