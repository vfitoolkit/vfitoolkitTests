function output=CoreInfHorzTPath_localsearch_noz(T,PricePath,ParamPath,n_d,n_a,d_grid,a_grid,Params,DiscountFactorParamNames,transpathoptionsbaseline,vfoptionsbaseline,n_d_GE,n_a_GE,d_grid_GE,a_grid_GE)
% Checks vfoptions.localsearch=1 WITHOUT z: no z and no e (stage L2), and no z with an iid e (stage L3), each without
% and with a decision variable. Written 2026-11-01, before the first run. The counterpart of CoreInfHorzTPath_localsearch_e
% (stage L1, z and e); read that header first, the sections and oracles here are its, with the same numbers.
%
% The reference is indexed by everything the search conditions on: [N_dref,N_a] without e, [N_dref,N_a,N_e] with e.
% Everything below is written in the (d,a,e) layout with NE=N_e (with e) or NE=1 (without), so one set of lines covers
% both families: an [N_dref,N_a] array IS an [N_dref,N_a,1] array.
% Without z and without e there is no uncertainty at all, so the stationary general eqm is degenerate (the interest
% rate is pinned at 1/beta-1 and the distribution is not unique); section (8) is therefore run with e only.

gpudev=gpuDevice;

ngi=5; % the value the rest of this bank uses for ngridinterp

N_a=prod(n_a);
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

if mod(N_a,2)==0
    error('This part needs an odd N_a: the window is 2n+1 points wide, so only an odd grid can be covered exactly by one n')
end
nwide=(N_a-1)/2; % the window is then 1,...,N_a at every state

n_z=0; z_grid=[]; pi_z=[]; % no z: the inputs stay positional, as empty

for fam=1:2
    for dcase=1:2
        vfoptions=struct();
        if fam==1
            famlab='noz';
            NE=1;
        else
            famlab='noz e';
            vfoptions.n_e=n_e;
            vfoptions.e_grid=e_grid;
            vfoptions.pi_e=pi_e;
            NE=prod(n_e);
        end
        if dcase==1
            lab=['localsearch ',famlab];
            n_dc=0; d_gridc=[]; d_gridvals=[];
            if fam==1
                ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
            else
                ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
            end
            ach=1; % the aprime channel of Policy without GI (and the lower coarse point channel with GI)
        else
            lab=['localsearch d ',famlab];
            n_dc=n_d; d_gridc=d_grid; d_gridvals=gpuArray(d_grid);
            if fam==1
                ReturnFn=@(d,aprime,a,r,w,sigma,eta,varphi) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,sigma,eta,varphi);
            else
                ReturnFn=@(d,aprime,a,e,r,w,sigma,eta,varphi) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,sigma,eta,varphi);
            end
            ach=2;
        end
        N_d=prod(n_dc);
        N_dref=max(N_d,1);

        [V_final,Policy_final]=ValueFnIter_InfHorz(n_dc,n_a,n_z,d_gridc,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);

        %% The unrestricted path, which the path-level checks are measured against
        [VPath0,PolicyPath0]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_dc,n_a,n_z, d_gridc,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
        VPath0=reshape(gather(VPath0),[N_a,NE,T]);
        V0=VPath0(:,:,1:T-1);
        PolicyPath0=reshape(gather(PolicyPath0),[size(PolicyPath0,1),N_a,NE,T]);
        Pol0=PolicyPath0(:,:,:,1:T-1);
        fprintf('%s: the unrestricted policy moves at most %i grid points from the current a index on this path \n',lab,max(abs(reshape(Pol0(ach,:,:,:),[N_a,NE,T-1])-(1:1:N_a)'),[],'all'))

        %% The ladder in nlocalsearch, along the path (localsearch=2 forces it in ValueFnOnTransPath_InfHorz)
        nladder=unique([1,2,5,10,25,nwide]);
        nladder=nladder(nladder<=nwide);
        monotone=1;
        Vprev=[]; Vn1=[]; Vnwide=[]; Polnwide=[];
        for n_c=1:length(nladder)
            vfoptionsLS=vfoptions;
            vfoptionsLS.localsearch=2;
            vfoptionsLS.nlocalsearch=nladder(n_c);
            [VPathn,PolicyPathn]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_dc,n_a,n_z, d_gridc,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptionsLS);
            VPathn=reshape(gather(VPathn),[N_a,NE,T]);
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
                PolicyPathn=reshape(gather(PolicyPathn),[size(PolicyPathn,1),N_a,NE,T]);
                Polnwide=PolicyPathn(:,:,:,1:T-1);
            end
        end
        fprintf('%s: (1) nlocalsearch=%i covers the grid, V vs the unrestricted path, this should be zero: %.3e \n',lab,nwide,max(abs(Vnwide-V0),[],'all'))
        fprintf('%s: (1) nlocalsearch=%i covers the grid, Policy vs the unrestricted path, this should be zero: %.3e \n',lab,nwide,max(abs(Polnwide-Pol0),[],'all'))
        fprintf('%s: (2) for scale, nlocalsearch=1 vs the unrestricted path, V, should NOT be zero: %.3e \n',lab,max(abs(Vn1-V0),[],'all'))
        fprintf('%s: (4) V is non-decreasing in nlocalsearch at every state and period, should be one: %i \n',lab,monotone)

        %% The single-step dispatcher, used by (3), (5), (6), (7) and (9)
        ParamsS=Params;
        ParamsS.r=PricePath.r(1);
        ParamsS.w=ParamPath.w(1);
        ReturnFnParamNames=ReturnFnParamNamesFn(ReturnFn,n_dc,n_a,n_z,0,vfoptions,ParamsS);
        a_gridG=gpuArray(a_grid);
        if fam==1
            VnextS=gpuArray(reshape(V_final,[N_a,1]));
            SS=@(Vn,ref,vfo) ValueFnIter_InfHorz_TPath_SingleStep_noz(Vn,n_dc,n_a,d_gridvals, a_gridG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref, vfo);
            lmlist=[]; % no exogenous state, so no lowmemory branch
        else
            e_gridvals=CreateGridvals(n_e,gpuArray(e_grid),1);
            pi_eG=gpuArray(pi_e);
            VnextS=gpuArray(reshape(V_final,[N_a,NE]));
            SS=@(Vn,ref,vfo) ValueFnIter_InfHorz_TPath_SingleStep_noz_e(Vn,n_dc,n_a,n_e,d_gridvals, a_gridG, e_gridvals, pi_eG, ReturnFn, ParamsS, DiscountFactorParamNames, ReturnFnParamNames, ref, vfo);
            lmlist=1; % loops over e (lowmemory=2 loops over z and e, so it does not exist without z)
        end

        vfoptionsS.exoticpreferences='None';
        vfoptionsS.lowmemory=0;
        vfoptionsS.gridinterplayer=0;
        vfoptionsS.divideandconquer=0;
        vfoptionsS.localsearch=0;
        vfoptionsS.localsearchseed=0;

        stateoffset=N_dref*(0:1:N_a-1)+N_dref*N_a*shiftdim((0:1:NE-1),-1); % [1,N_a,NE], into an [N_dref,N_a,NE] array

        [VS,PS]=SS(VnextS,[],vfoptionsS);
        vfoptionsW=vfoptionsS;
        vfoptionsW.localsearch=1;
        vfoptionsW.nlocalsearch=nwide;
        [VW,PW,refW]=SS(VnextS,[],vfoptionsW);
        refWg=gather(refW);
        fprintf('%s: (3) the reference has one entry per state the search conditions on, should be one: %i \n',lab,numel(refWg)==N_dref*N_a*NE && size(refWg,1)==N_dref)
        fprintf('%s: (3) the whole-grid window in one step vs the standard step, V, this should be zero: %.3e \n',lab,max(abs(gather(VW)-gather(VS)),[],'all'))

        %% (3) A reference pinned to the unrestricted optimum gives the unrestricted answer at every n
        pin3V=0; pin3P=0;
        for nls=[1,2,5]
            vfoptions3=vfoptionsS;
            vfoptions3.localsearch=1;
            vfoptions3.nlocalsearch=nls;
            [V3,P3]=SS(VnextS,refW,vfoptions3);
            pin3V=max(pin3V,max(abs(gather(V3)-gather(VW)),[],'all'));
            pin3P=max(pin3P,max(abs(gather(P3)-gather(PW)),[],'all'));
        end
        fprintf('%s: (3) the reference pinned to the unrestricted optimum, worst over nlocalsearch 1, 2 and 5, V, this should be zero: %.3e \n',lab,pin3V)
        fprintf('%s: (3) same, Policy, this should be zero: %.3e \n',lab,pin3P)
        PWg=reshape(gather(PW),[size(PW,1),N_a,NE]);
        if dcase==1
            dstarW=ones(1,N_a,NE);
        else
            dstarW=reshape(PWg(1,:,:),[1,N_a,NE]);
        end
        fprintf('%s: (3) Policy aprime vs the new reference read at the chosen d, this should be zero: %.3e \n',lab,max(abs(reshape(PWg(ach,:,:),[1,N_a,NE])-refWg(dstarW+stateoffset)),[],'all'))

        %% (5) The window, against a reference that VARIES with d (and with e)
        nls5=2;
        ref5=repmat((1:1:N_a),[N_dref,1,NE])+((1:1:N_dref)'-ceil(N_dref/2))+3*(shiftdim((1:1:NE),-1)-ceil(NE/2));
        ref5=min(max(ref5,1),N_a); % [N_dref,N_a,NE]
        vfoptions5=vfoptionsS;
        vfoptions5.localsearch=1;
        vfoptions5.nlocalsearch=nls5;
        [V5,P5,ref5New]=SS(VnextS,ref5,vfoptions5);
        P5g=reshape(gather(P5),[size(P5,1),N_a,NE]);
        if dcase==1
            dstar5=ones(1,N_a,NE);
        else
            dstar5=reshape(P5g(1,:,:),[1,N_a,NE]);
        end
        refatd=ref5(dstar5+stateoffset);
        dev5=abs(reshape(P5g(ach,:,:),[1,N_a,NE])-refatd);
        interior5=(refatd-nls5>=1 & refatd+nls5<=N_a);
        devall5=abs(reshape(gather(ref5New),[N_dref,N_a,NE])-ref5);
        interiorall5=(ref5-nls5>=1 & ref5+nls5<=N_a);
        fprintf('%s: (5) with a varying reference, |aprime-reference at the chosen d|<=2*nlocalsearch everywhere, should be one: %i \n',lab,max(dev5,[],'all')<=2*nls5)
        fprintf('%s: (5) and <=nlocalsearch wherever that reference is interior, should be one: %i \n',lab,max(dev5(interior5),[],'all')<=nls5)
        fprintf('%s: (5) the new reference for every d, <=nlocalsearch from the old wherever that is interior, should be one: %i \n',lab,max(devall5(interiorall5),[],'all')<=nls5)
        fprintf('%s: (5) the new reference has one entry per state, should be one: %i \n',lab,numel(gather(ref5New))==numel(ref5))
        if fam==2
            % The window is per e: shift the reference at ONE e only and nothing at the other e may change
            ref5b=ref5;
            ref5b(:,:,1)=min(ref5b(:,:,1)+7,N_a);
            [V5b,P5b]=SS(VnextS,ref5b,vfoptions5);
            V5g=gather(V5); V5bg=gather(V5b); P5bg=reshape(gather(P5b),[size(P5b,1),N_a,NE]);
            fprintf('%s: (5) moving the reference at the first e only leaves V at every other e unchanged, this should be zero: %.3e \n',lab,max(abs(V5bg(:,2:end)-V5g(:,2:end)),[],'all'))
            fprintf('%s: (5) same, Policy, this should be zero: %.3e \n',lab,max(abs(P5bg(:,:,2:end)-P5g(:,:,2:end)),[],'all'))
            fprintf('%s: (5) and it does change the answer at the first e, should NOT be zero: %.3e \n',lab,max(abs(P5bg(:,:,1)-P5g(:,:,1)),[],'all'))
        end

        %% (6) With the grid interpolation layer
        vfoptionsG=vfoptionsS;
        vfoptionsG.gridinterplayer=1;
        vfoptionsG.ngridinterp=ngi;
        [VG,PG]=SS(VnextS,[],vfoptionsG);
        PGg=gather(PG);
        refGI=max(min(refW,N_a-1),2);
        vfoptionsLG=vfoptionsG;
        vfoptionsLG.localsearch=1;
        vfoptionsLG.nlocalsearch=1;
        [VLG,PLG]=SS(VnextS,refGI,vfoptionsLG);
        fprintf('%s + GI: (6) nlocalsearch=1 centred where the GI raw centres, V, this should be zero: %.3e \n',lab,max(abs(gather(VLG)-gather(VG)),[],'all'))
        fprintf('%s + GI: (6) same, Policy in every channel, this should be zero: %.3e \n',lab,max(abs(gather(PLG)-PGg),[],'all'))
        vfoptionsLGw=vfoptionsG;
        vfoptionsLGw.localsearch=1;
        vfoptionsLGw.nlocalsearch=nwide;
        [VLGw,PLGw,refLGw]=SS(VnextS,[],vfoptionsLGw);
        fprintf('%s + GI: (6) the whole-grid window reproduces the standard GI raw, V, this should be zero: %.3e \n',lab,max(abs(gather(VLGw)-gather(VG)),[],'all'))
        fprintf('%s + GI: (6) same, Policy in every channel, this should be zero: %.3e \n',lab,max(abs(gather(PLGw)-PGg),[],'all'))
        fprintf('%s + GI: (6) its new reference is the unclamped coarse argmax, the no-GI wide reference, this should be zero: %.3e \n',lab,max(abs(gather(refLGw)-refWg),[],'all'))
        fprintf('%s + GI: (6) the GI step is weakly better than the coarse step at every state, should be one: %i \n',lab,min(gather(VG)-gather(VS),[],'all')>=-1e-12)
        vfoptions6t=vfoptionsG;
        vfoptions6t.localsearch=1;
        vfoptions6t.nlocalsearch=2;
        [~,P6t]=SS(VnextS,N_a*ones(N_dref,N_a,NE),vfoptions6t);
        P6tg=reshape(gather(P6t),[size(P6t,1),N_a,NE]);
        fprintf('%s + GI: (6) with the reference pinned to the top of the grid the coarse channel stays in [1,N_a-1], should be one: %i \n',lab,all(P6tg(ach,:,:)>=1 & P6tg(ach,:,:)<=N_a-1,'all'))

        %% (7) The seed sweeps, against the raws they replace
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
            [V7,P7,ref7]=SS(VnextS,[],vfoptions7);
            V7=gather(V7); P7=gather(P7); ref7=gather(ref7);
            fprintf('%s%s: (7) the seed sweep vs the standard sweep, V, this should be zero: %.3e \n',lab,labgi,max(abs(V7-Vstd),[],'all'))
            if dcase==2 && gi==0
                fprintf('%s%s: (7) for information, the seed and standard Policy differ at %i of %i states, which is exact ties broken differently \n',lab,labgi,sum(any(P7~=Pstd,1),'all'),N_a*NE)
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
                [V70,P70,ref70]=SS(VnextS,[],vfoptions7dc);
                V70=gather(V70); P70=gather(P70); ref70=gather(ref70);
                if dc==1
                    labdc=' + DC';
                    fprintf('%s%s%s: (7) the DC seed sweep vs the seed sweep, V, this should be zero: %.3e \n',lab,labgi,labdc,max(abs(V70-V7),[],'all'))
                    P70r=reshape(P70,[size(P70,1),N_a,NE]);
                    if dcase==1
                        dstar7=ones(1,N_a,NE);
                    else
                        dstar7=reshape(P70r(1,:,:),[1,N_a,NE]);
                    end
                    ref7atd=ref70(dstar7+stateoffset);
                    c7=reshape(P70r(ach,:,:),[1,N_a,NE]);
                    if gi==0
                        fprintf('%s%s%s: (7) the DC seed reference read at the chosen d vs that Policy aprime, this should be zero: %.3e \n',lab,labgi,labdc,max(abs(ref7atd-c7),[],'all'))
                    else
                        fprintf('%s%s%s: (7) the DC seed reference at the chosen d is the midpoint of the cell Policy sits in, should be one: %i \n',lab,labgi,labdc,all(c7==ref7atd | c7==ref7atd-1,'all'))
                        fprintf('%s%s%s: (7) the DC seed reference lies in [2,N_a-1] for every d, as a midpoint must, should be one: %i \n',lab,labgi,labdc,all(ref70>=2 & ref70<=N_a-1,'all'))
                    end
                else
                    labdc='';
                end
                for lm=lmlist
                    vfoptions7lm=vfoptions7dc;
                    vfoptions7lm.lowmemory=lm;
                    [V7lm,P7lm,ref7lm]=SS(VnextS,[],vfoptions7lm);
                    fprintf('%s%s%s: (7) the seed sweep under lowmemory=%i vs lowmemory=0, V, Policy and reference, this should be zero: %.3e \n',lab,labgi,labdc,lm,max([max(abs(gather(V7lm)-V70),[],'all'),max(abs(gather(P7lm)-P70),[],'all'),max(abs(gather(ref7lm)-ref70),[],'all')]))
                end
            end
            % THE DISCRIMINATOR: DC is exact, so only a missing level1n proves the seed route reached the DC raw
            vfoptionsNL=vfoptions7;
            vfoptionsNL.divideandconquer=1; % and no level1n, deliberately
            errNL=0;
            try
                SS(VnextS,[],vfoptionsNL);
            catch ME
                errNL=contains(ME.message,'level1n');
            end
            fprintf('%s%s: (7) the seed sweep with divideandconquer=1 and no level1n reaches the DC raw and fails on it, should be one: %i \n',lab,labgi,errNL)
        end
        % A RESTRICTED sweep ignores divideandconquer (and lowmemory): identical answer
        vfoptions7r=vfoptionsS;
        vfoptions7r.localsearch=1;
        vfoptions7r.nlocalsearch=2;
        [V7r0,P7r0,ref7r0]=SS(VnextS,ref5,vfoptions7r);
        vfoptions7r.divideandconquer=1; % no level1n: only an ignored option permits that
        if fam==2
            vfoptions7r.lowmemory=1;
        end
        [V7r1,P7r1,ref7r1]=SS(VnextS,ref5,vfoptions7r);
        fprintf('%s: (7) a restricted sweep under divideandconquer=1 (and lowmemory=1 with e) vs neither, V, Policy and reference, this should be zero: %.3e \n',lab,max([max(abs(gather(V7r1)-gather(V7r0)),[],'all'),max(abs(gather(P7r1)-gather(P7r0)),[],'all'),max(abs(gather(ref7r1)-gather(ref7r0)),[],'all')]))

        %% (9) The refusals on this family, each checked by its message
        err9=zeros(1,2);
        vfoptionsE=vfoptionsS; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
        try
            % a reference with a state too few: [N_dref,N_a-1,NE]
            SS(VnextS,ones(N_dref,N_a-1,NE),vfoptionsE);
        catch ME
            err9(1)=contains(ME.message,'one aprime index per state the search conditions on');
        end
        vfoptionsE=vfoptionsS; vfoptionsE.localsearch=2; vfoptionsE.nlocalsearch=1;
        try
            SS(VnextS,[],vfoptionsE);
        catch ME
            err9(2)=contains(ME.message,'must be 0 or 1 here');
        end
        fprintf('%s: (9) a reference of the wrong shape is refused by name, should be one: %i \n',lab,err9(1))
        fprintf('%s: (9) localsearch=2 reaching the single step is refused by name, should be one: %i \n',lab,err9(2))
    end
end

%% (8) Inside a general eqm transition path, no z with e, without and with d
% The model of CoreInfHorzTPath_localsearch_e section (8) without z: K is assets, N is efficiency labour (e, or d*e).
% The stationary solve is started from a GI solve, as in CoreInfHorzTPath_d_noz_e_nosemiz (K and N are step functions
% of the prices on coarse grids and the default start stalls).
% THE e HERE IS RISKIER THAN THE BANK'S (5 points, sd 0.4, against 3 points, sd 0.1). Without z, the bank's e is
% almost no risk: there is next to no precautionary saving, r sits at 1/beta-1 where asset demand is nearly flat, and
% the stationary solve stalled (CapitalMarket residual 1.7e-2 without d and 5.4e-2 with d on 2026-11-01), which made
% toleranceGEcondns=5*residual so loose that every check below passed without testing anything (Anderson with and
% without local search agreed to 0, as both stopped at the first evaluation). Hence also the guard: the stationary
% eqm must be solved to 1e-3, or the section says so as a failure rather than passing vacuously.
[e_gridGE,pi_eGE]=discretizeAR1_FarmerToda(0,0,0.4,5);
pi_eGE=pi_eGE(1,:)';
e_gridGE=exp(e_gridGE);
n_eGE=5;
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
for dcase=1:2
    vfoptions=struct();
    vfoptions.n_e=n_eGE;
    vfoptions.e_grid=e_gridGE;
    vfoptions.pi_e=pi_eGE;
    simoptions=struct();
    simoptions.n_e=n_eGE;
    simoptions.e_grid=e_gridGE;
    simoptions.pi_e=pi_eGE;
    FnsToEvaluateGE=struct();
    if dcase==1
        lab='localsearch noz e GE';
        n_dc=0; d_gridc=[];
        ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
        FnsToEvaluateGE.K=@(aprime,a,e) a;
        FnsToEvaluateGE.N=@(aprime,a,e) e;
    else
        lab='localsearch d noz e GE';
        n_dc=n_d_GE; d_gridc=d_grid_GE;
        ReturnFn=@(d,aprime,a,e,r,w,sigma,eta,varphi) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,sigma,eta,varphi);
        FnsToEvaluateGE.K=@(d,aprime,a,e) a;
        FnsToEvaluateGE.N=@(d,aprime,a,e) d*e;
    end
    GeneralEqmEqnsGE=struct();
    GeneralEqmEqnsGE.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
    GeneralEqmEqnsGE.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

    Params.r=0.04; Params.w=1; % initial guess
    ParamsStart=Params;
    % Both models get a start from a GI solve: without d as well, on this grid the default start stalled
    % (CapitalMarket residual 3.2e-3 on 2026-11-01, caught by the guard below); with d it got 3.1e-4 from the GI start
    vfoptionsGIs=vfoptions; vfoptionsGIs.gridinterplayer=1; vfoptionsGIs.ngridinterp=5;
    simoptionsGIs=simoptions; simoptionsGIs.gridinterplayer=1; simoptionsGIs.ngridinterp=5;
    p_eqmstart=HeteroAgentStationaryEqm_InfHorz(n_dc, n_a_GE, n_z, 0, pi_z, d_gridc, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptionsGIs, vfoptionsGIs);
    ParamsStart.r=p_eqmstart.r;
    ParamsStart.w=p_eqmstart.w;
    [p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_dc, n_a_GE, n_z, 0, pi_z, d_gridc, a_grid_GE, z_grid, ReturnFn, FnsToEvaluateGE, GeneralEqmEqnsGE, ParamsStart, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptions, vfoptions);
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
    fprintf('%s: (8) the stationary eqm is solved to 1e-3, so the checks below are not vacuous, should be one: %i \n',lab,GEfloor<=1e-3)
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
