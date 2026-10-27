function output=CoreFHorzTPath_localsearch(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params)
% vfoptions.localsearch=1 as a VALUE FUNCTION STEP, finite horizon, slowOLG (transpathoptions.fastOLG=0).
% No general eqm and no figures: these are direct calls to the single-step dispatcher
% ValueFnIter_FHorz_TPath_SingleStep, on all four corners -- without/with d, without/with the grid
% interpolation layer. The general eqm path is doPart(10).
%
% The FHorz counterpart of CoreInfHorzTPathTests doPart(4)/(5), minus the timing (settled on InfHorz).
% What the age dimension adds is that every check runs over every age at once, the terminal age j=N_j
% (no continuation value) included.
%
% THE STRONG CHECK IS (5): a restricted sweep at nlocalsearch=1, centred on the reference the seed sweep
% produced from the SAME next-period V, must reproduce the unrestricted V EXACTLY. The window then always
% contains the unrestricted argmax (also under GI, where the seed reference is the clamped midpoint and the
% two clamps coincide), so the restriction cannot bind. Because that identity would also hold if the
% window were silently ignored, (6) is the DISCRIMINATOR: the same restricted sweep centred on the wrong
% reference (the dispatcher's stay-put default) must give a DIFFERENT V.

N_a=prod(n_a);
N_z=prod(n_z);
nC=(N_a-1)/2; % the cap: 2*nC+1=N_a, so the window is the whole grid (n_a is odd here)

z_gridvals_J=repmat(gpuArray(z_grid),1,1,N_j);   % [N_z,l_z,N_j]
pi_z_J=repmat(gpuArray(pi_z),1,1,N_j);           % [N_z,N_z,N_j], slowOLG reads pi_z_J(:,:,jj)
a_grid=gpuArray(a_grid);
d_gridvals=gpuArray(d_grid);

Params.r=0.05;
Params.w=1;
DiscountFactorParamNames={'beta'};

cornername={'nod, GI0','nod, GI1','d, GI0','d, GI1'};
for cc=1:4
    hasd=(cc>=3);
    useGI=(cc==2 || cc==4);
    if hasd
        n_d_cc=n_d;
        N_dref=prod(n_d);
        ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        ReturnFnParamNames={'r','w','kappa_j','sigma','eta','varphi','agej','Jr','pension'};
    else
        n_d_cc=0;
        N_dref=1;
        ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
        ReturnFnParamNames={'r','w','kappa_j','sigma','agej','Jr','pension'};
    end
    lab=['localsearch FHorz step (',cornername{cc},')'];

    vfoptions=struct();
    vfoptions.exoticpreferences='None';
    vfoptions.divideandconquer=0;
    vfoptions.gridinterplayer=useGI;
    vfoptions.ngridinterp=5;
    vfoptions.lowmemory=0;
    vfoptions.verbose=0;
    vfoptions.localsearch=0;
    vfoptions.localsearchseed=0;

    % A non-trivial next-period V: one standard sweep from zeros
    [Vnext,~]=ValueFnIter_FHorz_TPath_SingleStep(zeros(N_a,N_z,N_j,'gpuArray'),n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,[],vfoptions);

    %% (1) The standard sweep, and the SEED sweep, which must agree with it
    [Vstd,Pstd]=ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,[],vfoptions);
    vfoptionsS=vfoptions;
    vfoptionsS.localsearchseed=1;
    [Vseed,Pseed,Ref]=ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,[],vfoptionsS);
    fprintf('%s: (1) seed sweep against standard sweep, V, this should be zero: %.3e \n',lab,gather(max(abs(Vseed(:)-Vstd(:)))))
    if hasd && ~useGI
        % LSseed_raw maxes over aprime per d then over d, the folded raw over (d,aprime) in one go:
        % same value, tie-breaking can differ. Reported, not checked.
        fprintf('%s: (1) seed sweep against standard sweep, Policy entries that differ (ties only): %i of %i \n',lab,gather(sum(Pseed(:)~=Pstd(:))),numel(Pstd))
    else
        % the seed sweep IS the standard raw here
        fprintf('%s: (1) seed sweep against standard sweep, Policy, this should be zero: %.3e \n',lab,gather(max(abs(Pseed(:)-Pstd(:)))))
    end
    fprintf('%s: (1) seed reference is [N_dref,N_a,N_z,N_j], should be one: %i \n',lab,isequal(size(Ref),[N_dref,N_a,N_z,N_j]))
    fprintf('%s: (1) seed reference is a coarse index in [1,N_a], should be one: %i \n',lab,gather(all(Ref(:)>=1) && all(Ref(:)<=N_a) && all(Ref(:)==round(Ref(:)))))
    if ~hasd && ~useGI
        % nothing folded and no interpolation layer, so the reference IS Policy's aprime
        fprintf('%s: (1) seed reference against the standard Policy, this should be zero: %.3e \n',lab,gather(max(abs(Ref(:)-Pstd(:)))))
    end

    %% (2) Restricted at the cap, where the window is the whole grid
    vfoptionsL=vfoptions;
    vfoptionsL.localsearch=1;
    vfoptionsL.nlocalsearch=nC;
    [Vcap,Pcap,RefCap]=ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,Ref,vfoptionsL);
    fprintf('%s: (2) restricted at the cap nlocalsearch=%i, V against the standard sweep, this should be zero: %.3e \n',lab,nC,gather(max(abs(Vcap(:)-Vstd(:)))))
    if hasd && useGI
        % LS1_GI1 takes layer 2 per d and then d; GI1_raw folds (d,L2). Same value, tie-breaking can differ.
        fprintf('%s: (2) restricted at the cap, Policy entries that differ from the standard sweep (ties only): %i of %i \n',lab,gather(sum(Pcap(:)~=Pstd(:))),numel(Pstd))
    else
        % same search order as the seed sweep, so identical including ties
        fprintf('%s: (2) restricted at the cap, Policy against the seed sweep, this should be zero: %.3e \n',lab,gather(max(abs(Pcap(:)-Pseed(:)))))
    end
    if ~useGI
        % Without GI the restricted argmax over the whole grid IS the seed's argmax, per d
        fprintf('%s: (2) restricted at the cap, new reference against the seed reference, this should be zero: %.3e \n',lab,gather(max(abs(RefCap(:)-Ref(:)))))
    end

    %% (3) The new reference has the reference's shape, and stays inside the window it was searched in
    vfoptionsL.nlocalsearch=1;
    [V1,P1,Ref1]=ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,Ref,vfoptionsL);
    fprintf('%s: (3) new reference is [N_dref,N_a,N_z,N_j], should be one: %i \n',lab,isequal(size(Ref1),[N_dref,N_a,N_z,N_j]))
    lo=min(max(Ref-1,1),N_a-2);
    fprintf('%s: (3) new reference lies inside its window [loweredge,loweredge+2], should be one: %i \n',lab,gather(all(Ref1(:)>=lo(:)) && all(Ref1(:)<=lo(:)+2)))

    %% (4) [placeholder number kept free, the InfHorz bank's (4) was timing]

    %% (5) THE STRONG CHECK: centred on the exact seed reference, nlocalsearch=1 is exact
    fprintf('%s: (5) restricted at nlocalsearch=1 around the exact seed reference, V against the standard sweep, this should be zero: %.3e \n',lab,gather(max(abs(V1(:)-Vstd(:)))))
    fprintf('%s: (5) restricted at nlocalsearch=1 around the exact seed reference, Policy entries that differ from the standard sweep (ties only): %i of %i \n',lab,gather(sum(P1(:)~=Pstd(:))),numel(Pstd))

    %% (6) THE DISCRIMINATOR: around the WRONG reference the window must bind
    % The dispatcher's default reference is the current a index (stay put). On this life-cycle model the
    % young save and the old dissave, so a window of three points around staying put cannot contain the
    % optimum everywhere: V has to move. Without this, (2) and (5) would pass even if the window were
    % ignored entirely.
    [Vw,Pw,Refw]=ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,[],vfoptionsL);
    fprintf('%s: (6) restricted around the stay-put reference, V against the standard sweep, should NOT be zero: %.3e \n',lab,gather(max(abs(Vw(:)-Vstd(:)))))
    if ~useGI
        % Without GI the restricted search set is a subset of the full one, so V can only fall.
        % (Under GI it is not: the polished sets do not nest, see the InfHorz memory notes.)
        fprintf('%s: (6) restricted around the stay-put reference, V is never above the standard sweep, this should be zero: %.3e \n',lab,gather(max(max(Vw(:)-Vstd(:)),0)))
    end
    lo=min(max(repmat((1:1:N_a),[N_dref,1,N_z,N_j])-1,1),N_a-2);
    fprintf('%s: (6) new reference lies inside the stay-put window, should be one: %i \n',lab,gather(all(Refw(:)>=lo(:)) && all(Refw(:)<=lo(:)+2)))

    %% (7) GI: Policy's coarse channel stays at or below N_a-1, which the agent distribution needs
    % It forms the upper grid point as P1+1. A value-fn check cannot see this, as it compares the decoded
    % fine point. The interior clamp is what guarantees it.
    if useGI
        if hasd
            crow=2; % (d, lower coarse, L2, L2flag)
        else
            crow=1; % (lower coarse, L2, L2flag)
        end
        fprintf('%s: (7) GI coarse channel within [1,N_a-1] at nlocalsearch=1 and around the stay-put reference, should be one: %i \n',lab,gather(all(P1(crow,:)>=1) && all(P1(crow,:)<=N_a-1) && all(Pw(crow,:)>=1) && all(Pw(crow,:)<=N_a-1)))
    end

    %% (8) The dispatcher's own refusals
    err8=0;
    vfoptionsE=vfoptionsL;
    vfoptionsE.localsearch=2;
    try
        ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,Ref,vfoptionsE);
    catch ME
        err8=contains(ME.message,'must be 0 or 1 here');
    end
    fprintf('%s: (8) localsearch=2 reaching the dispatcher is refused by name, should be one: %i \n',lab,err8)
    err8=0;
    try
        % a reference with one age too few
        ValueFnIter_FHorz_TPath_SingleStep(Vnext,n_d_cc,n_a,n_z,N_j,d_gridvals,a_grid,z_gridvals_J,pi_z_J,ReturnFn,Params,DiscountFactorParamNames,ReturnFnParamNames,Ref(:,:,:,1:N_j-1),vfoptionsL);
    catch ME
        err8=contains(ME.message,'aprimeReferencePolicy must hold one aprime index');
    end
    fprintf('%s: (8) a reference of the wrong shape is refused by name, should be one: %i \n',lab,err8)
end

output=1;

end
