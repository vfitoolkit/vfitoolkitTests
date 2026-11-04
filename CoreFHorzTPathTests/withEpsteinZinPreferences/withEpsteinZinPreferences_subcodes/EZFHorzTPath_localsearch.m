function output=EZFHorzTPath_localsearch(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params)
% vfoptions.localsearch=1 with Epstein-Zin preferences, as a VALUE FUNCTION STEP. The Epstein-Zin mirror of
% CoreFHorzTPath_localsearch (CoreFHorzTPathTests doPart(9)): no general eqm and no figures, direct calls to the
% single-step dispatchers, on every combination of
%   EZ case:  consumption units, utility units positive, utility units negative
%   layout:   slowOLG (transpathoptions.fastOLG=0) and fastOLG (=1)
%   shocks:   z, no z, z and e, e without z, and z with n_z=1, z and e with n_z=n_e=1
%   corner:   without/with d, without/with the grid interpolation layer, and d with n_d=1 without/with GI
% 3 x 12 x 6 = 216 cases, all with survival probabilities (sj(N_j)=0) and the warm-glow of bequests, so the
% terminal-age composition and the certainty equivalent of the continuation are both inside the window.
% The general eqm path is EZFHorzTPath_localsearch_shooting.
% The grids are built by ExogShockSetup_FHorz_TPath, exactly as TransitionPath_Case1_FHorz builds them,
% so each dispatcher sees the shapes it sees on a real path; V is shaped as that command shapes V_final.
% vfoptions is filled in here with everything the top-level commands would otherwise default (EZoneminusbeta,
% WarmGlowBequestsFnParamsNames), as the dispatchers are called directly.
%
% THE STRONG CHECK IS (5): a restricted sweep at nlocalsearch=1, centred on the reference the seed sweep
% produced from the SAME next-period V, must reproduce the unrestricted V EXACTLY. The local-search raws form
% the Epstein-Zin right-hand side exactly as the unrestricted raws do, so a window holding the argmax gives the
% same bits. (6) is the DISCRIMINATOR: centred on the wrong reference (the dispatcher's stay-put default) the
% window must bind and give a DIFFERENT V. (9) is the divide-and-conquer SEED route. (10) is lowmemory on the
% fastOLG seed sweeps.
% (11) and (12), added 2026-11-03 with the fix to the chosen-d lookups (which broke when n_d=1 made the per-d
% arrays vectors): n_d=1 against the no-d model with that d as a parameter, V, Policy and seed reference exactly;
% and lowmemory on the slowOLG seed sweeps (reached only by a direct call). Before the fix the n_d=1 corners with
% n_z=1 (ss=5), n_z=n_e=1 (ss=6) and (12) stopped the run with a dimension error.
% References are compared FLATTENED (N_dref-by-the rest): each raw keeps its own Policy dim order.

N_a=prod(n_a);
% The shocks as passed in; the single-point cases (ss=5,6) replace them inside the loop
n_z0=n_z; z_grid0=z_grid; pi_z0=pi_z;
n_e0=n_e; e_grid0=e_grid; pi_e0=pi_e;
nC=(N_a-1)/2; % the cap: 2*nC+1=N_a, so the window is the whole grid (n_a is odd here)
a_grid=gpuArray(a_grid);
Params.r=0.05;
Params.w=1;
Params.dfix=d_grid(ceil(n_d/2)); % for (11): the n_d=1 models' single d, as a parameter of the no-d model
DiscountFactorParamNames={'beta'};
level1n=floor(sqrt(N_a));

casename={'cons units','positive utils','negative utils'};
casefn={'cons','positiveUtils','negativeUtils'};
layoutname={'slowOLG','fastOLG'};
shockname={'z','noz','e','noz_e','z with n_z=1','e with n_z=1 and n_e=1'};
cornername={'nod, GI0','nod, GI1','d, GI0','d, GI1','d with n_d=1, GI0','d with n_d=1, GI1'};
for ezcase=1:3
    for ll=1:2
        fastOLG=ll-1;
        for ss=1:6
            % ss=5,6 are the shapes of ss=1,3 with single-point shocks: with n_d=1 too, every per-d array the raws
            % read at the chosen d is then a vector, the case the lookups must reshape (2026-11-03)
            sh=ss; % the shape
            n_z=n_z0; z_grid=z_grid0; pi_z=pi_z0; n_e=n_e0; e_grid=e_grid0; pi_e=pi_e0;
            if ss==5
                sh=1;
                n_z=1; z_grid=z_grid0(ceil(n_z0/2)); pi_z=1;
            elseif ss==6
                sh=3;
                n_z=1; z_grid=z_grid0(ceil(n_z0/2)); pi_z=1;
                n_e=1; e_grid=e_grid0(ceil(n_e0/2)); pi_e=1;
            end
            N_z=prod(n_z);
            N_e=prod(n_e);
            hasz=(sh==1 || sh==3);
            hase=(sh==3 || sh==4);
            % The exogenous shocks, set up exactly as TransitionPath_Case1_FHorz does
            transpathoptions=struct();
            transpathoptions.fastOLG=fastOLG;
            simoptions=struct();
            simoptions.fastOLG=1; % TransitionPath_Case1_FHorz's default
            if hase
                simoptions.n_e=n_e;
                simoptions.e_grid=e_grid;
                simoptions.pi_e=pi_e;
            end
            if hasz
                n_z_ss=n_z; z_grid_ss=z_grid; pi_z_ss=pi_z;
            else
                n_z_ss=0; z_grid_ss=[]; pi_z_ss=[];
            end
            [z_gridvals_J, pi_z_J, ~, e_gridvals_J, pi_e_J]=ExogShockSetup_FHorz_TPath(n_z_ss,z_grid_ss,pi_z_ss,N_a,N_j,T,Params,{},[],{},[],transpathoptions,simoptions,4);
            % V as TransitionPath_Case1_FHorz shapes V_final
            if fastOLG==0
                if sh==1
                    V0=zeros(N_a,N_z,N_j,'gpuArray');
                elseif sh==2
                    V0=zeros(N_a,N_j,'gpuArray');
                elseif sh==3
                    V0=zeros(N_a,N_z,N_e,N_j,'gpuArray');
                else
                    V0=zeros(N_a,N_e,N_j,'gpuArray');
                end
            else
                if sh==1
                    V0=zeros(N_a*N_j,N_z,'gpuArray');
                elseif sh==2
                    V0=zeros(N_a,N_j,'gpuArray');
                elseif sh==3
                    V0=zeros(N_a*N_j,N_z,N_e,'gpuArray');
                else
                    V0=zeros(N_a*N_j,N_e,'gpuArray');
                end
            end
            N_rest=N_a*N_j*max(N_z*hasz,1)*max(N_e*hase,1);

            for cc=1:6
                hasd=(cc>=3);
                useGI=(cc==2 || cc==4 || cc==6);
                lab=['EZ localsearch FHorz step (',casename{ezcase},', ',layoutname{ll},', ',shockname{ss},', ',cornername{cc},')'];
                if cc>=5
                    % a single d: the middle point of d_grid
                    n_d_cc=1;
                    N_dref=1;
                    d_gridvals=gpuArray(d_grid(ceil(n_d/2)));
                elseif hasd
                    n_d_cc=n_d;
                    N_dref=prod(n_d);
                    d_gridvals=gpuArray(d_grid);
                else
                    n_d_cc=0;
                    N_dref=1;
                    d_gridvals=[];
                end
                % The return fn for this model, from the bank's EZ set: the states, then the parameters
                statenames={};
                if hasd
                    statenames{end+1}='d';
                end
                statenames=[statenames,{'aprime','a'}];
                if hasz
                    statenames{end+1}='z';
                end
                if hase
                    statenames{end+1}='e';
                end
                ReturnFnParamNames={'r','w','kappa_j'};
                if ezcase>1
                    ReturnFnParamNames{end+1}='ezsigma';
                end
                if hasd
                    ReturnFnParamNames{end+1}='varphi';
                end
                ReturnFnParamNames=[ReturnFnParamNames,{'agej','Jr','pension'}];
                if hasd
                    dname='d';
                else
                    dname='nod';
                end
                if hasz
                    zname='z';
                else
                    zname='noz';
                end
                if hase
                    ename='e';
                else
                    ename='noe';
                end
                arglist=strjoin([statenames,ReturnFnParamNames],',');
                ReturnFn=str2func(['@(',arglist,') EZReturnFn_',casefn{ezcase},'_',dname,'_',zname,'_',ename,'_nosemiz(',arglist,')']);
                % The dispatcher for this layout and these shocks, with its own argument list
                if fastOLG==0
                    if sh==1
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep(V,nd,n_a,n_z,N_j,dg,a_grid,z_gridvals_J,pi_z_J,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    elseif sh==2
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_noz(V,nd,n_a,N_j,dg,a_grid,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    elseif sh==3
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_e(V,nd,n_a,n_z,n_e,N_j,dg,a_grid,z_gridvals_J,e_gridvals_J,pi_z_J,pi_e_J,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    else
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_noz_e(V,nd,n_a,n_e,N_j,dg,a_grid,e_gridvals_J,pi_e_J,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    end
                else
                    if sh==1
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_fastOLG(V,nd,n_a,n_z,N_j,dg,a_grid,z_gridvals_J,pi_z_J,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    elseif sh==2
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_fastOLG_noz(V,nd,n_a,N_j,dg,a_grid,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    elseif sh==3
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_fastOLG_e(V,nd,n_a,n_z,n_e,N_j,dg,a_grid,z_gridvals_J,e_gridvals_J,pi_z_J,pi_e_J,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    else
                        disp_fnX=@(V,ref,vfo,nd,dg,RF,RFN) ValueFnIter_FHorz_TPath_SingleStep_fastOLG_noz_e(V,nd,n_a,n_e,N_j,dg,a_grid,e_gridvals_J,pi_e_J,RF,Params,DiscountFactorParamNames,RFN,ref,vfo);
                    end
                end

                disp_fn=@(V,ref,vfo) disp_fnX(V,ref,vfo,n_d_cc,d_gridvals,ReturnFn,ReturnFnParamNames);


                vfoptions=struct();
                vfoptions.exoticpreferences='EpsteinZin';
                if ezcase==1
                    vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
                    vfoptions.WarmGlowBequestsFn=@(aprime,wg1,wg2) EZWarmGlowFn_cons(aprime,wg1,wg2);
                    vfoptions.WarmGlowBequestsFnParamsNames={'wg1','wg2'};
                elseif ezcase==2
                    vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
                    vfoptions.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_positiveUtils(aprime,wg1,wg2,wg3);
                    vfoptions.WarmGlowBequestsFnParamsNames={'wg1','wg2','wg3'};
                else
                    vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';
                    vfoptions.WarmGlowBequestsFn=@(aprime,wg1,wg2,wg3) EZWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);
                    vfoptions.WarmGlowBequestsFnParamsNames={'wg1','wg2','wg3'};
                end
                vfoptions.EZoneminusbeta=0; % the top-level commands' default
                vfoptions.survivalprobability='sj';
                vfoptions.divideandconquer=0;
                vfoptions.level1n=level1n;
                vfoptions.gridinterplayer=useGI;
                vfoptions.ngridinterp=5;
                vfoptions.lowmemory=0;
                vfoptions.verbose=0;
                vfoptions.EVpre=0; % as every transition path command sets it
                vfoptions.localsearch=0;
                vfoptions.localsearchseed=0;

                % A non-trivial next-period V: one standard sweep from zeros
                [Vnext,~]=disp_fn(V0,[],vfoptions);

                %% (1) The standard sweep, and the SEED sweep, which must agree with it
                [Vstd,Pstd]=disp_fn(Vnext,[],vfoptions);
                vfoptionsS=vfoptions;
                vfoptionsS.localsearchseed=1;
                [Vseed,Pseed,Ref]=disp_fn(Vnext,[],vfoptionsS);
                fprintf('%s: (1) seed sweep against standard sweep, V, this should be zero: %.3e \n',lab,gather(max(abs(Vseed(:)-Vstd(:)))))
                if hasd && ~useGI
                    % The seed raw maxes over aprime per d then over d, the folded raw over (d,aprime) in one go:
                    % same value, tie-breaking can differ. Reported, not checked.
                    fprintf('%s: (1) seed sweep against standard sweep, Policy entries that differ (ties only): %i of %i \n',lab,gather(sum(Pseed(:)~=Pstd(:))),numel(Pstd))
                else
                    % the seed sweep IS the standard raw here
                    fprintf('%s: (1) seed sweep against standard sweep, Policy, this should be zero: %.3e \n',lab,gather(max(abs(Pseed(:)-Pstd(:)))))
                end
                fprintf('%s: (1) seed reference holds N_dref-by-N_rest entries, should be one: %i \n',lab,size(Ref,1)==N_dref && numel(Ref)==N_dref*N_rest)
                fprintf('%s: (1) seed reference is a coarse index in [1,N_a], should be one: %i \n',lab,gather(all(Ref(:)>=1) && all(Ref(:)<=N_a) && all(Ref(:)==round(Ref(:)))))
                Ref=reshape(Ref,[N_dref,N_rest]);

                %% (2) Restricted at the cap, where the window is the whole grid
                vfoptionsL=vfoptions;
                vfoptionsL.localsearch=1;
                vfoptionsL.nlocalsearch=nC;
                [Vcap,Pcap,RefCap]=disp_fn(Vnext,Ref,vfoptionsL);
                fprintf('%s: (2) restricted at the cap nlocalsearch=%i, V against the standard sweep, this should be zero: %.3e \n',lab,nC,gather(max(abs(Vcap(:)-Vstd(:)))))
                if hasd && useGI
                    % per-d layer 2 then d, against GI1's folded (d,L2): same value, tie-breaking can differ
                    fprintf('%s: (2) restricted at the cap, Policy entries that differ from the standard sweep (ties only): %i of %i \n',lab,gather(sum(Pcap(:)~=Pstd(:))),numel(Pstd))
                else
                    fprintf('%s: (2) restricted at the cap, Policy against the seed sweep, this should be zero: %.3e \n',lab,gather(max(abs(Pcap(:)-Pseed(:)))))
                end
                if ~useGI
                    fprintf('%s: (2) restricted at the cap, new reference against the seed reference, this should be zero: %.3e \n',lab,gather(max(abs(RefCap(:)-Ref(:)))))
                end

                %% (3) The new reference stays inside the window it was searched in
                vfoptionsL.nlocalsearch=1;
                [V1,P1,Ref1]=disp_fn(Vnext,Ref,vfoptionsL);
                fprintf('%s: (3) new reference holds N_dref-by-N_rest entries, should be one: %i \n',lab,size(Ref1,1)==N_dref && numel(Ref1)==N_dref*N_rest)
                lo=min(max(Ref(:)-1,1),N_a-2);
                fprintf('%s: (3) new reference lies inside its window [loweredge,loweredge+2], should be one: %i \n',lab,gather(all(Ref1(:)>=lo) && all(Ref1(:)<=lo+2)))

                %% (5) THE STRONG CHECK: centred on the exact seed reference, nlocalsearch=1 is exact
                fprintf('%s: (5) restricted at nlocalsearch=1 around the exact seed reference, V against the standard sweep, this should be zero: %.3e \n',lab,gather(max(abs(V1(:)-Vstd(:)))))
                fprintf('%s: (5) restricted at nlocalsearch=1 around the exact seed reference, Policy entries that differ from the standard sweep (ties only): %i of %i \n',lab,gather(sum(P1(:)~=Pstd(:))),numel(Pstd))

                %% (6) THE DISCRIMINATOR: around the WRONG reference the window must bind
                [Vw,Pw,Refw]=disp_fn(Vnext,[],vfoptionsL);
                fprintf('%s: (6) restricted around the stay-put reference, V against the standard sweep, should NOT be zero: %.3e \n',lab,gather(max(abs(Vw(:)-Vstd(:)))))
                if ~useGI
                    % Without GI the restricted search set is a subset of the full one, and ^ezc7 is increasing
                    % in the right-hand side (ezc7>0 in all three cases here), so V can only fall
                    fprintf('%s: (6) restricted around the stay-put reference, V is never above the standard sweep, this should be zero: %.3e \n',lab,gather(max(max(Vw(:)-Vstd(:)),0)))
                end
                lo=reshape(min(max(repmat((1:1:N_a),[N_dref,1,N_rest/N_a])-1,1),N_a-2),[],1);
                fprintf('%s: (6) new reference lies inside the stay-put window, should be one: %i \n',lab,gather(all(Refw(:)>=lo) && all(Refw(:)<=lo+2)))

                %% (7) GI: Policy's coarse channel stays at or below N_a-1, which the agent distribution needs
                if useGI
                    if hasd
                        crow=2; % (d, lower coarse, L2, L2flag)
                    else
                        crow=1; % (lower coarse, L2, L2flag)
                    end
                    fprintf('%s: (7) GI coarse channel within [1,N_a-1] at nlocalsearch=1 and around the stay-put reference, should be one: %i \n',lab,gather(all(P1(crow,:)>=1) && all(P1(crow,:)<=N_a-1) && all(Pw(crow,:)>=1) && all(Pw(crow,:)<=N_a-1)))
                end

                %% (8) A reference of the wrong size is refused by name (the Epstein-Zin router's own check)
                err8=0;
                try
                    disp_fn(Vnext,Ref(:,1:end-1),vfoptionsL); % one entry too few per d
                catch ME
                    err8=contains(ME.message,'aprimeReferencePolicy must hold one aprime index');
                end
                fprintf('%s: (8) a reference of the wrong size is refused by name, should be one: %i \n',lab,err8)

                %% (9) The divide-and-conquer SEED route
                % A seed sweep with divideandconquer=1 goes through the EZ DC raws, which return the same object.
                vfoptionsD=vfoptionsS;
                vfoptionsD.divideandconquer=1;
                [Vdc,~,Refdc]=disp_fn(Vnext,[],vfoptionsD);
                fprintf('%s: (9) DC seed sweep against standard sweep, V, this should be zero: %.3e \n',lab,gather(max(abs(Vdc(:)-Vstd(:)))))
                fprintf('%s: (9) DC seed reference against the non-DC seed reference, entries that differ (ties only): %i of %i \n',lab,gather(sum(Refdc(:)~=Ref(:))),numel(Ref))
                % and it serves as a reference exactly as the non-DC one does. The restricted sweep IGNORES DC.
                vfoptionsLD=vfoptionsL;
                vfoptionsLD.divideandconquer=1;
                [V1dc,~,~]=disp_fn(Vnext,reshape(Refdc,[N_dref,N_rest]),vfoptionsLD);
                fprintf('%s: (9) restricted at nlocalsearch=1 around the DC seed reference, V against the standard sweep, this should be zero: %.3e \n',lab,gather(max(abs(V1dc(:)-Vstd(:)))))

                %% (10) lowmemory on the fastOLG seed sweeps (slowOLG transition paths refuse lowmemory>0)
                if fastOLG==1 && (hasz || hase)
                    vfoptionsM=vfoptionsS;
                    vfoptionsM.lowmemory=1;
                    [Vlm,~,Reflm]=disp_fn(Vnext,[],vfoptionsM);
                    fprintf('%s: (10) seed sweep with lowmemory=1 against lowmemory=0, V, this should be zero: %.3e \n',lab,gather(max(abs(Vlm(:)-Vseed(:)))))
                    fprintf('%s: (10) seed sweep with lowmemory=1 against lowmemory=0, reference, this should be zero: %.3e \n',lab,gather(max(abs(Reflm(:)-Ref(:)))))
                end

                %% (11) n_d=1 against the no-d model with that d as a parameter: the same problem, so the same answer.
                % A discriminator for the n_d=1 corners, not just a no-crash test: a lookup that read the wrong d or the
                % wrong state would give a different Policy or reference.
                if cc>=5
                    nodstates={'aprime','a'};
                    zn='noz'; en='noe';
                    if hasz
                        nodstates{end+1}='z'; zn='z';
                    end
                    if hase
                        nodstates{end+1}='e'; en='e';
                    end
                    fname=['EZReturnFn_',casefn{ezcase},'_d_',zn,'_',en,'_nosemiz'];
                    ReturnFnFix=str2func(['@(',strjoin([nodstates,ReturnFnParamNames,{'dfix'}],','),') ',fname,'(dfix,',strjoin([nodstates,ReturnFnParamNames],','),')']);
                    [Vn,Pn]=disp_fnX(Vnext,[],vfoptions,0,[],ReturnFnFix,[ReturnFnParamNames,{'dfix'}]);
                    [~,~,Refn]=disp_fnX(Vnext,[],vfoptionsS,0,[],ReturnFnFix,[ReturnFnParamNames,{'dfix'}]);
                    fprintf('%s: (11) n_d=1 against the no-d model with d as a parameter, V, this should be zero: %.3e \n',lab,gather(max(abs(Vstd(:)-Vn(:)))))
                    % Policy: with d it has a leading d row (always 1) where the no-d Policy has none, or folds d in (d+N_d*(aprime-1)=aprime)
                    Pd=reshape(Pstd,size(Pstd,1),[]);
                    Pn=reshape(Pn,size(Pn,1),[]);
                    if size(Pd,1)==size(Pn,1)+1
                        fprintf('%s: (11) n_d=1, the d row of Policy is all ones, should be one: %i \n',lab,gather(all(Pd(1,:)==1)))
                        Pd=Pd(2:end,:);
                    end
                    fprintf('%s: (11) n_d=1 against the no-d model with d as a parameter, Policy, this should be zero: %.3e \n',lab,gather(max(abs(Pd(:)-Pn(:)))))
                    fprintf('%s: (11) n_d=1 against the no-d model with d as a parameter, seed reference, this should be zero: %.3e \n',lab,gather(max(abs(Ref(:)-Refn(:)))))
                end

                %% (12) lowmemory on the slowOLG seed sweeps with d, without GI or DC (the LSseed raws' per-e and per-(z,e) loops)
                % Reached only by a direct call: TransitionPath_Case1_FHorz refuses lowmemory>0 with slowOLG, and the seed sweep is
                % only ever run by it. The loops are kept (as in the stationary raws) and must agree with the whole sweep.
                if fastOLG==0 && hasd && ~useGI && hase
                    if hasz
                        lowmemlist=[1,2];
                    else
                        lowmemlist=1;
                    end
                    for lm=lowmemlist
                        vfoptionsM=vfoptionsS;
                        vfoptionsM.lowmemory=lm;
                        [Vlm,Plm,Reflm]=disp_fn(Vnext,[],vfoptionsM);
                        fprintf('%s: (12) slowOLG seed sweep with lowmemory=%i against lowmemory=0, V, this should be zero: %.3e \n',lab,lm,gather(max(abs(Vlm(:)-Vseed(:)))))
                        fprintf('%s: (12) slowOLG seed sweep with lowmemory=%i against lowmemory=0, Policy, this should be zero: %.3e \n',lab,lm,gather(max(abs(Plm(:)-Pseed(:)))))
                        fprintf('%s: (12) slowOLG seed sweep with lowmemory=%i against lowmemory=0, reference, this should be zero: %.3e \n',lab,lm,gather(max(abs(Reflm(:)-Ref(:)))))
                    end
                end
            end
        end
    end
end

output=1;

end
