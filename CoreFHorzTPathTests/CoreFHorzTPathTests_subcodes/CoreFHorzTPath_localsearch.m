function output=CoreFHorzTPath_localsearch(T,n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params)
% vfoptions.localsearch=1 as a VALUE FUNCTION STEP. No general eqm and no figures: direct calls to the
% single-step dispatchers, on every combination of
%   layout:   slowOLG (transpathoptions.fastOLG=0) and fastOLG (=1)
%   shocks:   z, no z, z and e, e without z, and z with n_z=1, z and e with n_z=n_e=1
%   corner:   without/with d, without/with the grid interpolation layer, and d with n_d=1 without/with GI
% 12 x 6 = 72 cases. The general eqm path is doPart(10).
% The grids are built by ExogShockSetup_FHorz_TPath, exactly as TransitionPath_Case1_FHorz builds them,
% so each dispatcher sees the shapes it sees on a real path; V is shaped as that command shapes V_final.
%
% THE STRONG CHECK IS (5): a restricted sweep at nlocalsearch=1, centred on the reference the seed sweep
% produced from the SAME next-period V, must reproduce the unrestricted V EXACTLY. The window then always
% contains the unrestricted argmax (also under GI, where the seed reference is the clamped midpoint and the
% two clamps coincide). Because that identity would also hold if the window were silently ignored, (6) is
% the DISCRIMINATOR: the same restricted sweep centred on the wrong reference (the dispatcher's stay-put
% default) must give a DIFFERENT V.
% (9) is the divide-and-conquer SEED route: a DC seed sweep must give the unrestricted V and a reference
% that makes (5) hold again. (10) is lowmemory on the fastOLG seed sweeps.
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

layoutname={'slowOLG','fastOLG'};
shockname={'z','noz','e','noz_e','z with n_z=1','e with n_z=1 and n_e=1'};
cornername={'nod, GI0','nod, GI1','d, GI0','d, GI1','d with n_d=1, GI0','d with n_d=1, GI1'};
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
            lab=['localsearch FHorz step (',layoutname{ll},', ',shockname{ss},', ',cornername{cc},')'];
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
            % The return fn for this model, from the bank's set
            if hasd
                if sh==1
                    ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                elseif sh==2
                    ReturnFn=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                elseif sh==3
                    ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                else
                    ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
                end
                ReturnFnParamNames={'r','w','kappa_j','sigma','eta','varphi','agej','Jr','pension'};
            else
                if sh==1
                    ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
                elseif sh==2
                    ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
                elseif sh==3
                    ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
                else
                    ReturnFn=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
                end
                ReturnFnParamNames={'r','w','kappa_j','sigma','agej','Jr','pension'};
            end
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
            vfoptions.exoticpreferences='None';
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
                % Without GI the restricted search set is a subset of the full one, so V can only fall
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

            %% (8) The dispatcher's own refusals
            err8=0;
            vfoptionsE=vfoptionsL;
            vfoptionsE.localsearch=2;
            try
                disp_fn(Vnext,Ref,vfoptionsE);
            catch ME
                err8=contains(ME.message,'must be 0 or 1 here');
            end
            fprintf('%s: (8) localsearch=2 reaching the dispatcher is refused by name, should be one: %i \n',lab,err8)
            err8=0;
            try
                disp_fn(Vnext,Ref(:,1:end-1),vfoptionsL); % one entry too few per d
            catch ME
                err8=contains(ME.message,'aprimeReferencePolicy must hold one aprime index');
            end
            fprintf('%s: (8) a reference of the wrong size is refused by name, should be one: %i \n',lab,err8)

            %% (9) The divide-and-conquer SEED route
            % A seed sweep with divideandconquer=1 goes through the DC raws, which return the same object.
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
            if fastOLG==1
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
                fname=['ReturnFn_d_',zn,'_',en,'_nosemiz'];
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

output=1;

end
