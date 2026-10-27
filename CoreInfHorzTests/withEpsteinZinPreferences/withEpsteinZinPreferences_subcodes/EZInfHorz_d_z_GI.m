function output=EZInfHorz_d_z_GI(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% InfHorz Epstein-Zin, with d, with z, WITH the grid interpolation layer (postGI, the only EZ GI path).
%
% Raw reached: ValueFnIter_InfHorz_EpsteinZin_Refine_postGI_raw (new 2026-10-27). d is refined out of
% the RAW return function first (the Epstein-Zin right-hand side is increasing in F for every sign
% case, with a continuation independent of d, so the d that maximises F maximises the objective),
% and the rest is the nod raw on the refined matrices. lowmemory=1 IS a second path here (the
% refined return matrices are built one z at a time), unlike for the nod GI raw.
%
% Policy has four rows: (d, L1, L2, L2flag). Same checks as EZInfHorz_nod_z_GI, with the (d,aprime)
% policy: the identity uses the composite good x=c^varphi*(1-d)^(1-varphi) at the chosen d and the
% linear aprime, and the two-point interpolated transformed expectation.

N_d=prod(n_d);
N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_cons=@(d,aprime,a,z,r,w,varphi) EZReturnFn_cons_d_z(d,aprime,a,z,r,w,varphi);
ReturnFn_posU=@(d,aprime,a,z,r,w,varphi,ezsigma) EZReturnFn_positiveUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma);
ReturnFn_negU=@(d,aprime,a,z,r,w,varphi,ezsigma) EZReturnFn_negativeUtils_d_z(d,aprime,a,z,r,w,varphi,ezsigma);

% Shorthands for the identity checks (plain CPU doubles throughout)
beta=Params.beta;
r=Params.r;
w=Params.w;
varphi=Params.varphi;
ezgamma=Params.ezgamma;
ezphi=Params.ezphi;
ezrisk=Params.ezrisk;
ezsigma=Params.ezsigma;
a_gridcpu=gather(a_grid);
d_gridcpu=gather(d_grid);
z_gridcpu=gather(z_grid);
pi_zcpu=gather(pi_z);
avals=repmat(a_gridcpu,1,N_z);
zvals=repmat(z_gridcpu',N_a,1);
zmid=ceil(N_z/2);

vftol=10^(-10); % see the driver header on why not the default 1e-9
n2short=5; % vfoptions.ngridinterp
N_afine=N_a+(N_a-1)*n2short;
subsamp=(n_a_big-1)/(n_a-1); % a_grid_big(1:subsamp:end)==a_grid

for ezcase=1:3
    vfoptions=struct();
    vfoptions.exoticpreferences='EpsteinZin';
    vfoptions.tolerance=vftol;
    vfoptions.gridinterplayer=1;
    vfoptions.ngridinterp=n2short;
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptions.EZutils=0;
        vfoptions.EZriskaversion='ezgamma';
        vfoptions.EZeis='ezphi';
        ReturnFn=ReturnFn_cons;
        Paramsref=Params;
        Paramsref.ezsigma=Params.ezgamma; % vNM with the same risk aversion
        ReturnFnref=ReturnFn_negU;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=1;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_posU;
        Paramsref=Params;
        ReturnFnref=ReturnFn_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=0;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn=ReturnFn_negU;
        Paramsref=Params;
        ReturnFnref=ReturnFn_negU;
    end
    fprintf('\n--- [EZ %s] with d, with z, grid interpolation layer --- \n',casestr)

    %% The vNM reference solve, with the grid interpolation layer
    vfref=struct(); % deliberately no exoticpreferences field
    vfref.tolerance=vftol;
    vfref.gridinterplayer=1;
    vfref.ngridinterp=n2short;
    [Vref,Policyref]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFnref,Paramsref,DiscountFactorParamNames,[],vfref);
    if ezcase==1
        Vrefcomp=((1-ezgamma)*Vref).^(1/(1-ezgamma)); % into consumption units
    else
        Vrefcomp=Vref;
    end
    PolicyValsref=PolicyInd2Val_InfHorz(Policyref,n_d,n_a,n_z,d_grid,a_grid,vfref);

    %% The Epstein-Zin solve with the grid interpolation layer
    [V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
    PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
    fprintf('[EZ %s] (GI) PolicyInd2Val ran \n',casestr)

    % Shapes and index ranges
    fprintf('[EZ %s] (GI) size(V) should be [%d %d], is [%s] \n',casestr,N_a,N_z,num2str(size(V1)))
    fprintf('[EZ %s] (GI) size(Policy) should be [4 %d %d], is [%s] \n',casestr,N_a,N_z,num2str(size(Policy1)))
    dindex=gather(reshape(Policy1(1,:,:),[N_a,N_z]));
    L1=gather(reshape(Policy1(2,:,:),[N_a,N_z]));
    L2=gather(reshape(Policy1(3,:,:),[N_a,N_z]));
    L2flag=gather(reshape(Policy1(4,:,:),[N_a,N_z]));
    fprintf('[EZ %s] (GI) d index in range, this should be one: %d \n',casestr,all(dindex(:)>=1) && all(dindex(:)<=N_d))
    fprintf('[EZ %s] (GI) L1 in [1,N_a-1], this should be one: %d \n',casestr,all(L1(:)>=1) && all(L1(:)<=N_a-1))
    fprintf('[EZ %s] (GI) L2 in [1,ngridinterp+2], this should be one: %d \n',casestr,all(L2(:)>=1) && all(L2(:)<=n2short+2))
    fprintf('[EZ %s] (GI) L2flag in {1,2,3}, this should be one: %d \n',casestr,all(ismember(L2flag(:),[1,2,3])))
    fprintf('[EZ %s] (GI) L2flag is 2 everywhere (no -Inf refined return in this model), this should be one: %d \n',casestr,all(L2flag(:)==2))
    fprintf('[EZ %s] (GI) combined fine index in [1,%d], this should be one: %d \n',casestr,N_afine,all((n2short+1)*(L1(:)-1)+L2(:)>=1) && all((n2short+1)*(L1(:)-1)+L2(:)<=N_afine))
    fprintf('[EZ %s] (GI) V is finite everywhere, this should be one: %d \n',casestr,gather(all(isfinite(V1(:)))))
    if ezcase==3
        fprintf('[EZ %s] (GI) V is negative everywhere (negative-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)<0)))
    else
        fprintf('[EZ %s] (GI) V is positive everywhere (consumption units / positive-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)>0)))
    end
    fprintf('[EZ %s] (GI) the L2 channel is not degenerate, number of distinct L2 chosen should be WELL ABOVE zero: %d \n',casestr,numel(unique(L2(:))))

    %% The Epstein-Zin recursion with the two-point interpolated continuation, recomputed here from the formula at (d,aprime)
    Vg=gather(V1);
    wu=(L2-1)/(n2short+1); % weight on the upper node
    aprimevals=a_gridcpu(L1)+(a_gridcpu(L1+1)-a_gridcpu(L1)).*wu; % linear in the weight
    dvals=d_gridcpu(dindex);
    fprintf('[EZ %s] (GI) reconstructed aprime matches PolicyInd2Val, this should be zero: %2.8e \n',casestr,max(abs(aprimevals(:)-gather(reshape(PolicyVals1(2,:,:),[N_a*N_z,1])))))
    fprintf('[EZ %s] (GI) reconstructed d matches PolicyInd2Val, this should be zero: %2.8e \n',casestr,max(abs(dvals(:)-gather(reshape(PolicyVals1(1,:,:),[N_a*N_z,1])))))
    c=(1+r)*avals+w*dvals.*zvals-aprimevals;
    ok=(c>0)&(dvals<1);
    fprintf('[EZ %s] (GI) chosen (d,aprime) always feasible (c>0, d<1), this should be one: %d \n',casestr,all(ok(:)))
    x=zeros(N_a,N_z);
    x(ok)=(c(ok).^varphi).*((1-dvals(ok)).^(1-varphi));
    Fpol=-Inf(N_a,N_z);
    if ezcase==1
        Fpol(ok)=x(ok);
    elseif ezcase==2
        Fpol(ok)=((1+x(ok)).^(1-ezsigma)-1)/(1-ezsigma);
    else
        Fpol(ok)=(x(ok).^(1-ezsigma))/(1-ezsigma);
    end
    if ezcase==1
        TV=Vg.^(1-ezgamma);
        TValt=Vg.^(1-1/ezphi); % gamma=1/phi, the vNM collapse point
    elseif ezcase==2
        TV=Vg.^(1-ezrisk);
        TValt=Vg; % risk=0
    else
        TV=(-Vg).^(1+ezrisk);
        TValt=-Vg;
    end
    EVpow=zeros(N_a,N_z);
    EVpowalt=zeros(N_a,N_z);
    for z_c=1:N_z
        EVpow(:,z_c)=((1-wu(:,z_c)).*TV(L1(:,z_c),:)+wu(:,z_c).*TV(L1(:,z_c)+1,:))*pi_zcpu(z_c,:)';
        EVpowalt(:,z_c)=((1-wu(:,z_c)).*TValt(L1(:,z_c),:)+wu(:,z_c).*TValt(L1(:,z_c)+1,:))*pi_zcpu(z_c,:)';
    end
    if ezcase==1
        rho=1-1/ezphi;
        alpha=1-ezgamma;
        RHS=(Fpol.^rho+beta*EVpow.^(rho/alpha)).^(1/rho);
        RHSalt=(Fpol.^rho+beta*EVpowalt).^(1/rho);
    elseif ezcase==2
        RHS=Fpol+beta*EVpow.^(1/(1-ezrisk));
        RHSalt=Fpol+beta*EVpowalt;
    else
        RHS=Fpol-beta*EVpow.^(1/(1+ezrisk));
        RHSalt=Fpol-beta*EVpowalt;
    end
    resid=Vg-RHS;
    fprintf('[EZ %s] (GI) Epstein-Zin recursion with the interpolated continuation, recomputed from the formula at (d,aprime), this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(resid(:))))
    residalt=Vg-RHSalt;
    fprintf('[EZ %s] (GI) the SAME identity at the vNM-collapse risk aversion, this should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(residalt(:))))

    %% The chosen d maximises the return at the chosen aprime (the refinement is valid under Epstein-Zin: the objective is increasing in F)
    % F over all d at the chosen aprime, for every (a,z); the chosen d must attain the max
    Fall=-Inf(N_d,N_a,N_z);
    for d_c=1:N_d
        call=(1+r)*avals+w*d_gridcpu(d_c)*zvals-aprimevals;
        okall=(call>0)&(d_gridcpu(d_c)<1);
        xall=zeros(N_a,N_z);
        xall(okall)=(call(okall).^varphi).*((1-d_gridcpu(d_c)).^(1-varphi));
        Fd=-Inf(N_a,N_z);
        if ezcase==1
            Fd(okall)=xall(okall);
        elseif ezcase==2
            Fd(okall)=((1+xall(okall)).^(1-ezsigma)-1)/(1-ezsigma);
        else
            Fd(okall)=(xall(okall).^(1-ezsigma))/(1-ezsigma);
        end
        Fall(d_c,:,:)=reshape(Fd,[1,N_a,N_z]);
    end
    Fbest=reshape(max(Fall,[],1),[N_a,N_z]);
    fprintf('[EZ %s] (GI) chosen d maximises the return at the chosen aprime, this should be zero: %2.8e \n',casestr,max(abs(Fpol(:)-Fbest(:))))

    %% Howards improvement is only an accelerator (this is the two-point Howards block of the new raw)
    vfoptions_noH=vfoptions;
    vfoptions_noH.howards=0;
    [V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_noH);
    fprintf('[EZ %s] (GI) howards=0 (pure VFI), V, this should be zero to convergence tolerance: %.3e \n',casestr,max(abs(V1(:)-V1noH(:))))
    fprintf('[EZ %s] (GI) howards=0 (pure VFI), Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1noH(:))))

    %% lowmemory=1 builds the refined return matrices one z at a time; same arithmetic, so same answer
    vfoptions_lm=vfoptions;
    vfoptions_lm.lowmemory=1;
    [V1lm,Policy1lm]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_lm);
    fprintf('[EZ %s] (GI) lowmemory=1, V, this should be zero: %.3e \n',casestr,max(abs(V1(:)-V1lm(:))))
    fprintf('[EZ %s] (GI) lowmemory=1, Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1lm(:))))

    %% The maxaprimediff window: with d the default is ceil(n_a/5)=21; against a much wider one
    vfoptions_wide=vfoptions;
    vfoptions_wide.maxaprimediff=45;
    [V1w,Policy1w]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_wide);
    fprintf('[EZ %s] (GI) maxaprimediff=45 vs the default, V, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V1(:)-V1w(:))))
    fprintf('[EZ %s] (GI) maxaprimediff=45 vs the default, Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1w(:))))

    %% ngridinterp sensitivity: both approximate the same continuous problem
    vfoptions_ng=vfoptions;
    vfoptions_ng.ngridinterp=11;
    [V1ng,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_ng);
    fprintf('[EZ %s] (GI) ngridinterp=11 vs 5, V, this should be small (both approximate the same continuous problem): %2.8e \n',casestr,max(abs(V1(:)-V1ng(:))))

    %% GI on the coarse grid against the plain solve on the big grid
    vfoptions_plain=vfoptions;
    vfoptions_plain.gridinterplayer=0;
    vfoptions_plain=rmfield(vfoptions_plain,'ngridinterp');
    [V1plain,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_plain);
    [V1big,~]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_plain);
    V1bigcoarse=V1big(1:subsamp:end,:);
    dGI=max(abs(gather(V1(:))-gather(V1bigcoarse(:))));
    dplain=max(abs(gather(V1plain(:))-gather(V1bigcoarse(:))));
    fprintf('[EZ %s] |V(coarse grid, with GI) - V(big grid)|, DIAGNOSTIC: %2.8f; |V(coarse grid, plain) - V(big grid)|: %2.8f \n',casestr,dGI,dplain)
    fprintf('[EZ %s] GI on the coarse grid is closer to the big-grid solve than the plain coarse solve is, this should be one: %d \n',casestr,dGI<dplain)

    %% ValueFnFromPolicy_InfHorz with the grid interpolation layer
    VfromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    fprintf('[EZ %s] (GI) ValueFnFromPolicy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(VfromPolicy(:)-V1(:))))

    %% Anti-vacuity
    fprintf('[EZ %s] (GI) Epstein-Zin matters, max|V-Vref| should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(V1(:)-Vrefcomp(:))))
    fprintf('[EZ %s] (GI) Epstein-Zin matters, max|aprime-aprimeref| should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(PolicyVals1(2,:)-PolicyValsref(2,:))))
    fprintf('[EZ %s] (GI) labor supply is interior, share of states with 1<d index<%d should be WELL ABOVE zero: %2.4f \n',casestr,N_d,mean(dindex(:)>1 & dindex(:)<N_d))

    %% Figure (one per case: figure_c, 100+figure_c, 200+figure_c)
    fig=figure(figure_c+100*(ezcase-1)); %#ok<NASGU>
    subplot(1,3,1); plot(a_gridcpu,squeeze(gather(PolicyVals1(2,:,zmid))), a_gridcpu,squeeze(gather(PolicyValsref(2,:,zmid))), a_gridcpu,a_gridcpu,'k:')
    title(['[EZ ',casestr,'] GI: aprime policy (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference','45 degree','Location','northwest')
    subplot(1,3,2); plot(a_gridcpu,squeeze(gather(PolicyVals1(1,:,zmid))), a_gridcpu,squeeze(gather(PolicyValsref(1,:,zmid))))
    title(['[EZ ',casestr,'] GI: d policy (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference','Location','northeast')
    subplot(1,3,3); plot(a_gridcpu,gather(V1(:,zmid)), a_gridcpu,gather(Vrefcomp(:,zmid)))
    title(['[EZ ',casestr,'] GI: value fn (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference (same units)','Location','southeast')
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
