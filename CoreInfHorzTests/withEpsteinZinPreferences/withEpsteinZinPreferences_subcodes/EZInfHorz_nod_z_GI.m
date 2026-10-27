function output=EZInfHorz_nod_z_GI(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% InfHorz Epstein-Zin, without d, with z, WITH the grid interpolation layer (postGI, the only EZ GI path).
%
% Raw reached: ValueFnIter_InfHorz_EpsteinZin_postGI_nod_raw (new 2026-10-27, Tier 2 of the InfHorz EZ bank).
% The Epstein-Zin grid interpolation layer is always postGI with dense two-point Howards: the
% dispatcher refuses preGI=1, howardsgreedy>0 and postGIrepeat>0 (guards cross-test), and does
% not read howardssparse (as the quasi-hyperbolic path does not).
%
% Policy has three rows: (L1, L2, L2flag). L1 is the lower coarse-grid index, L2 the index into
% the ngridinterp+2 points spanning [L1,L1+1], and L2flag the -Inf-avoidance flag. The aprime
% value is linear in the weight w=(L2-1)/(ngridinterp+1), and so is the interpolated TRANSFORMED
% expectation E[(ezc4*V')^ezc5]: the raw interpolates that object (not V, and not the certainty
% equivalent), then applies the certainty-equivalent power pointwise. That is what the identity
% below rebuilds by hand, from the formula: for each (a,z)
%   EVpow = sum_z' pi(z,z') [ (1-w) T(V(L1,z')) + w T(V(L1+1,z')) ],  T(V)=V^(1-gamma) | V^(1-risk) | (-V)^(1+risk)
%   cons:  V = [ F^(1-1/phi) + beta*EVpow^((1-1/phi)/(1-gamma)) ]^(1/(1-1/phi))
%   pos:   V = F + beta*EVpow^(1/(1-risk))
%   neg:   V = F - beta*EVpow^(1/(1+risk))
% with F at the linear aprime. Holds to convergence tolerance.
%
% Three cases as in EZInfHorz_nod_z, with the GI-specific checks of QHDInfHorz_nod_z_GI: shapes of
% the three rows, the window (maxaprimediff) sensitivity, ngridinterp sensitivity, GI-on-the-coarse-
% grid against the plain solve on the big grid, ValueFnFromPolicy under GI, and the vNM collapse
% identities under GI (exact under GI too, because what is interpolated is the transformed
% expectation and the collapse is an affine map of it).

n_d=0; d_grid=[]; %#ok<NASGU>

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_cons=@(aprime,a,z,r,w) EZReturnFn_cons_nod_z(aprime,a,z,r,w);
ReturnFn_posU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_positiveUtils_nod_z(aprime,a,z,r,w,ezsigma);
ReturnFn_negU=@(aprime,a,z,r,w,ezsigma) EZReturnFn_negativeUtils_nod_z(aprime,a,z,r,w,ezsigma);

% Shorthands for the identity checks (plain CPU doubles throughout)
beta=Params.beta;
r=Params.r;
w=Params.w;
ezgamma=Params.ezgamma;
ezphi=Params.ezphi;
ezrisk=Params.ezrisk;
ezsigma=Params.ezsigma;
a_gridcpu=gather(a_grid);
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
    fprintf('\n--- [EZ %s] without d, with z, grid interpolation layer --- \n',casestr)

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
    fprintf('[EZ %s] (GI) size(Policy) should be [3 %d %d], is [%s] \n',casestr,N_a,N_z,num2str(size(Policy1)))
    L1=gather(reshape(Policy1(1,:,:),[N_a,N_z]));
    L2=gather(reshape(Policy1(2,:,:),[N_a,N_z]));
    L2flag=gather(reshape(Policy1(3,:,:),[N_a,N_z]));
    fprintf('[EZ %s] (GI) L1 in [1,N_a-1], this should be one: %d \n',casestr,all(L1(:)>=1) && all(L1(:)<=N_a-1))
    fprintf('[EZ %s] (GI) L2 in [1,ngridinterp+2], this should be one: %d \n',casestr,all(L2(:)>=1) && all(L2(:)<=n2short+2))
    fprintf('[EZ %s] (GI) L2flag in {1,2,3}, this should be one: %d \n',casestr,all(ismember(L2flag(:),[1,2,3])))
    fprintf('[EZ %s] (GI) L2flag is 2 everywhere (no -Inf return in this model), this should be one: %d \n',casestr,all(L2flag(:)==2))
    fprintf('[EZ %s] (GI) combined fine index in [1,%d], this should be one: %d \n',casestr,N_afine,all((n2short+1)*(L1(:)-1)+L2(:)>=1) && all((n2short+1)*(L1(:)-1)+L2(:)<=N_afine))
    fprintf('[EZ %s] (GI) V is finite everywhere, this should be one: %d \n',casestr,gather(all(isfinite(V1(:)))))
    if ezcase==3
        fprintf('[EZ %s] (GI) V is negative everywhere (negative-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)<0)))
    else
        fprintf('[EZ %s] (GI) V is positive everywhere (consumption units / positive-valued utility), this should be one: %d \n',casestr,gather(all(V1(:)>0)))
    end
    fprintf('[EZ %s] (GI) the L2 channel is not degenerate, number of distinct L2 chosen should be WELL ABOVE zero: %d \n',casestr,numel(unique(L2(:))))

    %% The Epstein-Zin recursion with the two-point interpolated continuation, recomputed here from the formula
    Vg=gather(V1);
    wu=(L2-1)/(n2short+1); % weight on the upper node
    aprimevals=a_gridcpu(L1)+(a_gridcpu(L1+1)-a_gridcpu(L1)).*wu; % linear in the weight
    fprintf('[EZ %s] (GI) reconstructed aprime matches PolicyInd2Val, this should be zero: %2.8e \n',casestr,max(abs(aprimevals(:)-gather(reshape(PolicyVals1(1,:,:),[N_a*N_z,1])))))
    c=(1+r)*avals+w*zvals-aprimevals;
    fprintf('[EZ %s] (GI) chosen aprime always affordable, this should be one: %d \n',casestr,all(c(:)>0))
    Fpol=-Inf(N_a,N_z);
    if ezcase==1
        Fpol(c>0)=c(c>0);
    elseif ezcase==2
        Fpol(c>0)=((1+c(c>0)).^(1-ezsigma)-1)/(1-ezsigma);
    else
        Fpol(c>0)=(c(c>0).^(1-ezsigma))/(1-ezsigma);
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
    fprintf('[EZ %s] (GI) Epstein-Zin recursion with the interpolated continuation, recomputed from the formula at Policy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(resid(:))))
    residalt=Vg-RHSalt;
    fprintf('[EZ %s] (GI) the SAME identity at the vNM-collapse risk aversion, this should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(residalt(:))))

    %% Howards improvement is only an accelerator (this is the two-point Howards block of the new raw)
    vfoptions_noH=vfoptions;
    vfoptions_noH.howards=0;
    [V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_noH);
    fprintf('[EZ %s] (GI) howards=0 (pure VFI), V, this should be zero to convergence tolerance: %.3e \n',casestr,max(abs(V1(:)-V1noH(:))))
    fprintf('[EZ %s] (GI) howards=0 (pure VFI), Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1noH(:))))

    %% lowmemory=1 is not a second path for the nod GI raw (it builds the return matrices once); with d it is (see fig 4)

    %% The maxaprimediff window: the default (5) against a much wider one; a window that binds would pin the policy at its edge
    vfoptions_wide=vfoptions;
    vfoptions_wide.maxaprimediff=25;
    [V1w,Policy1w]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_wide);
    fprintf('[EZ %s] (GI) maxaprimediff=25 vs the default 5, V, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V1(:)-V1w(:))))
    fprintf('[EZ %s] (GI) maxaprimediff=25 vs the default 5, Policy, this should be zero: %.3e \n',casestr,max(abs(Policy1(:)-Policy1w(:))))

    %% ngridinterp sensitivity: both approximate the same continuous problem
    vfoptions_ng=vfoptions;
    vfoptions_ng.ngridinterp=11;
    [V1ng,~]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_ng);
    fprintf('[EZ %s] (GI) ngridinterp=11 vs 5, V, this should be small (both approximate the same continuous problem): %2.8e \n',casestr,max(abs(V1(:)-V1ng(:))))

    %% GI on the coarse grid against the plain solve on the big grid: GI should sit closer to it than the plain coarse solve does
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

    %% ValueFnFromPolicy_InfHorz with the grid interpolation layer (the two-point lookup of the transformed expectation)
    VfromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    fprintf('[EZ %s] (GI) ValueFnFromPolicy, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(VfromPolicy(:)-V1(:))))

    %% Anti-vacuity: Epstein-Zin must matter under GI too
    fprintf('[EZ %s] (GI) Epstein-Zin matters, max|V-Vref| should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(V1(:)-Vrefcomp(:))))
    fprintf('[EZ %s] (GI) Epstein-Zin matters, max|aprime-aprimeref| should be WELL ABOVE zero: %2.8f \n',casestr,max(abs(PolicyVals1(:)-PolicyValsref(:))))

    %% The vNM collapse identities under the grid interpolation layer
    if ezcase==1
        Pi=Params;
        Pi.ezphi=1/Params.ezgamma;
        Pi.ezsigma=Params.ezgamma;
        [V1a,Policy1a]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Pi,DiscountFactorParamNames,[],vfoptions);
        [V1b,Policy1b]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_negU,Pi,DiscountFactorParamNames,[],vfref);
        V1bt=((1-Pi.ezgamma)*V1b).^(1/(1-Pi.ezgamma));
        fprintf('[EZ %s] (GI) gamma=1/phi collapse, V after transform (relative), this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V1a(:)-V1bt(:)))/max(abs(V1a(:))))
        fprintf('[EZ %s] (GI) gamma=1/phi collapse, raw Policy difference, DIAGNOSTIC (near-ties can flip a fine point): %.3e at %i of %i state-rows \n',casestr,max(abs(Policy1a(:)-Policy1b(:))),gather(sum(Policy1a(:)~=Policy1b(:))),3*N_a*N_z)
        V1alt=ValueFnFromPolicy_InfHorz(Policy1b,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_cons,Pi,DiscountFactorParamNames,vfoptions);
        fprintf('[EZ %s] (GI) gamma=1/phi collapse, EZ value of the vNM policy vs V_EZ (relative), this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V1a(:)-V1alt(:)))/max(abs(V1a(:))))
    else
        Pii=Params;
        Pii.ezrisk=0;
        [V2a,Policy2a]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Pii,DiscountFactorParamNames,[],vfoptions);
        [V2b,Policy2b]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Pii,DiscountFactorParamNames,[],vfref);
        fprintf('[EZ %s] (GI) EZriskaversion=0 vs vNM, V, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V2a(:)-V2b(:))))
        fprintf('[EZ %s] (GI) EZriskaversion=0 vs vNM, raw Policy difference, DIAGNOSTIC: %.3e at %i of %i state-rows \n',casestr,max(abs(Policy2a(:)-Policy2b(:))),gather(sum(Policy2a(:)~=Policy2b(:))),3*N_a*N_z)
        V2alt=ValueFnFromPolicy_InfHorz(Policy2b,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Pii,DiscountFactorParamNames,vfoptions);
        fprintf('[EZ %s] (GI) EZriskaversion=0 vs vNM, EZ value of the vNM policy vs V_EZ, this should be zero to convergence tolerance: %2.8e \n',casestr,max(abs(V2a(:)-V2alt(:))))
    end

    %% Figure (one per case: figure_c, 100+figure_c, 200+figure_c)
    fig=figure(figure_c+100*(ezcase-1)); %#ok<NASGU>
    subplot(1,2,1); plot(a_gridcpu,squeeze(gather(PolicyVals1(1,:,zmid))), a_gridcpu,squeeze(gather(PolicyValsref(1,:,zmid))), a_gridcpu,a_gridcpu,'k:')
    title(['[EZ ',casestr,'] GI: aprime policy (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference','45 degree','Location','northwest')
    subplot(1,2,2); plot(a_gridcpu,gather(V1(:,zmid)), a_gridcpu,gather(Vrefcomp(:,zmid)))
    title(['[EZ ',casestr,'] GI: value fn (mid z)']); xlabel('a'); legend('Epstein-Zin','vNM reference (same units)','Location','southeast')
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
