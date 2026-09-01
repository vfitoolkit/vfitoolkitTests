function output=EZRiskyAsset_CrossTests_degenerateu_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (degenerate u): RiskyAsset with no ACTUAL risky-return risk == the plain
% (non-riskyasset) EZ savings model. EPSTEIN-ZIN, nod1, with markov z. NEW test (no vNM analog):
% it ties this bank to the CoreFHorzTests EZ solvers through a completely different code path.
%
% Design (user-approved):
%   - Params.r=0 (stored/restored), riskyshare grid = {0} (n_d2=1, d2_grid=0), and the d3 savings
%     grid set equal to a_grid (n_d3=n_a). With riskyshare=0 and r=0 the aprimeFn gives
%     aprime = (1+r)*(1-riskyshare)*savings + (1+r+u)*riskyshare*savings = savings EXACTLY,
%     for every u, and savings is ON-GRID (d3_grid=a_grid). u itself is kept at its baseline
%     (n_u=3): it is genuinely irrelevant, which also checks that a degenerate-in-u lottery
%     drops out of the joint certainty-equivalent exactly.
%   - The plain side solves the SAME model with ValueFnIter_Case1_FHorz, n_d=0, aprime chosen
%     directly on a_grid, and the SAME EZ ReturnFn budget (the riskyasset budget has no (1+r)*a
%     term: c = w*kappa_j*z+a-aprime while working, pension+a-aprime after). This is achieved by
%     reusing the riskyasset ReturnFn with aprime passed in the savings slot -- no separate
%     plain-model ReturnFn is needed.
%   - Compare V exactly ('should be zero') and compare the riskyasset policy's savings VALUES
%     against the plain model's aprime VALUES (PolicyInd2Val_FHorz on both; 'should be zero').
% Run for the cons-units case and the negativeUtils case (per the approved design).
% TEST-FIRST: needs the (existing, never GPU-run) EZ riskyasset nod1 with-z raw and the plain
% EZ FHorz solver (GPU-verified in the main EZ bank) to agree.
%
% Bridge variants (user-approved extension): the same degenerate comparison repeated with
%   (+sj)              vfoptions.survivalprobability='sj' on both sides (driver provides the
%                      declining Params.sj=[linspace(1,0.6,N_j-1),0] and Params.oneminussj),
%   (+warm-glow)       vfoptions.WarmGlowBequestsFn on both sides with NO survivalprobability
%                      (both sides print the toolkit's terminal-only warm-glow warning), and
%   (+sj+warm-glow)    both together.
% The warm-glow fns are the bank's EZRiskyWarmGlowFn_{cons,negativeUtils} (same functional forms
% GPU-validated in the main EZ bank); as in the main bank, WarmGlowBequestsFnParamsNames is NOT
% set -- the toolkit reads (aprime,wg1,wg2[,wg3]) off the anonymous handle and pulls wg1/wg2/wg3
% from Params. KEY POINT of the warm-glow bridge: on the riskyasset side the warm-glow is
% evaluated at a2prime through the (d,u) lottery, which under the degenerate u and the on-grid
% d3_grid=a_grid collapses to exactly the plain side's WG(aprime), so the comparison is exact and
% INCLUDES cons-units (both sides share the same terminal additive-after-root convention). These
% legs are the first-ever exercise of the riskyasset EZ sj/warm-glow code paths.

% Store/restore r (Params is function-local, but mirror the approved design explicitly)
r_store=Params.r;
Params.r=0;

% Degenerate riskyasset decision grids: riskyshare={0}, savings grid = a_grid
n_d_deg=[1,n_a]; % [n_d2,n_d3]=[riskyshare,savings]
d_grid_deg=[0; a_grid];

% ReturnFns: the same EZRiskyReturnFn serves both sides (plain side passes aprime in the savings slot)
ReturnFn_risky_cons=@(savings,a,z,r,w,kappa_j,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_plain_cons=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) EZRiskyReturnFn_cons_nod1_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,agej,Jr,pension);
ReturnFn_risky_negU=@(savings,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);
ReturnFn_plain_negU=@(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension) EZRiskyReturnFn_negativeUtils_nod1_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,ezsigma,agej,Jr,pension);

% Riskyasset-side base options (u kept at baseline; the EZ fields are added per case below)
vfoptions_r_base=struct();
vfoptions_r_base.riskyasset=1;
vfoptions_r_base.refine_d=[0,1,1];
vfoptions_r_base.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions_r_base.n_u=vfoptionsbaseline.n_u; vfoptions_r_base.u_grid=vfoptionsbaseline.u_grid; vfoptions_r_base.pi_u=vfoptionsbaseline.pi_u;

%% The two EZ cases of the approved design: cons-units and negative utils
for ezcase=[1,3]
    vfoptions_r=vfoptions_r_base;
    vfoptions_r.exoticpreferences='EpsteinZin';
    vfoptions_p=struct(); % plain (non-riskyasset) side
    vfoptions_p.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptions_r.EZutils=0; vfoptions_r.EZriskaversion='ezgamma'; vfoptions_r.EZeis='ezphi';
        vfoptions_p.EZutils=0; vfoptions_p.EZriskaversion='ezgamma'; vfoptions_p.EZeis='ezphi';
        ReturnFn_risky=ReturnFn_risky_cons; ReturnFn_plain=ReturnFn_plain_cons;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptions_r.EZutils=1; vfoptions_r.EZpositiveutility=0; vfoptions_r.EZriskaversion='ezrisk';
        vfoptions_p.EZutils=1; vfoptions_p.EZpositiveutility=0; vfoptions_p.EZriskaversion='ezrisk';
        ReturnFn_risky=ReturnFn_risky_negU; ReturnFn_plain=ReturnFn_plain_negU;
    end

    % Riskyasset solve (degenerate risk: riskyshare=0 always, r=0, so aprime=savings on-grid)
    [V_r,Policy_r]=ValueFnIter_Case1_FHorz(n_d_deg,n_a,n_z,N_j,d_grid_deg,a_grid,z_grid,pi_z,ReturnFn_risky,Params,DiscountFactorParamNames,[],vfoptions_r);

    % Plain (non-riskyasset) EZ savings model
    [V_p,Policy_p]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_plain,Params,DiscountFactorParamNames,[],vfoptions_p);

    fprintf('CrossTest degenerateu (nod1) [EZ %s], V, this should be zero: %.3e \n',casestr,max(abs(V_r(:)-V_p(:))))

    % Compare savings VALUES: riskyasset PolicyValues row 2 is the savings(d3) value;
    % plain PolicyValues row 1 is the aprime value. Same grid (d3_grid=a_grid), so exact match.
    PolicyVals_r=PolicyInd2Val_FHorz(Policy_r,n_d_deg,n_a,n_z,N_j,d_grid_deg,a_grid,vfoptions_r);
    PolicyVals_p=PolicyInd2Val_FHorz(Policy_p,0,n_a,n_z,N_j,[],a_grid,vfoptions_p);
    savings_r=PolicyVals_r(2,:,:,:);
    aprime_p=PolicyVals_p(1,:,:,:);
    fprintf('CrossTest degenerateu (nod1) [EZ %s], savings-policy values, this should be zero: %.3e \n',casestr,max(abs(savings_r(:)-aprime_p(:))))

end

%% Bridge variants: +sj, +warm-glow (terminal-only), +sj+warm-glow; same two EZ cases each.
% Warm-glow fns follow the main-bank convention: anonymous handles, no WarmGlowBequestsFnParamsNames
% (the toolkit reads wg1/wg2/wg3 off the handle signature and pulls them from Params).
WGFn_cons=@(aprime,wg1,wg2) EZRiskyWarmGlowFn_cons(aprime,wg1,wg2);
WGFn_negU=@(aprime,wg1,wg2,wg3) EZRiskyWarmGlowFn_negativeUtils(aprime,wg1,wg2,wg3);

for bridgevariant=1:3
    if bridgevariant==1 % survival probabilities only
        variantstr='+sj';
    elseif bridgevariant==2 % warm-glow only (terminal-only default; expect the toolkit warning on BOTH sides)
        variantstr='+warm-glow';
    else % survival probabilities AND warm-glow (warm-glow weight 1-sj(j) at every age)
        variantstr='+sj+warm-glow';
    end

    for ezcase=[1,3]
        vfoptions_r=vfoptions_r_base;
        vfoptions_r.exoticpreferences='EpsteinZin';
        vfoptions_p=struct(); % plain (non-riskyasset) side
        vfoptions_p.exoticpreferences='EpsteinZin';
        if ezcase==1 % consumption-units (traditional Epstein-Zin)
            casestr='cons-units';
            vfoptions_r.EZutils=0; vfoptions_r.EZriskaversion='ezgamma'; vfoptions_r.EZeis='ezphi';
            vfoptions_p.EZutils=0; vfoptions_p.EZriskaversion='ezgamma'; vfoptions_p.EZeis='ezphi';
            ReturnFn_risky=ReturnFn_risky_cons; ReturnFn_plain=ReturnFn_plain_cons;
            WGFn=WGFn_cons;
        else % utility-units, negative-valued utility fn
            casestr='negative utils';
            vfoptions_r.EZutils=1; vfoptions_r.EZpositiveutility=0; vfoptions_r.EZriskaversion='ezrisk';
            vfoptions_p.EZutils=1; vfoptions_p.EZpositiveutility=0; vfoptions_p.EZriskaversion='ezrisk';
            ReturnFn_risky=ReturnFn_risky_negU; ReturnFn_plain=ReturnFn_plain_negU;
            WGFn=WGFn_negU;
        end
        if bridgevariant==1 || bridgevariant==3 % survival probabilities on both sides
            vfoptions_r.survivalprobability='sj';
            vfoptions_p.survivalprobability='sj';
        end
        if bridgevariant==2 || bridgevariant==3 % warm-glow on both sides
            vfoptions_r.WarmGlowBequestsFn=WGFn;
            vfoptions_p.WarmGlowBequestsFn=WGFn;
        end

        % Riskyasset solve (degenerate risk: riskyshare=0 always, r=0, so aprime=savings on-grid;
        % the warm-glow's a2prime argument likewise collapses to the on-grid savings value)
        [V_r,Policy_r]=ValueFnIter_Case1_FHorz(n_d_deg,n_a,n_z,N_j,d_grid_deg,a_grid,z_grid,pi_z,ReturnFn_risky,Params,DiscountFactorParamNames,[],vfoptions_r);

        % Plain (non-riskyasset) EZ savings model
        [V_p,Policy_p]=ValueFnIter_Case1_FHorz(0,n_a,n_z,N_j,[],a_grid,z_grid,pi_z,ReturnFn_plain,Params,DiscountFactorParamNames,[],vfoptions_p);

        fprintf('CrossTest degenerateu %s (nod1) [EZ %s], V, this should be zero: %.3e \n',variantstr,casestr,max(abs(V_r(:)-V_p(:))))

        % Compare savings VALUES: riskyasset PolicyValues row 2 is the savings(d3) value;
        % plain PolicyValues row 1 is the aprime value. Same grid (d3_grid=a_grid), so exact match.
        PolicyVals_r=PolicyInd2Val_FHorz(Policy_r,n_d_deg,n_a,n_z,N_j,d_grid_deg,a_grid,vfoptions_r);
        PolicyVals_p=PolicyInd2Val_FHorz(Policy_p,0,n_a,n_z,N_j,[],a_grid,vfoptions_p);
        savings_r=PolicyVals_r(2,:,:,:);
        aprime_p=PolicyVals_p(1,:,:,:);
        fprintf('CrossTest degenerateu %s (nod1) [EZ %s], savings-policy values, this should be zero: %.3e \n',variantstr,casestr,max(abs(savings_r(:)-aprime_p(:))))

    end
end

% Restore r
Params.r=r_store;

%%
output=struct();

end
