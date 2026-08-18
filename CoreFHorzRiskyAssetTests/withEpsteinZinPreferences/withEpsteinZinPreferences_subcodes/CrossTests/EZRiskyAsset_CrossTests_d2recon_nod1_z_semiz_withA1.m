function output=EZRiskyAsset_CrossTests_d2recon_nod1_z_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% REGRESSION GUARD (EPSTEIN-ZIN version): the "per-d4 d2 (riskyshare) reconstruction" bug class
% in the RiskyAsset nod1+WITH a1+semiz(+markov-z) raws, run for the three EZ cases.
% Mirrors CoreFHorzRiskyAsset_CrossTests_d2recon_nod1_z_semiz_withA1: whatever EZ raw ends up
% serving this shape must not read the riskyshare(d2) reconstruction from the d4=1 slice only
% (the vNM bug: Policy(1)=riskyshare was read from d2index_ford4_jj WITHOUT the chosen-dsemiz(d4)
% selection term). V is unaffected by that bug class; only the reconstructed riskyshare is wrong,
% and only where the optimal dsemiz differs from 1 AND the optimal riskyshare depends on dsemiz.
% Two Policy-sensitive oracles:
%   (1) ValueFnFromPolicy == V   -- a wrong riskyshare makes the stored policy suboptimal -> V gap.
%   (2) markov-z (routes the nod1 semiz withA1 raw) vs the SAME iid shock declared as e (routes
%       the nod1_noz_e semiz withA1 raw): V and the full Policy must match. Oracle (2) is strictly
%       more sensitive than (1) -- it flags a wrong policy index even if it is value-neutral.
% This reuses the z-as-e equivalence already asserted in EZ zase_nod1_semiz_withA1 (block B); the
% value added is the explicit ValueFnFromPolicy oracle, the trigger diagnostic, and this doc.
% TEST-FIRST: (a) EZ riskyasset semiz withA1 raws do not exist yet; (b) ValueFnFromPolicy under
% EZ riskyasset errors deliberately (no riskyasset branch in ValueFnFromPolicy_FHorz_EpsteinZin);
% (c) the e-side solve needs n_z=0 handling. All three must land before this guard runs through.

% Make the markov z iid so it can equivalently be declared as an e shock (needed for oracle 2).
% The d2<->d4 coupling that triggers the bug comes from semiz (untouched), not from z.
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

% ReturnFns: 2 shapes x 3 EZ cases
ReturnFn_z_cons=@(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_z_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_cons=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_cons_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_posU=@(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_z_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_posU=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_positiveUtils_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z_negU=@(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_z_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e_negU=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) EZRiskyReturnFn_negativeUtils_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,ezsigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);

% RiskyAsset + semiz base options; the markov-z solve routes the raw under guard.
vfoptions_base=struct();
vfoptions_base.riskyasset=1;
vfoptions_base.refine_d=[0,1,1,1]; % nod1: riskyshare(d2), savings(d3), dsemiz(d4)
vfoptions_base.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions_base.n_u=vfoptionsbaseline.n_u; vfoptions_base.u_grid=vfoptionsbaseline.u_grid; vfoptions_base.pi_u=vfoptionsbaseline.pi_u;
vfoptions_base.n_semiz=vfoptionsbaseline.n_semiz; vfoptions_base.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions_base.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

%% The three EZ cases
for ezcase=1:3
    vfoptions=vfoptions_base;
    vfoptions.exoticpreferences='EpsteinZin';
    if ezcase==1 % consumption-units (traditional Epstein-Zin)
        casestr='cons-units';
        vfoptions.EZutils=0;
        vfoptions.EZriskaversion='ezgamma';
        vfoptions.EZeis='ezphi';
        ReturnFn_z=ReturnFn_z_cons; ReturnFn_e=ReturnFn_e_cons;
    elseif ezcase==2 % utility-units, positive-valued utility fn
        casestr='positive utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=1;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn_z=ReturnFn_z_posU; ReturnFn_e=ReturnFn_e_posU;
    else % utility-units, negative-valued utility fn
        casestr='negative utils';
        vfoptions.EZutils=1;
        vfoptions.EZpositiveutility=0;
        vfoptions.EZriskaversion='ezrisk';
        ReturnFn_z=ReturnFn_z_negU; ReturnFn_e=ReturnFn_e_negU;
    end
    % Same options but the shock declared as e (routes the nod1_noz_e semiz withA1 raw).
    vfoptions_e=vfoptions;
    vfoptions_e.n_e=vfoptionsbaseline.n_e; vfoptions_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_e.pi_e=vfoptionsbaseline.pi_e;

    %% Oracle 1: ValueFnFromPolicy on the markov-z solve
    [V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
    V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,vfoptions);
    fprintf('d2recon guard (1) ValueFnFromPolicy [EZ %s], this should be zero: %2.8f \n',casestr,max(abs(V1fromPolicy(:)-V1(:))))

    %% Oracle 2: markov-z vs same iid shock as e -- Policy must match
    [V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_e);
    fprintf('d2recon guard (2) markov-z vs e [EZ %s], this should be zero: V %2.8f, Policy %2.8f \n',casestr,max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))))

    %% Trigger diagnostic: the bug can only manifest where the optimal dsemiz leaves 1.
    % Policy layout for nod1+withA1+semiz: (1)=riskyshare(d2), (2)=savings(d3), (3)=dsemiz(d4), (4)=a1prime.
    dsemizPolicy=Policy1(3,:,:,:,:);
    riskysharePolicy=Policy1(1,:,:,:,:);
    fprintf('d2recon guard trigger diagnostic [EZ %s]: dsemiz takes %d distinct value(s), fraction of states with dsemiz~=1 = %2.4f; riskyshare takes %d distinct value(s). \n',casestr,numel(unique(dsemizPolicy(:))),mean(dsemizPolicy(:)~=1),numel(unique(riskysharePolicy(:))))
    if isscalar(unique(dsemizPolicy(:)))
        fprintf('   WARNING: dsemiz is constant -> this calibration does NOT exercise the d2-reconstruction bug (guard is vacuous here). Strengthen the semiz/portfolio coupling.\n')
    end

end

output=struct();

end
