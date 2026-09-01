function output=CoreFHorzRiskyAsset_CrossTests_d2recon_nod1_z_semiz_withA1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% REGRESSION GUARD: the "per-d4 d2 (riskyshare) reconstruction" bug in
% ValueFnIter_FHorz_RiskyAssetSemiExo_nod1_raw.m (RiskyAsset, nod1, WITH a1, semiz + markov-z).
%
% The bug: Policy(1)=riskyshare (d2) was read from d2index_ford4_jj WITHOUT the chosen-dsemiz(d4)
% selection term, so it always used the d4=1 slice. V is unaffected (computed correctly); only the
% reconstructed riskyshare is wrong -- and only where the optimal dsemiz differs from 1 AND the
% optimal riskyshare depends on the dsemiz decision. Because V stays correct and every lowmemory
% level shares the same (wrong) reconstruction, neither the lowmemory checks nor any V-only compare
% can see it. This guard uses two Policy-sensitive oracles that can:
%   (1) ValueFnFromPolicy == V   -- a wrong riskyshare makes the stored policy suboptimal -> V gap.
%   (2) markov-z (routes nod1_raw) vs the SAME iid shock declared as e (routes the correct
%       nod1_noz_e raw): V and the full Policy must match. Oracle (2) is strictly more sensitive
%       than (1) -- it flags a wrong policy index even if it happens to be value-neutral (a d2 tie).
%
% This reuses the z-as-e equivalence already asserted in zase_nod1_semiz_withA1 (block B); the
% value added here is the explicit ValueFnFromPolicy oracle, the trigger diagnostic, and this doc.
%
% HOW TO VALIDATE THIS GUARD ACTUALLY BITES (mutation test):
%   In ValueFnIter_FHorz_RiskyAssetSemiExo_nod1_raw.m, delete the
%   "+N_d3*N_a1*N_bothz*shiftdim(maxindex-1,-1)" term from the six
%   "Policy(1,...)=shiftdim(d2index_ford4_jj(...),-1)" lines, then rerun: oracles (1) and (2) should
%   print NON-zero. Re-add the term -> back to zero. If they stay zero even after the mutation, this
%   calibration does not exercise the trigger (optimal riskyshare independent of the dsemiz choice):
%   see the diagnostic below and strengthen the semiz<->portfolio coupling (larger income gap across
%   semiz, dsemiz shifting the semiz transition more, wider u_grid) until dsemiz is reported varying.

% Make the markov z iid so it can equivalently be declared as an e shock (needed for oracle 2).
% The d2<->d4 coupling that triggers the bug comes from semiz (untouched), not from z.
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;

ReturnFn_z=@(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e=@(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz_withA1(savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,r_a1,agej,Jr,pension,uempbenefit,searcheffortcost);

% RiskyAsset + semiz options; the markov-z solve routes through nod1_raw (the raw under guard).
vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1,1]; % nod1: riskyshare(d2), savings(d3), dsemiz(d4)
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

% Same options but the shock declared as e (routes the correct nod1_noz_e raw).
vfoptions_e=vfoptions;
vfoptions_e.n_e=vfoptionsbaseline.n_e; vfoptions_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_e.pi_e=vfoptionsbaseline.pi_e;

n_semiz=vfoptions.n_semiz;

%% Oracle 1: ValueFnFromPolicy on the markov-z (nod1_raw) solve
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,vfoptions);
fprintf('d2recon guard (1) ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Oracle 2: markov-z (nod1_raw) vs same iid shock as e (correct nod1_noz_e raw) -- Policy must match
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_e);
fprintf('d2recon guard (2) markov-z vs e, this should be zero: V %.3e, Policy %.3e \n',max(abs(V1(:)-V2(:))),max(abs(Policy1(:)-Policy2(:))))

%% Trigger diagnostic: the bug can only manifest where the optimal dsemiz leaves 1.
% Policy layout for nod1+withA1+semiz: (1)=riskyshare(d2), (2)=savings(d3), (3)=dsemiz(d4), (4)=a1prime.
dsemizPolicy=Policy1(3,:,:,:,:);
riskysharePolicy=Policy1(1,:,:,:,:);
fprintf('d2recon guard trigger diagnostic: dsemiz takes %d distinct value(s), fraction of states with dsemiz~=1 = %2.4f; riskyshare takes %d distinct value(s). \n',numel(unique(dsemizPolicy(:))),mean(dsemizPolicy(:)~=1),numel(unique(riskysharePolicy(:))))
if numel(unique(dsemizPolicy(:)))==1
    fprintf('   WARNING: dsemiz is constant -> this calibration does NOT exercise the d2-reconstruction bug (guard is vacuous here). Strengthen the semiz/portfolio coupling.\n')
end

output=struct();

end
