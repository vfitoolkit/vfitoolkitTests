% Tests all FHorz toolkit features that consume FnsToEvaluate.
% Runs every consumer command + cross-validations + analytical-truth tests
% on three representative shock configurations.
%
% Consumers exercised:
%   EvalFnOnAgentDist_AggVars_FHorz_Case1
%   EvalFnOnAgentDist_AllStats_FHorz_Case1
%   EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1
%   LifeCycleProfiles_FHorz_Case1
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz
%   SimPanelValues_FHorz_Case1
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz (markov-only config 3)
%
% Configs:
%   d_z_e_nosemiz    -- with d, with z, with e, without semiz
%   nod1_z_noe_semiz -- without d1 (only d2), with z, without e, with semiz
%   d_z_noe_nosemiz  -- with d, with z, no e, no semiz (only config supporting AutoCorrTransProbs)

addpath('./CoreFHorzTests_subcodes/')
addpath('./CoreFHorzTests_Setup/')
addpath('./CoreFHorz_ReturnFns/')
addpath('./CoreFHorzTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorz_ReturnFns/Semiz_ReturnFns/')

%% Shared setup
CoreFHorz_setup

total_fail_count=0;

%% Config 1: with d, with z, with e, without semiz
out1=TestFnsToEvaluate_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
total_fail_count=total_fail_count+out1.fail_count;

%% Config 2: no d1, with z, no e, with semiz
out2=TestFnsToEvaluate_nod1_z_noe_semiz(n_d2_semiz,n_a,n_a_big,n_z,N_j,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
total_fail_count=total_fail_count+out2.fail_count;

%% Config 3: with d, with z, no e, no semiz (the only one supporting AutoCorrTransProbs_FHorz)
out3=TestFnsToEvaluate_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
total_fail_count=total_fail_count+out3.fail_count;

%% Final report
if total_fail_count==0
    fprintf('\n===== TestFnsToEvaluate: 0 tests outside tolerance. All clear. =====\n')
else
    warning('TestFnsToEvaluate: %d test(s) outside tolerance -- scroll up to find the failing line(s).', total_fail_count)
end
