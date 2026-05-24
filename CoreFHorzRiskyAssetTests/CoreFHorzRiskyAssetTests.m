% Implement core tests of the VFI Toolkit FHorz RiskyAsset commands.
% Phase 1: foundation + 4 sample subcodes (with/without d1, with/without semiz; all noz/noe).
%
% Riskyasset notes:
%   - aprimeFn(d2, d3, u, ...) — uses iid u shock, NOT current a
%   - vfoptions.refine_d is required (categorises d into d1, d2, d3, [d4 if semiz])
%       d1: in ReturnFn only         (labour supply h here)
%       d2: in aprimeFn only         (riskyshare here)
%       d3: in both                  (savings here)
%       d4: in ReturnFn only AND determines semiz transitions (semiz decision)
%   - divide-and-conquer and grid interpolation layer are NOT supported for riskyasset


addpath('./CoreFHorzRiskyAssetTests_subcodes/')
addpath('./CoreFHorzRiskyAssetTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzRiskyAssetTests_Setup/')
addpath('./CoreFHorzRiskyAsset_ReturnFns/')
addpath('./CoreFHorzRiskyAsset_ReturnFns/Semiz_ReturnFns/')


%% Setup
CoreFHorzRiskyAsset_setup


%% without d1, without z, without e, without semiz
figure_c=1;
output=CoreFHorzRiskyAsset_nod1_noz_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% should run cleanly

%% with d1, without z, without e, without semiz
figure_c=2;
output=CoreFHorzRiskyAsset_d1_noz_noe_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% should run cleanly

%% without d1, without z, without e, with semiz
figure_c=3;
output=CoreFHorzRiskyAsset_nod1_noz_noe_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% should run cleanly

%% with d1, without z, without e, with semiz
figure_c=4;
output=CoreFHorzRiskyAsset_d1_noz_noe_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% should run cleanly


%% Phase 1 done. Phase 2 (when ready): expand to full 16-subcode matrix (add z, e variants) + cross-tests.
