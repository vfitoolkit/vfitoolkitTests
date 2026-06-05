% PType tests for the Core FHorz commands.
%
% Tests 1-4 use only z (no e, no semiz) and rely on existing return functions
% from CoreFHorzTests/CoreFHorz_ReturnFns/.
%
% The ShockTests test at the end uses all 8 (z,e,semiz) combinations across 8
% PTypes and must be set up via Names_i + per-type structures, because the
% n_z, z_grid, pi_z and vfoptions pieces differ across types.

addpath('./CoreFHorzPTypeTests_subcodes/')
addpath('./CoreFHorzPTypeTests_subcodes/ShockTests/')
addpath('./CoreFHorzPTypeTests_Setup/')

% Reuse return functions and the existing setup from CoreFHorzTests
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/')
addpath('../CoreFHorzTests/CoreFHorz_ReturnFns/Semiz_ReturnFns/')
addpath('../CoreFHorzTests/CoreFHorzTests_Setup/')

%% Setup: builds on CoreFHorz_setup and adds PType-specific pieces
CoreFHorzPType_setup

%% 1. N_i (numeric) vs Names_i (cell) — should give identical output
output=CoreFHorzPType_NivsNames(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i);

%% 2. N_i vs no-PType with z having identity transitions (z encodes the type)
output=CoreFHorzPType_NivsZidentity(n_d,n_a,N_j,d_grid,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);

%% 3. N_i with all types identical ≡ a single no-PType solve
output=CoreFHorzPType_NiAllIdentical_vsNoPType(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i);

%% 4. N_i with two distinct types ≡ two no-PType solves stacked
output=CoreFHorzPType_NiTwoTypes_vsTwoSolves(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames);

%% 5. jequaloneDist accepted three ways: [n_a,n_z], [n_a,n_z,N_i], struct
output=CoreFHorzPType_jequaloneDist_3ways(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i);

%% 6. Per-type N_j (different lifespans) — PType as struct ≡ stacked solo solves
output=CoreFHorzPType_PerTypeNj(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames);

%% ShockTests: one PType solve with 8 different (z,e,semiz) combos
% Requires Names_i (per-type structures) because n_z/z_grid/pi_z and the e/semiz
% pieces of vfoptions differ across types.
output=CoreFHorzPType_ShockTests_8types(n_d,n_a,n_z,n_d_semiz,d_grid_semiz,n_d2_semiz,d2_grid_semiz,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,vfoptionsbaseline,simoptionsbaseline);

% All looks good!
