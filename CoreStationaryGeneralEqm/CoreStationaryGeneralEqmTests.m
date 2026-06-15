% Implement core Stationary General Eqm tests of the VFI Toolkit.
%
% Baseline model: Aiyagari-style incomplete markets with endogenous labor and
% a government. Three general eqm prices (r,Tr,tau_c) solve three general eqm
% eqns (CapitalMarket, GovBudget, ConsTax). The wage w is hardcoded from r via
% the firm FOC (not a price); the labor tax rate tau and spending G are fixed.
% See the setup file for details.
%
% For each horizon we:
%  (i)  solve with fminalgo=1, 5, 8, 4 and confirm all give the same answer
%        (fminalgo=4 is CMA-ES, stochastic/lower-accuracy, so a looser match)
%  (ii) re-solve using all three kinds of parameter constraint, each applied to
%        one of the prices: constrain0to1 on r, constrainpositive on Tr,
%        constrainAtoB on tau_c (and then all three at once); confirm each
%        reproduces the unconstrained answer
%
% This file runs: InfHorz, then FHorz, then the same again with permanent types
% (PType, N_i=2 types differing in sigma=2.2 and 1.8). Outputs: output1..output4
% for the non-PType tests, output1ptype..output4ptype for the PType ones.

%%
addpath('./CoreStationaryGeneralEqm_subcodes/')
addpath('./CoreStationaryGeneralEqm_Setup/')
addpath('./CoreStationaryGeneralEqm_ReturnFns/')

% Setup so that we use the same model in both halves
CoreStationaryGE_setup


%% ========================= InfHorz ==========================

%% (i)  fminalgo agreement
output1=CoreStationaryGE_InfHorz_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% (ii) parameter-constraint invariance
output2=CoreStationaryGE_InfHorz_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);


%% Worth doing a 'clear all' here, but not necessary.
% Mainly is so you can run the second half independent of the first half.


%% ========================== FHorz ===========================

%% (i)  fminalgo agreement
output3=CoreStationaryGE_FHorz_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% (ii) parameter-constraint invariance
output4=CoreStationaryGE_FHorz_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);


%% ===================== InfHorz with PType ===================
% Same tests again, but with N_i=2 permanent types differing in sigma (2.2 and 1.8).

%% (i)  fminalgo agreement
output1ptype=CoreStationaryGE_InfHorz_PType_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% (ii) parameter-constraint invariance
output2ptype=CoreStationaryGE_InfHorz_PType_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);


%% ====================== FHorz with PType ====================

%% (i)  fminalgo agreement
output3ptype=CoreStationaryGE_FHorz_PType_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% (ii) parameter-constraint invariance
output4ptype=CoreStationaryGE_FHorz_PType_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% Done!
