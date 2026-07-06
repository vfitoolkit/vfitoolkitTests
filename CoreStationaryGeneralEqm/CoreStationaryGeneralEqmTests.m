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

% fprintf('\n=== InfHorz: fminalgo agreement (r,Tr,tau_c) ===\n')
% fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm1.r,output1.p_eqm1.Tr,output1.p_eqm1.tau_c)
% fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm5.r,output1.p_eqm5.Tr,output1.p_eqm5.tau_c)
% fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm8.r,output1.p_eqm8.Tr,output1.p_eqm8.tau_c)
% fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',output1.p_eqm4.r,output1.p_eqm4.Tr,output1.p_eqm4.tau_c)
% d15=max(abs([output1.p_eqm1.r-output1.p_eqm5.r,output1.p_eqm1.Tr-output1.p_eqm5.Tr,output1.p_eqm1.tau_c-output1.p_eqm5.tau_c]));
% d18=max(abs([output1.p_eqm1.r-output1.p_eqm8.r,output1.p_eqm1.Tr-output1.p_eqm8.Tr,output1.p_eqm1.tau_c-output1.p_eqm8.tau_c]));
% d14=max(abs([output1.p_eqm1.r-output1.p_eqm4.r,output1.p_eqm1.Tr-output1.p_eqm4.Tr,output1.p_eqm1.tau_c-output1.p_eqm4.tau_c]));
% fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
% fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
% fprintf('fminalgo 1 vs 4, this should be near zero: %.8f \n',d14)


%% (ii) parameter-constraint invariance
output2=CoreStationaryGE_InfHorz_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

% fprintf('\n=== InfHorz: parameter-constraint invariance ===\n')
% fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqm0.r,output2.p_eqm0.Tr,output2.p_eqm0.tau_c)
% fprintf('constrain0to1 on r:   r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmA.r,output2.p_eqmA.Tr,output2.p_eqmA.tau_c)
% fprintf('constrainpositive Tr: r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmB.r,output2.p_eqmB.Tr,output2.p_eqmB.tau_c)
% fprintf('constrainAtoB tau_c:  r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmC.r,output2.p_eqmC.Tr,output2.p_eqmC.tau_c)
% fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',output2.p_eqmD.r,output2.p_eqmD.Tr,output2.p_eqmD.tau_c)
% dA=max(abs([output2.p_eqm0.r-output2.p_eqmA.r,output2.p_eqm0.Tr-output2.p_eqmA.Tr,output2.p_eqm0.tau_c-output2.p_eqmA.tau_c]));
% dB=max(abs([output2.p_eqm0.r-output2.p_eqmB.r,output2.p_eqm0.Tr-output2.p_eqmB.Tr,output2.p_eqm0.tau_c-output2.p_eqmB.tau_c]));
% dC=max(abs([output2.p_eqm0.r-output2.p_eqmC.r,output2.p_eqm0.Tr-output2.p_eqmC.Tr,output2.p_eqm0.tau_c-output2.p_eqmC.tau_c]));
% dD=max(abs([output2.p_eqm0.r-output2.p_eqmD.r,output2.p_eqm0.Tr-output2.p_eqmD.Tr,output2.p_eqm0.tau_c-output2.p_eqmD.tau_c]));
% fprintf('constrain0to1 on r, this should be near zero: %.8f \n',dA)
% fprintf('constrainpositive on Tr, this should be near zero: %.8f \n',dB)
% fprintf('constrainAtoB on tau_c, this should be near zero: %.8f \n',dC)
% fprintf('constrain all three, this should be near zero: %.8f \n',dD)





%% ========================== FHorz ===========================

%% (i)  fminalgo agreement
output3=CoreStationaryGE_FHorz_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

% fprintf('\n=== InfHorz: fminalgo agreement (r,Tr,tau_c) ===\n')
% fprintf('fminalgo=1: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm1.r,output3.p_eqm1.Tr,output3.p_eqm1.tau_c)
% fprintf('fminalgo=5: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm5.r,output3.p_eqm5.Tr,output3.p_eqm5.tau_c)
% fprintf('fminalgo=8: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm8.r,output3.p_eqm8.Tr,output3.p_eqm8.tau_c)
% fprintf('fminalgo=4: r=%.6f Tr=%.6f tau_c=%.6f \n',output3.p_eqm4.r,output3.p_eqm4.Tr,output3.p_eqm4.tau_c)
% d15=max(abs([output3.p_eqm1.r-output3.p_eqm5.r,output3.p_eqm1.Tr-output3.p_eqm5.Tr,output3.p_eqm1.tau_c-output3.p_eqm5.tau_c]));
% d18=max(abs([output3.p_eqm1.r-output3.p_eqm8.r,output3.p_eqm1.Tr-output3.p_eqm8.Tr,output3.p_eqm1.tau_c-output3.p_eqm8.tau_c]));
% d14=max(abs([output3.p_eqm1.r-output3.p_eqm4.r,output3.p_eqm1.Tr-output3.p_eqm4.Tr,output3.p_eqm1.tau_c-output3.p_eqm4.tau_c]));
% fprintf('fminalgo 1 vs 5, this should be near zero: %.8f \n',d15)
% fprintf('fminalgo 1 vs 8, this should be near zero: %.8f \n',d18)
% fprintf('fminalgo 1 vs 4, this should be near zero: %.8f \n',d14)

%% (ii) parameter-constraint invariance
output4=CoreStationaryGE_FHorz_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

% fprintf('\n=== FHorz: parameter-constraint invariance ===\n')
% fprintf('unconstrained:        r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqm0.r,output4.p_eqm0.Tr,output4.p_eqm0.tau_c)
% fprintf('constrain0to1 on r:   r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmA.r,output4.p_eqmA.Tr,output4.p_eqmA.tau_c)
% fprintf('constrainpositive Tr: r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmB.r,output4.p_eqmB.Tr,output4.p_eqmB.tau_c)
% fprintf('constrainAtoB tau_c:  r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmC.r,output4.p_eqmC.Tr,output4.p_eqmC.tau_c)
% fprintf('constrain all three:  r=%.6f Tr=%.6f tau_c=%.6f \n',output4.p_eqmD.r,output4.p_eqmD.Tr,output4.p_eqmD.tau_c)
% dA=max(abs([output4.p_eqm0.r-output4.p_eqmA.r,output4.p_eqm0.Tr-output4.p_eqmA.Tr,output4.p_eqm0.tau_c-output4.p_eqmA.tau_c]));
% dB=max(abs([output4.p_eqm0.r-output4.p_eqmB.r,output4.p_eqm0.Tr-output4.p_eqmB.Tr,output4.p_eqm0.tau_c-output4.p_eqmB.tau_c]));
% dC=max(abs([output4.p_eqm0.r-output4.p_eqmC.r,output4.p_eqm0.Tr-output4.p_eqmC.Tr,output4.p_eqm0.tau_c-output4.p_eqmC.tau_c]));
% dD=max(abs([output4.p_eqm0.r-output4.p_eqmD.r,output4.p_eqm0.Tr-output4.p_eqmD.Tr,output4.p_eqm0.tau_c-output4.p_eqmD.tau_c]));
% fprintf('constrain0to1 on r, this should be near zero: %.8f \n',dA)
% fprintf('constrainpositive on Tr, this should be near zero: %.8f \n',dB)
% fprintf('constrainAtoB on tau_c, this should be near zero: %.8f \n',dC)
% fprintf('constrain all three, this should be near zero: %.8f \n',dD)




%% ===================== InfHorz with PType ===================
% Same tests again, but with N_i=2 permanent types differing in sigma (2.2 and 1.8).

%% (i)  fminalgo agreement
output1ptype=CoreStationaryGE_InfHorz_PType_fminalgo(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% (ii) parameter-constraint invariance
output2ptype=CoreStationaryGE_InfHorz_PType_constraints(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);





%% ====================== FHorz with PType ====================
% Same tests again, but with N_i=2 permanent types differing in sigma (2.2 and 1.8).

%% (i)  fminalgo agreement
output3ptype=CoreStationaryGE_FHorz_PType_fminalgo(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);

%% (ii) parameter-constraint invariance
output4ptype=CoreStationaryGE_FHorz_PType_constraints(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,GEPriceParamNames,heteroagentoptionsbaseline,simoptionsbaseline,vfoptionsbaseline);



% I want to next redo all these tests, but using ptype and with GEbyptype.
% This remains to be built.
