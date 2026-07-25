% Implement core tests of the VFI Toolkit FHorz RiskyAsset commands.
% Full 16-subcode matrix: {nod1,d1} x {noz,z} x {noe,e} x {nosemiz,semiz}.
% Then rerun all 16 with a1 (a=[a1,a2]: a1=safe asset, a2=risky asset).
%
% Riskyasset notes:
%   - aprimeFn(d2, d3, u, ...) — uses iid u shock, NOT current a
%   - vfoptions.refine_d is required (categorises d into d1, d2, d3, [d4 if semiz])
%       d1: in ReturnFn only         (labour supply h here)
%       d2: in aprimeFn only         (riskyshare here)
%       d3: in both                  (savings here)
%       d4: in ReturnFn only AND determines semiz transitions (semiz decision)
%   - noa1 case: divide-and-conquer and grid interpolation layer do NOT apply (no a1 to refine).
%   - with-a1 case: DC, GI and DC+GI are all supported by the toolkit for RiskyAsset.
%
% Ordering mirrors the other core tests: the 8 nosemiz variants first, then the 8 semiz variants
% ({nod1,d1} x {noz,z} x {noe,e} within each). Then the whole 16 repeated with a1.

addpath('./CoreFHorzRiskyAssetTests_subcodes/')
addpath('./CoreFHorzRiskyAssetTests_Setup/')
addpath('./CoreFHorzRiskyAsset_ReturnFns/')

% Setup
CoreFHorzRiskyAsset_setup


%% without d1, without z, without e, without semiz
figure_c=1;
output=CoreFHorzRiskyAsset_nod1_noz_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, without e, without semiz
figure_c=2;
output=CoreFHorzRiskyAsset_d1_noz_noe_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, without e, without semiz
figure_c=3;
output=CoreFHorzRiskyAsset_nod1_z_noe_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, without e, without semiz
figure_c=4;
output=CoreFHorzRiskyAsset_d1_z_noe_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, without z, with e, without semiz
figure_c=5;
output=CoreFHorzRiskyAsset_nod1_noz_e_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, with e, without semiz
figure_c=6;
output=CoreFHorzRiskyAsset_d1_noz_e_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, with e, without semiz
figure_c=7;
output=CoreFHorzRiskyAsset_nod1_z_e_nosemiz(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, with e, without semiz
figure_c=8;
output=CoreFHorzRiskyAsset_d1_z_e_nosemiz(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good


%% Now repeat with semi-exogenous shock
addpath('./CoreFHorzRiskyAssetTests_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzRiskyAsset_ReturnFns/Semiz_ReturnFns/')

%% without d1, without z, without e, with semiz
figure_c=9;
output=CoreFHorzRiskyAsset_nod1_noz_noe_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, without e, with semiz
figure_c=10;
output=CoreFHorzRiskyAsset_d1_noz_noe_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, without e, with semiz
figure_c=11;
output=CoreFHorzRiskyAsset_nod1_z_noe_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, without e, with semiz
figure_c=12;
output=CoreFHorzRiskyAsset_d1_z_noe_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, without z, with e, with semiz
figure_c=13;
output=CoreFHorzRiskyAsset_nod1_noz_e_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, with e, with semiz
figure_c=14;
output=CoreFHorzRiskyAsset_d1_noz_e_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, with e, with semiz
figure_c=15;
output=CoreFHorzRiskyAsset_nod1_z_e_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, with e, with semiz
figure_c=16;
output=CoreFHorzRiskyAsset_d1_z_e_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good


%% Now repeat all 16 with a1 (a=[a1,a2]: a1=safe asset, a2=risky asset).
%% DC, GI and DC+GI are all supported by the toolkit for with-a1 RiskyAsset.
%% Uses n_a_withA1 / a_grid_withA1 (and big variants) from setup.

addpath('./CoreFHorzRiskyAssetTests_subcodes/WithA1_subcodes/')
addpath('./CoreFHorzRiskyAssetTests_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./CoreFHorzRiskyAsset_ReturnFns/WithA1_ReturnFns/')
addpath('./CoreFHorzRiskyAsset_ReturnFns/WithA1_ReturnFns/Semiz_ReturnFns/')

%% With a1

%% without d1, without z, without e, without semiz, with a1
figure_c=17;
output=CoreFHorzRiskyAsset_nod1_noz_noe_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, without e, without semiz, with a1
figure_c=18;
output=CoreFHorzRiskyAsset_d1_noz_noe_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, without e, without semiz, with a1
figure_c=19;
output=CoreFHorzRiskyAsset_nod1_z_noe_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, without e, without semiz, with a1
figure_c=20;
output=CoreFHorzRiskyAsset_d1_z_noe_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, without z, with e, without semiz, with a1
figure_c=21;
output=CoreFHorzRiskyAsset_nod1_noz_e_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, with e, without semiz, with a1
figure_c=22;
output=CoreFHorzRiskyAsset_d1_noz_e_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, with e, without semiz, with a1
figure_c=23;
output=CoreFHorzRiskyAsset_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, with e, without semiz, with a1
figure_c=24;
output=CoreFHorzRiskyAsset_d1_z_e_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good


%% With a1 and semiz

%% without d1, without z, without e, with semiz, with a1
figure_c=25;
output=CoreFHorzRiskyAsset_nod1_noz_noe_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, without e, with semiz, with a1
figure_c=26;
output=CoreFHorzRiskyAsset_d1_noz_noe_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, without e, with semiz, with a1
figure_c=27;
output=CoreFHorzRiskyAsset_nod1_z_noe_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, without e, with semiz, with a1
figure_c=28;
output=CoreFHorzRiskyAsset_d1_z_noe_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, without z, with e, with semiz, with a1
figure_c=29;
output=CoreFHorzRiskyAsset_nod1_noz_e_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, without z, with e, with semiz, with a1
figure_c=30;
output=CoreFHorzRiskyAsset_d1_noz_e_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% without d1, with z, with e, with semiz, with a1
figure_c=31;
output=CoreFHorzRiskyAsset_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good

%% with d1, with z, with e, with semiz, with a1
figure_c=32;
output=CoreFHorzRiskyAsset_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
% looks good





%% =====================================================================
%% CROSS TESTS: numerical equivalences (every printed diff should be ~0)
%%   ExpAssetu     : riskyasset == degenerate experienceassetu (aprimeFn ignores current a)
%%   zase          : z-as-e / iid-as-markov / z&e-collapse equivalences
%%   semizasz      : a JustAMarkov semiz == the equivalent z-markov
%%   plainvswithA1 : plain single-asset == withA1 with a degenerate n_a1=1
%% =====================================================================
addpath('./CoreFHorzRiskyAssetTests_subcodes/CrossTests/')

%% ExpAssetu family
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_d1_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_nod1_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_ExpAssetu_d1_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% z-as-e family
CoreFHorzRiskyAsset_CrossTests_zase_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_nod1_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_d1_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_nod1_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_zase_d1_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% semiz-as-z family (driver passes the semiz n_d; the z-side is derived internally)
CoreFHorzRiskyAsset_CrossTests_semizasz_nod1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_semizasz_d1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_semizasz_nod1_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_semizasz_d1_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% plain-vs-withA1 family (plain grids; the withA1 side is built internally)
CoreFHorzRiskyAsset_CrossTests_plainvswithA1_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_plainvswithA1_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_plainvswithA1_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
CoreFHorzRiskyAsset_CrossTests_plainvswithA1_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

% Regression guard for the per-d4 riskyshare(d2) reconstruction bug in the nod1+a1+semiz+markov-z raw
CoreFHorzRiskyAsset_CrossTests_d2recon_nod1_z_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
