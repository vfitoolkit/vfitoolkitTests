% Implement core tests of the VFI Toolkit FHorz RiskyAsset commands, under EPSTEIN-ZIN preferences.
% Full 16-subcode matrix: {nod1,d1} x {noz,z} x {noe,e} x {nosemiz,semiz}.
% Then rerun all 16 with a1 (a=[a1,a2]: a1=safe asset, a2=risky asset).
%
% NOTE: 2a1 (two standard endogenous assets alongside the risky asset) is NOT covered here;
% those tests WILL BE ADDED LATER.
%
% This is the EZ mirror of CoreFHorzRiskyAssetTests.m. Each subcode runs the model three times:
%   Case 1: consumption-units (traditional Epstein-Zin, vfoptions.EZutils=0)
%   Case 2: utility-units with positive-valued utility fn (EZutils=1, EZpositiveutility=1)
%   Case 3: utility-units with negative-valued utility fn (EZutils=1, EZpositiveutility=0)
% Each subcode then finishes with the EZ special tests:
%   (i)   consumption-units: gamma=1/phi collapses EZ to standard vNM (Policy identical; V related
%         by V_EZ=((1-ezgamma)*V_vNM).^(1/(1-ezgamma)))
%   (ii)  utility-units: EZriskaversion=0 collapses EZ to standard vNM (positive and negative cases)
%   (iii) the same collapses under the grid interpolation layer [withA1 subcodes only: noa1
%         riskyasset has no DC and no GI, as there is no a1 to refine]
%   (iv)  vfoptions.EZoneminusbeta=1 versus manually scaling the return fn
%
% The u-shock and the EZ certainty-equivalent: the u-expectation is taken INSIDE one joint CE
% over (u, semizprime, zprime, eprime). The collapse tests (i)/(ii) only hold exactly under this
% convention, so it is enforced by test. A consequence worth noting: unlike CoreFHorzEZTests
% figs 1-2, the noz_noe shapes here are NOT no-shock models (u always provides genuine risk),
% so no 'EZ without shocks' warning is expected anywhere in this bank.
%
% Riskyasset notes (as in the vNM bank):
%   - aprimeFn(d2, d3, u, ...) — uses iid u shock, NOT current a
%   - vfoptions.refine_d is required (categorises d into d1, d2, d3, [d4 if semiz])
%   - noa1 case: divide-and-conquer and grid interpolation layer do NOT apply.
%   - with-a1 case: DC, GI and DC+GI are all supported by the toolkit for RiskyAsset.
%
% TEST-FIRST STATE (see EZRiskyAsset_coverage_proposal.md). The following are EXPECTED to error
% until the corresponding toolkit features are implemented; they are here so the gap is covered
% and drives the implementation:
%  - EZ riskyasset noz shapes: no EZ riskyasset raw handles n_z=0 (figs 1,2,5,6 + semiz/withA1
%    counterparts error at the first solve).
%  - EZ riskyasset + DC / GI / DC+GI (withA1 figs 17-32): no EZ riskyasset DC/GI tiers exist yet.
%  - EZ riskyasset + semiz (figs 9-16, 25-32): only the nod1 semiz raw exists.
%  - ValueFnFromPolicy under EZ riskyasset: errors deliberately (no riskyasset branch in
%    ValueFnFromPolicy_FHorz_EpsteinZin yet).
%  - The 8 existing EZ riskyasset base raws carry the Change B (outer ezc1 removal) edits and two
%    pre-existing parse-error fixes; they have never been GPU-run, so even the 'existing' shapes
%    (figs 3,4,7,8 and withA1 basic solves) are unverified.
%
% Each subcode draws Case 1 as figure_c, Case 2 as 100+figure_c, and Case 3 as 200+figure_c.

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzRiskyAssetEZTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzRiskyAssetEZTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzRiskyAssetEZTestsdiary.txt

addpath('../CoreFHorzRiskyAssetTests_Setup/')
addpath('./EZRisky_ReturnFns/')
addpath('./withEpsteinZinPreferences_subcodes/Noa1_subcodes/')

% Setup so that use the same d,a,z,e,semiz,u in all the models that use them
CoreFHorzRiskyAsset_setup

% Epstein-Zin preference parameters (the EZ vfoptions themselves are set inside each subcode)
Params.ezgamma=3; % consumption-units risk aversion (vfoptions.EZriskaversion when EZutils=0; must be >1)
Params.ezphi=0.5; % consumption-units elasticity of intertemporal substitution (vfoptions.EZeis; cannot be 1)
Params.ezrisk=3; % utility-units additional risk aversion (vfoptions.EZriskaversion when EZutils=1)
Params.ezsigma=2; % curvature of the positiveUtils/negativeUtils utility fns (>1; the positiveUtils
% family uses a (1+x) shift so it stays strictly positive despite ezsigma>1)

% vfoptions.exoticpreferences='EpsteinZin';
% Case 1: vfoptions.EZutils=0; vfoptions.EZriskaversion='ezgamma'; vfoptions.EZeis='ezphi';
% Case 2: vfoptions.EZutils=1; vfoptions.EZpositiveutility=1; vfoptions.EZriskaversion='ezrisk';
% Case 3: vfoptions.EZutils=1; vfoptions.EZpositiveutility=0; vfoptions.EZriskaversion='ezrisk';

%% without d1, without z, without e, without semiz
figure_c=1;
output=EZRiskyAsset_nod1_noz_noe_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, without e, without semiz
figure_c=2;
output=EZRiskyAsset_d1_noz_noe_nosemiz_noa1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, without e, without semiz
figure_c=3;
output=EZRiskyAsset_nod1_z_noe_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, without e, without semiz
figure_c=4;
output=EZRiskyAsset_d1_z_noe_nosemiz_noa1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, without z, with e, without semiz
figure_c=5;
output=EZRiskyAsset_nod1_noz_e_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, with e, without semiz
figure_c=6;
output=EZRiskyAsset_d1_noz_e_nosemiz_noa1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, with e, without semiz
figure_c=7;
output=EZRiskyAsset_nod1_z_e_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, with e, without semiz
figure_c=8;
output=EZRiskyAsset_d1_z_e_nosemiz_noa1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)


%% Now repeat with semi-exogenous shock
addpath('./withEpsteinZinPreferences_subcodes/Noa1_subcodes/Semiz_subcodes/')
addpath('./EZRisky_ReturnFns/Semiz_ReturnFns/')

%% without d1, without z, without e, with semiz
figure_c=9;
output=EZRiskyAsset_nod1_noz_noe_semiz_noa1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, without e, with semiz
figure_c=10;
output=EZRiskyAsset_d1_noz_noe_semiz_noa1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, without e, with semiz
figure_c=11;
output=EZRiskyAsset_nod1_z_noe_semiz_noa1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, without e, with semiz
figure_c=12;
output=EZRiskyAsset_d1_z_noe_semiz_noa1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, without z, with e, with semiz
figure_c=13;
output=EZRiskyAsset_nod1_noz_e_semiz_noa1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, with e, with semiz
figure_c=14;
output=EZRiskyAsset_d1_noz_e_semiz_noa1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, with e, with semiz
figure_c=15;
output=EZRiskyAsset_nod1_z_e_semiz_noa1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, with e, with semiz
figure_c=16;
output=EZRiskyAsset_d1_z_e_semiz_noa1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)


%% Now repeat all 16 with a1 (a=[a1,a2]: a1=safe asset, a2=risky asset).
%% DC, GI and DC+GI are all supported by the toolkit for with-a1 RiskyAsset (test-first for EZ).
%% Uses n_a_withA1 / a_grid_withA1 (and big variants) from setup.

addpath('./withEpsteinZinPreferences_subcodes/WithA1_subcodes/')
addpath('./withEpsteinZinPreferences_subcodes/WithA1_subcodes/Semiz_subcodes/')
addpath('./EZRisky_ReturnFns/WithA1_ReturnFns/')
addpath('./EZRisky_ReturnFns/WithA1_ReturnFns/Semiz_ReturnFns/')

%% With a1

%% without d1, without z, without e, without semiz, with a1
figure_c=17;
output=EZRiskyAsset_nod1_noz_noe_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, without e, without semiz, with a1
figure_c=18;
output=EZRiskyAsset_d1_noz_noe_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, without e, without semiz, with a1
figure_c=19;
output=EZRiskyAsset_nod1_z_noe_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, without e, without semiz, with a1
figure_c=20;
output=EZRiskyAsset_d1_z_noe_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, without z, with e, without semiz, with a1
figure_c=21;
output=EZRiskyAsset_nod1_noz_e_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, with e, without semiz, with a1
figure_c=22;
output=EZRiskyAsset_d1_noz_e_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, with e, without semiz, with a1
figure_c=23;
output=EZRiskyAsset_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, with e, without semiz, with a1
figure_c=24;
output=EZRiskyAsset_d1_z_e_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)


%% With a1 and semiz

%% without d1, without z, without e, with semiz, with a1
figure_c=25;
output=EZRiskyAsset_nod1_noz_noe_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, without e, with semiz, with a1
figure_c=26;
output=EZRiskyAsset_d1_noz_noe_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, without e, with semiz, with a1
figure_c=27;
output=EZRiskyAsset_nod1_z_noe_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, without e, with semiz, with a1
figure_c=28;
output=EZRiskyAsset_d1_z_noe_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, without z, with e, with semiz, with a1
figure_c=29;
output=EZRiskyAsset_nod1_noz_e_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, without z, with e, with semiz, with a1
figure_c=30;
output=EZRiskyAsset_d1_noz_e_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% without d1, with z, with e, with semiz, with a1
figure_c=31;
output=EZRiskyAsset_nod1_z_e_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)

%% with d1, with z, with e, with semiz, with a1
figure_c=32;
output=EZRiskyAsset_d1_z_e_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_ConsUnits.png'],'Resolution',150)
exportgraphics(figure(100+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_PositiveUtils.png'],'Resolution',150)
exportgraphics(figure(200+figure_c),['../TestOutput/CoreFHorzRiskyAssetEZTests_Fig',num2str(figure_c),'_NegativeUtils.png'],'Resolution',150)


%% =====================================================================
%% CROSS TESTS: numerical equivalences (every printed diff should be ~0)
%%   zase          : z-as-e / iid-as-markov / z&e-collapse equivalences (nested-vs-joint CE
%%                   in the presence of the u shock — the riskyasset analogue of CrossTests3)
%%   semizasz      : a JustAMarkov semiz == the equivalent z-markov
%%   plainvswithA1 : plain single-asset == withA1 with a degenerate n_a1=1
%%   degenerateu   : riskyasset with degenerate risk (riskyshare grid {0}, Params.r=0,
%%                   d3 savings grid = a_grid so aprime=savings lands on-grid; u kept at
%%                   baseline — it is genuinely irrelevant, which also checks a
%%                   degenerate-in-u lottery drops out of the joint CE)
%%                   == the plain (non-riskyasset) EZ savings model — ties this bank to the
%%                   CoreFHorzTests EZ solvers through a completely different code path
%%   d2recon       : regression guard for the per-d4 riskyshare(d2) reconstruction bug path
%%
%% NOT included: the ExpAssetu family (riskyasset == degenerate experienceassetu). That
%% equivalence needs EZ-ExpAssetu solvers, which do not exist; restore it when they do
%% (recorded in EZRiskyAsset_coverage_proposal.md).
%% =====================================================================
addpath('./withEpsteinZinPreferences_subcodes/CrossTests/')

%% z-as-e family
EZRiskyAsset_CrossTests_zase_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_nod1_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_d1_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_nod1_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_zase_d1_semiz_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% semiz-as-z family (driver passes the semiz n_d; the z-side is derived internally)
EZRiskyAsset_CrossTests_semizasz_nod1(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_semizasz_d1(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_semizasz_nod1_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_semizasz_d1_withA1(n_d_withd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% plain-vs-withA1 family (plain grids; the withA1 side is built internally)
EZRiskyAsset_CrossTests_plainvswithA1_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_plainvswithA1_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_plainvswithA1_nod1_semiz(n_d_withoutd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withoutd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);
EZRiskyAsset_CrossTests_plainvswithA1_d1_semiz(n_d_withd1semiz,n_a,n_a_big,n_z,N_j,d_grid_withd1semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% degenerate-u family (riskyasset with no actual risk == plain EZ savings model)
EZRiskyAsset_CrossTests_degenerateu_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Regression guard for the per-d4 riskyshare(d2) reconstruction path in the nod1+a1+semiz+markov-z raw
EZRiskyAsset_CrossTests_d2recon_nod1_z_semiz_withA1(n_d_withoutd1semiz,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1semiz,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

diary off
