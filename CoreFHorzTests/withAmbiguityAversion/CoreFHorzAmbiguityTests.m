% Implement lots of tests of the core VFI Toolkit FHorz commands, under AMBIGUITY AVERSION
% (multiple priors over the exogenous shock transition probabilities; maxmin, so EV is the
% worst case over the priors)
% with/without d
% with z and/or e (ambiguity with no shocks at all is deliberately an error, so no noz+noe)
% with/without divide-and-conquer
% with/without grid interpolation
% with/without low memory (where appropriate; the z&e models also get lowmemory=2)
%
% This is the Ambiguity Aversion mirror of CoreFHorzTests.m. Ambiguity is over the transition
% probabilities ONLY: every prior shares the model shock grid. The regular pi_z/pi_e inputs are
% not used by the value fn solve (which uses vfoptions.ambiguity_pi_z/ambiguity_pi_e), but are
% used for the agent distribution; the toolkit warns if they are not equal to one of the priors.
%
% TEST-FIRST STATE: the main (single-asset) tier is complete and GPU-green (2026-09-01; see
% AmbiguityAversion_testbank_proposal.md). The with2A section at the bottom is test-first (see
% AmbiguityAversion_with2A_proposal.md): its plain-tier solves should pass day one (the AmbAverse
% plain raws are generic in n_a), but the DC2A/GI2A/DC2A_GI2A solves and GI ValueFnFromPolicy
% error until the 2A wave is written (the AmbAverse level-2 dispatchers error on ~isscalar(n_a)).

%% Diary of the command window output (figures are saved into the same folder as they are created)
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzAmbiguityTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzAmbiguityTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzAmbiguityTestsdiary.txt

addpath('../CoreFHorzTests_Setup/')
addpath('../CoreFHorz_ReturnFns/')

addpath('./withAmbiguityAversion_subcodes/')
addpath('./withAmbiguityAversion_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e,semiz in all the models that use them
CoreFHorz_setup

%% Ambiguity Aversion: the multiple priors
% Three priors, shared by every main subcode. Priors 2 and 3 are not ranked against each other,
% so the binding prior can switch across (a,z,j) and the min-over-priors genuinely bites.
vfoptionsbaseline.n_ambiguity=3;
% Priors for z: baseline, contamination toward the worst z', contamination toward uniform
ambiguity_pi_z=zeros(n_z,n_z,3);
ambiguity_pi_z(:,:,1)=pi_z;
ambiguity_pi_z(:,:,2)=0.9*pi_z+0.1*[ones(n_z,1),zeros(n_z,n_z-1)];
ambiguity_pi_z(:,:,3)=0.9*pi_z+0.1*ones(n_z,n_z)/n_z;
vfoptionsbaseline.ambiguity_pi_z=ambiguity_pi_z;
% Priors for e: same recipe
ambiguity_pi_e=zeros(n_e,3);
ambiguity_pi_e(:,1)=vfoptionsbaseline.pi_e;
ambiguity_pi_e(:,2)=0.9*vfoptionsbaseline.pi_e+0.1*[1;zeros(n_e-1,1)];
ambiguity_pi_e(:,3)=0.9*vfoptionsbaseline.pi_e+0.1*ones(n_e,1)/n_e;
vfoptionsbaseline.ambiguity_pi_e=ambiguity_pi_e;

% vfoptions.exoticpreferences='AmbiguityAversion' is set inside each subcode

%% The only functions worth testing are the value fn ones, as after you have Policy everything else is anyway ignoring the ambiguity aversion

%% Ambiguity aversion with no shocks at all is deliberately an error ('what is the point?')
vfoptionstemp.exoticpreferences='AmbiguityAversion';
vfoptionstemp.n_ambiguity=3;
ReturnFn_none=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
try
    [Vtemp,Policytemp]=ValueFnIter_Case1_FHorz(0,n_a,0,N_j,[],a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionstemp);
    fprintf('AmbiguityAversion with no shocks: FAIL, this should error and did not \n')
catch
    fprintf('AmbiguityAversion with no shocks errors as intended :) \n')
end
clear vfoptionstemp

%% without d, with z, without e
figure_c=1;
output=AmbFHorz_nod_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, with z, without e
figure_c=2;
output=AmbFHorz_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, without z, with e
figure_c=3;
output=AmbFHorz_nod_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, without z, with e
figure_c=4;
output=AmbFHorz_d_noz_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, with z, with e
figure_c=5;
output=AmbFHorz_nod_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, with z, with e
figure_c=6;
output=AmbFHorz_d_z_e_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% The cross tests (see the comments at the top of the cross-test subcodes for what they cover)
output=AmbFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=AmbFHorz_CrossTests_d_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% with2A: TWO standard endogenous states (triggers the DC2A/GI2A/DC2A_GI2A code paths), under Ambiguity Aversion
% Mirror of the QH bank's with2A section, reusing the same With2A ReturnFns as the exponential suite.
% TEST-FIRST: plain-tier solves should pass; DC2A/GI2A/DC2A_GI2A error until the 2A raws are written.

addpath('./withAmbiguityAversion_subcodes/With2A_subcodes/')
addpath('./withAmbiguityAversion_subcodes/With2A_subcodes/CrossTests/')
addpath('../CoreFHorz_ReturnFns/With2A_ReturnFns/')

% Redefine the asset grid to two endogenous states for this section
n_a_2A=[n_a,4];
n_a_2A_big=[n_a_big,4];
a2_grid_2A=[0;1;2;3];
a_grid_2A=[a_grid; a2_grid_2A];
a_grid_2A_big=[a_grid_big; a2_grid_2A];
n_a_notsobig=[501,4]; % to test Grid Interpolation without OOM in the d variants
a_grid_notsobig=[5*linspace(0,1,n_a_notsobig(1))'.^3; a2_grid_2A];
Params.phi1=3; % second endo-state preference params
Params.phi2=0.1;

%% without d, with z, without e (with2A)
figure_c=7;
output=AmbFHorz_nod_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, with z, without e (with2A)
figure_c=8;
output=AmbFHorz_d_z_noe_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, without z, with e (with2A)
figure_c=9;
output=AmbFHorz_nod_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, without z, with e (with2A)
figure_c=10;
output=AmbFHorz_d_noz_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% without d, with z, with e (with2A)
figure_c=11;
output=AmbFHorz_nod_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% with d, with z, with e (with2A)
figure_c=12;
output=AmbFHorz_d_z_e_nosemiz_with2A(n_d,n_a_2A,n_a_notsobig,n_z,N_j,d_grid,a_grid_2A,a_grid_notsobig,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
exportgraphics(figure(figure_c),['../TestOutput/CoreFHorzAmbiguityTests_Fig',num2str(figure_c),'.png'],'Resolution',150)

%% The with2A cross tests (see the comments at the top of the cross-test subcodes for what they cover)
output=AmbFHorz_CrossTests_nod_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

output=AmbFHorz_CrossTests_d_nosemiz_with2A(n_d,n_a_2A,n_a_2A_big,n_z,N_j,d_grid,a_grid_2A,a_grid_2A_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline);

%% Done
diary off
