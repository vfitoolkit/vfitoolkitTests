% Tests of the FHorz RiskyAsset commands under AMBIGUITY AVERSION.
% u is treated as AMBIGUITY, not as risk: the agent does not know the risky return distribution,
% so vfoptions.ambiguity_pi_u (the multiple priors over pi_u) is mandatory in this combination;
% the regular pi_u is only the true process (agent distribution etc.). ambiguity_pi_z and
% ambiguity_pi_e are likewise mandatory when z/e are present, exactly as in the standard-asset
% AmbiguityAversion family. The continuation is the worst case at each expectation stage
% (e-min, then z-min, then — after the aprime lottery, which is conditional on the prior — the
% u-min), and the d2 (riskyshare) refinement max comes after all the mins.
%
% Tiers: noa1 (base method only — no a1 to refine, so no DC/GI) and withA1 (all four methods).
% V/Policy/ValueFnFromPolicy only (exotic preferences change nothing downstream of Policy);
% no dist, no figures. See AmbiguityAversion_RiskyAsset_proposal.md in the toolkit repo.
%
% TEST-FIRST STATE: written before the toolkit side exists. Expect errors until the AA riskyasset
% wave (40 raws + dispatcher + setup + FromPolicy) is written; today Case1 errors upfront
% ('AmbiguityAversion preferences are not implemented for riskyasset').

%% Which parts to run
% One entry per part, in the order they appear below. Set an entry to zero to skip that part.
% That is for building and for rerunning: while one tier is being worked on there is no reason to
% rerun the ones that already pass, and a bank that dies partway (out-of-memory, most often) can
% be finished off by running just the parts that never got to run.
% doPart(1): fig 1
% doPart(2): cross-tests
%
% Skipping a part does NOT renumber the figures. figure_c is a literal in every block, so a given
% Fig number is the same test whatever doPart says, and a png from a previous run is never
% overwritten by a different test.
%
% Parts are independent: the setup, the addpaths and the grid/parameter preambles all sit OUTSIDE
% the if-blocks and so always run, and no part reads another part's output. Any subset can be run,
% in any combination. Anything added to this bank later must keep that true.
doPart=[1,1];

%% Diary of the command window output
if ~exist('../TestOutput','dir')
    mkdir('../TestOutput')
end
if exist('../TestOutput/CoreFHorzRiskyAssetAmbiguityTestsdiary.txt','file')
    delete('../TestOutput/CoreFHorzRiskyAssetAmbiguityTestsdiary.txt') % otherwise diary just appends to the previous run
end
diary ../TestOutput/CoreFHorzRiskyAssetAmbiguityTestsdiary.txt
fprintf('CoreFHorzRiskyAssetAmbiguityTests, run started %s, doPart=[%s] \n',char(datetime('now')),sprintf('%i',doPart))

addpath('../../SharedSubcodes/') % CoreSummary, shared by the Core banks
addpath('../CoreFHorzRiskyAssetTests_Setup/')
addpath('../CoreFHorzRiskyAsset_ReturnFns/')
addpath('../CoreFHorzRiskyAsset_ReturnFns/WithA1_ReturnFns/')

addpath('./withAmbiguityAversion_subcodes/Noa1_subcodes/')
addpath('./withAmbiguityAversion_subcodes/WithA1_subcodes/')
addpath('./withAmbiguityAversion_subcodes/CrossTests/')

% Setup so that use the same d,a,z,e,u in all the models that use them
CoreFHorzRiskyAsset_setup

%% Ambiguity Aversion: the multiple priors
% Three priors per shock, same recipe as the standard-asset AA bank: baseline, contamination
% toward the worst realisation, contamination toward uniform. Priors 2 and 3 are not ranked
% against each other, so the min-over-priors genuinely bites.
vfoptionsbaseline.n_ambiguity=3;
% Priors for u (the risky return distribution — the marquee ambiguity in this family)
ambiguity_pi_u=zeros(n_u,3);
ambiguity_pi_u(:,1)=pi_u;
ambiguity_pi_u(:,2)=0.9*pi_u+0.1*[1;zeros(n_u-1,1)];
ambiguity_pi_u(:,3)=0.9*pi_u+0.1*ones(n_u,1)/n_u;
vfoptionsbaseline.ambiguity_pi_u=ambiguity_pi_u;
% Priors for z
ambiguity_pi_z=zeros(n_z,n_z,3);
ambiguity_pi_z(:,:,1)=pi_z;
ambiguity_pi_z(:,:,2)=0.9*pi_z+0.1*[ones(n_z,1),zeros(n_z,n_z-1)];
ambiguity_pi_z(:,:,3)=0.9*pi_z+0.1*ones(n_z,n_z)/n_z;
vfoptionsbaseline.ambiguity_pi_z=ambiguity_pi_z;
% Priors for e
ambiguity_pi_e=zeros(n_e,3);
ambiguity_pi_e(:,1)=vfoptionsbaseline.pi_e;
ambiguity_pi_e(:,2)=0.9*vfoptionsbaseline.pi_e+0.1*[1;zeros(n_e-1,1)];
ambiguity_pi_e(:,3)=0.9*vfoptionsbaseline.pi_e+0.1*ones(n_e,1)/n_e;
vfoptionsbaseline.ambiguity_pi_e=ambiguity_pi_e;

%% RiskyAsset + AmbiguityAversion without ambiguity_pi_u is deliberately an error (u is ambiguity, not risk)
vfoptionstemp=struct();
vfoptionstemp.riskyasset=1;
vfoptionstemp.refine_d=[0,1,1];
vfoptionstemp.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptionstemp.n_u=vfoptionsbaseline.n_u;
vfoptionstemp.u_grid=vfoptionsbaseline.u_grid;
vfoptionstemp.pi_u=vfoptionsbaseline.pi_u;
vfoptionstemp.exoticpreferences='AmbiguityAversion';
vfoptionstemp.n_ambiguity=3;
ReturnFn_none=@(savings,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension);
try
    [Vtemp,Policytemp]=ValueFnIter_Case1_FHorz(n_d_withoutd1,n_a,0,N_j,d_grid_withoutd1,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionstemp);
    fprintf('RiskyAsset AmbiguityAversion without ambiguity_pi_u: FAIL, this should error and did not \n')
catch
    fprintf('RiskyAsset AmbiguityAversion without ambiguity_pi_u errors as intended :) \n')
end
clear vfoptionstemp

%% ===== doPart(1): fig 1 =====
if doPart(1)==1
    fprintf('\n===== doPart(1): fig 1 =====\n')
    %% Noa1 tier (base method only)
    figure_c=1;
    output=AmbRiskyAsset_nod1_noz_noe_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_noz_noe_nosemiz_noa1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_nod1_z_noe_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_z_noe_nosemiz_noa1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_nod1_noz_e_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,0,N_j,d_grid_withoutd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_noz_e_nosemiz_noa1(n_d_withd1,n_a,n_a_big,0,N_j,d_grid_withd1,a_grid,a_grid_big,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_nod1_z_e_nosemiz_noa1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_z_e_nosemiz_noa1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);

    %% WithA1 tier (all four methods)
    output=AmbRiskyAsset_nod1_noz_noe_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,0,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_noz_noe_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,0,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_nod1_z_noe_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_z_noe_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_nod1_noz_e_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,0,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_noz_e_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,0,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,[],[],Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_nod1_z_e_nosemiz_withA1(n_d_withoutd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withoutd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
    output=AmbRiskyAsset_d1_z_e_nosemiz_withA1(n_d_withd1,n_a_withA1,n_a_big_withA1,n_z,N_j,d_grid_withd1,a_grid_withA1,a_grid_big_withA1,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c);
end % doPart(1): fig 1

%% ===== doPart(2): cross-tests =====
if doPart(2)==1
    fprintf('\n===== doPart(2): cross-tests =====\n')
    %% The cross tests (see the comments at the top of the cross-test subcodes for what they cover)
    output=AmbRiskyAsset_CrossTests_nod1(n_d_withoutd1,n_a,n_a_big,n_z,N_j,d_grid_withoutd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,n_a_withA1,a_grid_withA1);

    output=AmbRiskyAsset_CrossTests_d1(n_d_withd1,n_a,n_a_big,n_z,N_j,d_grid_withd1,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,n_a_withA1,a_grid_withA1);
end % doPart(2): cross-tests

%% Done
%% One verdict for the whole run
% The bank prints a lot of checks; this reads the diary back and says plainly whether the run
% passed, how many checks it contained, and which parts they came from. See CoreSummary.m.
CoreSummary('../TestOutput/CoreFHorzRiskyAssetAmbiguityTestsdiary.txt')

diary off
