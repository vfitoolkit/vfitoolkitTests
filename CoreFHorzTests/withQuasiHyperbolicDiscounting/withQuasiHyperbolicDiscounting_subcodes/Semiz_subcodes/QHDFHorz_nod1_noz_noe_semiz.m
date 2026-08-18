function output=QHDFHorz_nod1_noz_noe_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)

% n_d=n_d2_semiz;
% d_grid=d2_grid_semiz;

% Setup vfoptions and simoptions
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.d_grid=d_grid;
% Do the current setup
n_z=0;
z_grid=[];
pi_z=[];
% zeros assets, mid points for any shocks
jequaloneDist=zeros(n_a_big,vfoptions.n_semiz,'gpuArray'); % Note: based on n_a_big, not n_a
jequaloneDist(1,ceil(vfoptions.n_semiz/2))=1;

ReturnFn=@(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d2,aprime,a,semiz) a;
FnsToEvaluate.earnings=@(d2,aprime,a,semiz,w,kappa_j) w*kappa_j*semiz;


%% Quasi-Hyperbolic Discounting
vfoptions.exoticpreferences='QuasiHyperbolic';
vfoptions.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
%% Naive
vfoptions.quasi_hyperbolic='Naive';

%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1,V1alt,Policy1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2,V2alt,Policy2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Divide-and-conquer (Valt, Naive), this should be zero: %2.8f \n',max(abs(V1alt(:)-V2alt(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B,V1Balt,Policy1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
fprintf('lowmemory=1 (Valt, Naive), this should be zero: %2.8f \n',max(abs(V1alt(:)-V1Balt(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B,V2Balt,Policy2Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
fprintf('lowmemory=1 (with DC) (Valt, Naive), this should be zero: %2.8f \n',max(abs(V2alt(:)-V2Balt(:))))
vfoptions2.lowmemory=0;

%%
% V from Policy
vfoptions1.Policyalt=Policy1alt; % Naive QH ValueFnFromPolicy requires the exp-discounter argmax
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy (Naive), this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))
fprintf('ValueFnFromPolicy (Valt, Naive), this should be zero: %2.8f \n',max(abs(V1altfromPolicy(:)-V1alt(:))))

clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B V1fromPolicy V2fromPolicy V1alt V2alt V1Balt V2Balt Policy1alt Policy2alt Policy1Balt Policy2Balt V1altfromPolicy V2altfromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3,V3alt,Policy3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4,V4alt,Policy4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer) (Valt, Naive), this should be zero: %2.8f \n',max(abs(V3alt(:)-V4alt(:))))

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B,V3Balt,Policy3Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
fprintf('lowmemory=1 (with GI) (Valt, Naive), this should be zero: %2.8f \n',max(abs(V3alt(:)-V3Balt(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B,V4Balt,Policy4Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
fprintf('lowmemory=1  (with DC+GI) (Valt, Naive), this should be zero: %2.8f \n',max(abs(V4alt(:)-V4Balt(:))))
vfoptions4.lowmemory=0;

%%
% V from Policy
vfoptions3.Policyalt=Policy3alt; % Naive QH ValueFnFromPolicy requires the exp-discounter argmax
[V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy (GI, Naive), this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))
fprintf('ValueFnFromPolicy (Valt, GI, Naive), this should be zero: %2.8f \n',max(abs(V3altfromPolicy(:)-V3alt(:))))

clear V3 V4 V3B V4B Policy3 Policy4 Policy3B Policy4B V3fromPolicy V4fromPolicy V3alt V4alt V3Balt V4Balt Policy3alt Policy4alt Policy3Balt Policy4Balt V3altfromPolicy V4altfromPolicy


%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

[V1b,Policy1b,V1balt,Policy1balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

[V3b,Policy3b,V3balt,Policy3balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3 V1balt V3balt Policy1balt Policy3balt

% This is also true if using divide-and-conquer
[V2b,Policy2b,V2balt,Policy2balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

[V4b,Policy4b,V4balt,Policy4balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; with divide-and-conquer) \n')
[AllStats2.assets.Mean,AllStats4.assets.Mean]
[AllStats2.earnings.Gini,AllStats4.earnings.Gini]
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean]
[AgeConditionalStats2.assets.StdDeviation; AgeConditionalStats4.assets.StdDeviation]

clear V2b V4b StationaryDist2 StationaryDist4 V2balt V4balt % Policy2b Policy4b 

%% Do some graphs of the age-conditional to see them
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean')
legend('1','2','3','4')
subplot(2,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats2.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation, 1:1:N_j,AgeConditionalStats4.assets.StdDeviation)
title('Assets Std Dev')
legend('1','2','3','4')
% ylim([0,0.01]) % If you want to make graph look nicer










%% Sophisticated
vfoptions.quasi_hyperbolic='Sophisticated';

%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1,V1alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2,V2alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Divide-and-conquer (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V1alt(:)-V2alt(:))))

% lowmemory
vfoptions1.lowmemory=1;
[V1B,Policy1B,V1Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1B(:))))
fprintf('lowmemory=1 (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V1alt(:)-V1Balt(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B,V2Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %2.8f \n',max(abs(Policy2(:)-Policy2B(:))))
fprintf('lowmemory=1 (with DC) (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V2alt(:)-V2Balt(:))))
vfoptions2.lowmemory=0;

%%
% V from Policy
[V1fromPolicy,V1altfromPolicy]=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy (Sophisticated), this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))
fprintf('ValueFnFromPolicy (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V1altfromPolicy(:)-V1alt(:))))

clear V1 V2 V1B V2B Policy1 Policy2 Policy1B Policy2B V1fromPolicy V2fromPolicy V1alt V2alt V1Balt V2Balt V1altfromPolicy V2altfromPolicy

%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3,V3alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);


% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4,V4alt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer) (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V3alt(:)-V4alt(:))))

% lowmemory
vfoptions3.lowmemory=1;
[V3B,Policy3B,V3Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3B(:))))
fprintf('lowmemory=1 (with GI) (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V3alt(:)-V3Balt(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B,V4Balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %2.8f \n',max(abs(Policy4(:)-Policy4B(:))))
fprintf('lowmemory=1  (with DC+GI) (Valt, Sophisticated), this should be zero: %2.8f \n',max(abs(V4alt(:)-V4Balt(:))))
vfoptions4.lowmemory=0;

%%
% V from Policy
[V3fromPolicy,V3altfromPolicy]=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy (GI, Sophisticated), this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))
fprintf('ValueFnFromPolicy (Valt, GI, Sophisticated), this should be zero: %2.8f \n',max(abs(V3altfromPolicy(:)-V3alt(:))))

clear V3 V4 V3B V4B Policy3 Policy4 Policy3B Policy4B V3fromPolicy V4fromPolicy V3alt V4alt V3Balt V4Balt V3altfromPolicy V4altfromPolicy


%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

[V1b,Policy1b,V1balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

[V3b,Policy3b,V3balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3 V1balt V3balt

% This is also true if using divide-and-conquer
[V2b,Policy2b,V2balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy2b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions2);
AllStats2=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);
AgeConditionalStats2=LifeCycleProfiles_FHorz_Case1(StationaryDist2,Policy2b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions2);

[V4b,Policy4b,V4balt]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy4b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions4);
AllStats4=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);
AgeConditionalStats4=LifeCycleProfiles_FHorz_Case1(StationaryDist4,Policy4b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions4);

fprintf('With/without grid interp, should get much the same moments (for big a_grid; with divide-and-conquer) \n')
[AllStats2.assets.Mean,AllStats4.assets.Mean]
[AllStats2.earnings.Gini,AllStats4.earnings.Gini]
[AgeConditionalStats2.earnings.Mean; AgeConditionalStats4.earnings.Mean]
[AgeConditionalStats2.assets.StdDeviation; AgeConditionalStats4.assets.StdDeviation]

clear V2b V4b StationaryDist2 StationaryDist4 V2balt V4balt % Policy2b Policy4b 

%% Do some graphs of the age-conditional to see them
fig=figure(100+figure_c); % Sophisticated (Naive is figure_c)
subplot(2,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats2.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean, 1:1:N_j,AgeConditionalStats4.earnings.Mean)
title('Earnings Mean')
legend('1','2','3','4')
subplot(2,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats2.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation, 1:1:N_j,AgeConditionalStats4.assets.StdDeviation)
title('Assets Std Dev')
legend('1','2','3','4')
% ylim([0,0.01]) % If you want to make graph look nicer







%% V_Jplus1: use Valt of period jstar as the terminal value function of a shorter model
% As in the other test banks: solve the model, then solve a shorter model that runs only periods
% 1,...,jstar-1 with Njs=jstar-1 and the age-dependent parameters trimmed to length Njs, and check
% we get the same answer for periods 1,...,jstar-1.
% IMPORTANT (quasi-hyperbolic specific): the continuation value that enters the Bellman equation is
% the STANDARD-discounted one -- V_std when Naive, Vunderbar when Sophisticated -- which is the
% third output Valt. It is NOT V (which is Vtilde/Vhat, the quasi-hyperbolic-discounted object).
% So it is Valt that gets fed back in as vfoptions.V_Jplus1. (Check: the raws use EVpre=Valt(:,:,jj+1)
% / EVpre=Vunderbar(:,:,jj+1) inside the backward loop.)
% Note: this block sits before the "Versus exponential discounting" section on purpose -- that
% section sets Params.beta0=1, which would collapse QH to exponential and stop these tests from
% exercising any of the quasi-hyperbolic machinery.
% Note: mewj is age-dependent, but is only used for the agent distribution, which is not computed
% here, so it is left alone.

%% V_Jplus1, Naive, without divide-and-conquer, without grid interpolation
jstar=round(3*N_j/4);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions1;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Naive';
[Vbase,Policybase,Valtbase,Policyaltbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
Policyaltbase=Policyaltbase(:,:,:,1:Njs);
[Vshort,Policyshort,Valtshort,Policyaltshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Naive), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive, Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive, Policyalt), this should be zero: %2.8f \n',jstar,max(abs(Policyaltbase(:)-Policyaltshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Naive), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive, Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Naive, with divide-and-conquer
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions2;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Naive';
[Vbase,Policybase,Valtbase,Policyaltbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
Policyaltbase=Policyaltbase(:,:,:,1:Njs);
[Vshort,Policyshort,Valtshort,Policyaltshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Naive (with DC)), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with DC)), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with DC), Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with DC), Policyalt), this should be zero: %2.8f \n',jstar,max(abs(Policyaltbase(:)-Policyaltshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Naive (with DC)), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive (with DC)), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive (with DC), Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Naive, with grid interpolation
jstar=round(N_j/2);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions3;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Naive';
[Vbase,Policybase,Valtbase,Policyaltbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
Policyaltbase=Policyaltbase(:,:,:,1:Njs);
[Vshort,Policyshort,Valtshort,Policyaltshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Naive (with GI)), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with GI)), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with GI), Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with GI), Policyalt), this should be zero: %2.8f \n',jstar,max(abs(Policyaltbase(:)-Policyaltshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Naive (with GI)), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive (with GI)), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive (with GI), Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Naive, with divide-and-conquer and grid interpolation
jstar=N_j;
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions4;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Naive';
[Vbase,Policybase,Valtbase,Policyaltbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
Policyaltbase=Policyaltbase(:,:,:,1:Njs);
[Vshort,Policyshort,Valtshort,Policyaltshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Naive (with DC+GI)), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with DC+GI)), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with DC+GI), Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
fprintf('V_Jplus1 (jstar=%i, Naive (with DC+GI), Policyalt), this should be zero: %2.8f \n',jstar,max(abs(Policyaltbase(:)-Policyaltshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Naive (with DC+GI)), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive (with DC+GI)), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Naive (with DC+GI), Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Sophisticated, without divide-and-conquer, without grid interpolation
jstar=round(3*N_j/4);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions1;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Sophisticated';
[Vbase,Policybase,Valtbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Sophisticated), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated, Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Sophisticated), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated, Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Sophisticated, with divide-and-conquer
jstar=round(2*N_j/3);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions2;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Sophisticated';
[Vbase,Policybase,Valtbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with DC)), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with DC)), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with DC), Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with DC)), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with DC)), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with DC), Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Sophisticated, with grid interpolation
jstar=round(N_j/2);
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions3;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Sophisticated';
[Vbase,Policybase,Valtbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with GI)), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with GI)), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with GI), Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with GI)), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with GI)), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with GI), Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

%% V_Jplus1, Sophisticated, with divide-and-conquer and grid interpolation
jstar=N_j;
Njs=jstar-1; % the shorter model runs periods 1,...,jstar-1
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
vfoptionsjs=vfoptions4;
vfoptionsjs.exoticpreferences='QuasiHyperbolic';
vfoptionsjs.quasi_hyperbolic='Sophisticated';
[Vbase,Policybase,Valtbase]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptionsjs);
vfoptionsjs.V_Jplus1=Valtbase(:,:,jstar); % Valt, not V
Vbase=Vbase(:,:,1:Njs);
Policybase=Policybase(:,:,:,1:Njs);
Valtbase=Valtbase(:,:,1:Njs);
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with DC+GI)), this should be zero: %2.8f \n',jstar,max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with DC+GI)), this should be zero: %2.8f \n',jstar,max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1 (jstar=%i, Sophisticated (with DC+GI), Valt), this should be zero: %2.8f \n',jstar,max(abs(Valtbase(:)-Valtshort(:))))
% lowmemory (the V_Jplus1 branch of each raw has its own lowmemory sub-branches)
vfoptionsjs.lowmemory=1;
[Vshort,Policyshort,Valtshort]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn,Paramsjs,DiscountFactorParamNames,[],vfoptionsjs);
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with DC+GI)), this should be zero: %2.8f \n',max(abs(Vbase(:)-Vshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with DC+GI)), this should be zero: %2.8f \n',max(abs(Policybase(:)-Policyshort(:))))
fprintf('V_Jplus1, lowmemory=1 (Sophisticated (with DC+GI), Valt), this should be zero: %2.8f \n',max(abs(Valtbase(:)-Valtshort(:))))
vfoptionsjs.lowmemory=0;

clear Vbase Policybase Valtbase Policyaltbase Vshort Policyshort Valtshort Policyaltshort


%% Versus exponential discounting
Params.beta0=1;

%% Basic (no need to test divide-and-conquer again, as we already know that gives the same as basic)
vfoptions1.exoticpreferences='None';
[V1a,Policy1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='QuasiHyperbolic';
vfoptions1.quasi_hyperbolic='Naive';
[V1b,Policy1b,V1balt,Policy1balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
vfoptions1.exoticpreferences='QuasiHyperbolic';
vfoptions1.quasi_hyperbolic='Sophisticated';
[V1c,Policy1c,V1calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(V1a(:)-V1b(:))))
fprintf('QH with beta0=1 (Valt, Naive): should give zero: %2.8f \n',max(abs(V1balt(:)-V1a(:))))
fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(V1a(:)-V1c(:))))
fprintf('QH with beta0=1 (Valt, Sophisticated): should give zero: %2.8f \n',max(abs(V1calt(:)-V1a(:))))
fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1b(:))))
fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(Policy1a(:)-Policy1c(:))))

%% Grid interpolation layer (no need to test divide-and-conquer again, as we already know that gives the same as basic)
vfoptions3.exoticpreferences='None';
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
vfoptions3.exoticpreferences='QuasiHyperbolic';
vfoptions3.quasi_hyperbolic='Naive';
[V3b,Policy3b,V3balt,Policy3balt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
vfoptions3.exoticpreferences='QuasiHyperbolic';
vfoptions3.quasi_hyperbolic='Sophisticated';
[V3c,Policy3c,V3calt]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(V3a(:)-V3b(:))))
fprintf('QH with beta0=1 (Valt, Naive): should give zero: %2.8f \n',max(abs(V3balt(:)-V3a(:))))
fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(V3a(:)-V3c(:))))
fprintf('QH with beta0=1 (Valt, Sophisticated): should give zero: %2.8f \n',max(abs(V3calt(:)-V3a(:))))
fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3b(:))))
fprintf('QH with beta0=1: should give zero: %2.8f \n',max(abs(Policy3a(:)-Policy3c(:))))


%% Since alternative preferences have no impact beyond the contents of Policy there is not point testing the EvalOnAgentDist functions.

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
