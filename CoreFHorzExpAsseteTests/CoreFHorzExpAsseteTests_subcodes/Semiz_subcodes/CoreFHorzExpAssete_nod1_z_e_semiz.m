function output=CoreFHorzExpAssete_nod1_z_e_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)

% Setup vfoptions and simoptions
vfoptions=struct();
simoptions=struct();
% semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
% e
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
% zeros assets, mid points for any shocks
jequaloneDist=zeros([n_a_big,vfoptions.n_semiz,n_z,vfoptions.n_e],'gpuArray'); % Note: based on n_a_big, not n_a
jequaloneDist(1,1,ceil(vfoptions.n_semiz/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

ReturnFn=@(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_nod1_z_e_semiz(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Setup some FnsToEvaluate
FnsToEvaluate.assets=@(d2,d3,a1prime,a1,a2,semiz,z,e) a1;
FnsToEvaluate.humancapital=@(d2,d3,a1prime,a1,a2,semiz,z,e) a2;
FnsToEvaluate.earnings=@(d2,d3,a1prime,a1,a2,semiz,z,e,w,kappa_j) w*kappa_j*d2*a2*semiz*z*e;

% Experience asset (e variant)
vfoptions.experienceassete=1;
simoptions.experienceassete=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.a_grid=a_grid;

%%
vfoptions1=vfoptions;
simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);

PolicyVals1=PolicyInd2Val_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_FHorz(Policy1,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions1);
fprintf('ValueFnFromPolicy, this should be zero: %.3e \n',max(abs(V1fromPolicy(:)-V1(:))))

% Solve with divide-and-conquer, should give same answer
vfoptions2=vfoptions;
vfoptions2.divideandconquer=1;
simoptions2=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);

fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Divide-and-conquer, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))

%
vfoptions1.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-V1B(:))))
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1B(:))))
vfoptions1.lowmemory=2;
[V1C,Policy1C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=2, this should be zero: %.3e \n',max(abs(V1(:)-V1C(:))))
fprintf('lowmemory=2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1C(:))))
vfoptions1.lowmemory=3;
[V1D,Policy1D]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
fprintf('lowmemory=3, this should be zero: %.3e \n',max(abs(V1(:)-V1D(:))))
fprintf('lowmemory=3, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy1D(:))))
vfoptions1.lowmemory=0;

vfoptions2.lowmemory=1;
[V2B,Policy2B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(V2(:)-V2B(:))))
fprintf('lowmemory=1 (with DC), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2B(:))))
vfoptions2.lowmemory=2;
[V2C,Policy2C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=2 (with DC), this should be zero: %.3e \n',max(abs(V2(:)-V2C(:))))
fprintf('lowmemory=2 (with DC), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2C(:))))
vfoptions2.lowmemory=3;
[V2D,Policy2D]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions2);
fprintf('lowmemory=3 (with DC), this should be zero: %.3e \n',max(abs(V2(:)-V2D(:))))
fprintf('lowmemory=3 (with DC), this should be zero: %.3e \n',max(abs(Policy2(:)-Policy2D(:))))
vfoptions2.lowmemory=0;

%%
clear V1 V2 V1B V2B V1C V2C V1D V2D Policy1 Policy2 Policy1B Policy2B Policy1C Policy2C Policy1D Policy2D PolicyVals1 V1fromPolicy
%% Solve with grid-interpolation
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);

PolicyVals3=PolicyInd2Val_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions3);

V3fromPolicy=ValueFnFromPolicy_FHorz(Policy3,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions3);
fprintf('ValueFnFromPolicy with grid interp, this should be zero: %.3e \n',max(abs(V3fromPolicy(:)-V3(:))))

% Solve with divide-and-conquer, should give same answer
vfoptions4=vfoptions;
vfoptions4.divideandconquer=1;
vfoptions4.gridinterplayer=1;
vfoptions4.ngridinterp=5;
simoptions4=simoptions;
simoptions4.gridinterplayer=vfoptions4.gridinterplayer;
simoptions4.ngridinterp=vfoptions4.ngridinterp;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);

fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %.3e \n',max(abs(V3(:)-V4(:))))
fprintf('Divide-and-conquer (with Grid Interp Layer), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy4(:))))

%
vfoptions3.lowmemory=1;
[V3B,Policy3B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3B(:))))
fprintf('lowmemory=1 (with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3B(:))))
vfoptions3.lowmemory=2;
[V3C,Policy3C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=2 (with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3C(:))))
fprintf('lowmemory=2 (with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3C(:))))
vfoptions3.lowmemory=3;
[V3D,Policy3D]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
fprintf('lowmemory=3 (with GI), this should be zero: %.3e \n',max(abs(V3(:)-V3D(:))))
fprintf('lowmemory=3 (with GI), this should be zero: %.3e \n',max(abs(Policy3(:)-Policy3D(:))))
vfoptions3.lowmemory=0;

vfoptions4.lowmemory=1;
[V4B,Policy4B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=1  (with DC+GI), this should be zero: %.3e \n',max(abs(V4(:)-V4B(:))))
fprintf('lowmemory=1  (with DC+GI), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4B(:))))
vfoptions4.lowmemory=2;
[V4C,Policy4C]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=2  (with DC+GI), this should be zero: %.3e \n',max(abs(V4(:)-V4C(:))))
fprintf('lowmemory=2  (with DC+GI), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4C(:))))
vfoptions4.lowmemory=3;
[V4D,Policy4D]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions4);
fprintf('lowmemory=3  (with DC+GI), this should be zero: %.3e \n',max(abs(V4(:)-V4D(:))))
fprintf('lowmemory=3  (with DC+GI), this should be zero: %.3e \n',max(abs(Policy4(:)-Policy4D(:))))
vfoptions4.lowmemory=0;

%%
clear V3 V4 V3B V4B V3C V4C V3D V4D Policy3 Policy4 Policy3B Policy4B Policy3C Policy3C Policy3D Policy4D PolicyVals3 V3fromPolicy
%% Use a really big a_grid, then the moments should be essentially the same with/without grid interpolation

simoptions1.a_grid=a_grid_big;
[V1b,Policy1b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy1b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions1);
AllStats1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);
AgeConditionalStats1=LifeCycleProfiles_FHorz_Case1(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions1);

simoptions3.a_grid=a_grid_big;
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions3);
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy3b,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions3);
AllStats3=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);
AgeConditionalStats3=LifeCycleProfiles_FHorz_Case1(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('With/without grid interp, should get much the same moments (for big a_grid) \n')
fprintf('StationaryDist with/without grid interp, this should be close to zero: %.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.earnings.Gini,AllStats3.earnings.Gini]
[AgeConditionalStats1.earnings.Mean; AgeConditionalStats3.earnings.Mean]
[AgeConditionalStats1.assets.StdDeviation; AgeConditionalStats3.assets.StdDeviation]

%% Sim panel and check it gives the same age conditional stats
% With and without grid interpolation layer
SimPanelValues1=SimPanelValues_FHorz_Case1(jequaloneDist,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z, simoptions1);
SimPanelValues3=SimPanelValues_FHorz_Case1(jequaloneDist,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z, simoptions3);

% Do two comparisons to age conditional stats
fprintf('Without grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats1.earnings.Mean; mean(SimPanelValues1.earnings,2)']
[AgeConditionalStats1.assets.Mean; mean(SimPanelValues1.assets,2)']
[AgeConditionalStats1.humancapital.Mean; mean(SimPanelValues1.humancapital,2)']
fprintf('With grid interp, sim panel data should give roughly the same age conditional stats \n')
[AgeConditionalStats3.earnings.Mean; mean(SimPanelValues3.earnings,2)']
[AgeConditionalStats3.assets.Mean; mean(SimPanelValues3.assets,2)']
[AgeConditionalStats3.humancapital.Mean; mean(SimPanelValues3.humancapital,2)']

%% Compare AllStats to the panel data simulation
% Note: mewj is uniform, so pooling all panel observations uses the same age-weighting as the stationary dist
fprintf('Without grid interp, sim panel data should give roughly the same AllStats \n')
[AllStats1.earnings.Mean, mean(SimPanelValues1.earnings(:))]
[AllStats1.assets.Mean, mean(SimPanelValues1.assets(:))]
[AllStats1.humancapital.Mean, mean(SimPanelValues1.humancapital(:))]
[AllStats1.earnings.StdDeviation, std(SimPanelValues1.earnings(:),1)]
[AllStats1.assets.StdDeviation, std(SimPanelValues1.assets(:),1)]
fprintf('With grid interp, sim panel data should give roughly the same AllStats \n')
[AllStats3.earnings.Mean, mean(SimPanelValues3.earnings(:))]
[AllStats3.assets.Mean, mean(SimPanelValues3.assets(:))]
[AllStats3.humancapital.Mean, mean(SimPanelValues3.humancapital(:))]
[AllStats3.earnings.StdDeviation, std(SimPanelValues3.earnings(:),1)]
[AllStats3.assets.StdDeviation, std(SimPanelValues3.assets(:),1)]

clear V1b V3b Policy1b Policy3b StationaryDist1 StationaryDist3

%% Do some graphs of the age-conditional to see them
fig=figure(figure_c);
subplot(3,1,1); plot(1:1:N_j,AgeConditionalStats1.earnings.Mean, 1:1:N_j,AgeConditionalStats3.earnings.Mean)
title('Earnings Mean (with semiz+z+e)')
legend('1','3')
subplot(3,1,2); plot(1:1:N_j,AgeConditionalStats1.assets.StdDeviation, 1:1:N_j,AgeConditionalStats3.assets.StdDeviation)
title('Assets Std Dev')
legend('1','3')
subplot(3,1,3); plot(1:1:N_j,AgeConditionalStats1.humancapital.Mean, 1:1:N_j,AgeConditionalStats3.humancapital.Mean)
title('Human capital mean')
legend('1','3')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
