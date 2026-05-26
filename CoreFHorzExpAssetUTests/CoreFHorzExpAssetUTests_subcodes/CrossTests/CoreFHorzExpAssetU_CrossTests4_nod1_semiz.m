function output=CoreFHorzExpAssetU_CrossTests4_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 4 (nod1, semiz): noa1 ExpAsset+SemiExo vs a1=1 degenerate ExpAsset+SemiExo.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid.
% Side A: noa1 (the path under test). Side B: a1+ExpAsset with n_a1=1 (a1 does nothing).
% Both sides use the same noa1 ReturnFn (wrapped on Side B to drop a1prime/a1).
% V, Policy d-channels, Dist should match bit-exact.
%
% n_d input = [n_d2; n_d3]. d_grid = [d2_grid; d3_grid].

l_d=length(n_d);
n_semiz=vfoptionsbaseline.n_semiz;

% Side B: a1=1 degenerate
n_a1_dummy=1;
a1_grid_dummy=0; % value irrelevant -- ignored by wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common semiz + ExpAsset setup
vfoptions=struct();
simoptions=struct();
vfoptions.experienceassetu=1;
simoptions.experienceassetu=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
vfoptions.n_u=vfoptionsbaseline.n_u;
vfoptions.u_grid=vfoptionsbaseline.u_grid;
vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions.n_u=vfoptions.n_u;
simoptions.u_grid=vfoptions.u_grid;
simoptions.pi_u=vfoptions.pi_u;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;

vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e;
vfoptions_withe.e_grid=vfoptionsbaseline.e_grid;
vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptions_withe=simoptions;
simoptions_withe.n_e=simoptionsbaseline.n_e;
simoptions_withe.e_grid=simoptionsbaseline.e_grid;
simoptions_withe.pi_e=simoptionsbaseline.pi_e;

simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;
simoptionsA_e=simoptions_withe; simoptionsA_e.a_grid=a_grid;
simoptionsB_e=simoptions_withe; simoptionsB_e.a_grid=a_grid_B;

%% (1) noz, noe
ReturnFn_A1=@(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_noa1_semiz(d2,d3,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B1=@(d2,d3,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_noa1_semiz(d2,d3,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

jequaloneDist_A1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_A1(1,ceil(n_semiz/2))=1;
jequaloneDist_B1=zeros([n_a_B,n_semiz],'gpuArray');
jequaloneDist_B1(1,1,ceil(n_semiz/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A1,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A1,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B1,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B1,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:);
FnsA.assets=@(d2,d3,a,semiz) a;
FnsB.assets=@(d2,d3,a1prime,a1,a2,semiz) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A1,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],[],simoptionsA);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B1,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],simoptionsB);
fprintf('Cross test 4 (noa1 vs a1=1, nod1 semiz noz noe): this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f, SimPanel.Mean %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean),abs(mean(SimPanel_A.assets,'all')-mean(SimPanel_B.assets,'all')))

%% (2) z, noe
ReturnFn_A2=@(d2,d3,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_noa1_semiz(d2,d3,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B2=@(d2,d3,a1prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_noa1_semiz(d2,d3,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

jequaloneDist_A2=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_A2(1,ceil(n_semiz/2),ceil(n_z/2))=1;
jequaloneDist_B2=zeros([n_a_B,n_semiz,n_z],'gpuArray');
jequaloneDist_B2(1,1,ceil(n_semiz/2),ceil(n_z/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A2,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A2,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B2,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B2,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:,:);
FnsA.assets=@(d2,d3,a,semiz,z) a;
FnsB.assets=@(d2,d3,a1prime,a1,a2,semiz,z) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A2,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptionsA);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B2,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,simoptionsB);
fprintf('Cross test 4 (noa1 vs a1=1, nod1 semiz z noe): this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f, SimPanel.Mean %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean),abs(mean(SimPanel_A.assets,'all')-mean(SimPanel_B.assets,'all')))

%% (3) noz, e
ReturnFn_A3=@(d2,d3,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_noa1_semiz(d2,d3,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B3=@(d2,d3,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_noa1_semiz(d2,d3,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

jequaloneDist_A3=zeros([n_a,n_semiz,vfoptions_withe.n_e],'gpuArray');
jequaloneDist_A3(1,ceil(n_semiz/2),ceil(vfoptions_withe.n_e/2))=1;
jequaloneDist_B3=zeros([n_a_B,n_semiz,vfoptions_withe.n_e],'gpuArray');
jequaloneDist_B3(1,1,ceil(n_semiz/2),ceil(vfoptions_withe.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A3,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A3,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA_e);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B3,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B3,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB_e);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:,:);
FnsA.assets=@(d2,d3,a,semiz,e) a;
FnsB.assets=@(d2,d3,a1prime,a1,a2,semiz,e) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA_e);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB_e);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A3,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],[],simoptionsA_e);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B3,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],simoptionsB_e);
fprintf('Cross test 4 (noa1 vs a1=1, nod1 semiz noz e): this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f, SimPanel.Mean %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean),abs(mean(SimPanel_A.assets,'all')-mean(SimPanel_B.assets,'all')))

%% (4) z, e
ReturnFn_A4=@(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_noa1_semiz(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B4=@(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_noa1_semiz(d2,d3,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

jequaloneDist_A4=zeros([n_a,n_semiz,n_z,vfoptions_withe.n_e],'gpuArray');
jequaloneDist_A4(1,ceil(n_semiz/2),ceil(n_z/2),ceil(vfoptions_withe.n_e/2))=1;
jequaloneDist_B4=zeros([n_a_B,n_semiz,n_z,vfoptions_withe.n_e],'gpuArray');
jequaloneDist_B4(1,1,ceil(n_semiz/2),ceil(n_z/2),ceil(vfoptions_withe.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A4,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A4,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA_e);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B4,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B4,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB_e);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:,:,:);
FnsA.assets=@(d2,d3,a,semiz,z,e) a;
FnsB.assets=@(d2,d3,a1prime,a1,a2,semiz,z,e) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA_e);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB_e);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A4,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptionsA_e);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B4,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,simoptionsB_e);
fprintf('Cross test 4 (noa1 vs a1=1, nod1 semiz z e): this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f, SimPanel.Mean %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean),abs(mean(SimPanel_A.assets,'all')-mean(SimPanel_B.assets,'all')))

output=struct();

end
