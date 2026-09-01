function output=CoreFHorzExpAssetsemiz_CrossTests4_d1_z_e_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 4 (d1, z, e): noa1 ExpAssetsemiz vs a1=1 degenerate ExpAssetsemiz.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid. n_a_big/a_grid_big are unused.
% Side A: noa1 (the path under test). Side B: a1+ExpAssetsemiz with n_a1=1 (a1 does nothing).
% Both sides use the same noa1 ReturnFn (wrapped on Side B to drop a1prime/a1).
% V, Policy d-channels, Dist should match bit-exact.
%
% n_d input = [n_d1, n_d2, n_d3]. d_grid = [d1_grid; d2_grid; d3_grid].

l_d=length(n_d);
n_semiz=vfoptionsbaseline.n_semiz;
n_e=vfoptionsbaseline.n_e;

% Side B: a1=1 degenerate
n_a1_dummy=1;
a1_grid_dummy=0; % value irrelevant -- ignored by wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common semiz + ExpAssetsemiz setup (with z, with e)
vfoptions=struct();
simoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
vfoptions.experienceassetsemiz=1;
simoptions.experienceassetsemiz=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
simoptions.z_grid=z_grid;

simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;

%% (semiz, z, e)
ReturnFn_A=@(d1,d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_d1_z_e_noa1(d1,d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_d1_z_e_noa1(d1,d2,d3,a2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);

jequaloneDist_A=zeros([n_a,n_semiz,n_z,n_e],'gpuArray');
jequaloneDist_A(1,ceil(n_semiz/2),ceil(n_z/2),ceil(n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,n_semiz,n_z,n_e],'gpuArray');
jequaloneDist_B(1,1,ceil(n_semiz/2),ceil(n_z/2),ceil(n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:,:,:);
FnsA.humancapital=@(d1,d2,d3,a,semiz,z,e) a;
FnsB.humancapital=@(d1,d2,d3,a1prime,a1,a2,semiz,z,e) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptionsA);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,simoptionsB);
fprintf('Cross test 4 (noa1 vs a1=1, d1 semiz z e): this should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.humancapital.Mean-AllStats_B.humancapital.Mean))
% Panels are independent draws: SimPanelIndexes_FHorz_* simulates inside a parfor,
% and rng(1) on the client does not reset the worker streams. So this is Monte-Carlo
% noise, not a disagreement -- it is checked as 'roughly equal', not as an exact zero.
fprintf('Cross test 4 (noa1 vs a1=1, d1 semiz z e): sim panel means should roughly match: %2.8f vs %2.8f \n',mean(SimPanel_A.humancapital,'all'),mean(SimPanel_B.humancapital,'all'))

output=struct();

end
