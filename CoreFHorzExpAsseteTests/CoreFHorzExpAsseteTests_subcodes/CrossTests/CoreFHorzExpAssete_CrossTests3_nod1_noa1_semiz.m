function output=CoreFHorzExpAssete_CrossTests3_nod1_noa1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 3 (nod1, semiz): noa1 ExpAssete+SemiExo vs a1=1 degenerate ExpAssete+SemiExo.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid.
% Side A: noa1 (the path under test). Side B: a1+ExpAssete with n_a1=1 (a1 does nothing).
% Both sides use the same noa1 ReturnFn (wrapped on Side B to drop a1prime/a1).
% V, Policy d-channels, Dist should match bit-exact.
% Follows CoreFHorzExpAsset_CrossTests4_nod1_semiz, adapted to experienceassete
% (e always present; aprimeFn(d2,a2,e,...)); uses the richest shock case: semiz, z and e.
% PENDING TOOLKIT SUPPORT: ExpAssete+SemiExo+noa1 raws do not exist yet (being built
% test-first), so Side A errors at the ValueFnIter call.
%
% n_d input = [n_d2,n_d3]. d_grid = [d2_grid; d3_grid].

l_d=length(n_d);
n_semiz=vfoptionsbaseline.n_semiz;

% Side B: a1=1 degenerate
n_a1_dummy=1;
a1_grid_dummy=0; % value irrelevant -- ignored by wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common semiz + ExpAssete setup (e always present)
vfoptions=struct();
simoptions=struct();
vfoptions.experienceassete=1;
simoptions.experienceassete=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;

%% semiz, z, e
ReturnFn_A=@(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_nod1_z_e_noa1_semiz(d2,d3,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B=@(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssete_nod1_z_e_noa1_semiz(d2,d3,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

jequaloneDist_A=zeros([n_a,n_semiz,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist_A(1,ceil(n_semiz/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,n_semiz,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist_B(1,1,ceil(n_semiz/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:,:,:);
FnsA.assets=@(d2,d3,a,semiz,z,e) a;
FnsB.assets=@(d2,d3,a1prime,a1,a2,semiz,z,e) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptionsA);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,simoptionsB);
fprintf('Cross test 3 (noa1 vs a1=1, nod1 semiz z e): this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))
% Panels are independent draws: SimPanelIndexes_FHorz_* simulates inside a parfor,
% and rng(1) on the client does not reset the worker streams. So this is Monte-Carlo
% noise, not a disagreement -- it is checked as 'roughly equal', not as an exact zero.
fprintf('Cross test 3 (noa1 vs a1=1, nod1 semiz z e): sim panel means should roughly match: %2.8f vs %2.8f \n',mean(SimPanel_A.assets,'all'),mean(SimPanel_B.assets,'all'))

output=struct();

end
