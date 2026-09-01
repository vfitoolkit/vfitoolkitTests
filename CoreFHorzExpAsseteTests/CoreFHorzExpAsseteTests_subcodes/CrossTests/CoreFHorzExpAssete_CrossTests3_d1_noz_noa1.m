function output=CoreFHorzExpAssete_CrossTests3_d1_noz_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 3 (d1): noa1 ExpAssete vs a1=1 degenerate ExpAssete.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid.
% Side A: noa1 (the path under test). Side B: a1+ExpAssete with n_a1=1 (a1 does nothing).
% Both sides use the same noa1 ReturnFn (wrapped on Side B to drop a1prime/a1).
% V, Policy d-channels, Dist should match bit-exact.
% Follows CoreFHorzExpAsset_CrossTests4_d1_nosemiz, adapted to experienceassete
% (e always present; aprimeFn(d2,a2,e,...)); uses the leanest shock case: e only (no ordinary z).
%
% n_d input = [n_d1,n_d2]. d_grid = [d1_grid; d2_grid].

l_d=length(n_d);

% Side B: a1=1 degenerate
n_a1_dummy=1;
a1_grid_dummy=0; % value irrelevant -- ignored by wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common ExpAssete setup (e always present)
vfoptions=struct();
simoptions=struct();
vfoptions.experienceassete=1;
simoptions.experienceassete=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;

%% noz, e
ReturnFn_A=@(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_noz_e_noa1(d1,d2,a,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_B=@(d1,d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_noz_e_noa1(d1,d2,a2,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

jequaloneDist_A=zeros([n_a,vfoptions.n_e],'gpuArray');
jequaloneDist_A(1,ceil(vfoptions.n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,vfoptions.n_e],'gpuArray');
jequaloneDist_B(1,1,ceil(vfoptions.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:);
FnsA.assets=@(d1,d2,a,e) a;
FnsB.assets=@(d1,d2,a1prime,a1,a2,e) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],[],simoptionsA);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],simoptionsB);
fprintf('Cross test 3 (noa1 vs a1=1, d1 noz e): this should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))
% Panels are independent draws: SimPanelIndexes_FHorz_* simulates inside a parfor,
% and rng(1) on the client does not reset the worker streams. So this is Monte-Carlo
% noise, not a disagreement -- it is checked as 'roughly equal', not as an exact zero.
fprintf('Cross test 3 (noa1 vs a1=1, d1 noz e): sim panel means should roughly match: %2.8f vs %2.8f \n',mean(SimPanel_A.assets,'all'),mean(SimPanel_B.assets,'all'))

output=struct();

end
