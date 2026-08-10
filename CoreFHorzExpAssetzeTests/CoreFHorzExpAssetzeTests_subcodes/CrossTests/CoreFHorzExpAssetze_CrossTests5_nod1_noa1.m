function output=CoreFHorzExpAssetze_CrossTests5_nod1_noa1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 5 (nod1, nosemiz): noa1 ExpAssetze vs a1=1 degenerate ExpAssetze.
% n_a is scalar (n_a_justexpasset); a_grid is the a2_grid.
% Side A: noa1 (the path under test). Side B: a1+ExpAssetze with n_a1=1 (a1 does nothing).
% Both sides use the same noa1 ReturnFn (wrapped on Side B to drop a1prime/a1).
% V, Policy d-channels, Dist should match bit-exact.
%
% n_d input = [n_d2] (nod1, single decision for expassetze). d_grid = d2_grid.
% (z and e are always present in the experienceassetze family, so just the one z+e case.)

l_d=length(n_d);

% Side B: a1=1 degenerate
n_a1_dummy=1;
a1_grid_dummy=0; % value irrelevant -- ignored by wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common ExpAssetze setup (z and e always present)
vfoptions=struct();
simoptions=struct();
vfoptions.experienceassetze=1;
simoptions.experienceassetze=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.d_grid=d_grid;
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
simoptions.z_grid=z_grid; % aprimeFn depends on z in the ze family

simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;

%% z, e
ReturnFn_A=@(d2,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetze_nod1_z_e_noa1(d2,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_B=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetze_nod1_z_e_noa1(d2,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);

jequaloneDist_A=zeros([n_a,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist_A(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,n_z,vfoptions.n_e],'gpuArray');
jequaloneDist_B(1,1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);

PolicyDpart_B=Policy_B(1:l_d,:,:,:,:,:);
FnsA.assets=@(d2,a2,z,e) a2;
FnsB.assets=@(d2,a1prime,a1,a2,z,e) a2;
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB);
rng(1); SimPanel_A=SimPanelValues_FHorz_Case1(jequaloneDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptionsA);
rng(1); SimPanel_B=SimPanelValues_FHorz_Case1(jequaloneDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,simoptionsB);
fprintf('Cross test 5 (noa1 vs a1=1, nod1 z e): this should be zero: V %2.8f, Policy %2.8f, Dist %2.8f, AllStats.Mean %2.8f, SimPanel.Mean %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(Policy_A(:)-PolicyDpart_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean),abs(mean(SimPanel_A.assets,'all')-mean(SimPanel_B.assets,'all')))

output=struct();

end
