function output=CoreFHorzExpAssetU_CrossTests_nod1_nosemiz_with2A1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A1): a DEGENERATE second standard asset a1_2 (single grid point {0}) must
% reduce the two-standard-asset model back to the plain single-standard-asset with-a1 model.
% Solve both (no z, no e, for simplicity) and check V and StationaryDist match exactly.
% Inputs: the with-a1 grids -- n_a=[n_a1, n_a2experience], a_grid=[a1_grid; a2_grid].

n_a1=n_a(1); n_a2=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_grid=a_grid(n_a1+1:end);
aprimeFn=vfoptionsbaseline.aprimeFn;

%% Model A: plain with-a1 (single standard asset + experience asset)
ReturnFn_A=@(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceassetu=1; vfoptionsA.n_u=vfoptionsbaseline.n_u; vfoptionsA.u_grid=vfoptionsbaseline.u_grid; vfoptionsA.pi_u=vfoptionsbaseline.pi_u; vfoptionsA.aprimeFn=aprimeFn;
simoptionsA=struct(); simoptionsA.experienceassetu=1; simoptionsA.n_u=vfoptionsbaseline.n_u; simoptionsA.u_grid=vfoptionsbaseline.u_grid; simoptionsA.pi_u=vfoptionsbaseline.pi_u; simoptionsA.aprimeFn=aprimeFn; simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid;
jequaloneDist_A=zeros([n_a1,n_a2],'gpuArray'); jequaloneDist_A(1,1)=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,1,Params,simoptionsA);

%% Model B: with2A1 but a1_2 degenerate (single grid point {0}) -> should equal Model A
n_a_B=[n_a1,1,n_a2];
a_grid_B=[a1_grid;0;a2_grid];
ReturnFn_B=@(d2,a1prime,a1_2prime,a1,a1_2,a2,r,r2,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz_with2A1(d2,a1prime,a1_2prime,a1,a1_2,a2,r,r2,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsB=struct(); vfoptionsB.experienceassetu=1; vfoptionsB.n_u=vfoptionsbaseline.n_u; vfoptionsB.u_grid=vfoptionsbaseline.u_grid; vfoptionsB.pi_u=vfoptionsbaseline.pi_u; vfoptionsB.aprimeFn=aprimeFn;
simoptionsB=struct(); simoptionsB.experienceassetu=1; simoptionsB.n_u=vfoptionsbaseline.n_u; simoptionsB.u_grid=vfoptionsbaseline.u_grid; simoptionsB.pi_u=vfoptionsbaseline.pi_u; simoptionsB.aprimeFn=aprimeFn; simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B;
jequaloneDist_B=zeros([n_a1,1,n_a2],'gpuArray'); jequaloneDist_B(1,1,1)=1;
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,1,Params,simoptionsB);

fprintf('Cross test (with2A1, nod1): degenerate a1_2 reduces to with-a1, this should be zero: V %2.8f, Dist %2.8f \n',max(abs(V_A(:)-V_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();
end
