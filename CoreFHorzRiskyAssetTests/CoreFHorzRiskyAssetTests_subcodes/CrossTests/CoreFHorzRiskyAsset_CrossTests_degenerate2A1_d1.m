function output=CoreFHorzRiskyAsset_CrossTests_degenerate2A1_d1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A1, d1): a DEGENERATE second standard asset a1_2 (single grid point {0})
% must reduce the two-standard-asset model back to the plain single-standard-asset with-a1 model.
% Side A goes through the ordinary RiskyAsset solvers, side B through the length(n_a1)>1 path
% (DC2A / GI2A / DC2A_GI2A). Uses the z, no-e configuration; check V and StationaryDist match exactly.
% Inputs: the with-a1 grids -- n_a=[n_a1, n_a2risky], a_grid=[a1_grid; a2_grid].
% n_a_big/a_grid_big are unused.
%
% n_d input = [n_d1,n_d2,n_d3]

n_a1=n_a(1); n_a2=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_grid=a_grid(n_a1+1:end);
Params.r_a1_2=0.08; % return on the (degenerate) second asset; irrelevant since a1_2 is always zero

%% Model A: plain with-a1 (single standard asset + riskyasset)
ReturnFn_A=@(h,savings,a1prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz_withA1(h,savings,a1prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,r_a1,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.riskyasset=1; vfoptionsA.refine_d=[1,1,1];
vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn; vfoptionsA.n_u=vfoptionsbaseline.n_u; vfoptionsA.u_grid=vfoptionsbaseline.u_grid; vfoptionsA.pi_u=vfoptionsbaseline.pi_u;
simoptionsA=vfoptionsA; simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid; simoptionsA.z_grid=z_grid;
jequaloneDist_A=zeros([n_a1,n_a2,n_z],'gpuArray'); jequaloneDist_A(1,1,ceil(n_z/2))=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model B: with2A1 but a1_2 degenerate (single grid point {0}) -> should equal Model A
n_a_B=[n_a1,1,n_a2];
a_grid_B=[a1_grid;0;a2_grid];
ReturnFn_B=@(h,savings,a1prime,a1_2prime,a1,a1_2,a2,z,r,w,kappa_j,sigma,eta,varphi,r_a1,r_a1_2,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz_with2A1(h,savings,a1prime,a1_2prime,a1,a1_2,a2,z,r,w,kappa_j,sigma,eta,varphi,r_a1,r_a1_2,agej,Jr,pension);
vfoptionsB=vfoptionsA;
simoptionsB=simoptionsA; simoptionsB.a_grid=a_grid_B;
jequaloneDist_B=zeros([n_a1,1,n_a2,n_z],'gpuArray'); jequaloneDist_B(1,1,1,ceil(n_z/2))=1;
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('Cross test (with2A1, d1): degenerate a1_2 reduces to with-a1, this should be zero: V %.3e, Dist %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();
end
