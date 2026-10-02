function output=CoreFHorzExpAssetze_CrossTests9_nod1_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, experienceassetze): SWAP SYMMETRY. Give the two experience assets the same
% grid and the same law of motion, and a return function symmetric in them (they enter only
% through their average). Then V, Policy and the agent distribution must be exactly invariant
% under swapping the two a2 dimensions.
% Unlike CrossTests7/8 this keeps BOTH assets alive, so it exercises the nested two-corner
% interpolation in the region where all four corners carry weight -- which the inert-asset tests,
% where one dimension always interpolates with probability one, never reach.
%
% NOTE the aprimeFn argument order for experienceassetZE: 'whicha' comes AFTER z and e.
% Expect V/Policy/dist symmetry at machine precision rather than exactly zero: the nested interp
% treats dim 1 as inner and dim 2 as outer, so the summation order differs between (i,j) and
% (j,i). The ExpAsset twin of this test lands at ~1e-14 for that reason.

n_a1=n_a(1);
a1_grid=a_grid(1:n_a1);
% A deliberately small, common a2 grid for both dimensions (the state space is n_a1*n_a2sym^2,
% and this family also carries z and e, so keep it tight).
n_a2sym=5;
a2sym_grid=linspace(0,10,n_a2sym)';

n_a_D=[n_a1,n_a2sym,n_a2sym];
a_grid_D=[a1_grid;a2sym_grid;a2sym_grid];
% Identical law of motion in both dimensions
aprimeFn_D=@(d2,a2_1,a2_2,z,e,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)*z*e+(1-phi2)*a2_1)+(whicha==2)*(phi1*(1-d2)*z*e+(1-phi2)*a2_2);
ReturnFn_D=@(d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetze_nod1_z_e_2A2sym(d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsD=struct(); vfoptionsD.experienceassetze=2; vfoptionsD.aprimeFn=aprimeFn_D;
vfoptionsD.n_e=vfoptionsbaseline.n_e; vfoptionsD.pi_e=vfoptionsbaseline.pi_e; vfoptionsD.e_grid=vfoptionsbaseline.e_grid;
simoptionsD=struct(); simoptionsD.experienceassetze=2; simoptionsD.aprimeFn=aprimeFn_D;
simoptionsD.d_grid=d_grid; simoptionsD.a_grid=a_grid_D; simoptionsD.z_grid=z_grid;
simoptionsD.n_e=simoptionsbaseline.n_e; simoptionsD.pi_e=simoptionsbaseline.pi_e; simoptionsD.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_D=zeros([n_a1,n_a2sym,n_a2sym,n_z,vfoptionsD.n_e],'gpuArray');
jequaloneDist_D(1,1,1,ceil(n_z/2),ceil(vfoptionsD.n_e/2))=1; % symmetric under the swap
[V_D,Policy_D]=ValueFnIter_Case1_FHorz(n_d,n_a_D,n_z,N_j,d_grid,a_grid_D,z_grid,pi_z,ReturnFn_D,Params,DiscountFactorParamNames,[],vfoptionsD);
StationaryDist_D=StationaryDist_FHorz_Case1(jequaloneDist_D,AgeWeightParamNames,Policy_D,n_d,n_a_D,n_z,N_j,pi_z,Params,simoptionsD);

%% Checks: swap the two a2 dimensions
% V_D is (a1,a2_1,a2_2,z,e,j); Policy_D is (daprime,a1,a2_1,a2_2,z,e,j); dist is like V.
Vswap=permute(V_D,[1,3,2,4,5,6]);
Pswap=permute(Policy_D,[1,2,4,3,5,6,7]);
Dswap=permute(StationaryDist_D,[1,3,2,4,5,6]);
fprintf('Cross test 9 (with2A2 ze, nod1): V symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(V_D(:)-Vswap(:))))
fprintf('Cross test 9 (with2A2 ze, nod1): Policy symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(Policy_D(:)-Pswap(:))))
fprintf('Cross test 9 (with2A2 ze, nod1): agent dist symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(StationaryDist_D(:)-Dswap(:))))

output=struct();

end
