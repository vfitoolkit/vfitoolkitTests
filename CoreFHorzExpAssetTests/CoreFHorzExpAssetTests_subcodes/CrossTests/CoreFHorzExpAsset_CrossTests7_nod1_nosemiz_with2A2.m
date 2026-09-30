function output=CoreFHorzExpAsset_CrossTests7_nod1_nosemiz_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2): SWAP SYMMETRY. Give the two experience assets the same grid and the same
% law of motion, and a return function that is symmetric in them (they enter only through their
% average). Then V, Policy and the agent distribution must be exactly invariant under swapping
% the two a2 dimensions.
% Unlike CrossTests5/6 this keeps BOTH assets alive, so it tests the nested two-corner
% interpolation in the region where all four corners carry weight -- which the inert-asset tests,
% where one dimension always interpolates with probability one, never reach.
% Inputs: the one-experience-asset with-a1 grids, n_a=[n_a1,n_a2_1], a_grid=[a1_grid;a2_1_grid].

n_a1=n_a(1);
a1_grid=a_grid(1:n_a1);
% Use a deliberately small, common a2 grid for both dimensions (the state space is n_a1*n_a2sym^2)
n_a2sym=6;
a2sym_grid=linspace(0,10,n_a2sym)';

n_a_D=[n_a1,n_a2sym,n_a2sym];
a_grid_D=[a1_grid;a2sym_grid;a2sym_grid];
% Identical law of motion in both dimensions
aprimeFn_D=@(d2,a2_1,a2_2,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)+(1-phi2)*a2_1)+(whicha==2)*(phi1*(1-d2)+(1-phi2)*a2_2);
ReturnFn_D=@(d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz_2A2sym(d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsD=struct(); vfoptionsD.experienceasset=2; vfoptionsD.aprimeFn=aprimeFn_D;
simoptionsD=struct(); simoptionsD.experienceasset=2; simoptionsD.aprimeFn=aprimeFn_D; simoptionsD.d_grid=d_grid; simoptionsD.a_grid=a_grid_D;
jequaloneDist_D=zeros([n_a1,n_a2sym,n_a2sym],'gpuArray'); jequaloneDist_D(1,1,1)=1; % symmetric under the swap
[V_D,Policy_D]=ValueFnIter_Case1_FHorz(n_d,n_a_D,0,N_j,d_grid,a_grid_D,[],[],ReturnFn_D,Params,DiscountFactorParamNames,[],vfoptionsD);
StationaryDist_D=StationaryDist_FHorz_Case1(jequaloneDist_D,AgeWeightParamNames,Policy_D,n_d,n_a_D,0,N_j,1,Params,simoptionsD);

%% Checks: V is [n_a1,n_a2sym,n_a2sym,N_j], Policy is [nP,n_a1,n_a2sym,n_a2sym,N_j]
Vswap=permute(V_D,[1,3,2,4]);
Pswap=permute(Policy_D,[1,2,4,3,5]);
Dswap=permute(StationaryDist_D,[1,3,2,4]);
fprintf('Cross test 7 (with2A2, nod1): V symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(V_D(:)-Vswap(:))))
fprintf('Cross test 7 (with2A2, nod1): Policy symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(Policy_D(:)-Pswap(:))))
fprintf('Cross test 7 (with2A2, nod1): agent dist symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(StationaryDist_D(:)-Dswap(:))))

output=struct();
end
