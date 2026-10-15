function output=CoreFHorzExpAssetsemiz_CrossTests8_d1_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, experienceassetsemiz): SWAP SYMMETRY. Give the two experience assets
% the same grid and the same law of motion, and a return function symmetric in them (they
% enter only through their average). Then V, Policy and the agent distribution must be
% invariant under swapping the two a2 dimensions.
% Unlike the inert-asset cross-tests this keeps BOTH assets alive, so it is the only test
% that reaches the region where all four corners of the a2 lattice carry weight -- the
% inert tests always interpolate one dimension with probability one.
% Expect machine precision rather than exact zero: the nested interp treats a2 dim 1 as
% inner and dim 2 as outer, so the summation order differs between (i,j) and (j,i).
% NOTE the aprimeFn argument order for experienceassetsemiz: the 'whicha' selector comes AFTER
% semiz, i.e. aprimeFn(d2,a2_1,a2_2,semiz,whicha,params). Both experience assets are driven by
% the semi-exogenous state.
%
% This family always carries semiz, and semiz is the FAST exogenous index, so V is
% (a1, a2_1[, a2_2], semiz, z, e, j) and the a2 slices sit before the shock subscripts.

n_a1=n_a(1);
a1_grid=a_grid(1:n_a1);

% A deliberately small, common a2 grid for both dimensions (the state space is
% n_a1*n_a2sym^2, and this family carries semiz, z and e as well, so keep it tight).
n_a2sym=5;
a2sym_grid=linspace(0,10,n_a2sym)';
n_a_D=[n_a1,n_a2sym,n_a2sym];
a_grid_D=[a1_grid;a2sym_grid;a2sym_grid];
% Identical law of motion in both dimensions
aprimeFn_D=@(d2,a2_1,a2_2,semiz,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)*semiz+(1-phi2)*a2_1)+(whicha==2)*(phi1*(1-d2)*semiz+(1-phi2)*a2_2);
ReturnFn_D=@(d1,d2,d3,a1prime,a1,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_CrossTest_2A2sym_d1(d1,d2,d3,a1prime,a1,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptionsD=struct(); vfoptionsD.experienceassetsemiz=2; vfoptionsD.aprimeFn=aprimeFn_D;
vfoptionsD.n_semiz=vfoptionsbaseline.n_semiz; vfoptionsD.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptionsD.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptionsD.n_e=vfoptionsbaseline.n_e; vfoptionsD.pi_e=vfoptionsbaseline.pi_e; vfoptionsD.e_grid=vfoptionsbaseline.e_grid;
simoptionsD=struct(); simoptionsD.experienceassetsemiz=2; simoptionsD.aprimeFn=vfoptionsD.aprimeFn;
simoptionsD.n_semiz=simoptionsbaseline.n_semiz; simoptionsD.semiz_grid=simoptionsbaseline.semiz_grid; simoptionsD.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsD.n_e=simoptionsbaseline.n_e; simoptionsD.pi_e=simoptionsbaseline.pi_e; simoptionsD.e_grid=simoptionsbaseline.e_grid;
simoptionsD.d_grid=d_grid; simoptionsD.a_grid=a_grid_D; simoptionsD.z_grid=z_grid;
jequaloneDist_D=zeros([n_a1,n_a2sym,n_a2sym,vfoptionsD.n_semiz,n_z,vfoptionsD.n_e],'gpuArray');
jequaloneDist_D(1,1,1,ceil(vfoptionsD.n_semiz/2),ceil(n_z/2),ceil(vfoptionsD.n_e/2))=1; % symmetric under the swap
[V_D,Policy_D]=ValueFnIter_Case1_FHorz(n_d,n_a_D,n_z,N_j,d_grid,a_grid_D,z_grid,pi_z,ReturnFn_D,Params,DiscountFactorParamNames,[],vfoptionsD);
StationaryDist_D=StationaryDist_FHorz_Case1(jequaloneDist_D,AgeWeightParamNames,Policy_D,n_d,n_a_D,n_z,N_j,pi_z,Params,simoptionsD);

%% Checks: swap the two a2 dimensions
% V_D is (a1,a2_1,a2_2,semiz,z,e,j); Policy_D is (daprime,a1,a2_1,a2_2,semiz,z,e,j); dist is like V.
Vswap=permute(V_D,[1,3,2,4,5,6,7]);
Pswap=permute(Policy_D,[1,2,4,3,5,6,7,8]);
Dswap=permute(StationaryDist_D,[1,3,2,4,5,6,7]);
fprintf('Cross test 8 (with2A2 semiz, d1): V symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(V_D(:)-Vswap(:))))
fprintf('Cross test 8 (with2A2 semiz, d1): Policy symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(Policy_D(:)-Pswap(:))))
fprintf('Cross test 8 (with2A2 semiz, d1): agent dist symmetric in the two experience assets, this should be zero: %.3e \n',max(abs(StationaryDist_D(:)-Dswap(:))))

output=struct();

end
