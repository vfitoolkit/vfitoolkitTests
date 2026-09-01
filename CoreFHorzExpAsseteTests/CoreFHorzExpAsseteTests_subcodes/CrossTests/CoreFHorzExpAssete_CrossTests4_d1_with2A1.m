function output=CoreFHorzExpAssete_CrossTests4_d1_with2A1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test 4 (d1, with2A1): a DEGENERATE second standard asset a1_2 (single grid point {0})
% must reduce the two-standard-asset model back to the plain single-standard-asset withA1 model.
% Follows CoreFHorzExpAsset_CrossTests_d1_nosemiz_with2A1, adapted to experienceassete
% (e always present; aprimeFn(d2,a2,e,...)); uses the richest shock case: z and e.
% Solve both and check V and StationaryDist match exactly.
% Inputs: the with-a1 grids -- n_a=[n_a1, n_a2experience], a_grid=[a1_grid; a2_grid].

n_a1=n_a(1); n_a2=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_grid=a_grid(n_a1+1:end);
aprimeFn=vfoptionsbaseline.aprimeFn;
Params.r2=0.08; % return on a1_2; irrelevant here as a1_2 is degenerate at 0, but must exist for the ReturnFn

%% Model A: plain withA1 (single standard asset + experience asset), z and e
ReturnFn_A=@(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_z_e(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceassete=1; vfoptionsA.aprimeFn=aprimeFn;
vfoptionsA.n_e=vfoptionsbaseline.n_e; vfoptionsA.e_grid=vfoptionsbaseline.e_grid; vfoptionsA.pi_e=vfoptionsbaseline.pi_e;
simoptionsA=struct(); simoptionsA.experienceassete=1; simoptionsA.aprimeFn=aprimeFn; simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid;
simoptionsA.n_e=simoptionsbaseline.n_e; simoptionsA.e_grid=simoptionsbaseline.e_grid; simoptionsA.pi_e=simoptionsbaseline.pi_e;
jequaloneDist_A=zeros([n_a1,n_a2,n_z,vfoptionsA.n_e],'gpuArray'); jequaloneDist_A(1,1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model B: with2A1 but a1_2 degenerate (single grid point {0}) -> should equal Model A
n_a_B=[n_a1,1,n_a2];
a_grid_B=[a1_grid;0;a2_grid];
ReturnFn_B=@(d1,d2,a1prime,a1_2prime,a1,a1_2,a2,z,e,r,r2,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssete_d1_z_e_with2A1(d1,d2,a1prime,a1_2prime,a1,a1_2,a2,z,e,r,r2,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsB=struct(); vfoptionsB.experienceassete=1; vfoptionsB.aprimeFn=aprimeFn;
vfoptionsB.n_e=vfoptionsbaseline.n_e; vfoptionsB.e_grid=vfoptionsbaseline.e_grid; vfoptionsB.pi_e=vfoptionsbaseline.pi_e;
simoptionsB=struct(); simoptionsB.experienceassete=1; simoptionsB.aprimeFn=aprimeFn; simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B;
simoptionsB.n_e=simoptionsbaseline.n_e; simoptionsB.e_grid=simoptionsbaseline.e_grid; simoptionsB.pi_e=simoptionsbaseline.pi_e;
jequaloneDist_B=zeros([n_a1,1,n_a2,n_z,vfoptionsB.n_e],'gpuArray'); jequaloneDist_B(1,1,1,ceil(n_z/2),ceil(vfoptionsB.n_e/2))=1;
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);

fprintf('Cross test 4 (with2A1, d1, z e): degenerate a1_2 reduces to withA1, this should be zero: V %.3e, Dist %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))))

output=struct();
end
