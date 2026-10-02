function output=CoreFHorzExpAssetze_CrossTests7_nod1_noa1_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, experienceassetze, NOA1): an INERT SECOND experience asset must reduce the
% two-experience-asset model back to the one-experience-asset model, when the experience assets
% are the ONLY endogenous states (no a1 at all).
%   Model A: experienceassetze=1, n_a=n_a2_1 (scalar)
%   Model B: experienceassetze=2, n_a=[n_a2_1,n_a2_2], a2_2'=a2_2 (identity) and pensionrate=0.
% This is the noa1 twin of CrossTests7. It matters separately because the noa1 raws are a distinct
% code path (no a1prime offset in the aprime index: the aprime index IS the a2 Kron index), and on
% the plain-ExpAsset side it was exactly this test that caught the first real defect.
%
% NOTE the aprimeFn argument order for experienceassetZE: 'whicha' comes AFTER z and e.
% Inputs: n_a is the scalar one-experience-asset grid size, a_grid is its a2 grid.

n_a2_1=n_a(1);
a2_1_grid=a_grid;

% The inert second experience asset (grid values irrelevant; needs >1 point to be a real check)
n_a2_2=3;
a2_2_grid=linspace(0,1,n_a2_2)';

%% Model A: one experience asset, no a1
ReturnFn_A=@(d2,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetze_nod1_z_e_noa1(d2,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceassetze=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptionsA.n_e=vfoptionsbaseline.n_e; vfoptionsA.pi_e=vfoptionsbaseline.pi_e; vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
simoptionsA=struct(); simoptionsA.experienceassetze=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn;
simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid; simoptionsA.z_grid=z_grid;
simoptionsA.n_e=simoptionsbaseline.n_e; simoptionsA.pi_e=simoptionsbaseline.pi_e; simoptionsA.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_A=zeros([n_a2_1,n_z,vfoptionsA.n_e],'gpuArray'); jequaloneDist_A(1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model B: two experience assets, the second one inert, still no a1
n_a_B=[n_a2_1,n_a2_2];
a_grid_B=[a2_1_grid;a2_2_grid];
aprimeFn_B=@(d2,a2_1,a2_2,z,e,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)*z*e+(1-phi2)*a2_1)+(whicha==2)*a2_2;
ReturnFn_B=@(d2,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate) ReturnFn_ExpAssetze_nod1_z_e_noa1_with2A2(d2,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate);
Params_B=Params; Params_B.pensionrate=0;
vfoptionsB=struct(); vfoptionsB.experienceassetze=2; vfoptionsB.aprimeFn=aprimeFn_B;
vfoptionsB.n_e=vfoptionsbaseline.n_e; vfoptionsB.pi_e=vfoptionsbaseline.pi_e; vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
simoptionsB=struct(); simoptionsB.experienceassetze=2; simoptionsB.aprimeFn=aprimeFn_B;
simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B; simoptionsB.z_grid=z_grid;
simoptionsB.n_e=simoptionsbaseline.n_e; simoptionsB.pi_e=simoptionsbaseline.pi_e; simoptionsB.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_B=zeros([n_a2_1,n_a2_2,n_z,vfoptionsB.n_e],'gpuArray');
jequaloneDist_B(1,1,ceil(n_z/2),ceil(vfoptionsB.n_e/2))=1;
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params_B,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params_B,simoptionsB);

%% Checks
% V_A is (a2_1,z,e,j); V_B is (a2_1,a2_2,z,e,j). Policy_B is (daprime,a2_1,a2_2,z,e,j).
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_2
    Vslice=V_B(:,kk,:,:,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_B(:,1,:,:,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_B(:,:,kk,:,:,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_B(:,1,:,:,:);
fprintf('Cross test 7 (with2A2 ze, nod1 noa1): V is flat in the inert a2_2, this should be zero: %.3e \n',flat)
fprintf('Cross test 7 (with2A2 ze, nod1 noa1): inert a2_2 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 7 (with2A2 ze, nod1 noa1): agent dist at the inert a2_2 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))
fprintf('Cross test 7 (with2A2 ze, nod1 noa1): no mass leaked off the inert a2_2 point, this should be zero: %.3e \n',abs(sum(StationaryDist_B(:))-sum(Dslice(:))))

output=struct();

end
