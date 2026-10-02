function output=CoreFHorzExpAssetze_CrossTests7_nod1_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, experienceassetze): an INERT SECOND experience asset must reduce the
% two-experience-asset model back to the one-experience-asset with-a1 model.
%   Model A: experienceassetze=1, n_a=[n_a1,n_a2_1]
%   Model B: experienceassetze=2, n_a=[n_a1,n_a2_1,n_a2_2], with a2_2'=a2_2 (the identity map,
%            which the aprimeFn interpolation reproduces exactly for on-grid values) and
%            pensionrate=0, so a2_2 never enters the return function.
% V_B must be exactly flat in the a2_2 dimension, and equal V_A at every a2_2 slice. Together
% with CrossTests8 (which makes the FIRST asset the inert one) this pins down the ordering of the
% two a2 dimensions: either on its own would still pass if the nested interpolation had them
% transposed. These are the machine-precision checks of the tier.
%
% NOTE the aprimeFn argument order for experienceassetZE: the 'whicha' selector comes AFTER z and
% e, i.e. aprimeFn(d2,a2_1,a2_2,z,e,whicha,params). That is NOT where it sits for the plain
% experience asset, where it follows the a2 inputs directly.
%
% Model A V is (a1,a2_1,z,e,j); Model B V is (a1,a2_1,a2_2,z,e,j), so the Model B slices carry two
% more subscripts than the plain-ExpAsset version of this test, because this family always has z,e.

n_a1=n_a(1); n_a2_1=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_1_grid=a_grid(n_a1+1:end);

% The inert second experience asset. Its grid values are irrelevant (a2_2'=a2_2); all that matters
% is that it has more than one point, so that 'flat in a2_2' is a real check.
n_a2_2=3;
a2_2_grid=linspace(0,1,n_a2_2)';

%% Model A: one experience asset
ReturnFn_A=@(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_ExpAssetze_nod1_z_e(d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceassetze=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptionsA.n_e=vfoptionsbaseline.n_e; vfoptionsA.pi_e=vfoptionsbaseline.pi_e; vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
simoptionsA=struct(); simoptionsA.experienceassetze=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn;
simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid; simoptionsA.z_grid=z_grid;
simoptionsA.n_e=simoptionsbaseline.n_e; simoptionsA.pi_e=simoptionsbaseline.pi_e; simoptionsA.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_A=zeros([n_a1,n_a2_1,n_z,vfoptionsA.n_e],'gpuArray'); jequaloneDist_A(1,1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model B: two experience assets, the second one inert
n_a_B=[n_a1,n_a2_1,n_a2_2];
a_grid_B=[a1_grid;a2_1_grid;a2_2_grid];
% whicha=1 -> the Model A law of motion; whicha=2 -> identity, so a2_2 never moves
aprimeFn_B=@(d2,a2_1,a2_2,z,e,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)*z*e+(1-phi2)*a2_1)+(whicha==2)*a2_2;
ReturnFn_B=@(d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate) ReturnFn_ExpAssetze_nod1_z_e_with2A2(d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate);
Params_B=Params; Params_B.pensionrate=0; % so a2_2 does not enter the return function
vfoptionsB=struct(); vfoptionsB.experienceassetze=2; vfoptionsB.aprimeFn=aprimeFn_B;
vfoptionsB.n_e=vfoptionsbaseline.n_e; vfoptionsB.pi_e=vfoptionsbaseline.pi_e; vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
simoptionsB=struct(); simoptionsB.experienceassetze=2; simoptionsB.aprimeFn=aprimeFn_B;
simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B; simoptionsB.z_grid=z_grid;
simoptionsB.n_e=simoptionsbaseline.n_e; simoptionsB.pi_e=simoptionsbaseline.pi_e; simoptionsB.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_B=zeros([n_a1,n_a2_1,n_a2_2,n_z,vfoptionsB.n_e],'gpuArray');
jequaloneDist_B(1,1,1,ceil(n_z/2),ceil(vfoptionsB.n_e/2))=1; % all mass at a2_2 index 1, where it stays
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params_B,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params_B,simoptionsB);

%% Checks
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_2
    Vslice=V_B(:,:,kk,:,:,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_B(:,:,1,:,:,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_B(:,:,:,kk,:,:,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_B(:,:,1,:,:,:);
fprintf('Cross test 7 (with2A2 ze, nod1): V is flat in the inert a2_2, this should be zero: %.3e \n',flat)
fprintf('Cross test 7 (with2A2 ze, nod1): inert a2_2 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 7 (with2A2 ze, nod1): agent dist at the inert a2_2 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))
fprintf('Cross test 7 (with2A2 ze, nod1): no mass leaked off the inert a2_2 point, this should be zero: %.3e \n',abs(sum(StationaryDist_B(:))-sum(Dslice(:))))

output=struct();

end
