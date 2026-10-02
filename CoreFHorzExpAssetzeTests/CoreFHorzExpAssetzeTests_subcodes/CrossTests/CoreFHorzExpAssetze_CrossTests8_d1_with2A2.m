function output=CoreFHorzExpAssetze_CrossTests8_d1_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, experienceassetze): the mirror image of CrossTests7 -- here the FIRST
% experience asset is the inert one and the SECOND carries the real law of motion.
%   Model A: experienceassetze=1, n_a=[n_a1,n_a2]
%   Model C: experienceassetze=2, n_a=[n_a1,n_a2_1inert,n_a2], a2_1'=a2_1 (identity), a2_2' is the
%            Model A law of motion, and the return function reads earnings off a2_2 (_2A2swap).
% V_C must be exactly flat in the a2_1 dimension and equal V_A at every a2_1 slice.
% CrossTests7 and CrossTests8 together pin down the ordering of the two a2 dimensions: either one
% on its own would still pass if the nested interpolation had them transposed.
%
% NOTE the aprimeFn argument order for experienceassetZE: 'whicha' comes AFTER z and e.

n_a1=n_a(1); n_a2=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_grid=a_grid(n_a1+1:end);

% The inert FIRST experience asset.
n_a2_1inert=3;
a2_1inert_grid=linspace(0,1,n_a2_1inert)';

%% Model A: one experience asset
ReturnFn_A=@(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssetze_d1_z_e(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceassetze=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptionsA.n_e=vfoptionsbaseline.n_e; vfoptionsA.pi_e=vfoptionsbaseline.pi_e; vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
simoptionsA=struct(); simoptionsA.experienceassetze=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn;
simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid; simoptionsA.z_grid=z_grid;
simoptionsA.n_e=simoptionsbaseline.n_e; simoptionsA.pi_e=simoptionsbaseline.pi_e; simoptionsA.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_A=zeros([n_a1,n_a2,n_z,vfoptionsA.n_e],'gpuArray'); jequaloneDist_A(1,1,ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model C: two experience assets, the FIRST one inert
n_a_C=[n_a1,n_a2_1inert,n_a2];
a_grid_C=[a1_grid;a2_1inert_grid;a2_grid];
% whicha=1 -> identity, so a2_1 never moves; whicha=2 -> the Model A law of motion
aprimeFn_C=@(d2,a2_1,a2_2,z,e,whicha,phi1,phi2) (whicha==1)*a2_1+(whicha==2)*(phi1*(1-d2)*z*e+(1-phi2)*a2_2);
ReturnFn_C=@(d1,d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_ExpAssetze_d1_z_e_2A2swap(d1,d2,a1prime,a1,a2_1,a2_2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsC=struct(); vfoptionsC.experienceassetze=2; vfoptionsC.aprimeFn=aprimeFn_C;
vfoptionsC.n_e=vfoptionsbaseline.n_e; vfoptionsC.pi_e=vfoptionsbaseline.pi_e; vfoptionsC.e_grid=vfoptionsbaseline.e_grid;
simoptionsC=struct(); simoptionsC.experienceassetze=2; simoptionsC.aprimeFn=aprimeFn_C;
simoptionsC.d_grid=d_grid; simoptionsC.a_grid=a_grid_C; simoptionsC.z_grid=z_grid;
simoptionsC.n_e=simoptionsbaseline.n_e; simoptionsC.pi_e=simoptionsbaseline.pi_e; simoptionsC.e_grid=simoptionsbaseline.e_grid;
jequaloneDist_C=zeros([n_a1,n_a2_1inert,n_a2,n_z,vfoptionsC.n_e],'gpuArray');
jequaloneDist_C(1,1,1,ceil(n_z/2),ceil(vfoptionsC.n_e/2))=1; % all mass at a2_1 index 1, where it stays
[V_C,Policy_C]=ValueFnIter_Case1_FHorz(n_d,n_a_C,n_z,N_j,d_grid,a_grid_C,z_grid,pi_z,ReturnFn_C,Params,DiscountFactorParamNames,[],vfoptionsC);
StationaryDist_C=StationaryDist_FHorz_Case1(jequaloneDist_C,AgeWeightParamNames,Policy_C,n_d,n_a_C,n_z,N_j,pi_z,Params,simoptionsC);

%% Checks
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_1inert
    Vslice=V_C(:,kk,:,:,:,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_C(:,1,:,:,:,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_C(:,:,kk,:,:,:,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_C(:,1,:,:,:,:);
fprintf('Cross test 8 (with2A2 ze, d1): V is flat in the inert a2_1, this should be zero: %.3e \n',flat)
fprintf('Cross test 8 (with2A2 ze, d1): inert a2_1 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 8 (with2A2 ze, d1): agent dist at the inert a2_1 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))
fprintf('Cross test 8 (with2A2 ze, d1): no mass leaked off the inert a2_1 point, this should be zero: %.3e \n',abs(sum(StationaryDist_C(:))-sum(Dslice(:))))

output=struct();

end
