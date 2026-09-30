function output=CoreFHorzExpAsset_CrossTests6_d1_nosemiz_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2): the mirror image of CrossTests5 -- here the FIRST experience asset is the
% inert one and the SECOND carries the real law of motion.
%   Model A: experienceasset=1, n_a=[n_a1,n_a2_1]
%   Model C: experienceasset=2, n_a=[n_a1,n_a2_1inert,n_a2_2], a2_1'=a2_1 (identity), a2_2' is the
%            Model A law of motion, and the return function reads earnings off a2_2.
% V_C must be exactly flat in the a2_1 dimension and equal V_A at every a2_1 slice.
% CrossTests5 and CrossTests6 together pin down the ordering of the two a2 dimensions: either one
% on its own would still pass if the nested interpolation had them transposed.
% Inputs: the one-experience-asset with-a1 grids, n_a=[n_a1,n_a2_1], a_grid=[a1_grid;a2_1_grid].

n_a1=n_a(1); n_a2=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_grid=a_grid(n_a1+1:end);

% The inert FIRST experience asset.
n_a2_1inert=3;
a2_1inert_grid=linspace(0,1,n_a2_1inert)';

%% Model A: one experience asset
ReturnFn_A=@(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceasset=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptionsA=struct(); simoptionsA.experienceasset=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn; simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid;
jequaloneDist_A=zeros([n_a1,n_a2],'gpuArray'); jequaloneDist_A(1,1)=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,1,Params,simoptionsA);

%% Model C: two experience assets, the FIRST one inert
n_a_C=[n_a1,n_a2_1inert,n_a2];
a_grid_C=[a1_grid;a2_1inert_grid;a2_grid];
% whicha=1 -> identity, so a2_1 never moves; whicha=2 -> the Model A law of motion
aprimeFn_C=@(d2,a2_1,a2_2,whicha,phi1,phi2) (whicha==1)*a2_1+(whicha==2)*(phi1*(1-d2)+(1-phi2)*a2_2);
ReturnFn_C=@(d1,d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz_2A2swap(d1,d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsC=struct(); vfoptionsC.experienceasset=2; vfoptionsC.aprimeFn=aprimeFn_C;
simoptionsC=struct(); simoptionsC.experienceasset=2; simoptionsC.aprimeFn=aprimeFn_C; simoptionsC.d_grid=d_grid; simoptionsC.a_grid=a_grid_C;
jequaloneDist_C=zeros([n_a1,n_a2_1inert,n_a2],'gpuArray'); jequaloneDist_C(1,1,1)=1; % all mass at a2_1 index 1, where it stays
[V_C,Policy_C]=ValueFnIter_Case1_FHorz(n_d,n_a_C,0,N_j,d_grid,a_grid_C,[],[],ReturnFn_C,Params,DiscountFactorParamNames,[],vfoptionsC);
StationaryDist_C=StationaryDist_FHorz_Case1(jequaloneDist_C,AgeWeightParamNames,Policy_C,n_d,n_a_C,0,N_j,1,Params,simoptionsC);

%% Checks
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_1inert
    Vslice=V_C(:,kk,:,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_C(:,1,:,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_C(:,:,kk,:,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_C(:,1,:,:);
fprintf('Cross test 6 (with2A2, d1): V is flat in the inert a2_1, this should be zero: %.3e \n',flat)
fprintf('Cross test 6 (with2A2, d1): inert a2_1 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 6 (with2A2, d1): agent dist at the inert a2_1 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))
fprintf('Cross test 6 (with2A2, d1): no mass leaked off the inert a2_1 point, this should be zero: %.3e \n',abs(sum(StationaryDist_C(:))-sum(Dslice(:))))

output=struct();
end
