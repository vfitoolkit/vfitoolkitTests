function output=CoreFHorzExpAsset_CrossTests5_d1_noa1_nosemiz_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, noa1): an INERT SECOND experience asset must reduce the two-experience-asset
% model back to the one-experience-asset noa1 model. Same idea as the with-a1 CrossTests5, but
% with no standard endogenous asset at all, so the two experience assets are the only endogenous
% states (N_a1==0, a separate code path in every raw).
% Inputs: the noa1 grids, n_a=n_a_justexpasset (scalar), a_grid=a2_grid.

n_a2_1=n_a(1);
a2_1_grid=a_grid;

% The inert second experience asset.
n_a2_2=3;
a2_2_grid=linspace(0,1,n_a2_2)';

%% Model A: one experience asset, no a1
ReturnFn_A=@(d1,d2,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_noa1_nosemiz(d1,d2,a,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceasset=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptionsA=struct(); simoptionsA.experienceasset=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn; simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid;
jequaloneDist_A=zeros(n_a2_1,1,'gpuArray'); jequaloneDist_A(1)=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,1,Params,simoptionsA);

%% Model B: two experience assets, the second one inert, still no a1
n_a_B=[n_a2_1,n_a2_2];
a_grid_B=[a2_1_grid;a2_2_grid];
aprimeFn_B=@(d2,a2_1,a2_2,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)+(1-phi2)*a2_1)+(whicha==2)*a2_2;
ReturnFn_B=@(d1,d2,a2_1,a2_2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,pensionrate) ReturnFn_d1_noz_noe_noa1_nosemiz_with2A2(d1,d2,a2_1,a2_2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension,pensionrate);
Params_B=Params; Params_B.pensionrate=0; % so a2_2 does not enter the return function
vfoptionsB=struct(); vfoptionsB.experienceasset=2; vfoptionsB.aprimeFn=aprimeFn_B;
simoptionsB=struct(); simoptionsB.experienceasset=2; simoptionsB.aprimeFn=aprimeFn_B; simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B;
jequaloneDist_B=zeros([n_a2_1,n_a2_2],'gpuArray'); jequaloneDist_B(1,1)=1;
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params_B,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,1,Params_B,simoptionsB);

%% Checks
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_2
    Vslice=V_B(:,kk,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_B(:,1,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_B(:,:,kk,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_B(:,1,:);
fprintf('Cross test 5 (with2A2, noa1, d1): V is flat in the inert a2_2, this should be zero: %.3e \n',flat)
fprintf('Cross test 5 (with2A2, noa1, d1): inert a2_2 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 5 (with2A2, noa1, d1): agent dist at the inert a2_2 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))

output=struct();
end
