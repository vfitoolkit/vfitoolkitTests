function output=CoreFHorzExpAsset_CrossTests5_nod1_nosemiz_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2): an INERT SECOND experience asset must reduce the two-experience-asset
% model back to the one-experience-asset with-a1 model.
%   Model A: experienceasset=1, n_a=[n_a1,n_a2_1]
%   Model B: experienceasset=2, n_a=[n_a1,n_a2_1,n_a2_2], with a2_2'=a2_2 (the identity map, which
%            the aprimeFn interpolation reproduces exactly for on-grid values) and pensionrate=0
%            so that a2_2 never enters the return function.
% V_B must be exactly flat in the a2_2 dimension, and must equal V_A at every a2_2 slice.
% Together with CrossTests6 (which makes the FIRST asset the inert one) this pins down the
% ordering of the two a2 dimensions: CrossTests5 alone would still pass if they were transposed.
% These are the machine-precision checks of the tier; the figures are consistency checks only.
% Inputs: the one-experience-asset with-a1 grids, n_a=[n_a1,n_a2_1], a_grid=[a1_grid;a2_1_grid].

n_a1=n_a(1); n_a2_1=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_1_grid=a_grid(n_a1+1:end);

% The inert second experience asset. Its grid values are irrelevant (a2_2'=a2_2), only that it
% has more than one point, so that 'flat in a2_2' is a real check.
n_a2_2=3;
a2_2_grid=linspace(0,1,n_a2_2)';

%% Model A: one experience asset
ReturnFn_A=@(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(d2,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension);
vfoptionsA=struct(); vfoptionsA.experienceasset=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
simoptionsA=struct(); simoptionsA.experienceasset=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn; simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid;
jequaloneDist_A=zeros([n_a1,n_a2_1],'gpuArray'); jequaloneDist_A(1,1)=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,1,Params,simoptionsA);

%% Model B: two experience assets, the second one inert
n_a_B=[n_a1,n_a2_1,n_a2_2];
a_grid_B=[a1_grid;a2_1_grid;a2_2_grid];
% whicha=1 -> the Model A law of motion; whicha=2 -> identity, so a2_2 never moves
aprimeFn_B=@(d2,a2_1,a2_2,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)+(1-phi2)*a2_1)+(whicha==2)*a2_2;
ReturnFn_B=@(d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate) ReturnFn_nod1_noz_noe_nosemiz_with2A2(d2,a1prime,a1,a2_1,a2_2,r,w,kappa_j,sigma,agej,Jr,pension,pensionrate);
Params_B=Params; Params_B.pensionrate=0; % so a2_2 does not enter the return function
vfoptionsB=struct(); vfoptionsB.experienceasset=2; vfoptionsB.aprimeFn=aprimeFn_B;
simoptionsB=struct(); simoptionsB.experienceasset=2; simoptionsB.aprimeFn=aprimeFn_B; simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B;
jequaloneDist_B=zeros([n_a1,n_a2_1,n_a2_2],'gpuArray'); jequaloneDist_B(1,1,1)=1; % all mass at a2_2 index 1, where it stays
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params_B,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,1,Params_B,simoptionsB);

%% Checks
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_2
    Vslice=V_B(:,:,kk,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_B(:,:,1,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_B(:,:,:,kk,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_B(:,:,1,:);
fprintf('Cross test 5 (with2A2, nod1): V is flat in the inert a2_2, this should be zero: %.3e \n',flat)
fprintf('Cross test 5 (with2A2, nod1): inert a2_2 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 5 (with2A2, nod1): agent dist at the inert a2_2 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))
fprintf('Cross test 5 (with2A2, nod1): no mass leaked off the inert a2_2 point, this should be zero: %.3e \n',abs(sum(StationaryDist_B(:))-sum(Dslice(:))))

output=struct();
end
