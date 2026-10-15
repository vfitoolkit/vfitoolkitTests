function output=CoreFHorzExpAssetsemiz_CrossTests6_nod1_with2A2(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-test (with2A2, experienceassetsemiz): an INERT experience asset must reduce the
% two-experience-asset model back to the one-experience-asset with-a1 model.
%   Model A: experienceassetsemiz=1, n_a=[n_a1,n_a2_1]
%   Model B: experienceassetsemiz=2 with a2_2 held by the identity map a2'=a2, and
%            pensionrate=0 so a2_2 never enters the return function.
% V_B must be exactly flat in a2_2 and equal V_A at every slice of it. Run together with
% its sibling cross-test (which makes the OTHER asset the inert one), these pin down the
% ORDERING of the two a2 dimensions: either alone would still pass if the nested
% interpolation had them transposed.
% NOTE the aprimeFn argument order for experienceassetsemiz: the 'whicha' selector comes AFTER
% semiz, i.e. aprimeFn(d2,a2_1,a2_2,semiz,whicha,params). Both experience assets are driven by
% the semi-exogenous state.
%
% This family always carries semiz, and semiz is the FAST exogenous index, so V is
% (a1, a2_1[, a2_2], semiz, z, e, j) and the a2 slices sit before the shock subscripts.

n_a1=n_a(1); n_a2_1=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_1_grid=a_grid(n_a1+1:end);

% The inert experience asset. Its grid values are irrelevant (a2'=a2); what matters is that
% it has more than one point, so that 'flat in it' is a real check.
n_a2_inert=3;
a2_inert_grid=linspace(0,1,n_a2_inert)';

%% Model A: one experience asset
ReturnFn_A=@(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_ExpAssetsemiz_nod1_z_e(d2,d3,a1prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfoptionsA=struct(); vfoptionsA.experienceassetsemiz=1; vfoptionsA.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptionsA.n_semiz=vfoptionsbaseline.n_semiz; vfoptionsA.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptionsA.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptionsA.n_e=vfoptionsbaseline.n_e; vfoptionsA.pi_e=vfoptionsbaseline.pi_e; vfoptionsA.e_grid=vfoptionsbaseline.e_grid;
simoptionsA=struct(); simoptionsA.experienceassetsemiz=1; simoptionsA.aprimeFn=vfoptionsA.aprimeFn;
simoptionsA.n_semiz=simoptionsbaseline.n_semiz; simoptionsA.semiz_grid=simoptionsbaseline.semiz_grid; simoptionsA.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsA.n_e=simoptionsbaseline.n_e; simoptionsA.pi_e=simoptionsbaseline.pi_e; simoptionsA.e_grid=simoptionsbaseline.e_grid;
simoptionsA.d_grid=d_grid; simoptionsA.a_grid=a_grid; simoptionsA.z_grid=z_grid;
jequaloneDist_A=zeros([n_a1,n_a2_1,vfoptionsA.n_semiz,n_z,vfoptionsA.n_e],'gpuArray');
jequaloneDist_A(1,1,ceil(vfoptionsA.n_semiz/2),ceil(n_z/2),ceil(vfoptionsA.n_e/2))=1;
[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

%% Model B: two experience assets, one of them inert
n_a_B=[n_a1,n_a2_1,n_a2_inert];
a_grid_B=[a1_grid;a2_1_grid;a2_inert_grid];
aprimeFn_B=@(d2,a2_1,a2_2,semiz,whicha,phi1,phi2) (whicha==1)*(phi1*(1-d2)*semiz+(1-phi2)*a2_1)+(whicha==2)*a2_2;
ReturnFn_B=@(d2,d3,a1prime,a1,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost,pensionrate) ReturnFn_ExpAssetsemiz_nod1_z_e_with2A2(d2,d3,a1prime,a1,a2_1,a2_2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost,pensionrate);
Params_B=Params; Params_B.pensionrate=0; % so the inert asset does not enter the return function
vfoptionsB=struct(); vfoptionsB.experienceassetsemiz=2; vfoptionsB.aprimeFn=aprimeFn_B;
vfoptionsB.n_semiz=vfoptionsbaseline.n_semiz; vfoptionsB.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptionsB.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfoptionsB.n_e=vfoptionsbaseline.n_e; vfoptionsB.pi_e=vfoptionsbaseline.pi_e; vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
simoptionsB=struct(); simoptionsB.experienceassetsemiz=2; simoptionsB.aprimeFn=vfoptionsB.aprimeFn;
simoptionsB.n_semiz=simoptionsbaseline.n_semiz; simoptionsB.semiz_grid=simoptionsbaseline.semiz_grid; simoptionsB.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptionsB.n_e=simoptionsbaseline.n_e; simoptionsB.pi_e=simoptionsbaseline.pi_e; simoptionsB.e_grid=simoptionsbaseline.e_grid;
simoptionsB.d_grid=d_grid; simoptionsB.a_grid=a_grid_B; simoptionsB.z_grid=z_grid;
jequaloneDist_B=zeros([n_a1,n_a2_1,n_a2_inert,vfoptionsB.n_semiz,n_z,vfoptionsB.n_e],'gpuArray');
jequaloneDist_B(1,1,1,ceil(vfoptionsB.n_semiz/2),ceil(n_z/2),ceil(vfoptionsB.n_e/2))=1; % all mass at the inert point, where it stays
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params_B,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params_B,simoptionsB);

%% Checks
flat=0; matchV=0; matchP=0;
for kk=1:n_a2_inert
    Vslice=V_B(:,:,kk,:,:,:,:);
    flat=max(flat,max(abs(Vslice(:)-reshape(V_B(:,:,1,:,:,:,:),size(Vslice(:))))));
    matchV=max(matchV,max(abs(Vslice(:)-V_A(:))));
    Pslice=Policy_B(:,:,:,kk,:,:,:,:);
    matchP=max(matchP,max(abs(Pslice(:)-Policy_A(:))));
end
Dslice=StationaryDist_B(:,:,1,:,:,:,:);
fprintf('Cross test 6 (with2A2 semiz, nod1): V is flat in the inert a2_2, this should be zero: %.3e \n',flat)
fprintf('Cross test 6 (with2A2 semiz, nod1): inert a2_2 reduces to one experience asset, this should be zero: V %.3e, Policy %.3e \n',matchV,matchP)
fprintf('Cross test 6 (with2A2 semiz, nod1): agent dist at the inert a2_2 point, this should be zero: %.3e \n',max(abs(Dslice(:)-StationaryDist_A(:))))
fprintf('Cross test 6 (with2A2 semiz, nod1): no mass leaked off the inert a2_2 point, this should be zero: %.3e \n',abs(sum(StationaryDist_B(:))-sum(Dslice(:))))

output=struct();

end
