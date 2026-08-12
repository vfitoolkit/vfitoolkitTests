function output=QHDFHorz_CrossTests_nod1_semiz_with2A(n_d_semiz,n_d2_semiz,n_a,n_a_big,n_z,N_j,d_grid_semiz,d2_grid_semiz,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross tests for the two-endo SemiExo (nod1) code paths.
% Test 1: d2-driven semiz mirrors age-alternating pi_z_J
% Test 2: single-point d2 with semiz reduces to ordinary markov z
% Test 3: collapsing z to 1pt in (semiz+z+e) matches (semiz+e)
% Test 4: collapsing e to 1pt in (semiz+z+e) matches (semiz+z)

n_d=0; % no decision in the z-only comparison models
d_grid=[];

% Common: semiz_grid is [0;1] (from baseline)
semiz_grid=vfoptionsbaseline.semiz_grid;

%% Quasi-Hyperbolic Discounting
% vfoptions.exoticpreferences='QuasiHyperbolic';
%% Naive
% vfoptions.quasi_hyperbolic='Naive';

%% Test 1: d2-driven semiz vs age-alternating pi_z_J
% Use 2-state z that lives on the same grid as semiz ([0;1]).
% pi_z_J alternates between pi_odd (odd jj) and pi_even (even jj).
% In the semiz version, n_d2=2 with d2_grid=[1;2]. pi_semiz(:,:,1)=pi_odd
% and pi_semiz(:,:,2)=pi_even (constant in jj). ReturnFn forces optimal
% d2=1 at odd ages, d2=2 at even ages, so the realised semiz transitions
% match pi_z_J age-by-age.

n_z2=2;
z_grid2=semiz_grid; % [0;1]
pi_odd=[0.8,0.2; 0.3,0.7];
pi_even=[0.5,0.5; 0.1,0.9];
pi_z_J=zeros(n_z2,n_z2,N_j,'gpuArray');
for jj=1:N_j
    if mod(jj,2)==1
        pi_z_J(:,:,jj)=pi_odd;
    else
        pi_z_J(:,:,jj)=pi_even;
    end
end

% z-only model (no semiz, no d)
ReturnFn_z_test2=@(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
vfopts_z=struct();
vfopts_z.exoticpreferences='QuasiHyperbolic';
vfopts_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts_z.quasi_hyperbolic='Naive';
[Vz,Polz]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z2,N_j,d_grid,a_grid,z_grid2,pi_z_J,ReturnFn_z_test2,Params,DiscountFactorParamNames,[],vfopts_z);

% semiz model (n_d2=2, d2_grid=[1;2]); ReturnFn forces d2 by age parity
n_d2_t2=2;
d2_grid_t2=[1;2];
vfopts_sz=struct();
vfopts_sz.n_semiz=n_z2;
vfopts_sz.semiz_grid=z_grid2;
vfopts_sz.SemiExoStateFn=@(n,nprime,dsemiz,p_odd_00,p_odd_01,p_odd_10,p_odd_11,p_even_00,p_even_01,p_even_10,p_even_11) CoreFHorzTwoEndoSetup_SemiExoStateFn_TwoMarkovs(n,nprime,dsemiz,p_odd_00,p_odd_01,p_odd_10,p_odd_11,p_even_00,p_even_01,p_even_10,p_even_11);
Params.p_odd_00=pi_odd(1,1); Params.p_odd_01=pi_odd(1,2); Params.p_odd_10=pi_odd(2,1); Params.p_odd_11=pi_odd(2,2);
Params.p_even_00=pi_even(1,1); Params.p_even_01=pi_even(1,2); Params.p_even_10=pi_even(2,1); Params.p_even_11=pi_even(2,2);
vfopts_sz.exoticpreferences='QuasiHyperbolic';
vfopts_sz.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts_sz.quasi_hyperbolic='Naive';

ReturnFn_sz_test2=@(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod1_noz_noe_semiz_test2_with2A(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
[Vsz,Polsz]=ValueFnIter_Case1_FHorz(n_d2_t2,n_a,0,N_j,d2_grid_t2,a_grid,[],[],ReturnFn_sz_test2,Params,DiscountFactorParamNames,[],vfopts_sz);

% V should match: both are (a1,a2,semiz/z,N_j) of shape [n_a,n_z2,N_j]
fprintf('Test 1 (d2-driven semiz vs alternating pi_z_J), should be zero: %2.8f \n',max(abs(Vz(:)-Vsz(:))))


%% Test 2: collapse z to 1pt in (semiz+z+e) -> matches (semiz+e)
% Run nod1_z_e_semiz with n_z=1, z_grid=1, pi_z=1; compare to nod1_noz_e_semiz.
% Consumption: w*kappa_j*semiz*z*e + uempbenefit*(1-semiz). With z=1, this is
% w*kappa_j*semiz*e + uempbenefit*(1-semiz), matching nod1_noz_e_semiz.

vfopts5=struct();
vfopts5.n_semiz=vfoptionsbaseline.n_semiz;
vfopts5.semiz_grid=vfoptionsbaseline.semiz_grid;
vfopts5.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfopts5.n_e=vfoptionsbaseline.n_e;
vfopts5.e_grid=vfoptionsbaseline.e_grid;
vfopts5.pi_e=vfoptionsbaseline.pi_e;

ReturnFn_5_zfull=@(d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_z_e_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);
ReturnFn_5_noz=@(d2,a1prime,a2prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_noz_e_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);

% z collapsed to 1pt
vfopts5.exoticpreferences='QuasiHyperbolic';
vfopts5.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts5.quasi_hyperbolic='Naive';
[V5a_zfull,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,1,N_j,d2_grid_semiz,a_grid,1,1,ReturnFn_5_zfull,Params,DiscountFactorParamNames,[],vfopts5);
% No z at all
vfopts5.exoticpreferences='QuasiHyperbolic';
vfopts5.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts5.quasi_hyperbolic='Naive';
[V5a_noz,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,0,N_j,d2_grid_semiz,a_grid,[],[],ReturnFn_5_noz,Params,DiscountFactorParamNames,[],vfopts5);

V5a_zfull=squeeze(V5a_zfull);
fprintf('Test 2 (collapse z=1pt vs no z), should be zero: %2.8f \n',max(abs(V5a_zfull(:)-V5a_noz(:))))


%% Test 3: collapse e to 1pt in (semiz+z+e) -> matches (semiz+z)
vfopts5b=struct();
vfopts5b.n_semiz=vfoptionsbaseline.n_semiz;
vfopts5b.semiz_grid=vfoptionsbaseline.semiz_grid;
vfopts5b.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfopts5b.exoticpreferences='QuasiHyperbolic';
vfopts5b.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts5b.quasi_hyperbolic='Naive';

ReturnFn_5_noe=@(d2,a1prime,a2prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_z_noe_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);

% e collapsed to 1pt (n_e=1, e_grid=1, pi_e=1)
vfopts5b_efull=vfopts5b;
vfopts5b_efull.n_e=1; vfopts5b_efull.e_grid=1; vfopts5b_efull.pi_e=1;
[V5b_efull,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn_5_zfull,Params,DiscountFactorParamNames,[],vfopts5b_efull);

% No e at all (semiz+z only)
[V5b_noe,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn_5_noe,Params,DiscountFactorParamNames,[],vfopts5b);

V5b_efull=squeeze(V5b_efull);
fprintf('Test 3 (collapse e=1pt vs no e), should be zero: %2.8f \n',max(abs(V5b_efull(:)-V5b_noe(:))))



%% Sophisticated

%% Test 1: d2-driven semiz vs age-alternating pi_z_J
n_z2=2;
z_grid2=semiz_grid; % [0;1]
pi_odd=[0.8,0.2; 0.3,0.7];
pi_even=[0.5,0.5; 0.1,0.9];
pi_z_J=zeros(n_z2,n_z2,N_j,'gpuArray');
for jj=1:N_j
    if mod(jj,2)==1
        pi_z_J(:,:,jj)=pi_odd;
    else
        pi_z_J(:,:,jj)=pi_even;
    end
end

% z-only model (no semiz, no d)
ReturnFn_z_test2=@(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
vfopts_z=struct();
vfopts_z.exoticpreferences='QuasiHyperbolic';
vfopts_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts_z.quasi_hyperbolic='Sophisticated';
[Vz,Polz]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z2,N_j,d_grid,a_grid,z_grid2,pi_z_J,ReturnFn_z_test2,Params,DiscountFactorParamNames,[],vfopts_z);

% semiz model (n_d2=2, d2_grid=[1;2]); ReturnFn forces d2 by age parity
n_d2_t2=2;
d2_grid_t2=[1;2];
vfopts_sz=struct();
vfopts_sz.n_semiz=n_z2;
vfopts_sz.semiz_grid=z_grid2;
vfopts_sz.SemiExoStateFn=@(n,nprime,dsemiz,p_odd_00,p_odd_01,p_odd_10,p_odd_11,p_even_00,p_even_01,p_even_10,p_even_11) CoreFHorzTwoEndoSetup_SemiExoStateFn_TwoMarkovs(n,nprime,dsemiz,p_odd_00,p_odd_01,p_odd_10,p_odd_11,p_even_00,p_even_01,p_even_10,p_even_11);
Params.p_odd_00=pi_odd(1,1); Params.p_odd_01=pi_odd(1,2); Params.p_odd_10=pi_odd(2,1); Params.p_odd_11=pi_odd(2,2);
Params.p_even_00=pi_even(1,1); Params.p_even_01=pi_even(1,2); Params.p_even_10=pi_even(2,1); Params.p_even_11=pi_even(2,2);
vfopts_sz.exoticpreferences='QuasiHyperbolic';
vfopts_sz.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts_sz.quasi_hyperbolic='Sophisticated';

ReturnFn_sz_test2=@(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod1_noz_noe_semiz_test2_with2A(d2,a1prime,a2prime,a1,a2,semiz,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
[Vsz,Polsz]=ValueFnIter_Case1_FHorz(n_d2_t2,n_a,0,N_j,d2_grid_t2,a_grid,[],[],ReturnFn_sz_test2,Params,DiscountFactorParamNames,[],vfopts_sz);

% V should match: both are (a1,a2,semiz/z,N_j) of shape [n_a,n_z2,N_j]
fprintf('Test 1 (d2-driven semiz vs alternating pi_z_J), should be zero: %2.8f \n',max(abs(Vz(:)-Vsz(:))))


%% Test 2: collapse z to 1pt in (semiz+z+e) -> matches (semiz+e)
vfopts5=struct();
vfopts5.n_semiz=vfoptionsbaseline.n_semiz;
vfopts5.semiz_grid=vfoptionsbaseline.semiz_grid;
vfopts5.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfopts5.n_e=vfoptionsbaseline.n_e;
vfopts5.e_grid=vfoptionsbaseline.e_grid;
vfopts5.pi_e=vfoptionsbaseline.pi_e;

ReturnFn_5_zfull=@(d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_z_e_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);
ReturnFn_5_noz=@(d2,a1prime,a2prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_noz_e_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);

% z collapsed to 1pt
vfopts5.exoticpreferences='QuasiHyperbolic';
vfopts5.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts5.quasi_hyperbolic='Sophisticated';
[V5a_zfull,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,1,N_j,d2_grid_semiz,a_grid,1,1,ReturnFn_5_zfull,Params,DiscountFactorParamNames,[],vfopts5);
% No z at all
vfopts5.exoticpreferences='QuasiHyperbolic';
vfopts5.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts5.quasi_hyperbolic='Sophisticated';
[V5a_noz,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,0,N_j,d2_grid_semiz,a_grid,[],[],ReturnFn_5_noz,Params,DiscountFactorParamNames,[],vfopts5);

V5a_zfull=squeeze(V5a_zfull);
fprintf('Test 2 (collapse z=1pt vs no z), should be zero: %2.8f \n',max(abs(V5a_zfull(:)-V5a_noz(:))))


%% Test 3: collapse e to 1pt in (semiz+z+e) -> matches (semiz+z)
vfopts5b=struct();
vfopts5b.n_semiz=vfoptionsbaseline.n_semiz;
vfopts5b.semiz_grid=vfoptionsbaseline.semiz_grid;
vfopts5b.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
vfopts5b.exoticpreferences='QuasiHyperbolic';
vfopts5b.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfopts5b.quasi_hyperbolic='Sophisticated';

ReturnFn_5_noe=@(d2,a1prime,a2prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2) ReturnFn_nod1_z_noe_semiz_with2A(d2,a1prime,a2prime,a1,a2,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost,phi1,phi2);

% e collapsed to 1pt (n_e=1, e_grid=1, pi_e=1)
vfopts5b_efull=vfopts5b;
vfopts5b_efull.n_e=1; vfopts5b_efull.e_grid=1; vfopts5b_efull.pi_e=1;
[V5b_efull,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn_5_zfull,Params,DiscountFactorParamNames,[],vfopts5b_efull);

% No e at all (semiz+z only)
[V5b_noe,~]=ValueFnIter_Case1_FHorz(n_d2_semiz,n_a,n_z,N_j,d2_grid_semiz,a_grid,z_grid,pi_z,ReturnFn_5_noe,Params,DiscountFactorParamNames,[],vfopts5b);

V5b_efull=squeeze(V5b_efull);
fprintf('Test 3 (collapse e=1pt vs no e), should be zero: %2.8f \n',max(abs(V5b_efull(:)-V5b_noe(:))))


%%
output=struct(); % not currently used

end
