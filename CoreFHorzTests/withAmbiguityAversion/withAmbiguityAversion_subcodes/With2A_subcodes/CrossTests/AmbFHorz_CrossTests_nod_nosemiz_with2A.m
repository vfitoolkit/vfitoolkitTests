function output=AmbFHorz_CrossTests_nod_nosemiz_with2A(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% The with2A ambiguity-aversion cross tests. The main-tier cross tests already gate the pi-plumbing,
% so these have one job: pin each 2A ambiguity raw family to something known-good. Run in a markov-z
% and an iid-e flavour, all against exoticpreferences='None' — exact zeros, no tolerances.
% Test 1-2A: three identical priors = standard preferences, AT ALL FOUR TIERS (plain/DC2A/GI2A/DC2A+GI2A)
%           — pins each ambiguity 2A raw to its own exponential donor at the same tier.
% Test 3a/3b-2A: an unambiguously worse pi binds, either prior slot (plain tier; the maxmin logic was
%           proven at the main tier, this confirms it survives the 2A state space).
% Test 4-2A: age-varying n_ambiguity via V_Jplus1, AT ALL FOUR TIERS — day-one runtime coverage of the
%           DC2A/GI2A/DC2A_GI2A raws' V_Jplus1 branches.

n_d=0;
d_grid=[];

ReturnFn_z=@(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_z_noe_nosemiz_with2A(a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);
ReturnFn_e=@(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2) ReturnFn_nod_noz_e_nosemiz_with2A(a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension,phi1,phi2);

% The three baseline priors (prior 1 is the regular pi_z/pi_e)
ambiguity_pi_z=vfoptionsbaseline.ambiguity_pi_z; % [n_z,n_z,3]
ambiguity_pi_e=vfoptionsbaseline.ambiguity_pi_e; % [n_e,3]
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

half=N_j/2; % for cross test 4 (N_j=20)

%% ==================== Ambiguity on a markov z ====================

% Standard preferences (vNM) at each of the four tiers; used by cross tests 1-2A and 4-2A
vfoptions_none1=struct();
[Vstd1,Policystd1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none1);
vfoptions_none2=struct();
vfoptions_none2.divideandconquer=1;
[Vstd2,Policystd2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none2);
vfoptions_none3=struct();
vfoptions_none3.gridinterplayer=1;
vfoptions_none3.ngridinterp=5;
[Vstd3,Policystd3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none3);
vfoptions_none4=struct();
vfoptions_none4.divideandconquer=1;
vfoptions_none4.gridinterplayer=1;
vfoptions_none4.ngridinterp=5;
[Vstd4,Policystd4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none4);

%% Cross test 1-2A (z): three identical priors = standard preferences, at all four tiers
vfoptions_t11=struct();
vfoptions_t11.exoticpreferences='AmbiguityAversion';
vfoptions_t11.n_ambiguity=3;
vfoptions_t11.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V11,Policy11]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t11);
fprintf('Cross test 1-2A (z, plain): three identical priors, this should be zero: %.3e \n',max(abs(V11(:)-Vstd1(:))))
fprintf('Cross test 1-2A (z, plain): three identical priors, this should be zero: %.3e \n',max(abs(Policy11(:)-Policystd1(:))))
vfoptions_t12=struct();
vfoptions_t12.exoticpreferences='AmbiguityAversion';
vfoptions_t12.n_ambiguity=3;
vfoptions_t12.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
vfoptions_t12.divideandconquer=1;
[V12,Policy12]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t12);
fprintf('Cross test 1-2A (z, DC): three identical priors, this should be zero: %.3e \n',max(abs(V12(:)-Vstd2(:))))
fprintf('Cross test 1-2A (z, DC): three identical priors, this should be zero: %.3e \n',max(abs(Policy12(:)-Policystd2(:))))
vfoptions_t13=struct();
vfoptions_t13.exoticpreferences='AmbiguityAversion';
vfoptions_t13.n_ambiguity=3;
vfoptions_t13.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
vfoptions_t13.gridinterplayer=1;
vfoptions_t13.ngridinterp=5;
[V13,Policy13]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t13);
fprintf('Cross test 1-2A (z, GI): three identical priors, this should be zero: %.3e \n',max(abs(V13(:)-Vstd3(:))))
fprintf('Cross test 1-2A (z, GI): three identical priors, this should be zero: %.3e \n',max(abs(Policy13(:)-Policystd3(:))))
vfoptions_t14=struct();
vfoptions_t14.exoticpreferences='AmbiguityAversion';
vfoptions_t14.n_ambiguity=3;
vfoptions_t14.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
vfoptions_t14.divideandconquer=1;
vfoptions_t14.gridinterplayer=1;
vfoptions_t14.ngridinterp=5;
[V14,Policy14]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t14);
fprintf('Cross test 1-2A (z, DC+GI): three identical priors, this should be zero: %.3e \n',max(abs(V14(:)-Vstd4(:))))
fprintf('Cross test 1-2A (z, DC+GI): three identical priors, this should be zero: %.3e \n',max(abs(Policy14(:)-Policystd4(:))))
clear V11 V12 V13 V14 Policy11 Policy12 Policy13 Policy14

%% Cross tests 3a/3b-2A (z): an unambiguously worse pi binds, in either prior slot (plain tier)
% The worse pi shifts half of each row's best-column probability into the worst column
shiftz=0.5*pi_z(:,n_z);
pi_z_worse=pi_z;
pi_z_worse(:,n_z)=pi_z(:,n_z)-shiftz;
pi_z_worse(:,1)=pi_z(:,1)+shiftz;
% Standard preferences under the worse pi
[Vstdw,Policystdw]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_worse,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none1);
% 3a: worse pi in prior slot 1
vfoptions_t3=struct();
vfoptions_t3.exoticpreferences='AmbiguityAversion';
vfoptions_t3.n_ambiguity=2;
vfoptions_t3.ambiguity_pi_z=cat(3,pi_z_worse,pi_z);
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t3);
fprintf('Cross test 3a-2A (z): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(V3a(:)-Vstdw(:))))
fprintf('Cross test 3a-2A (z): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(Policy3a(:)-Policystdw(:))))
% 3b: worse pi in prior slot 2 (ordering of the priors is irrelevant)
vfoptions_t3.ambiguity_pi_z=cat(3,pi_z,pi_z_worse);
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t3);
fprintf('Cross test 3b-2A (z): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(V3b(:)-Vstdw(:))))
fprintf('Cross test 3b-2A (z): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(Policy3b(:)-Policystdw(:))))
clear V3a V3b Vstdw Policy3a Policy3b Policystdw

%% Cross test 4-2A (z): age-varying n_ambiguity via V_Jplus1, at all four tiers
% Three priors in the first half, a single (baseline) prior in the second half; the second half is
% then just vNM, and the first half must equal a short solve seeded with the vNM V at half+1.
% Slice jj of pi_z_J is the transition from period jj to jj+1, so the three priors sit in slices 1..half
ambiguity_pi_z_J=repmat(reshape(pi_z,[n_z,n_z,1,1]),[1,1,N_j,3]);
ambiguity_pi_z_J(:,:,1:half,2)=repmat(ambiguity_pi_z(:,:,2),[1,1,half]);
ambiguity_pi_z_J(:,:,1:half,3)=repmat(ambiguity_pi_z(:,:,3),[1,1,half]);
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
% (mewj is age-dependent but is only used for the agent distribution, which is not computed here)
% --- plain tier ---
vfoptions_t4a1=struct();
vfoptions_t4a1.exoticpreferences='AmbiguityAversion';
vfoptions_t4a1.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a1.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a1,Policy4a1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a1);
V4a_r=reshape(V4a1,[],N_j);
P4a_r=reshape(Policy4a1,[],N_j);
Vstd_r=reshape(Vstd1,[],N_j);
Pstd_r=reshape(Policystd1,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, plain): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, plain): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c1=struct();
vfoptions_t4c1.exoticpreferences='AmbiguityAversion';
vfoptions_t4c1.n_ambiguity=3;
vfoptions_t4c1.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c1.V_Jplus1=Vstd1(:,:,:,half+1);
[V4c1,Policy4c1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c1);
temp=V4a_r(:,1:half)-reshape(V4c1,[],Njs);
fprintf('Cross test 4-2A (z, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c1,[],Njs);
fprintf('Cross test 4-2A (z, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a1 V4c1 Policy4a1 Policy4c1
% --- DC tier ---
vfoptions_t4a2=struct();
vfoptions_t4a2.exoticpreferences='AmbiguityAversion';
vfoptions_t4a2.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a2.ambiguity_pi_z_J=ambiguity_pi_z_J;
vfoptions_t4a2.divideandconquer=1;
[V4a2,Policy4a2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a2);
V4a_r=reshape(V4a2,[],N_j);
P4a_r=reshape(Policy4a2,[],N_j);
Vstd_r=reshape(Vstd2,[],N_j);
Pstd_r=reshape(Policystd2,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c2=struct();
vfoptions_t4c2.exoticpreferences='AmbiguityAversion';
vfoptions_t4c2.n_ambiguity=3;
vfoptions_t4c2.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c2.divideandconquer=1;
vfoptions_t4c2.V_Jplus1=Vstd2(:,:,:,half+1);
[V4c2,Policy4c2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c2);
temp=V4a_r(:,1:half)-reshape(V4c2,[],Njs);
fprintf('Cross test 4-2A (z, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c2,[],Njs);
fprintf('Cross test 4-2A (z, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a2 V4c2 Policy4a2 Policy4c2
% --- GI tier ---
vfoptions_t4a3=struct();
vfoptions_t4a3.exoticpreferences='AmbiguityAversion';
vfoptions_t4a3.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a3.ambiguity_pi_z_J=ambiguity_pi_z_J;
vfoptions_t4a3.gridinterplayer=1;
vfoptions_t4a3.ngridinterp=5;
[V4a3,Policy4a3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a3);
V4a_r=reshape(V4a3,[],N_j);
P4a_r=reshape(Policy4a3,[],N_j);
Vstd_r=reshape(Vstd3,[],N_j);
Pstd_r=reshape(Policystd3,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c3=struct();
vfoptions_t4c3.exoticpreferences='AmbiguityAversion';
vfoptions_t4c3.n_ambiguity=3;
vfoptions_t4c3.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c3.gridinterplayer=1;
vfoptions_t4c3.ngridinterp=5;
vfoptions_t4c3.V_Jplus1=Vstd3(:,:,:,half+1);
[V4c3,Policy4c3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c3);
temp=V4a_r(:,1:half)-reshape(V4c3,[],Njs);
fprintf('Cross test 4-2A (z, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c3,[],Njs);
fprintf('Cross test 4-2A (z, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a3 V4c3 Policy4a3 Policy4c3
% --- DC+GI tier ---
vfoptions_t4a4=struct();
vfoptions_t4a4.exoticpreferences='AmbiguityAversion';
vfoptions_t4a4.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a4.ambiguity_pi_z_J=ambiguity_pi_z_J;
vfoptions_t4a4.divideandconquer=1;
vfoptions_t4a4.gridinterplayer=1;
vfoptions_t4a4.ngridinterp=5;
[V4a4,Policy4a4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a4);
V4a_r=reshape(V4a4,[],N_j);
P4a_r=reshape(Policy4a4,[],N_j);
Vstd_r=reshape(Vstd4,[],N_j);
Pstd_r=reshape(Policystd4,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (z, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c4=struct();
vfoptions_t4c4.exoticpreferences='AmbiguityAversion';
vfoptions_t4c4.n_ambiguity=3;
vfoptions_t4c4.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c4.divideandconquer=1;
vfoptions_t4c4.gridinterplayer=1;
vfoptions_t4c4.ngridinterp=5;
vfoptions_t4c4.V_Jplus1=Vstd4(:,:,:,half+1);
[V4c4,Policy4c4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c4);
temp=V4a_r(:,1:half)-reshape(V4c4,[],Njs);
fprintf('Cross test 4-2A (z, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c4,[],Njs);
fprintf('Cross test 4-2A (z, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a4 V4c4 Policy4a4 Policy4c4
clear Vstd1 Vstd2 Vstd3 Vstd4 Policystd1 Policystd2 Policystd3 Policystd4

%% ==================== Ambiguity on an iid e ====================

% Standard preferences (vNM) at each of the four tiers; used by cross tests 1-2A and 4-2A
vfoptions_none1e=struct();
vfoptions_none1e.n_e=n_e;
vfoptions_none1e.e_grid=e_grid;
vfoptions_none1e.pi_e=pi_e;
[Vstd1e,Policystd1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none1e);
vfoptions_none2e=struct();
vfoptions_none2e.n_e=n_e;
vfoptions_none2e.e_grid=e_grid;
vfoptions_none2e.pi_e=pi_e;
vfoptions_none2e.divideandconquer=1;
[Vstd2e,Policystd2e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none2e);
vfoptions_none3e=struct();
vfoptions_none3e.n_e=n_e;
vfoptions_none3e.e_grid=e_grid;
vfoptions_none3e.pi_e=pi_e;
vfoptions_none3e.gridinterplayer=1;
vfoptions_none3e.ngridinterp=5;
[Vstd3e,Policystd3e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none3e);
vfoptions_none4e=struct();
vfoptions_none4e.n_e=n_e;
vfoptions_none4e.e_grid=e_grid;
vfoptions_none4e.pi_e=pi_e;
vfoptions_none4e.divideandconquer=1;
vfoptions_none4e.gridinterplayer=1;
vfoptions_none4e.ngridinterp=5;
[Vstd4e,Policystd4e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none4e);

%% Cross test 1-2A (e): three identical priors = standard preferences, at all four tiers
vfoptions_t11e=struct();
vfoptions_t11e.exoticpreferences='AmbiguityAversion';
vfoptions_t11e.n_ambiguity=3;
vfoptions_t11e.ambiguity_pi_e=[pi_e,pi_e,pi_e];
vfoptions_t11e.n_e=n_e;
vfoptions_t11e.e_grid=e_grid;
vfoptions_t11e.pi_e=pi_e;
[V11e,Policy11e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t11e);
fprintf('Cross test 1-2A (e, plain): three identical priors, this should be zero: %.3e \n',max(abs(V11e(:)-Vstd1e(:))))
fprintf('Cross test 1-2A (e, plain): three identical priors, this should be zero: %.3e \n',max(abs(Policy11e(:)-Policystd1e(:))))
vfoptions_t12e=struct();
vfoptions_t12e.exoticpreferences='AmbiguityAversion';
vfoptions_t12e.n_ambiguity=3;
vfoptions_t12e.ambiguity_pi_e=[pi_e,pi_e,pi_e];
vfoptions_t12e.n_e=n_e;
vfoptions_t12e.e_grid=e_grid;
vfoptions_t12e.pi_e=pi_e;
vfoptions_t12e.divideandconquer=1;
[V12e,Policy12e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t12e);
fprintf('Cross test 1-2A (e, DC): three identical priors, this should be zero: %.3e \n',max(abs(V12e(:)-Vstd2e(:))))
fprintf('Cross test 1-2A (e, DC): three identical priors, this should be zero: %.3e \n',max(abs(Policy12e(:)-Policystd2e(:))))
vfoptions_t13e=struct();
vfoptions_t13e.exoticpreferences='AmbiguityAversion';
vfoptions_t13e.n_ambiguity=3;
vfoptions_t13e.ambiguity_pi_e=[pi_e,pi_e,pi_e];
vfoptions_t13e.n_e=n_e;
vfoptions_t13e.e_grid=e_grid;
vfoptions_t13e.pi_e=pi_e;
vfoptions_t13e.gridinterplayer=1;
vfoptions_t13e.ngridinterp=5;
[V13e,Policy13e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t13e);
fprintf('Cross test 1-2A (e, GI): three identical priors, this should be zero: %.3e \n',max(abs(V13e(:)-Vstd3e(:))))
fprintf('Cross test 1-2A (e, GI): three identical priors, this should be zero: %.3e \n',max(abs(Policy13e(:)-Policystd3e(:))))
vfoptions_t14e=struct();
vfoptions_t14e.exoticpreferences='AmbiguityAversion';
vfoptions_t14e.n_ambiguity=3;
vfoptions_t14e.ambiguity_pi_e=[pi_e,pi_e,pi_e];
vfoptions_t14e.n_e=n_e;
vfoptions_t14e.e_grid=e_grid;
vfoptions_t14e.pi_e=pi_e;
vfoptions_t14e.divideandconquer=1;
vfoptions_t14e.gridinterplayer=1;
vfoptions_t14e.ngridinterp=5;
[V14e,Policy14e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t14e);
fprintf('Cross test 1-2A (e, DC+GI): three identical priors, this should be zero: %.3e \n',max(abs(V14e(:)-Vstd4e(:))))
fprintf('Cross test 1-2A (e, DC+GI): three identical priors, this should be zero: %.3e \n',max(abs(Policy14e(:)-Policystd4e(:))))
clear V11e V12e V13e V14e Policy11e Policy12e Policy13e Policy14e

%% Cross tests 3a/3b-2A (e): an unambiguously worse pi binds, in either prior slot (plain tier)
shifte=0.5*pi_e(n_e);
pi_e_worse=pi_e;
pi_e_worse(n_e)=pi_e(n_e)-shifte;
pi_e_worse(1)=pi_e(1)+shifte;
% Standard preferences under the worse pi_e
vfoptions_noneew=vfoptions_none1e;
vfoptions_noneew.pi_e=pi_e_worse;
[Vstdwe,Policystdwe]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_noneew);
% 3a: worse pi in prior slot 1
vfoptions_t3e=struct();
vfoptions_t3e.exoticpreferences='AmbiguityAversion';
vfoptions_t3e.n_ambiguity=2;
vfoptions_t3e.ambiguity_pi_e=[pi_e_worse,pi_e];
vfoptions_t3e.n_e=n_e;
vfoptions_t3e.e_grid=e_grid;
vfoptions_t3e.pi_e=pi_e;
[V3ae,Policy3ae]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t3e);
fprintf('Cross test 3a-2A (e): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(V3ae(:)-Vstdwe(:))))
fprintf('Cross test 3a-2A (e): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(Policy3ae(:)-Policystdwe(:))))
% 3b: worse pi in prior slot 2 (ordering of the priors is irrelevant)
vfoptions_t3e.ambiguity_pi_e=[pi_e,pi_e_worse];
[V3be,Policy3be]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t3e);
fprintf('Cross test 3b-2A (e): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(V3be(:)-Vstdwe(:))))
fprintf('Cross test 3b-2A (e): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(Policy3be(:)-Policystdwe(:))))
clear V3ae V3be Vstdwe Policy3ae Policy3be Policystdwe

%% Cross test 4-2A (e): age-varying n_ambiguity via V_Jplus1, at all four tiers
% Three priors in the first half, a single (baseline) prior in the second half; the second half is
% then just vNM, and the first half must equal a short solve seeded with the vNM V at half+1.
% Column jj of pi_e_J is the distribution of the e realized in period jj; age-jj expectations read
% column jj+1, so the three priors sit in columns 2..half+1
ambiguity_pi_e_J=repmat(pi_e,[1,N_j,3]);
ambiguity_pi_e_J(:,2:half+1,2)=repmat(ambiguity_pi_e(:,2),[1,half]);
ambiguity_pi_e_J(:,2:half+1,3)=repmat(ambiguity_pi_e(:,3),[1,half]);
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
% (mewj is age-dependent but is only used for the agent distribution, which is not computed here)
% --- plain tier ---
vfoptions_t4a1e=struct();
vfoptions_t4a1e.exoticpreferences='AmbiguityAversion';
vfoptions_t4a1e.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a1e.ambiguity_pi_e_J=ambiguity_pi_e_J;
vfoptions_t4a1e.n_e=n_e;
vfoptions_t4a1e.e_grid=e_grid;
vfoptions_t4a1e.pi_e=pi_e;
[V4a1e,Policy4a1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a1e);
V4a_r=reshape(V4a1e,[],N_j);
P4a_r=reshape(Policy4a1e,[],N_j);
Vstd_r=reshape(Vstd1e,[],N_j);
Pstd_r=reshape(Policystd1e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, plain): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, plain): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c1e=struct();
vfoptions_t4c1e.exoticpreferences='AmbiguityAversion';
vfoptions_t4c1e.n_ambiguity=3;
vfoptions_t4c1e.ambiguity_pi_e=ambiguity_pi_e;
vfoptions_t4c1e.n_e=n_e;
vfoptions_t4c1e.e_grid=e_grid;
vfoptions_t4c1e.pi_e=pi_e;
vfoptions_t4c1e.V_Jplus1=Vstd1e(:,:,:,half+1);
[V4c1e,Policy4c1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c1e);
temp=V4a_r(:,1:half)-reshape(V4c1e,[],Njs);
fprintf('Cross test 4-2A (e, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c1e,[],Njs);
fprintf('Cross test 4-2A (e, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a1e V4c1e Policy4a1e Policy4c1e
% --- DC tier ---
vfoptions_t4a2e=struct();
vfoptions_t4a2e.exoticpreferences='AmbiguityAversion';
vfoptions_t4a2e.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a2e.ambiguity_pi_e_J=ambiguity_pi_e_J;
vfoptions_t4a2e.n_e=n_e;
vfoptions_t4a2e.e_grid=e_grid;
vfoptions_t4a2e.pi_e=pi_e;
vfoptions_t4a2e.divideandconquer=1;
[V4a2e,Policy4a2e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a2e);
V4a_r=reshape(V4a2e,[],N_j);
P4a_r=reshape(Policy4a2e,[],N_j);
Vstd_r=reshape(Vstd2e,[],N_j);
Pstd_r=reshape(Policystd2e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c2e=struct();
vfoptions_t4c2e.exoticpreferences='AmbiguityAversion';
vfoptions_t4c2e.n_ambiguity=3;
vfoptions_t4c2e.ambiguity_pi_e=ambiguity_pi_e;
vfoptions_t4c2e.n_e=n_e;
vfoptions_t4c2e.e_grid=e_grid;
vfoptions_t4c2e.pi_e=pi_e;
vfoptions_t4c2e.divideandconquer=1;
vfoptions_t4c2e.V_Jplus1=Vstd2e(:,:,:,half+1);
[V4c2e,Policy4c2e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c2e);
temp=V4a_r(:,1:half)-reshape(V4c2e,[],Njs);
fprintf('Cross test 4-2A (e, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c2e,[],Njs);
fprintf('Cross test 4-2A (e, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a2e V4c2e Policy4a2e Policy4c2e
% --- GI tier ---
vfoptions_t4a3e=struct();
vfoptions_t4a3e.exoticpreferences='AmbiguityAversion';
vfoptions_t4a3e.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a3e.ambiguity_pi_e_J=ambiguity_pi_e_J;
vfoptions_t4a3e.n_e=n_e;
vfoptions_t4a3e.e_grid=e_grid;
vfoptions_t4a3e.pi_e=pi_e;
vfoptions_t4a3e.gridinterplayer=1;
vfoptions_t4a3e.ngridinterp=5;
[V4a3e,Policy4a3e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a3e);
V4a_r=reshape(V4a3e,[],N_j);
P4a_r=reshape(Policy4a3e,[],N_j);
Vstd_r=reshape(Vstd3e,[],N_j);
Pstd_r=reshape(Policystd3e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c3e=struct();
vfoptions_t4c3e.exoticpreferences='AmbiguityAversion';
vfoptions_t4c3e.n_ambiguity=3;
vfoptions_t4c3e.ambiguity_pi_e=ambiguity_pi_e;
vfoptions_t4c3e.n_e=n_e;
vfoptions_t4c3e.e_grid=e_grid;
vfoptions_t4c3e.pi_e=pi_e;
vfoptions_t4c3e.gridinterplayer=1;
vfoptions_t4c3e.ngridinterp=5;
vfoptions_t4c3e.V_Jplus1=Vstd3e(:,:,:,half+1);
[V4c3e,Policy4c3e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c3e);
temp=V4a_r(:,1:half)-reshape(V4c3e,[],Njs);
fprintf('Cross test 4-2A (e, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c3e,[],Njs);
fprintf('Cross test 4-2A (e, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a3e V4c3e Policy4a3e Policy4c3e
% --- DC+GI tier ---
vfoptions_t4a4e=struct();
vfoptions_t4a4e.exoticpreferences='AmbiguityAversion';
vfoptions_t4a4e.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
vfoptions_t4a4e.ambiguity_pi_e_J=ambiguity_pi_e_J;
vfoptions_t4a4e.n_e=n_e;
vfoptions_t4a4e.e_grid=e_grid;
vfoptions_t4a4e.pi_e=pi_e;
vfoptions_t4a4e.divideandconquer=1;
vfoptions_t4a4e.gridinterplayer=1;
vfoptions_t4a4e.ngridinterp=5;
[V4a4e,Policy4a4e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a4e);
V4a_r=reshape(V4a4e,[],N_j);
P4a_r=reshape(Policy4a4e,[],N_j);
Vstd_r=reshape(Vstd4e,[],N_j);
Pstd_r=reshape(Policystd4e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4-2A (e, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c4e=struct();
vfoptions_t4c4e.exoticpreferences='AmbiguityAversion';
vfoptions_t4c4e.n_ambiguity=3;
vfoptions_t4c4e.ambiguity_pi_e=ambiguity_pi_e;
vfoptions_t4c4e.n_e=n_e;
vfoptions_t4c4e.e_grid=e_grid;
vfoptions_t4c4e.pi_e=pi_e;
vfoptions_t4c4e.divideandconquer=1;
vfoptions_t4c4e.gridinterplayer=1;
vfoptions_t4c4e.ngridinterp=5;
vfoptions_t4c4e.V_Jplus1=Vstd4e(:,:,:,half+1);
[V4c4e,Policy4c4e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c4e);
temp=V4a_r(:,1:half)-reshape(V4c4e,[],Njs);
fprintf('Cross test 4-2A (e, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c4e,[],Njs);
fprintf('Cross test 4-2A (e, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a4e V4c4e Policy4a4e Policy4c4e
clear Vstd1e Vstd2e Vstd3e Vstd4e Policystd1e Policystd2e Policystd3e Policystd4e

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
