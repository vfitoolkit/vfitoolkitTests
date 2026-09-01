function output=AmbFHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% The five ambiguity-aversion cross tests, each run twice: once with the ambiguity on a markov z,
% once with it on an iid e. Every comparison is against exoticpreferences='None' and is an exact zero.
% Test 1: three identical priors = standard preferences (run in both flat and _J input forms)
% Test 2: duplicated priors: 3 pi vs 9 pi (five copies of prior 1, two copies each of priors 2 and 3)
% Test 3a/3b: an unambiguously worse pi binds, in either prior slot (prior ordering is irrelevant)
% Test 4: age-varying n_ambiguity via V_Jplus1 (a single prior in the second half is just vNM there)
% Plus: the pi-consistency warning (the regular pi should be one of the priors) stays silent/fires as appropriate.

n_d=0;
d_grid=[];

ReturnFn_z=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_e=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);

% The three baseline priors (prior 1 is the regular pi_z/pi_e)
ambiguity_pi_z=vfoptionsbaseline.ambiguity_pi_z; % [n_z,n_z,3]
ambiguity_pi_e=vfoptionsbaseline.ambiguity_pi_e; % [n_e,3]
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

half=N_j/2; % for cross test 4 (N_j=20)

%% ==================== Ambiguity on a markov z ====================

% Standard preferences (vNM) under the baseline pi_z; used by cross tests 1 and 4
vfoptions_none=struct();
[Vstd,Policystd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none);

%% Cross test 1 (z): three identical priors = standard preferences
vfoptions_t1=struct();
vfoptions_t1.exoticpreferences='AmbiguityAversion';
vfoptions_t1.n_ambiguity=3;
vfoptions_t1.ambiguity_pi_z=cat(3,pi_z,pi_z,pi_z);
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t1);
fprintf('Cross test 1 (z): three identical priors, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z): three identical priors, this should be zero: %.3e \n',max(abs(Policy1(:)-Policystd(:))))
% Twin run: the same priors input in the pre-built _J form
vfoptions_t1J=struct();
vfoptions_t1J.exoticpreferences='AmbiguityAversion';
vfoptions_t1J.n_ambiguity=3;
vfoptions_t1J.ambiguity_pi_z_J=repmat(reshape(pi_z,[n_z,n_z,1,1]),[1,1,N_j,3]);
[V1J,Policy1J]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t1J);
fprintf('Cross test 1 (z): flat vs _J input form, this should be zero: %.3e \n',max(abs(V1J(:)-V1(:))))
fprintf('Cross test 1 (z): flat vs _J input form, this should be zero: %.3e \n',max(abs(Policy1J(:)-Policy1(:))))

%% Cross test 2 (z): duplicated priors, 3 pi vs 9 pi
vfoptions_t2=struct();
vfoptions_t2.exoticpreferences='AmbiguityAversion';
vfoptions_t2.n_ambiguity=3;
vfoptions_t2.ambiguity_pi_z=ambiguity_pi_z; % the three distinct baseline priors
[V2a,Policy2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t2);
vfoptions_t2b=struct();
vfoptions_t2b.exoticpreferences='AmbiguityAversion';
vfoptions_t2b.n_ambiguity=9;
vfoptions_t2b.ambiguity_pi_z=ambiguity_pi_z(:,:,[1,2,3,1,2,3,1,1,1]); % five copies of prior 1, two copies each of priors 2 and 3
[V2b,Policy2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t2b);
fprintf('Cross test 2 (z): 3 pi vs 9 pi, this should be zero: %.3e \n',max(abs(V2a(:)-V2b(:))))
fprintf('Cross test 2 (z): 3 pi vs 9 pi, this should be zero: %.3e \n',max(abs(Policy2a(:)-Policy2b(:))))

%% Cross tests 3a/3b (z): an unambiguously worse pi binds, in either prior slot
% The worse pi shifts half of each row's best-column probability into the worst column; since V is
% increasing in z this pi has (weakly) lower EV everywhere, so the min always selects it.
shiftz=0.5*pi_z(:,n_z);
pi_z_worse=pi_z;
pi_z_worse(:,n_z)=pi_z(:,n_z)-shiftz;
pi_z_worse(:,1)=pi_z(:,1)+shiftz;
% Standard preferences under the worse pi
[Vstdw,Policystdw]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_worse,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none);
% 3a: worse pi in prior slot 1
vfoptions_t3=struct();
vfoptions_t3.exoticpreferences='AmbiguityAversion';
vfoptions_t3.n_ambiguity=2;
vfoptions_t3.ambiguity_pi_z=cat(3,pi_z_worse,pi_z);
[V3a,Policy3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t3);
fprintf('Cross test 3a (z): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(V3a(:)-Vstdw(:))))
fprintf('Cross test 3a (z): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(Policy3a(:)-Policystdw(:))))
% 3b: worse pi in prior slot 2 (ordering of the priors is irrelevant)
vfoptions_t3.ambiguity_pi_z=cat(3,pi_z,pi_z_worse);
[V3b,Policy3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t3);
fprintf('Cross test 3b (z): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(V3b(:)-Vstdw(:))))
fprintf('Cross test 3b (z): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(Policy3b(:)-Policystdw(:))))

%% Cross test 4 (z): age-varying n_ambiguity via V_Jplus1
% (a) full horizon: three priors in the first half, a single (baseline) prior in the second half
vfoptions_t4a=struct();
vfoptions_t4a.exoticpreferences='AmbiguityAversion';
vfoptions_t4a.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
% Slice jj of pi_z_J is the transition from period jj to jj+1, so the three priors sit in slices 1..half
ambiguity_pi_z_J=repmat(reshape(pi_z,[n_z,n_z,1,1]),[1,1,N_j,3]);
ambiguity_pi_z_J(:,:,1:half,2)=repmat(ambiguity_pi_z(:,:,2),[1,1,half]);
ambiguity_pi_z_J(:,:,1:half,3)=repmat(ambiguity_pi_z(:,:,3),[1,1,half]);
vfoptions_t4a.ambiguity_pi_z_J=ambiguity_pi_z_J;
[V4a,Policy4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a);
% The second half of (a) is a single-prior model, so it must equal the standard solve there
V4a_r=reshape(V4a,[],N_j);
P4a_r=reshape(Policy4a,[],N_j);
Vstd_r=reshape(Vstd,[],N_j);
Pstd_r=reshape(Policystd,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
% (c) the first half of (a) again, as a short model with V_Jplus1 from the standard solve
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
% (mewj is age-dependent but is only used for the agent distribution, which is not computed here)
vfoptions_t4c=struct();
vfoptions_t4c.exoticpreferences='AmbiguityAversion';
vfoptions_t4c.n_ambiguity=3;
vfoptions_t4c.ambiguity_pi_z=ambiguity_pi_z;
vfoptions_t4c.V_Jplus1=Vstd(:,:,half+1);
[V4c,Policy4c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c);
temp=V4a_r(:,1:half)-reshape(V4c,[],Njs);
fprintf('Cross test 4 (z): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c,[],Njs);
fprintf('Cross test 4 (z): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))

% --- Companion coverage: repeat cross test 4 (z) at the DC, GI and DC+GI tiers, giving the
% DC1/GI1/DC1_GI1 raws' V_Jplus1 branches runtime coverage (they are otherwise only reached in-loop) ---
% DC tier
vfoptions_none2=struct();
vfoptions_none2.divideandconquer=1;
[Vstd2,Policystd2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none2);
vfoptions_t4a2=vfoptions_t4a;
vfoptions_t4a2.divideandconquer=1;
[V4a2,Policy4a2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a2);
V4a_r=reshape(V4a2,[],N_j);
P4a_r=reshape(Policy4a2,[],N_j);
Vstd_r=reshape(Vstd2,[],N_j);
Pstd_r=reshape(Policystd2,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c2=vfoptions_t4c;
vfoptions_t4c2.divideandconquer=1;
vfoptions_t4c2.V_Jplus1=Vstd2(:,:,half+1);
[V4c2,Policy4c2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c2);
temp=V4a_r(:,1:half)-reshape(V4c2,[],Njs);
fprintf('Cross test 4 (z, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c2,[],Njs);
fprintf('Cross test 4 (z, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a2 V4c2 Policy4a2 Policy4c2 Vstd2 Policystd2
% GI tier
vfoptions_none3=struct();
vfoptions_none3.gridinterplayer=1;
vfoptions_none3.ngridinterp=5;
[Vstd3,Policystd3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none3);
vfoptions_t4a3=vfoptions_t4a;
vfoptions_t4a3.gridinterplayer=1;
vfoptions_t4a3.ngridinterp=5;
[V4a3,Policy4a3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a3);
V4a_r=reshape(V4a3,[],N_j);
P4a_r=reshape(Policy4a3,[],N_j);
Vstd_r=reshape(Vstd3,[],N_j);
Pstd_r=reshape(Policystd3,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c3=vfoptions_t4c;
vfoptions_t4c3.gridinterplayer=1;
vfoptions_t4c3.ngridinterp=5;
vfoptions_t4c3.V_Jplus1=Vstd3(:,:,half+1);
[V4c3,Policy4c3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c3);
temp=V4a_r(:,1:half)-reshape(V4c3,[],Njs);
fprintf('Cross test 4 (z, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c3,[],Njs);
fprintf('Cross test 4 (z, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a3 V4c3 Policy4a3 Policy4c3 Vstd3 Policystd3
% DC+GI tier
vfoptions_none4=struct();
vfoptions_none4.divideandconquer=1;
vfoptions_none4.gridinterplayer=1;
vfoptions_none4.ngridinterp=5;
[Vstd4,Policystd4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_none4);
vfoptions_t4a4=vfoptions_t4a;
vfoptions_t4a4.divideandconquer=1;
vfoptions_t4a4.gridinterplayer=1;
vfoptions_t4a4.ngridinterp=5;
[V4a4,Policy4a4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t4a4);
V4a_r=reshape(V4a4,[],N_j);
P4a_r=reshape(Policy4a4,[],N_j);
Vstd_r=reshape(Vstd4,[],N_j);
Pstd_r=reshape(Policystd4,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (z, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c4=vfoptions_t4c;
vfoptions_t4c4.divideandconquer=1;
vfoptions_t4c4.gridinterplayer=1;
vfoptions_t4c4.ngridinterp=5;
vfoptions_t4c4.V_Jplus1=Vstd4(:,:,half+1);
[V4c4,Policy4c4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c4);
temp=V4a_r(:,1:half)-reshape(V4c4,[],Njs);
fprintf('Cross test 4 (z, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c4,[],Njs);
fprintf('Cross test 4 (z, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a4 V4c4 Policy4a4 Policy4c4 Vstd4 Policystd4


%% pi-consistency warning (z): the regular pi_z should be one of the priors
lastwarn(''); % reset
[Vw,Policyw]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t2);
warnmsg=lastwarn();
if isempty(warnmsg)
    fprintf('pi-consistency (z): no warning when pi_z is one of the priors :) \n')
else
    fprintf('pi-consistency (z): FAIL, got an unexpected warning: %s \n',warnmsg)
end
lastwarn(''); % reset
pi_z_off=0.5*pi_z+0.5*ones(n_z,n_z)/n_z; % not equal to any of the priors
[Vw,Policyw]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_off,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_t2);
warnmsg=lastwarn();
if isempty(warnmsg)
    fprintf('pi-consistency (z): FAIL, expected a warning (pi_z is not one of the priors) and got none \n')
else
    fprintf('pi-consistency (z): warning fires when pi_z is not one of the priors :) \n')
end

clear V1 V1J V2a V2b V3a V3b V4a V4c Vw Vstd Vstdw Policy1 Policy1J Policy2a Policy2b Policy3a Policy3b Policy4a Policy4c Policyw Policystd Policystdw V4a_r P4a_r Vstd_r Pstd_r

%% ==================== Ambiguity on an iid e ====================
% Same five tests with the ambiguity on the iid e (a separate code path: ambiguity_pi_e, iid expectations)

% Standard preferences (vNM) under the baseline pi_e; used by cross tests 1 and 4
vfoptions_nonee=struct();
vfoptions_nonee.n_e=n_e;
vfoptions_nonee.e_grid=e_grid;
vfoptions_nonee.pi_e=pi_e;
[Vstde,Policystde]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_nonee);

%% Cross test 1 (e): three identical priors = standard preferences
vfoptions_t1e=vfoptions_nonee;
vfoptions_t1e.exoticpreferences='AmbiguityAversion';
vfoptions_t1e.n_ambiguity=3;
vfoptions_t1e.ambiguity_pi_e=[pi_e,pi_e,pi_e];
[V1e,Policy1e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t1e);
fprintf('Cross test 1 (e): three identical priors, this should be zero: %.3e \n',max(abs(V1e(:)-Vstde(:))))
fprintf('Cross test 1 (e): three identical priors, this should be zero: %.3e \n',max(abs(Policy1e(:)-Policystde(:))))
% Twin run: the same priors input in the pre-built _J form
vfoptions_t1Je=vfoptions_nonee;
vfoptions_t1Je.exoticpreferences='AmbiguityAversion';
vfoptions_t1Je.n_ambiguity=3;
vfoptions_t1Je.ambiguity_pi_e_J=repmat(pi_e,[1,N_j,3]);
[V1Je,Policy1Je]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t1Je);
fprintf('Cross test 1 (e): flat vs _J input form, this should be zero: %.3e \n',max(abs(V1Je(:)-V1e(:))))
fprintf('Cross test 1 (e): flat vs _J input form, this should be zero: %.3e \n',max(abs(Policy1Je(:)-Policy1e(:))))

%% Cross test 2 (e): duplicated priors, 3 pi vs 9 pi
vfoptions_t2e=vfoptions_nonee;
vfoptions_t2e.exoticpreferences='AmbiguityAversion';
vfoptions_t2e.n_ambiguity=3;
vfoptions_t2e.ambiguity_pi_e=ambiguity_pi_e; % the three distinct baseline priors
[V2ae,Policy2ae]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t2e);
vfoptions_t2be=vfoptions_nonee;
vfoptions_t2be.exoticpreferences='AmbiguityAversion';
vfoptions_t2be.n_ambiguity=9;
vfoptions_t2be.ambiguity_pi_e=ambiguity_pi_e(:,[1,2,3,1,2,3,1,1,1]); % five copies of prior 1, two copies each of priors 2 and 3
[V2be,Policy2be]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t2be);
fprintf('Cross test 2 (e): 3 pi vs 9 pi, this should be zero: %.3e \n',max(abs(V2ae(:)-V2be(:))))
fprintf('Cross test 2 (e): 3 pi vs 9 pi, this should be zero: %.3e \n',max(abs(Policy2ae(:)-Policy2be(:))))

%% Cross tests 3a/3b (e): an unambiguously worse pi binds, in either prior slot
shifte=0.5*pi_e(n_e);
pi_e_worse=pi_e;
pi_e_worse(n_e)=pi_e(n_e)-shifte;
pi_e_worse(1)=pi_e(1)+shifte;
% Standard preferences under the worse pi_e
vfoptions_noneew=vfoptions_nonee;
vfoptions_noneew.pi_e=pi_e_worse;
[Vstdwe,Policystdwe]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_noneew);
% 3a: worse pi in prior slot 1
vfoptions_t3e=vfoptions_nonee;
vfoptions_t3e.exoticpreferences='AmbiguityAversion';
vfoptions_t3e.n_ambiguity=2;
vfoptions_t3e.ambiguity_pi_e=[pi_e_worse,pi_e];
[V3ae,Policy3ae]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t3e);
fprintf('Cross test 3a (e): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(V3ae(:)-Vstdwe(:))))
fprintf('Cross test 3a (e): worse pi (slot 1) binds, this should be zero: %.3e \n',max(abs(Policy3ae(:)-Policystdwe(:))))
% 3b: worse pi in prior slot 2 (ordering of the priors is irrelevant)
vfoptions_t3e.ambiguity_pi_e=[pi_e,pi_e_worse];
[V3be,Policy3be]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t3e);
fprintf('Cross test 3b (e): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(V3be(:)-Vstdwe(:))))
fprintf('Cross test 3b (e): worse pi (slot 2) binds, this should be zero: %.3e \n',max(abs(Policy3be(:)-Policystdwe(:))))

%% Cross test 4 (e): age-varying n_ambiguity via V_Jplus1
% (a) full horizon: three priors in the first half, a single (baseline) prior in the second half
vfoptions_t4ae=vfoptions_nonee;
vfoptions_t4ae.exoticpreferences='AmbiguityAversion';
vfoptions_t4ae.n_ambiguity=[3*ones(1,half),ones(1,N_j-half)];
% Column jj of pi_e_J is the distribution of the e realized in period jj; age-jj expectations read
% column jj+1, so the three priors sit in columns 2..half+1
ambiguity_pi_e_J=repmat(pi_e,[1,N_j,3]);
ambiguity_pi_e_J(:,2:half+1,2)=repmat(ambiguity_pi_e(:,2),[1,half]);
ambiguity_pi_e_J(:,2:half+1,3)=repmat(ambiguity_pi_e(:,3),[1,half]);
vfoptions_t4ae.ambiguity_pi_e_J=ambiguity_pi_e_J;
[V4ae,Policy4ae]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4ae);
% The second half of (a) is a single-prior model, so it must equal the standard solve there
V4ae_r=reshape(V4ae,[],N_j);
P4ae_r=reshape(Policy4ae,[],N_j);
Vstde_r=reshape(Vstde,[],N_j);
Pstde_r=reshape(Policystde,[],N_j);
temp=V4ae_r(:,half+1:N_j)-Vstde_r(:,half+1:N_j);
fprintf('Cross test 4 (e): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4ae_r(:,half+1:N_j)-Pstde_r(:,half+1:N_j);
fprintf('Cross test 4 (e): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
% (c) the first half of (a) again, as a short model with V_Jplus1 from the standard solve
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
% (mewj is age-dependent but is only used for the agent distribution, which is not computed here)
vfoptions_t4ce=vfoptions_nonee;
vfoptions_t4ce.exoticpreferences='AmbiguityAversion';
vfoptions_t4ce.n_ambiguity=3;
vfoptions_t4ce.ambiguity_pi_e=ambiguity_pi_e;
vfoptions_t4ce.V_Jplus1=Vstde(:,:,half+1);
[V4ce,Policy4ce]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4ce);
temp=V4ae_r(:,1:half)-reshape(V4ce,[],Njs);
fprintf('Cross test 4 (e): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4ae_r(:,1:half)-reshape(Policy4ce,[],Njs);
fprintf('Cross test 4 (e): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))

% --- Companion coverage: repeat cross test 4 (e) at the DC, GI and DC+GI tiers, giving the
% DC1/GI1/DC1_GI1 raws' V_Jplus1 branches runtime coverage (they are otherwise only reached in-loop) ---
% DC tier
vfoptions_none2e=vfoptions_nonee;
vfoptions_none2e.divideandconquer=1;
[Vstd2e,Policystd2e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none2e);
vfoptions_t4a2e=vfoptions_t4ae;
vfoptions_t4a2e.divideandconquer=1;
[V4a2e,Policy4a2e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a2e);
V4a_r=reshape(V4a2e,[],N_j);
P4a_r=reshape(Policy4a2e,[],N_j);
Vstd_r=reshape(Vstd2e,[],N_j);
Pstd_r=reshape(Policystd2e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (e, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (e, DC): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c2e=vfoptions_t4ce;
vfoptions_t4c2e.divideandconquer=1;
vfoptions_t4c2e.V_Jplus1=Vstd2e(:,:,half+1);
[V4c2e,Policy4c2e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c2e);
temp=V4a_r(:,1:half)-reshape(V4c2e,[],Njs);
fprintf('Cross test 4 (e, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c2e,[],Njs);
fprintf('Cross test 4 (e, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a2e V4c2e Policy4a2e Policy4c2e Vstd2e Policystd2e
% GI tier
vfoptions_none3e=vfoptions_nonee;
vfoptions_none3e.gridinterplayer=1;
vfoptions_none3e.ngridinterp=5;
[Vstd3e,Policystd3e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none3e);
vfoptions_t4a3e=vfoptions_t4ae;
vfoptions_t4a3e.gridinterplayer=1;
vfoptions_t4a3e.ngridinterp=5;
[V4a3e,Policy4a3e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a3e);
V4a_r=reshape(V4a3e,[],N_j);
P4a_r=reshape(Policy4a3e,[],N_j);
Vstd_r=reshape(Vstd3e,[],N_j);
Pstd_r=reshape(Policystd3e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (e, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (e, GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c3e=vfoptions_t4ce;
vfoptions_t4c3e.gridinterplayer=1;
vfoptions_t4c3e.ngridinterp=5;
vfoptions_t4c3e.V_Jplus1=Vstd3e(:,:,half+1);
[V4c3e,Policy4c3e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c3e);
temp=V4a_r(:,1:half)-reshape(V4c3e,[],Njs);
fprintf('Cross test 4 (e, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c3e,[],Njs);
fprintf('Cross test 4 (e, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a3e V4c3e Policy4a3e Policy4c3e Vstd3e Policystd3e
% DC+GI tier
vfoptions_none4e=vfoptions_nonee;
vfoptions_none4e.divideandconquer=1;
vfoptions_none4e.gridinterplayer=1;
vfoptions_none4e.ngridinterp=5;
[Vstd4e,Policystd4e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_none4e);
vfoptions_t4a4e=vfoptions_t4ae;
vfoptions_t4a4e.divideandconquer=1;
vfoptions_t4a4e.gridinterplayer=1;
vfoptions_t4a4e.ngridinterp=5;
[V4a4e,Policy4a4e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t4a4e);
V4a_r=reshape(V4a4e,[],N_j);
P4a_r=reshape(Policy4a4e,[],N_j);
Vstd_r=reshape(Vstd4e,[],N_j);
Pstd_r=reshape(Policystd4e,[],N_j);
temp=V4a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 4 (e, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 4 (e, DC+GI): second half (single prior) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfoptions_t4c4e=vfoptions_t4ce;
vfoptions_t4c4e.divideandconquer=1;
vfoptions_t4c4e.gridinterplayer=1;
vfoptions_t4c4e.ngridinterp=5;
vfoptions_t4c4e.V_Jplus1=Vstd4e(:,:,half+1);
[V4c4e,Policy4c4e]=ValueFnIter_Case1_FHorz(n_d,n_a,0,Njs,d_grid,a_grid,[],[],ReturnFn_e,Paramsjs,DiscountFactorParamNames,[],vfoptions_t4c4e);
temp=V4a_r(:,1:half)-reshape(V4c4e,[],Njs);
fprintf('Cross test 4 (e, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P4a_r(:,1:half)-reshape(Policy4c4e,[],Njs);
fprintf('Cross test 4 (e, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V4a4e V4c4e Policy4a4e Policy4c4e Vstd4e Policystd4e


%% pi-consistency warning (e): the regular pi_e should be one of the priors
lastwarn(''); % reset
[Vwe,Policywe]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t2e);
warnmsg=lastwarn();
if isempty(warnmsg)
    fprintf('pi-consistency (e): no warning when pi_e is one of the priors :) \n')
else
    fprintf('pi-consistency (e): FAIL, got an unexpected warning: %s \n',warnmsg)
end
lastwarn(''); % reset
vfoptions_t2eoff=vfoptions_t2e;
vfoptions_t2eoff.pi_e=0.5*pi_e+0.5*ones(n_e,1)/n_e; % not equal to any of the priors
[Vwe,Policywe]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_t2eoff);
warnmsg=lastwarn();
if isempty(warnmsg)
    fprintf('pi-consistency (e): FAIL, expected a warning (pi_e is not one of the priors) and got none \n')
else
    fprintf('pi-consistency (e): warning fires when pi_e is not one of the priors :) \n')
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
