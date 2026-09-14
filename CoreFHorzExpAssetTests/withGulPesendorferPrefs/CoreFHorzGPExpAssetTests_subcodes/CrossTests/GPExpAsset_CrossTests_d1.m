function output=GPExpAsset_CrossTests_d1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% The Gul-Pesendorfer x ExperienceAsset cross tests (with d1): the decisions are d1
% (labour), d2 (which drives the experience asset a2prime=aprimeFn(d2,a2)), and the
% a1prime savings choice. The lambdaGP=0 limiting case switches the temptation machinery
% OFF entirely, so limiting cases alone never exercise it; hence the hand-rolled
% brute-force recursion in test 4, which assembles the ExpAsset a2prime expectation by
% hand from aprimeFn (lower grid index plus probability-of-lower, matching
% CreateExperienceAssetFnMatrix incl. its off-grid clamps and the skipinterp/endpoint
% conventions of the solver raws). The hand-rolled GI recursion lives in the nod1 cross
% file only (as in the GP core bank, where the GI hand-roll is in the nod file).
% Test 1: lambdaGP=0 (no temptation) equals standard ExpAsset. All four shock combos
%         (z_noe, noz_e, z_e, noz_noe) at all four tiers (plain, DC1, GI1, DC1+GI1), so
%         every raw gets a limiting-case anchor.
% Test 2: constant temptation (lambdaGP=0, shiftGP=0.7) equals standard: a constant v
%         cancels exactly against the most-tempting term.
% Test 3: shift invariance: v and v+0.7 give identical V and Policy.
% Test 4: hand-rolled brute-force GP-ExpAsset recursion (small grid, plain MATLAB loops,
%         most-tempting term a max over the FULL joint (d1,d2,a1prime) choice set, rows
%         ordered d1 fastest then d2 then a1prime, matching the toolkit's kron ordering)
%         vs the toolkit at the plain and DC1 tiers. The most tempting choice (the
%         consumption splurge: low a1prime) is far from the optimum, so a DC
%         implementation that wrongly took the temptation max over its restricted window
%         would show up here.
% Test 5: age-dependent lambdaGP (temptation in the first half of life only): second half
%         equals standard, first half equals a V_Jplus1 short solve. At all four tiers
%         (this also gives the new raws' V_Jplus1 branches runtime coverage).
% Test 5b: age-dependent pi_z_J (the KLM2021 configuration): lambdaGP=0 equals standard
%         under an age-varying transition at all four tiers (catches a dropped ,jj page
%         index in the GP raws' EV), plus the test-5 half-life split at the plain tier
%         under the age-varying transition at baseline lambdaGP.
% Test 6: sign check: temptation weakly lowers welfare, V_GP <= V_standard everywhere.

n_d1=n_d(1);
n_d2=n_d(2);
d1_grid=d_grid(1:n_d1);
d2_grid=d_grid(n_d1+1:end);
n_a1=n_a(1);
n_a2=n_a(2);
a1_grid=a_grid(1:n_a1);
a2_grid=a_grid(n_a1+1:end);

ReturnFn_z=@(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_e=@(d1,d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_e_nosemiz(d1,d2,a1prime,a1,a2,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_ze=@(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_z_e_nosemiz(d1,d2,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);
ReturnFn_none=@(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension) ReturnFn_d1_noz_noe_nosemiz(d1,d2,a1prime,a1,a2,r,w,kappa_j,sigma,varphi,eta,agej,Jr,pension);

TemptationFn_z=@(d1,d2,a1prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d1_z_noe_nosemiz(d1,d2,a1prime,a1,a2,z,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_e=@(d1,d2,a1prime,a1,a2,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d1_noz_e_nosemiz(d1,d2,a1prime,a1,a2,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_ze=@(d1,d2,a1prime,a1,a2,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d1_z_e_nosemiz(d1,d2,a1prime,a1,a2,z,e,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);
TemptationFn_none=@(d1,d2,a1prime,a1,a2,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension) GPTemptationFn_d1_noz_noe_nosemiz(d1,d2,a1prime,a1,a2,lambdaGP,shiftGP,r,w,kappa_j,sigma,agej,Jr,pension);

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

half=N_j/2; % for cross test 5 (N_j=20)

% Every solve (standard and GP) is an experienceasset solve
vfEA=struct();
vfEA.experienceasset=1;
vfEA.aprimeFn=vfoptionsbaseline.aprimeFn;

% Params with the temptation switched off (used by tests 1 and 2)
Params0=Params;
Params0.lambdaGP=0;

%% ==================== Cross test 1: lambdaGP=0 equals standard ====================
% --- combo z_noe ---
vfstd1=vfEA;
vfstd2=vfEA; vfstd2.divideandconquer=1;
vfstd3=vfEA; vfstd3.gridinterplayer=1; vfstd3.ngridinterp=5;
vfstd4=vfEA; vfstd4.divideandconquer=1; vfstd4.gridinterplayer=1; vfstd4.ngridinterp=5;
vfgp1=vfstd1; vfgp1.exoticpreferences='GulPesendorfer'; vfgp1.temptationFn=TemptationFn_z;
vfgp2=vfstd2; vfgp2.exoticpreferences='GulPesendorfer'; vfgp2.temptationFn=TemptationFn_z;
vfgp3=vfstd3; vfgp3.exoticpreferences='GulPesendorfer'; vfgp3.temptationFn=TemptationFn_z;
vfgp4=vfstd4; vfgp4.exoticpreferences='GulPesendorfer'; vfgp4.temptationFn=TemptationFn_z;
[Vstd1,Pstd1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd1);
[Vstd2,Pstd2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd2);
[Vstd3,Pstd3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd3);
[Vstd4,Pstd4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd4);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 1 (z_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 1 (z_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 1 (z_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 1 (z_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 1 (z_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 1 (z_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 1 (z_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 1 (z_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V1 V2 V3 V4 P1 P2 P3 P4
% (Vstd1..4/Pstd1..4 for the z_noe combo are kept: tests 2, 5 and 6 reuse them)

% --- combo noz_e ---
vfstd1e=vfEA; vfstd1e.n_e=n_e; vfstd1e.e_grid=e_grid; vfstd1e.pi_e=pi_e;
vfstd2e=vfstd1e; vfstd2e.divideandconquer=1;
vfstd3e=vfstd1e; vfstd3e.gridinterplayer=1; vfstd3e.ngridinterp=5;
vfstd4e=vfstd2e; vfstd4e.gridinterplayer=1; vfstd4e.ngridinterp=5;
vfgp1e=vfstd1e; vfgp1e.exoticpreferences='GulPesendorfer'; vfgp1e.temptationFn=TemptationFn_e;
vfgp2e=vfstd2e; vfgp2e.exoticpreferences='GulPesendorfer'; vfgp2e.temptationFn=TemptationFn_e;
vfgp3e=vfstd3e; vfgp3e.exoticpreferences='GulPesendorfer'; vfgp3e.temptationFn=TemptationFn_e;
vfgp4e=vfstd4e; vfgp4e.exoticpreferences='GulPesendorfer'; vfgp4e.temptationFn=TemptationFn_e;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd1e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp1e);
fprintf('Cross test 1 (noz_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd2e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp2e);
fprintf('Cross test 1 (noz_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd3e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp3e);
fprintf('Cross test 1 (noz_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfstd4e);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params0,DiscountFactorParamNames,[],vfgp4e);
fprintf('Cross test 1 (noz_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

% --- combo z_e ---
vfstd1ze=vfEA; vfstd1ze.n_e=n_e; vfstd1ze.e_grid=e_grid; vfstd1ze.pi_e=pi_e;
vfstd2ze=vfstd1ze; vfstd2ze.divideandconquer=1;
vfstd3ze=vfstd1ze; vfstd3ze.gridinterplayer=1; vfstd3ze.ngridinterp=5;
vfstd4ze=vfstd2ze; vfstd4ze.gridinterplayer=1; vfstd4ze.ngridinterp=5;
vfgp1ze=vfstd1ze; vfgp1ze.exoticpreferences='GulPesendorfer'; vfgp1ze.temptationFn=TemptationFn_ze;
vfgp2ze=vfstd2ze; vfgp2ze.exoticpreferences='GulPesendorfer'; vfgp2ze.temptationFn=TemptationFn_ze;
vfgp3ze=vfstd3ze; vfgp3ze.exoticpreferences='GulPesendorfer'; vfgp3ze.temptationFn=TemptationFn_ze;
vfgp4ze=vfstd4ze; vfgp4ze.exoticpreferences='GulPesendorfer'; vfgp4ze.temptationFn=TemptationFn_ze;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd1ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp1ze);
fprintf('Cross test 1 (z_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd2ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp2ze);
fprintf('Cross test 1 (z_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd3ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp3ze);
fprintf('Cross test 1 (z_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfstd4ze);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params0,DiscountFactorParamNames,[],vfgp4ze);
fprintf('Cross test 1 (z_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (z_e, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

% --- combo noz_noe ---
vfstd1n=vfEA;
vfstd2n=vfEA; vfstd2n.divideandconquer=1;
vfstd3n=vfEA; vfstd3n.gridinterplayer=1; vfstd3n.ngridinterp=5;
vfstd4n=vfEA; vfstd4n.divideandconquer=1; vfstd4n.gridinterplayer=1; vfstd4n.ngridinterp=5;
vfgp1n=vfstd1n; vfgp1n.exoticpreferences='GulPesendorfer'; vfgp1n.temptationFn=TemptationFn_none;
vfgp2n=vfstd2n; vfgp2n.exoticpreferences='GulPesendorfer'; vfgp2n.temptationFn=TemptationFn_none;
vfgp3n=vfstd3n; vfgp3n.exoticpreferences='GulPesendorfer'; vfgp3n.temptationFn=TemptationFn_none;
vfgp4n=vfstd4n; vfgp4n.exoticpreferences='GulPesendorfer'; vfgp4n.temptationFn=TemptationFn_none;
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd1n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp1n);
fprintf('Cross test 1 (noz_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, plain): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd2n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp2n);
fprintf('Cross test 1 (noz_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, DC): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd3n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp3n);
fprintf('Cross test 1 (noz_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
[Vstd,Pstd]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfstd4n);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params0,DiscountFactorParamNames,[],vfgp4n);
fprintf('Cross test 1 (noz_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd(:))))
fprintf('Cross test 1 (noz_noe, DC+GI): lambda=0 equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd(:))))
clear V1 P1 Vstd Pstd

%% ==================== Cross test 2: constant temptation equals standard ====================
% v = 0.7 on the feasible set (lambdaGP=0, shiftGP=0.7): the constant cancels exactly against
% the most-tempting term, so this equals standard preferences DESPITE the temptation machinery
% being fully engaged (nonzero MostTempting). Catches a dropped or double-counted subtraction.
Params2=Params;
Params2.lambdaGP=0;
Params2.shiftGP=0.7;
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 2 (z_noe, plain): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V1(:)-Vstd1(:))))
fprintf('Cross test 2 (z_noe, plain): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P1(:)-Pstd1(:))))
[V2,P2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 2 (z_noe, DC): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V2(:)-Vstd2(:))))
fprintf('Cross test 2 (z_noe, DC): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P2(:)-Pstd2(:))))
[V3,P3]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 2 (z_noe, GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V3(:)-Vstd3(:))))
fprintf('Cross test 2 (z_noe, GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P3(:)-Pstd3(:))))
[V4,P4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params2,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 2 (z_noe, DC+GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(V4(:)-Vstd4(:))))
fprintf('Cross test 2 (z_noe, DC+GI): constant temptation equals standard, this should be zero: %.3e \n',max(abs(P4(:)-Pstd4(:))))
clear V1 V2 V3 V4 P1 P2 P3 P4

%% ==================== Cross test 3: shift invariance ====================
% v and v+0.7 (lambdaGP=0.1, shiftGP 0 vs 0.7) give identical V and Policy
Params3a=Params; % lambdaGP=0.1, shiftGP=0 (the baseline)
Params3b=Params;
Params3b.shiftGP=0.7;
[V1a,P1a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp1);
[V1b,P1b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 3 (z_noe, plain): shift invariance, this should be zero: %.3e \n',max(abs(V1a(:)-V1b(:))))
fprintf('Cross test 3 (z_noe, plain): shift invariance, this should be zero: %.3e \n',max(abs(P1a(:)-P1b(:))))
[V2a,P2a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp2);
[V2b,P2b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 3 (z_noe, DC): shift invariance, this should be zero: %.3e \n',max(abs(V2a(:)-V2b(:))))
fprintf('Cross test 3 (z_noe, DC): shift invariance, this should be zero: %.3e \n',max(abs(P2a(:)-P2b(:))))
[V3a,P3a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp3);
[V3b,P3b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 3 (z_noe, GI): shift invariance, this should be zero: %.3e \n',max(abs(V3a(:)-V3b(:))))
fprintf('Cross test 3 (z_noe, GI): shift invariance, this should be zero: %.3e \n',max(abs(P3a(:)-P3b(:))))
[V4a,P4a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3a,DiscountFactorParamNames,[],vfgp4);
[V4b,P4b]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params3b,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 3 (z_noe, DC+GI): shift invariance, this should be zero: %.3e \n',max(abs(V4a(:)-V4b(:))))
fprintf('Cross test 3 (z_noe, DC+GI): shift invariance, this should be zero: %.3e \n',max(abs(P4a(:)-P4b(:))))
clear V2a V2b V3a V3b V4a V4b P2a P2b P3a P3b P4a P4b
% (V1a/P1a, the plain-tier GP solve at the baseline lambdaGP=0.1, is kept for test 6)

%% ==================== Cross test 4: hand-rolled brute force (plain and DC) ====================
% Small grid, plain MATLAB loops on the CPU, calling the very same ReturnFn/TemptationFn anons.
% V_j(a1,a2,z) = max_{d1,d2,a1'} [ u + v + beta*E V_{j+1} ] - max_{d1,d2,a1'} v, with BOTH
% maxes over the full joint (d1,d2,a1prime) choice set (rows ordered d1 fastest then d2 then
% a1prime, matching the toolkit's kron ordering). The continuation E V_{j+1} is assembled by
% hand: a2prime=aprimeFn(d2,a2) is converted to (lower a2 grid index, probability of lower)
% exactly as CreateExperienceAssetFnMatrix does it (incl. the off-grid clamps), the
% a2prime-mix uses the raws' endpoint conventions (prob 1 -> lower exactly, prob 0 -> upper
% exactly, and skipinterp when the two nodes are equal), and only then is the z' expectation
% taken. Operation order matches the raws: max((U+T)+beta*EV) with the most-tempting term
% subtracted AFTER the max. Toolkit Policy rows (d1, d2, a1prime) are recombined into the
% joint index d1 + n_d1*(d2-1) + n_d1*n_d2*(a1prime-1) for the comparison.
% A reduced a2 grid (5 points) keeps the CPU loops affordable; the toolkit solves it is
% compared against use the same small grids.
d1_grid_cpu=gather(d1_grid);
d2_grid_cpu=gather(d2_grid);
n_a1_bf=21;
a1_grid_bf=5*linspace(0,1,n_a1_bf)'.^3;
n_a2_bf=5;
a2_grid_bf=linspace(0,10,n_a2_bf)';
a_grid_bf=[a1_grid_bf;a2_grid_bf];
n_a_bf=[n_a1_bf,n_a2_bf];
pi_z_cpu=gather(pi_z);
z_grid_cpu=gather(z_grid);
a2_griddiff=a2_grid_bf(2:end)-a2_grid_bf(1:end-1);
% a2prime lower index and prob-of-lower, per (d2,a2): exactly the CreateExperienceAssetFnMatrix rules
a2p_lo=zeros(n_d2,n_a2_bf);
a2p_p=zeros(n_d2,n_a2_bf);
for d2_c=1:n_d2
    for a2_c=1:n_a2_bf
        a2pval=Params.phi1*(1-d2_grid_cpu(d2_c))+(1-Params.phi2)*a2_grid_bf(a2_c);
        lo=find(a2_grid_bf<=a2pval,1,'last');
        lo=max(lo,1);
        p=1-(a2pval-a2_grid_bf(lo))/a2_griddiff(min(lo,n_a2_bf-1));
        if a2pval>=a2_grid_bf(end)
            lo=n_a2_bf-1; p=0;
        elseif a2pval<=a2_grid_bf(1)
            lo=1; p=1;
        end
        a2p_lo(d2_c,a2_c)=lo;
        a2p_p(d2_c,a2_c)=p;
    end
end
U_j=zeros(n_d1,n_d2,n_a1_bf,n_a1_bf,n_a2_bf,n_z);
T_j=zeros(n_d1,n_d2,n_a1_bf,n_a1_bf,n_a2_bf,n_z);
V_bf=zeros(n_a1_bf,n_a2_bf,n_z,N_j);
P_bf=zeros(n_a1_bf,n_a2_bf,n_z,N_j); % joint (d1,d2,a1prime) index, d1 fastest
for jj=N_j:-1:1
    for z_c=1:n_z
        for a2_c=1:n_a2_bf
            for a1_c=1:n_a1_bf
                for ap_c=1:n_a1_bf
                    for d2_c=1:n_d2
                        for d1_c=1:n_d1
                            U_j(d1_c,d2_c,ap_c,a1_c,a2_c,z_c)=ReturnFn_z(d1_grid_cpu(d1_c),d2_grid_cpu(d2_c),a1_grid_bf(ap_c),a1_grid_bf(a1_c),a2_grid_bf(a2_c),z_grid_cpu(z_c),Params.r,Params.w,Params.kappa_j(jj),Params.sigma,Params.varphi,Params.eta,Params.agej(jj),Params.Jr,Params.pension);
                            T_j(d1_c,d2_c,ap_c,a1_c,a2_c,z_c)=TemptationFn_z(d1_grid_cpu(d1_c),d2_grid_cpu(d2_c),a1_grid_bf(ap_c),a1_grid_bf(a1_c),a2_grid_bf(a2_c),z_grid_cpu(z_c),Params.lambdaGP,Params.shiftGP,Params.r,Params.w,Params.kappa_j(jj),Params.sigma,Params.agej(jj),Params.Jr,Params.pension);
                        end
                    end
                end
            end
        end
    end
    % Continuation: EVterm(d2,a1prime,a2,z) = E[V_{j+1}(a1',aprimeFn(d2,a2),z')|z] (d1 does not affect it)
    EVterm=zeros(n_d2,n_a1_bf,n_a2_bf,n_z);
    if jj<N_j
        Vnext=V_bf(:,:,:,jj+1); % (a1',a2',z')
        for d2_c=1:n_d2
            for a2_c=1:n_a2_bf
                lo=a2p_lo(d2_c,a2_c);
                p=a2p_p(d2_c,a2_c);
                Vl=squeeze(Vnext(:,lo,:));   % (a1',z')
                Vu=squeeze(Vnext(:,lo+1,:));
                if p==1
                    EVmix=Vl;
                elseif p==0
                    EVmix=Vu;
                else
                    EVmix=p*Vl+(1-p)*Vu;
                    sk=(Vl==Vu); EVmix(sk)=Vu(sk); % skipinterp, as in the raws
                end
                EVterm(d2_c,:,a2_c,:)=reshape(EVmix*pi_z_cpu',[1,n_a1_bf,1,n_z]); % E over z'
            end
        end
    end
    MostT=max(reshape(T_j,[n_d1*n_d2*n_a1_bf,n_a1_bf,n_a2_bf,n_z]),[],1); % full joint (d1,d2,a1') most-tempting
    MostT(MostT==-Inf)=0; % a state where EVERY choice is infeasible leaves this -Inf: V there must be -Inf (as in the toolkit raws, which zero a -Inf MostTempting), not -Inf-(-Inf)=NaN
    entireRHS=U_j+T_j+Params.beta*reshape(EVterm,[1,n_d2,n_a1_bf,1,n_a2_bf,n_z]);
    [Vtmp,Ptmp]=max(reshape(entireRHS,[n_d1*n_d2*n_a1_bf,n_a1_bf,n_a2_bf,n_z]),[],1);
    V_bf(:,:,:,jj)=shiftdim(Vtmp,1)-shiftdim(MostT,1);
    P_bf(:,:,:,jj)=shiftdim(Ptmp,1);
end
[Vtk1,Ptk1]=ValueFnIter_Case1_FHorz(n_d,n_a_bf,n_z,N_j,d_grid,a_grid_bf,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp1);
Ptk1joint=squeeze(gather(Ptk1(1,:,:,:,:)))+n_d1*(squeeze(gather(Ptk1(2,:,:,:,:)))-1)+n_d1*n_d2*(squeeze(gather(Ptk1(3,:,:,:,:)))-1);
fprintf('Cross test 4 (z_noe, plain): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(gather(Vtk1(:))-V_bf(:))))
fprintf('Cross test 4 (z_noe, plain): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(Ptk1joint(:)-P_bf(:))))
[Vtk2,Ptk2]=ValueFnIter_Case1_FHorz(n_d,n_a_bf,n_z,N_j,d_grid,a_grid_bf,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfgp2);
Ptk2joint=squeeze(gather(Ptk2(1,:,:,:,:)))+n_d1*(squeeze(gather(Ptk2(2,:,:,:,:)))-1)+n_d1*n_d2*(squeeze(gather(Ptk2(3,:,:,:,:)))-1);
fprintf('Cross test 4 (z_noe, DC): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(gather(Vtk2(:))-V_bf(:))))
fprintf('Cross test 4 (z_noe, DC): toolkit vs hand-rolled brute force, this should be zero: %.3e \n',max(abs(Ptk2joint(:)-P_bf(:))))
clear U_j T_j V_bf P_bf Vtk1 Vtk2 Ptk1 Ptk2 Ptk1joint Ptk2joint EVterm

%% ==================== Cross test 5: age-dependent lambdaGP via V_Jplus1 ====================
% Temptation in the first half of life only: lambdaGP_j=0.1 for j<=half, 0 after. The second
% half must equal the standard solve; the first half must equal a short (N_j=half) GP solve
% with constant lambdaGP and V_Jplus1 taken from the standard solve. At all four tiers, so the
% new raws' V_Jplus1 branches get runtime coverage (they are otherwise only reached in-loop).
Params5=Params;
Params5.lambdaGP=[0.1*ones(1,half),zeros(1,N_j-half)];
Njs=half;
Paramsjs=Params;
Paramsjs.agej=Params.agej(1:Njs);
Paramsjs.kappa_j=Params.kappa_j(1:Njs);
Paramsjs.lambdaGP=0.1;
% plain tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp1);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd1,[],N_j);
Pstd_r=reshape(Pstd1,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, plain): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, plain): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp1c=vfgp1;
vfgp1c.V_Jplus1=Vstd1(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp1c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, plain): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp2);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd2,[],N_j);
Pstd_r=reshape(Pstd2,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp2c=vfgp2;
vfgp2c.V_Jplus1=Vstd2(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp2c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp3);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd3,[],N_j);
Pstd_r=reshape(Pstd3,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp3c=vfgp3;
vfgp3c.V_Jplus1=Vstd3(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp3c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
% DC+GI tier
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp4);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstd4,[],N_j);
Pstd_r=reshape(Pstd4,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC+GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5 (z_noe, DC+GI): second half (lambda=0) vs standard, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp4c=vfgp4;
vfgp4c.V_Jplus1=Vstd4(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp4c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5 (z_noe, DC+GI): first half vs V_Jplus1 solve, this should be zero: %.3e \n',max(abs(temp(:))))
clear V5a V5c P5a P5c V5a_r P5a_r Vstd_r Pstd_r

%% ==================== Cross test 5b: age-dependent pi_z_J (the KLM2021 configuration) ====================
% pi_z_J slice j is the transition from period j to period j+1. lambdaGP=0 must equal the
% standard ExpAsset solve under the same age-varying transition, at all four tiers: this
% catches a dropped ,jj page index in the GP raws' EV (the classic V_Jplus1/copy-paste drift
% class). Then the test-5 half-life split is repeated at the plain tier under the age-varying
% transition, which exercises GP at baseline lambdaGP jointly with age-dependent pi_z_J.
pi_z_J=gather(pi_z).*ones(1,1,N_j);
pi_z_J(:,:,1:2:N_j)=0.5*pi_z_J(:,:,1:2:N_j)+0.5*eye(n_z); % make it genuinely age-dependent
[Vstdp,Pstdp]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd1);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp1);
fprintf('Cross test 5b (z_noe, plain): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(V1(:)-Vstdp(:))))
fprintf('Cross test 5b (z_noe, plain): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(P1(:)-Pstdp(:))))
[Vstdp2,Pstdp2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd2);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp2);
fprintf('Cross test 5b (z_noe, DC): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(V1(:)-Vstdp2(:))))
fprintf('Cross test 5b (z_noe, DC): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(P1(:)-Pstdp2(:))))
[Vstdp2,Pstdp2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd3);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp3);
fprintf('Cross test 5b (z_noe, GI): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(V1(:)-Vstdp2(:))))
fprintf('Cross test 5b (z_noe, GI): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(P1(:)-Pstdp2(:))))
[Vstdp2,Pstdp2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params,DiscountFactorParamNames,[],vfstd4);
[V1,P1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params0,DiscountFactorParamNames,[],vfgp4);
fprintf('Cross test 5b (z_noe, DC+GI): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(V1(:)-Vstdp2(:))))
fprintf('Cross test 5b (z_noe, DC+GI): lambda=0 equals standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(P1(:)-Pstdp2(:))))
clear V1 P1 Vstdp2 Pstdp2
% half-life split at baseline lambdaGP under the age-varying transition (plain tier)
[V5a,P5a]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z_J,ReturnFn_z,Params5,DiscountFactorParamNames,[],vfgp1);
V5a_r=reshape(V5a,[],N_j);
P5a_r=reshape(P5a,[],N_j);
Vstd_r=reshape(Vstdp,[],N_j);
Pstd_r=reshape(Pstdp,[],N_j);
temp=V5a_r(:,half+1:N_j)-Vstd_r(:,half+1:N_j);
fprintf('Cross test 5b (z_noe, plain): second half (lambda=0) vs standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,half+1:N_j)-Pstd_r(:,half+1:N_j);
fprintf('Cross test 5b (z_noe, plain): second half (lambda=0) vs standard, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(temp(:))))
vfgp1c=vfgp1;
vfgp1c.V_Jplus1=Vstdp(:,:,:,half+1);
[V5c,P5c]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,Njs,d_grid,a_grid,z_grid,pi_z_J(:,:,1:Njs),ReturnFn_z,Paramsjs,DiscountFactorParamNames,[],vfgp1c);
temp=V5a_r(:,1:half)-reshape(V5c,[],Njs);
fprintf('Cross test 5b (z_noe, plain): first half vs V_Jplus1 solve, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(temp(:))))
temp=P5a_r(:,1:half)-reshape(P5c,[],Njs);
fprintf('Cross test 5b (z_noe, plain): first half vs V_Jplus1 solve, age-dependent pi_z_J, this should be zero: %.3e \n',max(abs(temp(:))))
clear V5a V5c P5a P5c V5a_r P5a_r Vstd_r Pstd_r Vstdp Pstdp

%% ==================== Cross test 6: temptation weakly lowers welfare ====================
% Self-control cost is nonnegative, so V_GP <= V_standard everywhere (lambdaGP=0.1 baseline)
temp=max(V1a(:)-Vstd1(:));
fprintf('Cross test 6 (z_noe, plain): max(V_GP - V_standard), this should be weakly negative: %.3e \n',temp)

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
