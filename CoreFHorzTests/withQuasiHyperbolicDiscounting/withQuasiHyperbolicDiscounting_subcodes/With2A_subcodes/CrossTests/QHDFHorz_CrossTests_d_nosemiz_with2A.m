function output=QHDFHorz_CrossTests_d_nosemiz_with2A(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% V AND POLICY ONLY. The cross tests below compare V and Policy, and nothing downstream of them.
% They used to also compare the agent distributions (56 zero-check lines across the ten
% QH cross-test subcodes, removed 2026-09-23). Those comparisons did not test quasi-hyperbolic
% discounting. StationaryDist is built FROM Policy by code that never looks at
% vfoptions.exoticpreferences, so the dist comparison on each of those lines was the baseline
% z-vs-e / semiz-vs-z distribution cross-test being re-run at exotic-preference prices: with the
% Policy comparison on the line immediately above it already printing an exact zero, the only thing
% the dist line could still catch was a disagreement between two StationaryDist branches, and
% CoreFHorzTests.m makes exactly that comparison, on exactly these shapes, already.
%
% This does NOT extend to the with/without-grid-interpolation moment comparisons in the figure
% subcodes ('should get much the same moments (for big a_grid)'). Those ARE a convergence check on
% the QH solver itself and must stay: grid interpolation returns policies off the coarse grid, so
% the tiers never agree elementwise and V/Policy cannot carry the check - the moment level is the
% only place it can be made.

% For crosstests, set up z to just be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in same place in earnings

% Setup vfoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;

ReturnFn_none=@(d,a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2) ReturnFn_d_noz_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2);
ReturnFn_z=@(d,a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2) ReturnFn_d_z_noe_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2);
ReturnFn_e=@(d,a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2) ReturnFn_d_noz_e_nosemiz_with2A(d,a1prime,a2prime,a1,a2,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2);
ReturnFn_ze=@(d,a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2) ReturnFn_d_z_e_nosemiz_with2A(d,a1prime,a2prime,a1,a2,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,phi1,phi2);

% Setup some FnsToEvaluate
% FnsToEvaluate_z.assets=@(d,aprime,a,z) a;
% FnsToEvaluate_z.earnings=@(d,aprime,a,z,w,kappa_j) w*kappa_j*z*d;
% FnsToEvaluate_e.assets=@(d,aprime,a,e) a;
% FnsToEvaluate_e.earnings=@(d,aprime,a,e,w,kappa_j) w*kappa_j*e*d;


%% Quasi-Hyperbolic Discounting
% vfoptions.exoticpreferences='QuasiHyperbolic';
%% Naive
% vfoptions.quasi_hyperbolic='Naive';


%% Solving with just a single points for z with value 1 and prob 1 gives us same as no shocks

vfoptions_z.exoticpreferences='QuasiHyperbolic';
vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_z.quasi_hyperbolic='Naive';
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_z);

vfoptions_z.exoticpreferences='QuasiHyperbolic';
vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_z.quasi_hyperbolic='Naive';
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy0z(:))))

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid

vfoptions_z.exoticpreferences='QuasiHyperbolic';
vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_z.quasi_hyperbolic='Naive';
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);

vfoptions_e=vfoptions;
vfoptions_e.exoticpreferences='QuasiHyperbolic';
vfoptions_e.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_e.quasi_hyperbolic='Naive';
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_e);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
vfoptions_ze1=vfoptions;
vfoptions_ze1.exoticpreferences='QuasiHyperbolic';
vfoptions_ze1.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_ze1.quasi_hyperbolic='Naive';
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze1);
V3=squeeze(V3);
Policy3=squeeze(Policy3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))

% Second, make e just 1
vfoptions_ze2.n_e=1;
vfoptions_ze2.e_grid=1;
vfoptions_ze2.pi_e=1;
vfoptions_ze2.exoticpreferences='QuasiHyperbolic';
vfoptions_ze2.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_ze2.quasi_hyperbolic='Naive';
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
V4=squeeze(V4);
Policy4=squeeze(Policy4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))



%% Sophisticated
vfoptions.quasi_hyperbolic='Sophisticated';

%% Solving with just a single points for z with value 1 and prob 1 gives us same as no shocks

vfoptions_z.exoticpreferences='QuasiHyperbolic';
vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_z.quasi_hyperbolic='Sophisticated';
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions_z);

vfoptions_z.exoticpreferences='QuasiHyperbolic';
vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_z.quasi_hyperbolic='Sophisticated';
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy0(:)-Policy0z(:))))

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid

vfoptions_z.exoticpreferences='QuasiHyperbolic';
vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_z.quasi_hyperbolic='Sophisticated';
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);

vfoptions_e=vfoptions;
vfoptions_e.exoticpreferences='QuasiHyperbolic';
vfoptions_e.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_e.quasi_hyperbolic='Sophisticated';
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptions_e);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy2(:))))

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
vfoptions_ze1=vfoptions;
vfoptions_ze1.exoticpreferences='QuasiHyperbolic';
vfoptions_ze1.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_ze1.quasi_hyperbolic='Sophisticated';
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze1);
V3=squeeze(V3);
Policy3=squeeze(Policy3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy3(:))))

% Second, make e just 1
vfoptions_ze2.n_e=1;
vfoptions_ze2.e_grid=1;
vfoptions_ze2.pi_e=1;
vfoptions_ze2.exoticpreferences='QuasiHyperbolic';
vfoptions_ze2.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
vfoptions_ze2.quasi_hyperbolic='Sophisticated';
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptions_ze2);
V4=squeeze(V4);
Policy4=squeeze(Policy4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(Policy1(:)-Policy4(:))))



%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
