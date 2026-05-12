function output=CoreFHorz_CrossTests_nod1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)

% n_d=n_d2_semiz;
% d_grid=d2_grid_semiz;

vfoptions.divideandconquer=1;

% Setup semiz
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d_grid;
% For convenience
n_semiz=vfoptionsbaseline.n_semiz;

% For crosstests, set up z to just be a copy of e
n_z=vfoptionsbaseline.n_e;
pi_z=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_grid=vfoptionsbaseline.e_grid;
% NOTE: z & e appear in same place in earnings

% Setup vfoptions and simoptions
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;
% Note: all of these have semiz
ReturnFn_none=@(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_z=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_e=@(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_ze=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost)...
    ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

% Setup some FnsToEvaluate
% FnsToEvaluate_z.assets=@(aprime,a,z) a;
% FnsToEvaluate_z.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
% FnsToEvaluate_e.assets=@(aprime,a,e) a;
% FnsToEvaluate_e.earnings=@(aprime,a,e,w,kappa_j) w*kappa_j*e;


%% Solving with just a single points for z with value 1 and prob 1 gives us same as no shocks (both with semiz)
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(n_semiz/2))=1; % no assets

% optionsA: just semiz (no e)
vfoptionsA.n_semiz=vfoptions.n_semiz;
vfoptionsA.semiz_grid=vfoptions.semiz_grid;
vfoptionsA.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptionsA.n_semiz=simoptions.n_semiz;
simoptionsA.semiz_grid=simoptions.semiz_grid;
simoptionsA.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptionsA.d_grid=simoptions.d_grid;

% Use vfoptionsA, which has semiz but nothing else
[V0,Policy0]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0,n_d,n_a,0,N_j,[],Params,simoptionsA);

% Use vfoptionsA, which has semiz but nothing else
[V0z,Policy0z]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy0z,n_d,n_a,1,N_j,1,Params,simoptionsA);

fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

clear V0 Policy0 V0z Policy0z StationaryDist0 StationaryDist0z

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid
% zeros assets, mid points for any shocks
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1; % no assets, midpoint shock

% Use vfoptionsA, which has semiz but nothing else
[V1,Policy1]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptionsA);
StationaryDist1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy1,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);

% Use semiz and e
vfoptionsB=vfoptions;
simoptionsB=simoptions;
[V2,Policy2]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy2,n_d,n_a,0,N_j,[],Params,simoptionsB);

fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist2(:))))

clear V2 Policy2 StationaryDist2

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1
% Use semiz and e (vfoptionsB)
[V3,Policy3]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptionsB);
jequaloneDist3=zeros([n_a,n_semiz,1,vfoptionsB.n_e],'gpuArray');
jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptionsB.n_e/2))=1; % no assets, midpoint shock
StationaryDist3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,Policy3,n_d,n_a,1,N_j,1,Params,simoptionsB);
V3=squeeze(V3);
Policy3=squeeze(Policy3);
StationaryDist3=squeeze(StationaryDist3);

fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

% Second, make e just 1 (with semiz)
vfoptionsC=vfoptionsA; % semiz
vfoptionsC.n_e=1; % and e=1 as single point
vfoptionsC.e_grid=1;
vfoptionsC.pi_e=1;
simoptionsC=simoptionsB;
simoptionsC.n_e=1;
simoptionsC.e_grid=1;
simoptionsC.pi_e=1;
[V4,Policy4]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptionsC);
jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray');
jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1; % no assets, midpoint shock
StationaryDist4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,Policy4,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsC);
V4=squeeze(V4);
Policy4=squeeze(Policy4);
StationaryDist4=squeeze(StationaryDist4);

fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end