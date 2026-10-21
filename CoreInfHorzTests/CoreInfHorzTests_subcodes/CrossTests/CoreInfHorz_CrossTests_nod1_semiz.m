function output=CoreInfHorz_CrossTests_nod1_semiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for InfHorz with a semi-exogenous state semiz, without d1 (the only decision variable is d2, which
% determines the semiz transitions). The idea (as in the FHorz cross-tests) is to run things that should give the
% same answer, disguised in different ways, and check that they do.
% With e: a markov z that is just an iid in disguise gives the same answer as the e, and the z-and-e code with either
% one of them a single trivial point gives the same answer again (as in the FHorz cross-tests).

% Setup semiz
vfoptions=struct();
simoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d_grid;

% Note: all of these have semiz
ReturnFn_none=@(d2,aprime,a,semiz,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,sigma,uempbenefit,searcheffortcost);
ReturnFn_z=@(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost);
ReturnFn_e=@(d2,aprime,a,semiz,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,sigma,uempbenefit,searcheffortcost);
ReturnFn_ze=@(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost);
% ReturnFn_none and ReturnFn_z are the same formula when z=1 (likewise ReturnFn_ze when z=1 or e=1); z and e appear
% in the same place in earnings, so ReturnFn_z with z on the e grid is ReturnFn_e

%% Solving with just a single point for z with value 1 and prob 1 gives us the same as no z (both with semiz)
[V0,Policy0]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0=StationaryDist_InfHorz(Policy0,n_d,n_a,0,[],simoptions,Params,[]);

[V0z,Policy0z]=ValueFnIter_InfHorz(n_d,n_a,1,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist0z=StationaryDist_InfHorz(Policy0z,n_d,n_a,1,1,simoptions,Params,[]);

fprintf('Cross test (semiz, trivial z == noz), V should be zero:             %.3e \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test (semiz, trivial z == noz), Policy should be zero:        %.3e \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test (semiz, trivial z == noz), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

clear V0 Policy0 V0z Policy0z StationaryDist0 StationaryDist0z

%% Solve using a markov which is just an iid in disguise (with semiz). Should give the same result as the iid
% z is a copy of e: the e grid, with every row of pi_z equal to pi_e
n_zcopy=vfoptionsbaseline.n_e;
z_gridcopy=vfoptionsbaseline.e_grid;
pi_zcopy=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);

% semiz and z
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_zcopy,d_grid,a_grid,z_gridcopy,pi_zcopy,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist1=StationaryDist_InfHorz(Policy1,n_d,n_a,n_zcopy,pi_zcopy,simoptions,Params,[]);

% semiz and e
vfoptionsB=vfoptions;
vfoptionsB.n_e=vfoptionsbaseline.n_e;
vfoptionsB.e_grid=vfoptionsbaseline.e_grid;
vfoptionsB.pi_e=vfoptionsbaseline.pi_e;
simoptionsB=simoptions;
simoptionsB.n_e=simoptionsbaseline.n_e;
simoptionsB.e_grid=simoptionsbaseline.e_grid;
simoptionsB.pi_e=simoptionsbaseline.pi_e;
[V2,Policy2]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_e,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist2=StationaryDist_InfHorz(Policy2,n_d,n_a,0,[],simoptionsB,Params,[]);

fprintf('Cross test (semiz, z as e), V should be zero:             %.3e \n',max(abs(V1(:)-V2(:))))
fprintf('Cross test (semiz, z as e), Policy should be zero:        %.3e \n',max(abs(Policy1(:)-Policy2(:))))
fprintf('Cross test (semiz, z as e), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist1(:)-StationaryDist2(:))))

clear V2 Policy2 StationaryDist2

%% The z-and-e code (with semiz), with the 'other' shock a single point with value 1 and prob 1, gives the same answer again
% First, z is the trivial point (and e is the e)
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,1,d_grid,a_grid,1,1,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptionsB);
StationaryDist3=StationaryDist_InfHorz(Policy3,n_d,n_a,1,1,simoptionsB,Params,[]);

fprintf('Cross test (semiz, z and e, trivial z), V should be zero:             %.3e \n',max(abs(V1(:)-V3(:))))
fprintf('Cross test (semiz, z and e, trivial z), Policy should be zero:        %.3e \n',max(abs(Policy1(:)-Policy3(:))))
fprintf('Cross test (semiz, z and e, trivial z), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))

clear V3 Policy3 StationaryDist3

% Second, e is the trivial point (and z is the copy of e)
vfoptionsC=vfoptions;
vfoptionsC.n_e=1;
vfoptionsC.e_grid=1;
vfoptionsC.pi_e=1;
simoptionsC=simoptions;
simoptionsC.n_e=1;
simoptionsC.e_grid=1;
simoptionsC.pi_e=1;
[V4,Policy4]=ValueFnIter_InfHorz(n_d,n_a,n_zcopy,d_grid,a_grid,z_gridcopy,pi_zcopy,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfoptionsC);
StationaryDist4=StationaryDist_InfHorz(Policy4,n_d,n_a,n_zcopy,pi_zcopy,simoptionsC,Params,[]);

fprintf('Cross test (semiz, z and e, trivial e), V should be zero:             %.3e \n',max(abs(V1(:)-V4(:))))
fprintf('Cross test (semiz, z and e, trivial e), Policy should be zero:        %.3e \n',max(abs(Policy1(:)-Policy4(:))))
fprintf('Cross test (semiz, z and e, trivial e), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist1(:)-StationaryDist4(:))))

clear V1 Policy1 StationaryDist1 V4 Policy4 StationaryDist4

%%
output=struct();

end
