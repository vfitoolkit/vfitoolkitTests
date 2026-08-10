function output=CoreInfHorz_CrossTests_nod_nosemiz(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for InfHorz, without d.
% The idea (as in the FHorz cross-tests) is to run things that should give the
% same answer, disguised in different ways, and check that they do.
% This bank has no e, so the tests use only d and z.

n_d=0; d_grid=[];

ReturnFn_none=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
ReturnFn_z=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
% At z=1 these two give identical consumption c=(1+r)a+w-aprime, hence identical F

FnsToEvaluate_none.assets=@(aprime,a) a;
FnsToEvaluate_none.earnings=@(aprime,a,w) w;
FnsToEvaluate_z.assets=@(aprime,a,z) a;
FnsToEvaluate_z.earnings=@(aprime,a,z,w) w*z;

%% Cross test A: a single trivial z point (value 1, prob 1) reproduces the noz code path
[V0,Policy0]=ValueFnIter_InfHorz(n_d,n_a,0,d_grid,a_grid,[],[],ReturnFn_none,Params,DiscountFactorParamNames,[],struct());
StationaryDist0=StationaryDist_InfHorz(Policy0,n_d,n_a,0,[],struct(),Params,[]);

[V0z,Policy0z]=ValueFnIter_InfHorz(n_d,n_a,1,d_grid,a_grid,1,1,ReturnFn_z,Params,DiscountFactorParamNames,[],struct());
StationaryDist0z=StationaryDist_InfHorz(Policy0z,n_d,n_a,1,1,struct(),Params,[]);

fprintf('Cross test A (trivial z == noz), V should be zero:             %2.8f \n',max(abs(V0(:)-V0z(:))))
fprintf('Cross test A (trivial z == noz), Policy should be zero:        %2.8f \n',max(abs(Policy0(:)-Policy0z(:))))
fprintf('Cross test A (trivial z == noz), StationaryDist should be zero:%2.8f \n',max(abs(StationaryDist0(:)-StationaryDist0z(:))))

%% Cross test B: iid disguised as a markov (identical rows) => stationary marginal over z equals that row
probz=pi_z(1,:); % a valid probability vector (a row of the markov)
pi_z_iid=repmat(probz,n_z,1); % identical rows => z is iid

[~,Policyiid]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z_iid,ReturnFn_z,Params,DiscountFactorParamNames,[],struct());
StationaryDistiid=StationaryDist_InfHorz(Policyiid,n_d,n_a,n_z,pi_z_iid,struct(),Params,[]);
marginalz=sum(StationaryDistiid,1); % sum over a => marginal over z
fprintf('Cross test B (iid as markov), stationary marginal over z should equal the pi_z row: %2.8f \n',max(abs(marginalz(:)-probz(:))))

% For nod, earnings=w*z is exogenous, so its mean is analytically w*E[z]
AggVarsiid=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDistiid,Policyiid,FnsToEvaluate_z,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
fprintf('Cross test B (iid as markov), earnings mean should equal w*E[z]: %2.8f \n',abs(AggVarsiid.earnings.Mean-Params.w*sum(z_grid(:).*probz(:))))

%% Cross test C: relabelling (permuting) the z states leaves all aggregate statistics unchanged
% Baseline (original z ordering)
[~,Policyb]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],struct());
StationaryDistb=StationaryDist_InfHorz(Policyb,n_d,n_a,n_z,pi_z,struct(),Params,[]);
AllStatsb=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDistb,Policyb,FnsToEvaluate_z,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());

% Permuted (reverse the z ordering; apply same permutation to grid and to pi_z rows+cols)
perm=n_z:-1:1;
z_grid_perm=z_grid(perm);
pi_z_perm=pi_z(perm,perm);
[~,Policyp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid_perm,pi_z_perm,ReturnFn_z,Params,DiscountFactorParamNames,[],struct());
StationaryDistp=StationaryDist_InfHorz(Policyp,n_d,n_a,n_z,pi_z_perm,struct(),Params,[]);
AllStatsp=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDistp,Policyp,FnsToEvaluate_z,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid_perm,struct());

fprintf('Cross test C (z permutation invariance), assets Mean diff should be zero:   %2.8f \n',abs(AllStatsb.assets.Mean-AllStatsp.assets.Mean))
fprintf('Cross test C (z permutation invariance), assets StdDev diff should be zero:  %2.8f \n',abs(AllStatsb.assets.StdDeviation-AllStatsp.assets.StdDeviation))
fprintf('Cross test C (z permutation invariance), earnings Mean diff should be zero:  %2.8f \n',abs(AllStatsb.earnings.Mean-AllStatsp.earnings.Mean))
fprintf('Cross test C (z permutation invariance), earnings Gini diff should be zero:  %2.8f \n',abs(AllStatsb.earnings.Gini-AllStatsp.earnings.Gini))

%%
output=struct();

end
