function output=CoreInfHorz_CrossTests_nod_nosemiz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for InfHorz with an iid e (and no z), without d.
% The idea (as in the FHorz cross-tests) is to run things that should give the
% same answer, disguised in different ways, and check that they do.

n_d=0; d_grid=[];
n_z=0; z_grid=[]; pi_z=[];

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

ReturnFn_none=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
ReturnFn_e=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
ReturnFn_z=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
% ReturnFn_e and ReturnFn_z are the same formula, c=(1+r)a+w*shock-aprime, so at the same shock
% grid and the same shock process they must give identical F

FnsToEvaluate_none.assets=@(aprime,a) a;
FnsToEvaluate_none.earnings=@(aprime,a,w) w;
FnsToEvaluate_e.assets=@(aprime,a,e) a;
FnsToEvaluate_e.earnings=@(aprime,a,e,w) w*e;
FnsToEvaluate_z.assets=@(aprime,a,z) a;
FnsToEvaluate_z.earnings=@(aprime,a,z,w) w*z;

%% Cross test D: an iid e, disguised as a markov z with identical rows, gives the same answer
% This is the headline test of the e code paths: an iid shock declared as e, versus the very same
% iid shock declared as a markov z whose rows are all identical. Same model, two code paths.
vfo_e=struct(); vfo_e.n_e=n_e; vfo_e.e_grid=e_grid; vfo_e.pi_e=pi_e;
simo_e=struct(); simo_e.n_e=n_e; simo_e.e_grid=e_grid; simo_e.pi_e=pi_e;

[Ve,Policye]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params,DiscountFactorParamNames,[],vfo_e);
StationaryDiste=StationaryDist_InfHorz(Policye,n_d,n_a,n_z,pi_z,simo_e,Params,[]);
AggVarse=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDiste,Policye,FnsToEvaluate_e,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simo_e);

% The same iid shock, but declared as a markov z (identical rows => iid)
n_z_asz=n_e;
z_grid_asz=e_grid;
pi_z_asz=repmat(pi_e',n_e,1); % every row is pi_e => z is iid with distribution pi_e

[Vz,Policyz]=ValueFnIter_InfHorz(n_d,n_a,n_z_asz,d_grid,a_grid,z_grid_asz,pi_z_asz,ReturnFn_z,Params,DiscountFactorParamNames,[],struct());
StationaryDistz=StationaryDist_InfHorz(Policyz,n_d,n_a,n_z_asz,pi_z_asz,struct(),Params,[]);
AggVarsz=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDistz,Policyz,FnsToEvaluate_z,Params,[],n_d,n_a,n_z_asz,d_grid,a_grid,z_grid_asz,struct());

fprintf('Cross test D (iid e == same iid disguised as markov z), V should be zero:              %2.8f \n',max(abs(Ve(:)-Vz(:))))
fprintf('Cross test D (iid e == same iid disguised as markov z), Policy should be zero:         %2.8f \n',max(abs(Policye(:)-Policyz(:))))
fprintf('Cross test D (iid e == same iid disguised as markov z), StationaryDist should be zero: %2.8f \n',max(abs(StationaryDiste(:)-StationaryDistz(:))))
fprintf('Cross test D (iid e == same iid disguised as markov z), assets Mean should be zero:    %2.8f \n',abs(AggVarse.assets.Mean-AggVarsz.assets.Mean))
fprintf('Cross test D (iid e == same iid disguised as markov z), earnings Mean should be zero:  %2.8f \n',abs(AggVarse.earnings.Mean-AggVarsz.earnings.Mean))

%% Cross test E: a single trivial e point (value 1, prob 1) reproduces the noe code path
vfo_e1=struct(); vfo_e1.n_e=1; vfo_e1.e_grid=1; vfo_e1.pi_e=1;
simo_e1=struct(); simo_e1.n_e=1; simo_e1.e_grid=1; simo_e1.pi_e=1;

[V0,Policy0]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_none,Params,DiscountFactorParamNames,[],struct());
StationaryDist0=StationaryDist_InfHorz(Policy0,n_d,n_a,n_z,pi_z,struct(),Params,[]);

[V0e,Policy0e]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params,DiscountFactorParamNames,[],vfo_e1);
StationaryDist0e=StationaryDist_InfHorz(Policy0e,n_d,n_a,n_z,pi_z,simo_e1,Params,[]);

fprintf('Cross test E (trivial e == noe), V should be zero:             %2.8f \n',max(abs(V0(:)-V0e(:))))
fprintf('Cross test E (trivial e == noe), Policy should be zero:        %2.8f \n',max(abs(Policy0(:)-Policy0e(:))))
fprintf('Cross test E (trivial e == noe), StationaryDist should be zero:%2.8f \n',max(abs(StationaryDist0(:)-StationaryDist0e(:))))

%% Cross test F: relabelling (permuting) the e states leaves all aggregate statistics unchanged
AllStatse=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDiste,Policye,FnsToEvaluate_e,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simo_e);

perm=n_e:-1:1;
vfo_ep=struct(); vfo_ep.n_e=n_e; vfo_ep.e_grid=e_grid(perm); vfo_ep.pi_e=pi_e(perm);
simo_ep=struct(); simo_ep.n_e=n_e; simo_ep.e_grid=e_grid(perm); simo_ep.pi_e=pi_e(perm);

[~,Policyp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params,DiscountFactorParamNames,[],vfo_ep);
StationaryDistp=StationaryDist_InfHorz(Policyp,n_d,n_a,n_z,pi_z,simo_ep,Params,[]);
AllStatsp=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDistp,Policyp,FnsToEvaluate_e,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simo_ep);

fprintf('Cross test F (e permutation invariance), assets Mean diff should be zero:   %2.8f \n',abs(AllStatse.assets.Mean-AllStatsp.assets.Mean))
fprintf('Cross test F (e permutation invariance), assets StdDev diff should be zero: %2.8f \n',abs(AllStatse.assets.StdDeviation-AllStatsp.assets.StdDeviation))
fprintf('Cross test F (e permutation invariance), earnings Mean diff should be zero: %2.8f \n',abs(AllStatse.earnings.Mean-AllStatsp.earnings.Mean))
fprintf('Cross test F (e permutation invariance), earnings Gini diff should be zero: %2.8f \n',abs(AllStatse.earnings.Gini-AllStatsp.earnings.Gini))

%%
output=struct();

end
