function output=CoreInfHorz_CrossTests_d_nosemiz_e(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% Cross-tests for InfHorz with an iid e, with d: D, E and F without z, G with both z and e.
% The idea (as in the FHorz cross-tests) is to run things that should give the
% same answer, disguised in different ways, and check that they do.

% Cross tests D, E and F have no z; cross test G (at the end) has both z and e, so keep the z inputs for it
n_z_withz=n_z; z_grid_withz=z_grid; pi_z_withz=pi_z;
n_z=0; z_grid=[]; pi_z=[];

n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;

ReturnFn_none=@(d,aprime,a,r,w,sigma,eta,varphi) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,sigma,eta,varphi);
ReturnFn_e=@(d,aprime,a,e,r,w,sigma,eta,varphi) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,sigma,eta,varphi);
ReturnFn_z=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
% ReturnFn_e and ReturnFn_z are the same formula, c=(1+r)a+w*d*shock-aprime, so at the same shock
% grid and the same shock process they must give identical F

FnsToEvaluate_none.assets=@(d,aprime,a) a;
FnsToEvaluate_none.earnings=@(d,aprime,a,w) w*d;
FnsToEvaluate_e.assets=@(d,aprime,a,e) a;
FnsToEvaluate_e.earnings=@(d,aprime,a,e,w) w*d*e;
FnsToEvaluate_z.assets=@(d,aprime,a,z) a;
FnsToEvaluate_z.earnings=@(d,aprime,a,z,w) w*d*z;

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

fprintf('Cross test D (iid e == same iid disguised as markov z), V should be zero:              %.3e \n',max(abs(Ve(:)-Vz(:))))
fprintf('Cross test D (iid e == same iid disguised as markov z), Policy should be zero:         %.3e \n',max(abs(Policye(:)-Policyz(:))))
fprintf('Cross test D (iid e == same iid disguised as markov z), StationaryDist should be zero: %.3e \n',max(abs(StationaryDiste(:)-StationaryDistz(:))))
fprintf('Cross test D (iid e == same iid disguised as markov z), assets Mean should be zero:    %.3e \n',abs(AggVarse.assets.Mean-AggVarsz.assets.Mean))
fprintf('Cross test D (iid e == same iid disguised as markov z), earnings Mean should be zero:  %.3e \n',abs(AggVarse.earnings.Mean-AggVarsz.earnings.Mean))

%% Cross test E: a single trivial e point (value 1, prob 1) reproduces the noe code path
vfo_e1=struct(); vfo_e1.n_e=1; vfo_e1.e_grid=1; vfo_e1.pi_e=1;
simo_e1=struct(); simo_e1.n_e=1; simo_e1.e_grid=1; simo_e1.pi_e=1;

[V0,Policy0]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_none,Params,DiscountFactorParamNames,[],struct());
StationaryDist0=StationaryDist_InfHorz(Policy0,n_d,n_a,n_z,pi_z,struct(),Params,[]);

[V0e,Policy0e]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params,DiscountFactorParamNames,[],vfo_e1);
StationaryDist0e=StationaryDist_InfHorz(Policy0e,n_d,n_a,n_z,pi_z,simo_e1,Params,[]);

fprintf('Cross test E (trivial e == noe), V should be zero:             %.3e \n',max(abs(V0(:)-V0e(:))))
fprintf('Cross test E (trivial e == noe), Policy should be zero:        %.3e \n',max(abs(Policy0(:)-Policy0e(:))))
fprintf('Cross test E (trivial e == noe), StationaryDist should be zero:%.3e \n',max(abs(StationaryDist0(:)-StationaryDist0e(:))))

%% Cross test F: relabelling (permuting) the e states leaves all aggregate statistics unchanged
AllStatse=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDiste,Policye,FnsToEvaluate_e,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simo_e);

perm=n_e:-1:1;
vfo_ep=struct(); vfo_ep.n_e=n_e; vfo_ep.e_grid=e_grid(perm); vfo_ep.pi_e=pi_e(perm);
simo_ep=struct(); simo_ep.n_e=n_e; simo_ep.e_grid=e_grid(perm); simo_ep.pi_e=pi_e(perm);

[~,Policyp]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params,DiscountFactorParamNames,[],vfo_ep);
StationaryDistp=StationaryDist_InfHorz(Policyp,n_d,n_a,n_z,pi_z,simo_ep,Params,[]);
AllStatsp=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDistp,Policyp,FnsToEvaluate_e,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simo_ep);

fprintf('Cross test F (e permutation invariance), assets Mean diff should be zero:   %.3e \n',abs(AllStatse.assets.Mean-AllStatsp.assets.Mean))
fprintf('Cross test F (e permutation invariance), assets StdDev diff should be zero: %.3e \n',abs(AllStatse.assets.StdDeviation-AllStatsp.assets.StdDeviation))
fprintf('Cross test F (e permutation invariance), earnings Mean diff should be zero: %.3e \n',abs(AllStatse.earnings.Mean-AllStatsp.earnings.Mean))
fprintf('Cross test F (e permutation invariance), earnings Gini diff should be zero: %.3e \n',abs(AllStatse.earnings.Gini-AllStatsp.earnings.Gini))

%% Cross test G: with both z and e, (z,e) disguised as a two-dimensional markov z gives the same answer
% The markov z and the iid e together are the same thing as a two-dimensional markov shock (z,e) whose transition
% matrix is kron(ones(n_e,1)*pi_e',pi_z): row/column (z,e) with z varying first, and probability pi_z(z,zprime)*pi_e(eprime).
ReturnFn_ze=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
FnsToEvaluate_ze.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate_ze.earnings=@(d,aprime,a,z,e,w) w*d*z*e;

[Vze,Policyze]=ValueFnIter_InfHorz(n_d,n_a,n_z_withz,d_grid,a_grid,z_grid_withz,pi_z_withz,ReturnFn_ze,Params,DiscountFactorParamNames,[],vfo_e);
StationaryDistze=StationaryDist_InfHorz(Policyze,n_d,n_a,n_z_withz,pi_z_withz,simo_e,Params,[]);
AggVarsze=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDistze,Policyze,FnsToEvaluate_ze,Params,[],n_d,n_a,n_z_withz,d_grid,a_grid,z_grid_withz,simo_e);

% The same shocks, but declared as a two-dimensional markov z (the ReturnFn and FnsToEvaluate are unchanged, the e input is just the second z)
n_z_asz2=[n_z_withz,n_e];
z_grid_asz2=[z_grid_withz;e_grid];
pi_z_asz2=kron(ones(n_e,1)*pi_e',pi_z_withz);

[Vz2,Policyz2]=ValueFnIter_InfHorz(n_d,n_a,n_z_asz2,d_grid,a_grid,z_grid_asz2,pi_z_asz2,ReturnFn_ze,Params,DiscountFactorParamNames,[],struct());
StationaryDistz2=StationaryDist_InfHorz(Policyz2,n_d,n_a,n_z_asz2,pi_z_asz2,struct(),Params,[]);
AggVarsz2=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDistz2,Policyz2,FnsToEvaluate_ze,Params,[],n_d,n_a,n_z_asz2,d_grid,a_grid,z_grid_asz2,struct());

fprintf('Cross test G ((z,e) == two-dim markov z), V should be zero:              %.3e \n',max(abs(Vze(:)-Vz2(:))))
fprintf('Cross test G ((z,e) == two-dim markov z), Policy should be zero:         %.3e \n',max(abs(Policyze(:)-Policyz2(:))))
fprintf('Cross test G ((z,e) == two-dim markov z), StationaryDist should be zero: %.3e \n',max(abs(StationaryDistze(:)-StationaryDistz2(:))))
fprintf('Cross test G ((z,e) == two-dim markov z), assets Mean should be zero:    %.3e \n',abs(AggVarsze.assets.Mean-AggVarsz2.assets.Mean))
fprintf('Cross test G ((z,e) == two-dim markov z), earnings Mean should be zero:  %.3e \n',abs(AggVarsze.earnings.Mean-AggVarsz2.earnings.Mean))

%%
output=struct();

end
