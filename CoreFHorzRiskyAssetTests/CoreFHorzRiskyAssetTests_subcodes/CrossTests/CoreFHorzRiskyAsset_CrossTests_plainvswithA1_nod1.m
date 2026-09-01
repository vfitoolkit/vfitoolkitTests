function output=CoreFHorzRiskyAsset_CrossTests_plainvswithA1_nod1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (plain vs degenerate withA1): RiskyAsset, nod1, nosemiz.
% Side A: plain single-asset riskyasset (the risky asset is the only endogenous state).
% Side B: withA1 riskyasset with n_a1=1 (a degenerate safe asset that does nothing).
% Side B wraps the PLAIN ReturnFn into the withA1 argument slot (drops a1prime/a1), so the
% two models are identical; V, Policy (d-part), Dist should match bit-exact.
% n_a is scalar (the risky asset); a_grid is its grid.
%
% Analog of CoreFHorzExpAsset_CrossTests4 (noa1 vs a1=1), adapted to riskyasset.

l_d=length(n_d);

% Side B: degenerate a1=1
n_a1_dummy=1; a1_grid_dummy=0; % value irrelevant -- ignored by the wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common riskyasset options
vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.refine_d=[0,1,1];
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions=struct();
simoptions.riskyasset=1; simoptions.refine_d=vfoptions.refine_d; simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.n_u=vfoptions.n_u; simoptions.u_grid=vfoptions.u_grid; simoptions.pi_u=vfoptions.pi_u;
simoptions.d_grid=d_grid;
simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;

vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e; vfoptions_withe.e_grid=vfoptionsbaseline.e_grid; vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptionsA_e=simoptionsA; simoptionsA_e.n_e=simoptionsbaseline.n_e; simoptionsA_e.e_grid=simoptionsbaseline.e_grid; simoptionsA_e.pi_e=simoptionsbaseline.pi_e;
simoptionsB_e=simoptionsB; simoptionsB_e.n_e=simoptionsbaseline.n_e; simoptionsB_e.e_grid=simoptionsbaseline.e_grid; simoptionsB_e.pi_e=simoptionsbaseline.pi_e;

%% (1) noz, noe
ReturnFn_A=@(savings,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_B=@(savings,a1prime,a1,a2,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_noe_nosemiz(savings,a2,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_A=zeros(n_a,1,'gpuArray'); jequaloneDist_A(1)=1;
jequaloneDist_B=zeros(n_a_B,'gpuArray'); jequaloneDist_B(1,1)=1;
FnsA.assets=@(riskyshare,savings,a) a; FnsB.assets=@(riskyshare,savings,a1prime,a1,a2) a2;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA);
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB);
PdA=Policy_A(1:l_d,:); PdB=Policy_B(1:l_d,:);
fprintf('CrossTest plainvswithA1 (nod1 noz noe): should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(PdA(:)-PdB(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))

%% (2) z, noe
ReturnFn_A=@(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(savings,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_B=@(savings,a1prime,a1,a2,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_noe_nosemiz(savings,a2,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_A=zeros([n_a,n_z],'gpuArray'); jequaloneDist_A(1,ceil(n_z/2))=1;
jequaloneDist_B=zeros([n_a_B,n_z],'gpuArray'); jequaloneDist_B(1,1,ceil(n_z/2))=1;
FnsA.assets=@(riskyshare,savings,a,z) a; FnsB.assets=@(riskyshare,savings,a1prime,a1,a2,z) a2;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA);
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB);
PdA=Policy_A(1:l_d,:); PdB=Policy_B(1:l_d,:);
fprintf('CrossTest plainvswithA1 (nod1 z noe): should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(PdA(:)-PdB(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))

%% (3) noz, e
ReturnFn_A=@(savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_e_nosemiz(savings,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_B=@(savings,a1prime,a1,a2,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_noz_e_nosemiz(savings,a2,e,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_A=zeros([n_a,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_A(1,ceil(vfoptions_withe.n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_B(1,1,ceil(vfoptions_withe.n_e/2))=1;
FnsA.assets=@(riskyshare,savings,a,e) a; FnsB.assets=@(riskyshare,savings,a1prime,a1,a2,e) a2;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA_e);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA_e);
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB_e);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB_e);
PdA=Policy_A(1:l_d,:); PdB=Policy_B(1:l_d,:);
fprintf('CrossTest plainvswithA1 (nod1 noz e): should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(PdA(:)-PdB(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))

%% (4) z, e
ReturnFn_A=@(savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(savings,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
ReturnFn_B=@(savings,a1prime,a1,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod1_z_e_nosemiz(savings,a2,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_A=zeros([n_a,n_z,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_A(1,ceil(n_z/2),ceil(vfoptions_withe.n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,n_z,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_B(1,1,ceil(n_z/2),ceil(vfoptions_withe.n_e/2))=1;
FnsA.assets=@(riskyshare,savings,a,z,e) a; FnsB.assets=@(riskyshare,savings,a1prime,a1,a2,z,e) a2;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA_e);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptionsA_e);
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,pi_z,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,n_z,N_j,pi_z,Params,simoptionsB_e);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,n_z,N_j,d_grid,a_grid_B,z_grid,simoptionsB_e);
PdA=Policy_A(1:l_d,:); PdB=Policy_B(1:l_d,:);
fprintf('CrossTest plainvswithA1 (nod1 z e): should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(PdA(:)-PdB(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))

%%
output=struct();

end
