function output=CoreFHorzRiskyAsset_CrossTests_plainvswithA1_d1_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (plain vs degenerate withA1): RiskyAsset, d1, +semiz.
% As CoreFHorzRiskyAsset_CrossTests_plainvswithA1_nod1_semiz but with a d1 (labour supply h).
%   n_d=[h,riskyshare,savings,dsemiz].
% Side A: plain single-asset riskyasset + genuine employment semiz.
% Side B: withA1 riskyasset with n_a1=1 + same semiz. Side B wraps the PLAIN *_semiz ReturnFn into the
% withA1 slot (drops a1prime/a1, passes a2 as the asset), so the two models are identical.
% Only 2 shock blocks: (noz,noe,+semiz) and (noz,e,+semiz).

l_d=length(n_d);

% Side B: degenerate a1=1
n_a1_dummy=1; a1_grid_dummy=0; % value irrelevant -- ignored by the wrapped ReturnFn
n_a_B=[n_a1_dummy, n_a];
a_grid_B=[a1_grid_dummy; a_grid];

% Common riskyasset options (+semiz)
vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.refine_d=[1,1,1,1];
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
vfoptions.n_semiz=vfoptionsbaseline.n_semiz; vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid; vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions=struct();
simoptions.riskyasset=1; simoptions.refine_d=vfoptions.refine_d; simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.n_u=vfoptions.n_u; simoptions.u_grid=vfoptions.u_grid; simoptions.pi_u=vfoptions.pi_u;
simoptions.n_semiz=simoptionsbaseline.n_semiz; simoptions.semiz_grid=simoptionsbaseline.semiz_grid; simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d_grid;
simoptionsA=simoptions; simoptionsA.a_grid=a_grid;
simoptionsB=simoptions; simoptionsB.a_grid=a_grid_B;

vfoptions_withe=vfoptions;
vfoptions_withe.n_e=vfoptionsbaseline.n_e; vfoptions_withe.e_grid=vfoptionsbaseline.e_grid; vfoptions_withe.pi_e=vfoptionsbaseline.pi_e;
simoptionsA_e=simoptionsA; simoptionsA_e.n_e=simoptionsbaseline.n_e; simoptionsA_e.e_grid=simoptionsbaseline.e_grid; simoptionsA_e.pi_e=simoptionsbaseline.pi_e;
simoptionsB_e=simoptionsB; simoptionsB_e.n_e=simoptionsbaseline.n_e; simoptionsB_e.e_grid=simoptionsbaseline.e_grid; simoptionsB_e.pi_e=simoptionsbaseline.pi_e;

%% (1) noz, noe, +semiz
ReturnFn_A=@(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B=@(h,savings,dsemiz,a1prime,a1,a2,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(h,savings,dsemiz,a2,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
jequaloneDist_A=zeros([n_a,vfoptions.n_semiz],'gpuArray'); jequaloneDist_A(1,ceil(vfoptions.n_semiz/2))=1;
jequaloneDist_B=zeros([n_a_B,vfoptions.n_semiz],'gpuArray'); jequaloneDist_B(1,1,ceil(vfoptions.n_semiz/2))=1;
FnsA.assets=@(h,riskyshare,savings,dsemiz,a,semiz) a; FnsB.assets=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz) a2;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA);
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB);
PdA=Policy_A(1:l_d,:); PdB=Policy_B(1:l_d,:);
fprintf('CrossTest plainvswithA1 (d1 noz noe, +semiz): should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(PdA(:)-PdB(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))

%% (2) noz, e, +semiz
ReturnFn_A=@(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
ReturnFn_B=@(h,savings,dsemiz,a1prime,a1,a2,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(h,savings,dsemiz,a2,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);
jequaloneDist_A=zeros([n_a,vfoptions.n_semiz,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_A(1,ceil(vfoptions.n_semiz/2),ceil(vfoptions_withe.n_e/2))=1;
jequaloneDist_B=zeros([n_a_B,vfoptions.n_semiz,vfoptions_withe.n_e],'gpuArray'); jequaloneDist_B(1,1,ceil(vfoptions.n_semiz/2),ceil(vfoptions_withe.n_e/2))=1;
FnsA.assets=@(h,riskyshare,savings,dsemiz,a,semiz,e) a; FnsB.assets=@(h,riskyshare,savings,dsemiz,a1prime,a1,a2,semiz,e) a2;

[V_A,Policy_A]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_A,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_A=StationaryDist_FHorz_Case1(jequaloneDist_A,AgeWeightParamNames,Policy_A,n_d,n_a,0,N_j,[],Params,simoptionsA_e);
AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_A,Policy_A,FnsA,Params,[],n_d,n_a,0,N_j,d_grid,a_grid,[],simoptionsA_e);
[V_B,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],[],ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions_withe);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a_B,0,N_j,[],Params,simoptionsB_e);
AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsB,Params,[],n_d,n_a_B,0,N_j,d_grid,a_grid_B,[],simoptionsB_e);
PdA=Policy_A(1:l_d,:); PdB=Policy_B(1:l_d,:);
fprintf('CrossTest plainvswithA1 (d1 noz e, +semiz): should be zero: V %.3e, Policy %.3e, Dist %.3e, AllStats.Mean %.3e \n',max(abs(V_A(:)-V_B(:))),max(abs(PdA(:)-PdB(:))),max(abs(StationaryDist_A(:)-StationaryDist_B(:))),abs(AllStats_A.assets.Mean-AllStats_B.assets.Mean))

%%
output=struct();

end
