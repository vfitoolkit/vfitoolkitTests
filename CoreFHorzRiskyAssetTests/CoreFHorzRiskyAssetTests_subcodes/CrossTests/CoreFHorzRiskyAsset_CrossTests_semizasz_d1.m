function output=CoreFHorzRiskyAsset_CrossTests_semizasz_d1(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% CrossTest (semiz-as-z): RiskyAsset, d1, single asset.
% A semiz whose transitions do NOT depend on the decision (JustAMarkov) is just a markov.
% With searcheffortcost=0 (so dsemiz is inert) and uempbenefit=0, the semiz model reproduces
% the z-markov model exactly. V and Dist match bit-exact; Policy matches on [h,riskyshare,savings]
% (the inert dsemiz policy row is dropped, as its argmax is arbitrary).
%
% Driver passes the SEMIZ n_d=[n_d1,n_d2,n_d3,n_d4]=[h,riskyshare,savings,dsemiz]; d_grid=[d1;d2;d3;d4].
% The z-side uses just the first three decisions.

% z-side decisions (drop the semiz decision dsemiz)
n_d_z=n_d(1:3);
d_grid_z=d_grid(1:(n_d(1)+n_d(2)+n_d(3)));

% Make dsemiz inert and unemployment payoff zero, so semiz==z-markov
Params.searcheffortcost=0;
Params.uempbenefit=0;

% z-markov equal to the JustAMarkov semiz: grid [0;1], transitions from probfindjob/problosejob
z1=0; z2=1;
Params.z1=z1; Params.z2=z2;
z_grid=[z1;z2];
pi_z=[1-Params.probfindjob, Params.probfindjob; Params.problosejob, 1-Params.problosejob];
n_z=2;

SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzRiskyAssetSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);

% Common riskyasset options
vfoptions=struct();
vfoptions.riskyasset=1;
vfoptions.aprimeFn=vfoptionsbaseline.aprimeFn;
vfoptions.n_u=vfoptionsbaseline.n_u; vfoptions.u_grid=vfoptionsbaseline.u_grid; vfoptions.pi_u=vfoptionsbaseline.pi_u;
simoptions=struct();
simoptions.riskyasset=1; simoptions.aprimeFn=vfoptions.aprimeFn;
simoptions.n_u=vfoptions.n_u; simoptions.u_grid=vfoptions.u_grid; simoptions.pi_u=vfoptions.pi_u;
simoptions.a_grid=a_grid;

% z-side options (3 decisions, refine_d=[1,1,1])
vfoptions_z=vfoptions; vfoptions_z.refine_d=[1,1,1];
simoptions_z=simoptions; simoptions_z.refine_d=vfoptions_z.refine_d; simoptions_z.d_grid=d_grid_z;

% semiz-side options (4 decisions, refine_d=[1,1,1,1], JustAMarkov semiz)
vfoptions_s=vfoptions; vfoptions_s.refine_d=[1,1,1,1];
vfoptions_s.n_semiz=2; vfoptions_s.semiz_grid=[z1;z2]; vfoptions_s.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions_s=simoptions; simoptions_s.refine_d=vfoptions_s.refine_d; simoptions_s.d_grid=d_grid;
simoptions_s.n_semiz=vfoptions_s.n_semiz; simoptions_s.semiz_grid=vfoptions_s.semiz_grid; simoptions_s.SemiExoStateFn=vfoptions_s.SemiExoStateFn;

%% (1) no e
jequaloneDist_z=zeros([n_a,n_z],'gpuArray'); jequaloneDist_z(1,1)=1;
jequaloneDist_s=zeros([n_a,vfoptions_s.n_semiz],'gpuArray'); jequaloneDist_s(1,1)=1;

ReturnFn_z=@(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_z_noe_nosemiz(h,savings,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_s=@(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(h,savings,dsemiz,a,semiz,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);

[V_z,Policy_z]=ValueFnIter_Case1_FHorz(n_d_z,n_a,n_z,N_j,d_grid_z,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z);
StationaryDist_z=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy_z,n_d_z,n_a,n_z,N_j,pi_z,Params,simoptions_z);

[V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_s,Params,DiscountFactorParamNames,[],vfoptions_s);
StationaryDist_s=StationaryDist_FHorz_Case1(jequaloneDist_s,AgeWeightParamNames,Policy_s,n_d,n_a,0,N_j,[],Params,simoptions_s);

Pd_z=Policy_z(1:3,:); Pd_s=Policy_s(1:3,:);
fprintf('CrossTest semizasz (d1 noe): should be zero: V %2.8f, Policy(h,riskyshare,savings) %2.8f, Dist %2.8f \n',max(abs(V_z(:)-V_s(:))),max(abs(Pd_z(:)-Pd_s(:))),max(abs(StationaryDist_z(:)-StationaryDist_s(:))))

%% (2) with e
vfoptions_z_e=vfoptions_z; vfoptions_z_e.n_e=vfoptionsbaseline.n_e; vfoptions_z_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_z_e.pi_e=vfoptionsbaseline.pi_e;
simoptions_z_e=simoptions_z; simoptions_z_e.n_e=simoptionsbaseline.n_e; simoptions_z_e.e_grid=simoptionsbaseline.e_grid; simoptions_z_e.pi_e=simoptionsbaseline.pi_e;
vfoptions_s_e=vfoptions_s; vfoptions_s_e.n_e=vfoptionsbaseline.n_e; vfoptions_s_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_s_e.pi_e=vfoptionsbaseline.pi_e;
simoptions_s_e=simoptions_s; simoptions_s_e.n_e=simoptionsbaseline.n_e; simoptions_s_e.e_grid=simoptionsbaseline.e_grid; simoptions_s_e.pi_e=simoptionsbaseline.pi_e;

jequaloneDist_z=zeros([n_a,n_z,vfoptions_z_e.n_e],'gpuArray'); jequaloneDist_z(1,1,ceil(vfoptions_z_e.n_e/2))=1;
jequaloneDist_s=zeros([n_a,vfoptions_s.n_semiz,vfoptions_s_e.n_e],'gpuArray'); jequaloneDist_s(1,1,ceil(vfoptions_s_e.n_e/2))=1;

ReturnFn_z=@(h,savings,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d1_z_e_nosemiz(h,savings,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_s=@(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(h,savings,dsemiz,a,semiz,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension,uempbenefit,searcheffortcost);

[V_z,Policy_z]=ValueFnIter_Case1_FHorz(n_d_z,n_a,n_z,N_j,d_grid_z,a_grid,z_grid,pi_z,ReturnFn_z,Params,DiscountFactorParamNames,[],vfoptions_z_e);
StationaryDist_z=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy_z,n_d_z,n_a,n_z,N_j,pi_z,Params,simoptions_z_e);

[V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],ReturnFn_s,Params,DiscountFactorParamNames,[],vfoptions_s_e);
StationaryDist_s=StationaryDist_FHorz_Case1(jequaloneDist_s,AgeWeightParamNames,Policy_s,n_d,n_a,0,N_j,[],Params,simoptions_s_e);

Pd_z=Policy_z(1:3,:); Pd_s=Policy_s(1:3,:);
fprintf('CrossTest semizasz (d1 e): should be zero: V %2.8f, Policy(h,riskyshare,savings) %2.8f, Dist %2.8f \n',max(abs(V_z(:)-V_s(:))),max(abs(Pd_z(:)-Pd_s(:))),max(abs(StationaryDist_z(:)-StationaryDist_s(:))))

%%
output=struct();

end
