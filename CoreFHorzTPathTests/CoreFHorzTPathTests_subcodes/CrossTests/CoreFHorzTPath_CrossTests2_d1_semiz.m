function output=CoreFHorzTPath_CrossTests2_d1_semiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% TPath version of CoreFHorz_CrossTests2_d1_semiz: set up a semi-exo state whose
% transition does NOT depend on the decision (SemiExoStateFn_JustAMarkov), so it is
% really just a markov z. Solve as semiz and as a plain markov z and compare the whole
% paths (V, Policy, AgentDist) which should be identical.

% n_d=n_d_semiz;  d_grid=d_grid_semiz;

% For models without semiz we still need
n_d1=n_d(1);
d1_grid=d_grid(1:n_d(1));

Params.uempbenefit=0; % Need this to make return fns the same
Params.searcheffortcost=0; % d2 does nothing in these CrossTests2 examples anyway (d2=0 is a valid choice), so the return fns are the same.

% WARNING: THE z and semiz IN HERE ARE PRETTY HARDCODED!
n_z=2;
z_grid=[0.6;1.4];
pi_z=[1-Params.probfindjob, Params.probfindjob;...
    Params.problosejob, 1-Params.problosejob];
Params.z1=z_grid(1);
Params.z2=z_grid(2);

% The semiz is just a duplicate of z (transition independent of the decision)
SemiExoStateFn_JustAMarkov=@(n,nprime,dsemiz,probfindjob,problosejob,z1,z2) CoreFHorzSetup_SemiExoStateFn_JustAMarkov(n,nprime,dsemiz,probfindjob,problosejob,z1,z2);
vfoptions.n_semiz=n_z;
vfoptions.semiz_grid=z_grid;
vfoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.n_semiz=n_z;
simoptions.semiz_grid=z_grid;
simoptions.SemiExoStateFn=SemiExoStateFn_JustAMarkov;
simoptions.d_grid=d_grid;
% For convenience
n_semiz=n_z;

% Setup vfoptions and simoptions (also carry e for the second block)
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

ReturnFn_semiz=@(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_semize=@(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);

ReturnFn_z=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)...
    ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_ze=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension)...
    ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

%% Solving for model with one markov
jequaloneDist1=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2))=1; % no assets

% First, just use z (without semiz)
vfoptions1A.divideandconquer=0;
simoptions1A=struct();
V_final1A=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final1A=ones([2,n_a,n_z,N_j],'gpuArray');
[VPath1A,PolicyPath1A]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final1A, Policy_final1A, Params, n_d1, n_a, n_z, N_j, d1_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, ReturnFn_z, transpathoptionsbaseline, vfoptions1A);
AgentDist_initial1A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,PolicyPath1A(:,:,:,:,1),n_d1,n_a,n_z,N_j,pi_z,Params,simoptions1A);
AgentDistPath1A=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial1A, jequaloneDist1, PricePath, ParamPath, PolicyPath1A, AgeWeightParamNames,n_d1,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions1A);

% Second, use semiz (without z)
vfoptions1B.divideandconquer=0;
vfoptions1B.n_semiz=vfoptions.n_semiz;
vfoptions1B.semiz_grid=vfoptions.semiz_grid;
vfoptions1B.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions1B.n_semiz=simoptions.n_semiz;
simoptions1B.semiz_grid=simoptions.semiz_grid;
simoptions1B.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptions1B.d_grid=simoptions.d_grid;
V_final1B=zeros([n_a,n_semiz,N_j],'gpuArray');
Policy_final1B=ones([3,n_a,n_semiz,N_j],'gpuArray');
[VPath1B,PolicyPath1B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final1B, Policy_final1B, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[],[], DiscountFactorParamNames, ReturnFn_semiz, transpathoptionsbaseline, vfoptions1B);
AgentDist_initial1B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,PolicyPath1B(:,:,:,:,1),n_d,n_a,0,N_j,[],Params,simoptions1B);
AgentDistPath1B=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial1B, jequaloneDist1, PricePath, ParamPath, PolicyPath1B, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptions1B);

PolicyPath1Bshort=PolicyPath1B([1,3],:,:,:,:); % remove the d2 policy (keep d1 and aprime; not in PolicyPath1A)

fprintf('Cross test: semiz as z, this should be zero: %.3e \n',max(abs(VPath1A(:)-VPath1B(:))))
fprintf('Cross test: semiz as z, this should be zero: %.3e \n',max(abs(PolicyPath1A(:)-PolicyPath1Bshort(:))))
fprintf('Cross test: semiz as z, this should be zero: %.3e \n',max(abs(AgentDistPath1A(:)-AgentDistPath1B(:))))

clear VPath1A PolicyPath1A AgentDistPath1A VPath1B PolicyPath1B PolicyPath1Bshort AgentDistPath1B

%% Solving for model with one markov and one e
jequaloneDist1=zeros([n_a,n_semiz,vfoptions.n_e],'gpuArray');
jequaloneDist1(1,ceil(n_semiz/2),ceil(vfoptions.n_e/2))=1; % no assets

% First, just use z and e (without semiz)
vfoptions2A.n_e=vfoptions.n_e;
vfoptions2A.e_grid=vfoptions.e_grid;
vfoptions2A.pi_e=vfoptions.pi_e;
simoptions2A.n_e=simoptions.n_e;
simoptions2A.e_grid=simoptions.e_grid;
simoptions2A.pi_e=simoptions.pi_e;
V_final2A=zeros([n_a,n_z,vfoptions.n_e,N_j],'gpuArray');
Policy_final2A=ones([2,n_a,n_z,vfoptions.n_e,N_j],'gpuArray');
[VPath2A,PolicyPath2A]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final2A, Policy_final2A, Params, n_d1, n_a, n_z, N_j, d1_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, ReturnFn_ze, transpathoptionsbaseline, vfoptions2A);
AgentDist_initial2A=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,PolicyPath2A(:,:,:,:,:,1),n_d1,n_a,n_z,N_j,pi_z,Params,simoptions2A);
AgentDistPath2A=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial2A, jequaloneDist1, PricePath, ParamPath, PolicyPath2A, AgeWeightParamNames,n_d1,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions2A);

% Second, use semiz and e (without z)
vfoptions2B.n_semiz=vfoptions.n_semiz;
vfoptions2B.semiz_grid=vfoptions.semiz_grid;
vfoptions2B.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptions2B.n_semiz=simoptions.n_semiz;
simoptions2B.semiz_grid=simoptions.semiz_grid;
simoptions2B.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptions2B.d_grid=simoptions.d_grid;
vfoptions2B.n_e=vfoptions.n_e;
vfoptions2B.e_grid=vfoptions.e_grid;
vfoptions2B.pi_e=vfoptions.pi_e;
simoptions2B.n_e=simoptions.n_e;
simoptions2B.e_grid=simoptions.e_grid;
simoptions2B.pi_e=simoptions.pi_e;
V_final2B=zeros([n_a,n_semiz,vfoptions.n_e,N_j],'gpuArray');
Policy_final2B=ones([3,n_a,n_semiz,vfoptions.n_e,N_j],'gpuArray');
[VPath2B,PolicyPath2B]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final2B, Policy_final2B, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[],[], DiscountFactorParamNames, ReturnFn_semize, transpathoptionsbaseline, vfoptions2B);
AgentDist_initial2B=StationaryDist_FHorz_Case1(jequaloneDist1,AgeWeightParamNames,PolicyPath2B(:,:,:,:,:,1),n_d,n_a,0,N_j,[],Params,simoptions2B);
AgentDistPath2B=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial2B, jequaloneDist1, PricePath, ParamPath, PolicyPath2B, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptions2B);

PolicyPath2Bshort=PolicyPath2B([1,3],:,:,:,:,:); % remove the d2 policy (keep d1 and aprime; not in PolicyPath2A)

fprintf('Cross test: semiz as z (with e), this should be zero: %.3e \n',max(abs(VPath2A(:)-VPath2B(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %.3e \n',max(abs(PolicyPath2A(:)-PolicyPath2Bshort(:))))
fprintf('Cross test: semiz as z (with e), this should be zero: %.3e \n',max(abs(AgentDistPath2A(:)-AgentDistPath2B(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
