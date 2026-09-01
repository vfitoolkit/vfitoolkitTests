function output=CoreFHorzTPath_CrossTests_d1_semiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% TPath version of CoreFHorz_CrossTests_d1_semiz: each ValueFnIter+StationaryDist
% is replaced by ValueFnOnTransPath_Case1_FHorz + AgentDistOnTransPath_Case1_FHorz,
% and we compare the whole paths (V, Policy, AgentDist) which should be identical.

% n_d=n_d_semiz;  d_grid=d_grid_semiz; (d1 plus the binary d2 that drives the semi-exo state)

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

% Setup vfoptions and simoptions (semiz plus e)
vfoptions.n_e=vfoptionsbaseline.n_e;
vfoptions.e_grid=vfoptionsbaseline.e_grid;
vfoptions.pi_e=vfoptionsbaseline.pi_e;
simoptions.n_e=simoptionsbaseline.n_e;
simoptions.e_grid=simoptionsbaseline.e_grid;
simoptions.pi_e=simoptionsbaseline.pi_e;

ReturnFn_none=@(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_z=@(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_e=@(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
ReturnFn_ze=@(d1,d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost)...
    ReturnFn_d1_z_e_semiz(d1,d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);

% optionsA: just semiz (no e)
vfoptionsA.n_semiz=vfoptions.n_semiz;
vfoptionsA.semiz_grid=vfoptions.semiz_grid;
vfoptionsA.SemiExoStateFn=vfoptions.SemiExoStateFn;
simoptionsA.n_semiz=simoptions.n_semiz;
simoptionsA.semiz_grid=simoptions.semiz_grid;
simoptionsA.SemiExoStateFn=simoptions.SemiExoStateFn;
simoptionsA.d_grid=simoptions.d_grid;
% optionsB: semiz and e
vfoptionsB=vfoptions;
simoptionsB=simoptions;

%% Solving with just a single point for z with value 1 and prob 1 gives us same as no shocks (both with semiz)
jequaloneDist_none=zeros([n_a,n_semiz],'gpuArray');
jequaloneDist_none(1,ceil(n_semiz/2))=1; % no assets

% semiz, no z
V_final0=zeros([n_a,n_semiz,N_j],'gpuArray');
Policy_final0=ones([3,n_a,n_semiz,N_j],'gpuArray');
[VPath0,PolicyPath0]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final0, Policy_final0, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[],[], DiscountFactorParamNames, ReturnFn_none, transpathoptionsbaseline, vfoptionsA);
AgentDist_initial0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,PolicyPath0(:,:,:,:,1),n_d,n_a,0,N_j,[],Params,simoptionsA);
AgentDistPath0=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial0, jequaloneDist_none, PricePath, ParamPath, PolicyPath0, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptionsA);

% semiz, single point z=1
V_final0z=zeros([n_a,n_semiz,1,N_j],'gpuArray');
Policy_final0z=ones([3,n_a,n_semiz,1,N_j],'gpuArray');
[VPath0z,PolicyPath0z]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final0z, Policy_final0z, Params, n_d, n_a, 1, N_j, d_grid, a_grid,1,1, DiscountFactorParamNames, ReturnFn_z, transpathoptionsbaseline, vfoptionsA);
AgentDist_initial0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,PolicyPath0z(:,:,:,:,:,1),n_d,n_a,1,N_j,1,Params,simoptionsA);
AgentDistPath0z=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial0z, jequaloneDist_none, PricePath, ParamPath, PolicyPath0z, AgeWeightParamNames,n_d,n_a,1,N_j,1, T,Params, transpathoptionsbaseline, simoptionsA);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(VPath0(:)-VPath0z(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(PolicyPath0(:)-PolicyPath0z(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(AgentDistPath0(:)-AgentDistPath0z(:))))

clear VPath0 PolicyPath0 VPath0z PolicyPath0z AgentDistPath0 AgentDistPath0z

%% Solve using a markov which is just an iid in disguise. Should give same result as the iid
% zeros assets, mid points for any shocks
jequaloneDist_z=zeros([n_a,n_semiz,n_z],'gpuArray');
jequaloneDist_z(1,ceil(n_semiz/2),ceil(n_z/2))=1; % no assets, midpoint shock

% semiz with markov z (iid in disguise)
V_final1=zeros([n_a,n_semiz,n_z,N_j],'gpuArray');
Policy_final1=ones([3,n_a,n_semiz,n_z,N_j],'gpuArray');
[VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final1, Policy_final1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, ReturnFn_z, transpathoptionsbaseline, vfoptionsA);
AgentDist_initial1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,PolicyPath1(:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptionsA);
AgentDistPath1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial1, jequaloneDist_z, PricePath, ParamPath, PolicyPath1, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptionsA);

% semiz with iid e
V_final2=zeros([n_a,n_semiz,vfoptionsB.n_e,N_j],'gpuArray');
Policy_final2=ones([3,n_a,n_semiz,vfoptionsB.n_e,N_j],'gpuArray');
[VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final2, Policy_final2, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[],[], DiscountFactorParamNames, ReturnFn_e, transpathoptionsbaseline, vfoptionsB);
AgentDist_initial2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,PolicyPath2(:,:,:,:,:,1),n_d,n_a,0,N_j,[],Params,simoptionsB);
AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial2, jequaloneDist_z, PricePath, ParamPath, PolicyPath2, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptionsB);

fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(VPath1(:)-VPath2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))
fprintf('Cross test: z as e, this should be zero: %.3e \n',max(abs(AgentDistPath1(:)-AgentDistPath2(:))))

clear VPath2 PolicyPath2 AgentDistPath2

%% Now use code with z and e, but just set the 'other' to be a single point with value 1 and prob 1
% So it should again give same answer

% First, make z just 1 (semiz and e, with z=1 single point)
V_final3=zeros([n_a,n_semiz,1,vfoptionsB.n_e,N_j],'gpuArray');
Policy_final3=ones([3,n_a,n_semiz,1,vfoptionsB.n_e,N_j],'gpuArray');
[VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final3, Policy_final3, Params, n_d, n_a, 1, N_j, d_grid, a_grid,1,1, DiscountFactorParamNames, ReturnFn_ze, transpathoptionsbaseline, vfoptionsB);
jequaloneDist3=zeros([n_a,n_semiz,1,vfoptionsB.n_e],'gpuArray');
jequaloneDist3(1,ceil(n_semiz/2),1,ceil(vfoptionsB.n_e/2))=1; % no assets, midpoint shock
AgentDist_initial3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,PolicyPath3(:,:,:,:,:,:,1),n_d,n_a,1,N_j,1,Params,simoptionsB);
AgentDistPath3=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial3, jequaloneDist3, PricePath, ParamPath, PolicyPath3, AgeWeightParamNames,n_d,n_a,1,N_j,1, T,Params, transpathoptionsbaseline, simoptionsB);
VPath3=squeeze(VPath3);
PolicyPath3=squeeze(PolicyPath3);
AgentDistPath3=squeeze(AgentDistPath3);

fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(VPath1(:)-VPath3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath3(:))))
fprintf('Cross test: z and e 1, this should be zero: %.3e \n',max(abs(AgentDistPath1(:)-AgentDistPath3(:))))

% Second, make e just 1 (semiz with markov z, e=1 single point)
vfoptionsC=vfoptionsA; % semiz
vfoptionsC.n_e=1; % and e=1 as single point
vfoptionsC.e_grid=1;
vfoptionsC.pi_e=1;
simoptionsC=simoptionsB;
simoptionsC.n_e=1;
simoptionsC.e_grid=1;
simoptionsC.pi_e=1;
V_final4=zeros([n_a,n_semiz,n_z,1,N_j],'gpuArray');
Policy_final4=ones([3,n_a,n_semiz,n_z,1,N_j],'gpuArray');
[VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final4, Policy_final4, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, ReturnFn_ze, transpathoptionsbaseline, vfoptionsC);
jequaloneDist4=zeros([n_a,n_semiz,n_z,1],'gpuArray');
jequaloneDist4(1,ceil(n_semiz/2),ceil(n_z/2),1)=1; % no assets, midpoint shock
AgentDist_initial4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,PolicyPath4(:,:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptionsC);
AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial4, jequaloneDist4, PricePath, ParamPath, PolicyPath4, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptionsC);
VPath4=squeeze(VPath4);
PolicyPath4=squeeze(PolicyPath4);
AgentDistPath4=squeeze(AgentDistPath4);

fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(VPath1(:)-VPath4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(PolicyPath1(:)-PolicyPath4(:))))
fprintf('Cross test: z and e 2, this should be zero: %.3e \n',max(abs(AgentDistPath1(:)-AgentDistPath4(:))))


%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
