function output=QHDFHorzTPath_CrossTests_d_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)

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

ReturnFn_none=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_z=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_e=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_ze=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

% V_final / Policy_final stubs (zeros — cross-test equality holds for any matching boundary)
V_final_none=zeros([n_a,N_j],'gpuArray');
Policy_final_none=zeros([2,n_a,N_j],'gpuArray');
V_final_none2=zeros([n_a,1,N_j],'gpuArray');
Policy_final_none2=zeros([2,n_a,1,N_j],'gpuArray');
V_final_z=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final_z=zeros([2,n_a,n_z,N_j],'gpuArray');
V_final_e=zeros([n_a,vfoptions.n_e,N_j],'gpuArray');
Policy_final_e=zeros([2,n_a,vfoptions.n_e,N_j],'gpuArray');
V_final_ze1=zeros([n_a,1,vfoptions.n_e,N_j],'gpuArray');
Policy_final_ze1=zeros([2,n_a,1,vfoptions.n_e,N_j],'gpuArray');
V_final_ze2=zeros([n_a,n_z,1,N_j],'gpuArray');
Policy_final_ze2=zeros([2,n_a,n_z,1,N_j],'gpuArray');

for qhregime={'Naive','Sophisticated'}
    qh=qhregime{1};
    fprintf('\n--- QH-TPath CrossTests (d, nosemiz): %s ---\n',qh)

    %% Solving with just a single point for z with value 1 and prob 1 gives us same as no shocks
    jequaloneDist_none=zeros(n_a,1,'gpuArray');
    jequaloneDist_none(1)=1; % no assets

    vfoptions_z=struct();
    simoptions_z=struct();
    vfoptions_z.exoticpreferences='QuasiHyperbolic';
    vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptions_z.quasi_hyperbolic=qh;
    [VPath0,PolicyPath0]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_none, Policy_final_none, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[], [], DiscountFactorParamNames, ReturnFn_none, transpathoptionsbaseline, vfoptions_z);
    AgentDist_initial0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,PolicyPath0(:,:,:,1),n_d,n_a,0,N_j,[],Params,simoptions_z);
    AgentDistPath0=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial0, jequaloneDist_none, PricePath, ParamPath, PolicyPath0, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptions_z);

    vfoptions_z=struct();
    simoptions_z=struct();
    vfoptions_z.exoticpreferences='QuasiHyperbolic';
    vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptions_z.quasi_hyperbolic=qh;
    [VPath0z,PolicyPath0z]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_none2, Policy_final_none2, Params, n_d, n_a, 1, N_j, d_grid, a_grid,1,1, DiscountFactorParamNames, ReturnFn_z, transpathoptionsbaseline, vfoptions_z);
    AgentDist_initial0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,PolicyPath0z(:,:,:,:,1),n_d,n_a,1,N_j,1,Params,simoptions_z);
    AgentDistPath0z=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial0z, jequaloneDist_none, PricePath, ParamPath, PolicyPath0z, AgeWeightParamNames,n_d,n_a,1,N_j,1, T,Params, transpathoptionsbaseline, simoptions_z);

    fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(VPath0(:)-VPath0z(:))))
    fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(PolicyPath0(:)-PolicyPath0z(:))))
    fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(AgentDistPath0(:)-AgentDistPath0z(:))))

    %% Solve using a markov which is just an iid in disguise. Should give same result as the iid
    jequaloneDist_z=zeros(n_a,n_z,'gpuArray');
    jequaloneDist_z(1,ceil(n_z/2))=1; % no assets, midpoint shock

    vfoptions_z=struct();
    simoptions_z=struct();
    vfoptions_z.exoticpreferences='QuasiHyperbolic';
    vfoptions_z.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptions_z.quasi_hyperbolic=qh;
    [VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_z, Policy_final_z, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, ReturnFn_z, transpathoptionsbaseline, vfoptions_z);
    AgentDist_initial1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,PolicyPath1(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions_z);
    AgentDistPath1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial1, jequaloneDist_z, PricePath, ParamPath, PolicyPath1, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions_z);

    vfoptions_e=vfoptions;
    simoptions_e=simoptions;
    vfoptions_e.exoticpreferences='QuasiHyperbolic';
    vfoptions_e.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptions_e.quasi_hyperbolic=qh;
    [VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_e, Policy_final_e, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[],[], DiscountFactorParamNames, ReturnFn_e, transpathoptionsbaseline, vfoptions_e);
    AgentDist_initial2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,PolicyPath2(:,:,:,:,1),n_d,n_a,0,N_j,[],Params,simoptions_e);
    AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial2, jequaloneDist_z, PricePath, ParamPath, PolicyPath2, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptions_e);

    fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath2(:))))
    fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath2(:))))
    fprintf('Cross test: z as e, this should be zero: %2.8f \n',max(abs(AgentDistPath1(:)-AgentDistPath2(:))))

    %% Use code with z and e, but set the 'other' to a single point with value 1 and prob 1
    % First, make z just 1
    vfoptions_ze1=vfoptions;
    simoptions_ze1=simoptions;
    vfoptions_ze1.exoticpreferences='QuasiHyperbolic';
    vfoptions_ze1.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptions_ze1.quasi_hyperbolic=qh;
    [VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_ze1, Policy_final_ze1, Params, n_d, n_a, 1, N_j, d_grid, a_grid,1,1, DiscountFactorParamNames, ReturnFn_ze, transpathoptionsbaseline, vfoptions_ze1);
    jequaloneDist3=zeros(n_a,1,vfoptions_ze1.n_e,'gpuArray');
    jequaloneDist3(1,1,ceil(vfoptions_ze1.n_e/2))=1; % no assets, midpoint shock
    AgentDist_initial3=StationaryDist_FHorz_Case1(jequaloneDist3,AgeWeightParamNames,PolicyPath3(:,:,:,:,:,1),n_d,n_a,1,N_j,1,Params,simoptions_ze1);
    AgentDistPath3=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial3, jequaloneDist3, PricePath, ParamPath, PolicyPath3, AgeWeightParamNames,n_d,n_a,1,N_j,1, T,Params, transpathoptionsbaseline, simoptions_ze1);
    VPath3=squeeze(VPath3);
    PolicyPath3=squeeze(PolicyPath3);
    AgentDistPath3=squeeze(AgentDistPath3);

    fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath3(:))))
    fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath3(:))))
    fprintf('Cross test: z and e 1, this should be zero: %2.8f \n',max(abs(AgentDistPath1(:)-AgentDistPath3(:))))

    % Second, make e just 1
    vfoptions_ze2.n_e=1;
    vfoptions_ze2.e_grid=1;
    vfoptions_ze2.pi_e=1;
    vfoptions_ze2.exoticpreferences='QuasiHyperbolic';
    vfoptions_ze2.QHadditionaldiscount=vfoptionsbaseline.QHadditionaldiscount;
    vfoptions_ze2.quasi_hyperbolic=qh;
    simoptions_ze2=vfoptions_ze2;
    simoptions_ze2=rmfield(simoptions_ze2,{'exoticpreferences','QHadditionaldiscount','quasi_hyperbolic'});
    [VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final_ze2, Policy_final_ze2, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, ReturnFn_ze, transpathoptionsbaseline, vfoptions_ze2);
    jequaloneDist4=zeros(n_a,n_z,1,'gpuArray');
    jequaloneDist4(1,ceil(n_z/2),1)=1; % no assets, midpoint shock
    AgentDist_initial4=StationaryDist_FHorz_Case1(jequaloneDist4,AgeWeightParamNames,PolicyPath4(:,:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions_ze2);
    AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial4, jequaloneDist4, PricePath, ParamPath, PolicyPath4, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions_ze2);
    VPath4=squeeze(VPath4);
    PolicyPath4=squeeze(PolicyPath4);
    AgentDistPath4=squeeze(AgentDistPath4);

    fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(VPath1(:)-VPath4(:))))
    fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(PolicyPath1(:)-PolicyPath4(:))))
    fprintf('Cross test: z and e 2, this should be zero: %2.8f \n',max(abs(AgentDistPath1(:)-AgentDistPath4(:))))

    clear vfoptions_ze2
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
