function output=EZFHorzTPath_CrossTests_d_nosemiz(T,PricePath,ParamPath,n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline,simoptionsbaseline)
% Epstein-Zin mirror of QHDFHorzTPath_CrossTests_d_nosemiz: the same model written with different shock
% set-ups must give the same value fn, policy and agent distribution along the (moving) transition path,
% for each of the three Epstein-Zin cases:
%   (1) no shocks against a single-point z (value 1, probability 1)
%   (2) a markov z that is an iid in disguise against the same shock as e
%   (3) z and e with z a single point, against e alone (here: against (2)'s z)
%   (4) z and e with e a single point, against z alone
%   (5) markov z plus iid e against the same two shocks as one joint markov, n_z=[n_z,n_e] and
%       pi_z=kron(ones(n_e,1)*pi_e',pi_z). Specific to Epstein-Zin: for vNM it is trivially true, for EZ
%       it is the sharpest test that the codes take ONE certainty equivalent over (zprime,eprime) jointly,
%       not a CE over e nested inside a CE over z (equal only at the gamma=1/phi knife-edge). On the path
%       each period's certainty equivalent is taken by a transition-path raw, which (5) checks against its
%       own joint-markov twin.
% Under Epstein-Zin a terminal value fn of zeros is not 'nothing after the path' but a catastrophe (a
% continuation of exactly zero is a value, see EZ_feasibility_conventions_proposal.md), so unlike the
% QH cross tests every configuration starts from its own stationary Epstein-Zin solve.
% (1) has no shocks, so the Epstein-Zin no-shocks warning is expected output there.

% with d: n_d and d_grid as passed in

% For the crosstests, z is set up as a copy of e (z and e appear in the same place in earnings)
n_zc=vfoptionsbaseline.n_e;
pi_zc=repmat(vfoptionsbaseline.pi_e',vfoptionsbaseline.n_e,1);
z_gridc=vfoptionsbaseline.e_grid;
n_e=vfoptionsbaseline.n_e;

ReturnFn_cons={@(d,aprime,a,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,varphi,agej,Jr,pension), ...
    @(d,aprime,a,z,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,varphi,agej,Jr,pension), ...
    @(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,varphi,agej,Jr,pension), ...
    @(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension) EZReturnFn_cons_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,varphi,agej,Jr,pension)};
ReturnFn_posU={@(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_positiveUtils_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension), ...
    @(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_positiveUtils_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension), ...
    @(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_positiveUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension), ...
    @(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_positiveUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)};
ReturnFn_negU={@(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension), ...
    @(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension), ...
    @(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension), ...
    @(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension) EZReturnFn_negativeUtils_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,ezsigma,varphi,agej,Jr,pension)};
% index into those: 1 none, 2 z, 3 e, 4 z and e

vfoptions_e=struct();
vfoptions_e.n_e=n_e; vfoptions_e.e_grid=vfoptionsbaseline.e_grid; vfoptions_e.pi_e=vfoptionsbaseline.pi_e;
simoptions_e=struct();
simoptions_e.n_e=n_e; simoptions_e.e_grid=simoptionsbaseline.e_grid; simoptions_e.pi_e=simoptionsbaseline.pi_e;
vfoptions_e1=struct(); vfoptions_e1.n_e=1; vfoptions_e1.e_grid=1; vfoptions_e1.pi_e=1;
simoptions_e1=vfoptions_e1;

jequaloneDist_none=zeros(n_a,1,'gpuArray'); jequaloneDist_none(1)=1; % no assets
jequaloneDist_z=zeros(n_a,n_zc,'gpuArray'); jequaloneDist_z(1,ceil(n_zc/2))=1; % no assets, midpoint shock
jequaloneDist_ze1=zeros(n_a,1,n_e,'gpuArray'); jequaloneDist_ze1(1,1,ceil(n_e/2))=1;
jequaloneDist_ze2=zeros(n_a,n_zc,1,'gpuArray'); jequaloneDist_ze2(1,ceil(n_zc/2),1)=1;

% (5): z from the setup plus e, against the two as one joint markov (z first, e second)
n_zB=[n_z,n_e];
z_gridB=[z_grid; vfoptionsbaseline.e_grid];
pi_zB=kron(ones(n_e,1)*vfoptionsbaseline.pi_e',pi_z);
jequaloneDist_5=zeros(n_a,n_z,n_e,'gpuArray'); jequaloneDist_5(1,ceil(n_z/2),ceil(n_e/2))=1;

casename={'cons units','positive utils','negative utils'};
for cc=1:3
    vfEZ=struct();
    vfEZ.exoticpreferences='EpsteinZin';
    if cc==1
        vfEZ.EZutils=0; vfEZ.EZriskaversion='ezgamma'; vfEZ.EZeis='ezphi';
        RF=ReturnFn_cons;
    elseif cc==2
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=1; vfEZ.EZriskaversion='ezrisk';
        RF=ReturnFn_posU;
    else
        vfEZ.EZutils=1; vfEZ.EZpositiveutility=0; vfEZ.EZriskaversion='ezrisk';
        RF=ReturnFn_negU;
    end
    fprintf('\n--- EZ-TPath CrossTests (d, nosemiz): %s ---\n',casename{cc})
    vfEZ_e=vfEZ; vfEZ_e.n_e=vfoptions_e.n_e; vfEZ_e.e_grid=vfoptions_e.e_grid; vfEZ_e.pi_e=vfoptions_e.pi_e;
    vfEZ_e1=vfEZ; vfEZ_e1.n_e=1; vfEZ_e1.e_grid=1; vfEZ_e1.pi_e=1;

    %% (1) no shocks against a single-point z
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],RF{1},Params,DiscountFactorParamNames,[],vfEZ);
    [VPath0,PolicyPath0]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[], [], DiscountFactorParamNames, RF{1}, transpathoptionsbaseline, vfEZ);
    AgentDist_initial0=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy_s,n_d,n_a,0,N_j,[],Params,struct());
    AgentDistPath0=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial0, jequaloneDist_none, PricePath, ParamPath, PolicyPath0, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, struct());

    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,RF{2},Params,DiscountFactorParamNames,[],vfEZ);
    [VPath0z,PolicyPath0z]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, 1, N_j, d_grid, a_grid,1,1, DiscountFactorParamNames, RF{2}, transpathoptionsbaseline, vfEZ);
    AgentDist_initial0z=StationaryDist_FHorz_Case1(jequaloneDist_none,AgeWeightParamNames,Policy_s,n_d,n_a,1,N_j,1,Params,struct());
    AgentDistPath0z=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial0z, jequaloneDist_none, PricePath, ParamPath, PolicyPath0z, AgeWeightParamNames,n_d,n_a,1,N_j,1, T,Params, transpathoptionsbaseline, struct());

    fprintf('Cross test (1) %s: no shocks against single-point z, this should be zero, V: %.3e \n',casename{cc},max(abs(VPath0(:)-VPath0z(:))))
    fprintf('Cross test (1) %s: no shocks against single-point z, this should be zero, Policy: %.3e \n',casename{cc},max(abs(PolicyPath0(:)-PolicyPath0z(:))))
    fprintf('Cross test (1) %s: no shocks against single-point z, this should be zero, AgentDist: %.3e \n',casename{cc},max(abs(AgentDistPath0(:)-AgentDistPath0z(:))))
    clear VPath0 PolicyPath0 AgentDistPath0 VPath0z PolicyPath0z AgentDistPath0z

    %% (2) a markov z that is an iid in disguise, against the same shock as e
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,RF{2},Params,DiscountFactorParamNames,[],vfEZ);
    [VPath1,PolicyPath1]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, n_zc, N_j, d_grid, a_grid,z_gridc,pi_zc, DiscountFactorParamNames, RF{2}, transpathoptionsbaseline, vfEZ);
    AgentDist_initial1=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy_s,n_d,n_a,n_zc,N_j,pi_zc,Params,struct());
    AgentDistPath1=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial1, jequaloneDist_z, PricePath, ParamPath, PolicyPath1, AgeWeightParamNames,n_d,n_a,n_zc,N_j,pi_zc, T,Params, transpathoptionsbaseline, struct());

    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,0,N_j,d_grid,a_grid,[],[],RF{3},Params,DiscountFactorParamNames,[],vfEZ_e);
    [VPath2,PolicyPath2]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, 0, N_j, d_grid, a_grid,[],[], DiscountFactorParamNames, RF{3}, transpathoptionsbaseline, vfEZ_e);
    AgentDist_initial2=StationaryDist_FHorz_Case1(jequaloneDist_z,AgeWeightParamNames,Policy_s,n_d,n_a,0,N_j,[],Params,simoptions_e);
    AgentDistPath2=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial2, jequaloneDist_z, PricePath, ParamPath, PolicyPath2, AgeWeightParamNames,n_d,n_a,0,N_j,[], T,Params, transpathoptionsbaseline, simoptions_e);

    fprintf('Cross test (2) %s: iid markov z against e, this should be zero, V: %.3e \n',casename{cc},max(abs(VPath1(:)-VPath2(:))))
    fprintf('Cross test (2) %s: iid markov z against e, this should be zero, Policy: %.3e \n',casename{cc},max(abs(PolicyPath1(:)-PolicyPath2(:))))
    fprintf('Cross test (2) %s: iid markov z against e, this should be zero, AgentDist: %.3e \n',casename{cc},max(abs(AgentDistPath1(:)-AgentDistPath2(:))))
    clear VPath2 PolicyPath2 AgentDistPath2

    %% (3) z and e with z a single point, against (2)'s z
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,1,N_j,d_grid,a_grid,1,1,RF{4},Params,DiscountFactorParamNames,[],vfEZ_e);
    [VPath3,PolicyPath3]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, 1, N_j, d_grid, a_grid,1,1, DiscountFactorParamNames, RF{4}, transpathoptionsbaseline, vfEZ_e);
    AgentDist_initial3=StationaryDist_FHorz_Case1(jequaloneDist_ze1,AgeWeightParamNames,Policy_s,n_d,n_a,1,N_j,1,Params,simoptions_e);
    AgentDistPath3=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial3, jequaloneDist_ze1, PricePath, ParamPath, PolicyPath3, AgeWeightParamNames,n_d,n_a,1,N_j,1, T,Params, transpathoptionsbaseline, simoptions_e);

    fprintf('Cross test (3) %s: z and e with z a single point, this should be zero, V: %.3e \n',casename{cc},max(abs(VPath1(:)-VPath3(:))))
    fprintf('Cross test (3) %s: z and e with z a single point, this should be zero, Policy: %.3e \n',casename{cc},max(abs(PolicyPath1(:)-PolicyPath3(:))))
    fprintf('Cross test (3) %s: z and e with z a single point, this should be zero, AgentDist: %.3e \n',casename{cc},max(abs(AgentDistPath1(:)-AgentDistPath3(:))))
    clear VPath3 PolicyPath3 AgentDistPath3

    %% (4) z and e with e a single point, against (2)'s z
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zc,N_j,d_grid,a_grid,z_gridc,pi_zc,RF{4},Params,DiscountFactorParamNames,[],vfEZ_e1);
    [VPath4,PolicyPath4]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, n_zc, N_j, d_grid, a_grid,z_gridc,pi_zc, DiscountFactorParamNames, RF{4}, transpathoptionsbaseline, vfEZ_e1);
    AgentDist_initial4=StationaryDist_FHorz_Case1(jequaloneDist_ze2,AgeWeightParamNames,Policy_s,n_d,n_a,n_zc,N_j,pi_zc,Params,simoptions_e1);
    AgentDistPath4=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial4, jequaloneDist_ze2, PricePath, ParamPath, PolicyPath4, AgeWeightParamNames,n_d,n_a,n_zc,N_j,pi_zc, T,Params, transpathoptionsbaseline, simoptions_e1);

    fprintf('Cross test (4) %s: z and e with e a single point, this should be zero, V: %.3e \n',casename{cc},max(abs(VPath1(:)-VPath4(:))))
    fprintf('Cross test (4) %s: z and e with e a single point, this should be zero, Policy: %.3e \n',casename{cc},max(abs(PolicyPath1(:)-PolicyPath4(:))))
    fprintf('Cross test (4) %s: z and e with e a single point, this should be zero, AgentDist: %.3e \n',casename{cc},max(abs(AgentDistPath1(:)-AgentDistPath4(:))))
    clear VPath1 PolicyPath1 AgentDistPath1 VPath4 PolicyPath4 AgentDistPath4

    %% (5) markov z plus iid e, against the two as one joint markov
    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,RF{4},Params,DiscountFactorParamNames,[],vfEZ_e);
    [VPathA,PolicyPathA]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid,pi_z, DiscountFactorParamNames, RF{4}, transpathoptionsbaseline, vfEZ_e);
    AgentDist_initialA=StationaryDist_FHorz_Case1(jequaloneDist_5,AgeWeightParamNames,Policy_s,n_d,n_a,n_z,N_j,pi_z,Params,simoptions_e);
    AgentDistPathA=AgentDistOnTransPath_Case1_FHorz(AgentDist_initialA, jequaloneDist_5, PricePath, ParamPath, PolicyPathA, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptionsbaseline, simoptions_e);

    [V_s,Policy_s]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zB,N_j,d_grid,a_grid,z_gridB,pi_zB,RF{4},Params,DiscountFactorParamNames,[],vfEZ);
    [VPathB,PolicyPathB]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_s, Policy_s, Params, n_d, n_a, n_zB, N_j, d_grid, a_grid,z_gridB,pi_zB, DiscountFactorParamNames, RF{4}, transpathoptionsbaseline, vfEZ);
    AgentDist_initialB=StationaryDist_FHorz_Case1(jequaloneDist_5,AgeWeightParamNames,Policy_s,n_d,n_a,n_zB,N_j,pi_zB,Params,struct());
    AgentDistPathB=AgentDistOnTransPath_Case1_FHorz(AgentDist_initialB, jequaloneDist_5, PricePath, ParamPath, PolicyPathB, AgeWeightParamNames,n_d,n_a,n_zB,N_j,pi_zB, T,Params, transpathoptionsbaseline, struct());

    fprintf('Cross test (5) %s: (z,e) against the joint markov, this should be zero, V: %.3e \n',casename{cc},max(abs(VPathA(:)-VPathB(:))))
    fprintf('Cross test (5) %s: (z,e) against the joint markov, this should be zero, Policy: %.3e \n',casename{cc},max(abs(PolicyPathA(:)-PolicyPathB(:))))
    fprintf('Cross test (5) %s: (z,e) against the joint markov, this should be zero, AgentDist: %.3e \n',casename{cc},max(abs(AgentDistPathA(:)-AgentDistPathB(:))))
    clear VPathA PolicyPathA AgentDistPathA VPathB PolicyPathB AgentDistPathB V_s Policy_s
end

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
