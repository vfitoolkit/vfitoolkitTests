function output=CoreFHorzTPathPType_GEptype_Bequests(T,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i,Names_i,figure_c)
% Test 7: 'general eqm by ptype' (GEptype): a 'bequests left equal bequests
% received' general eqm condition that holds conditional on permanent type,
% with two types that differ by a fixed effect fe in their earnings.
%
% Bequests stay within type, so the types fully decouple and the PType+GEptype
% solve must equal two independent one-type GE solves: checked for the initial
% and final stationary eqms (HeteroAgentStationaryEqm_Case1_FHorz_PType with
% heteroagentoptions.GEptype) and for the GE transition path
% (TransitionPath_Case1_FHorz_PType with transpathoptions.GEptype). Also
% checks directly that bequests left equal bequests received (per type) at
% the stationary eqms and along the converged path.
%
% The transition is run twice, written out in full: first with fastOLG=1,
% then with fastOLG=0.
%
% Only z is used (no d, no e, no semiz). All checks here are 'close to zero'
% (solver tolerances), not machine zero.

n_d=0;
d_grid=[];

%% Model pieces specific to this test
% Fixed effect: entry ii is the fe of PType ii
Params.fe=[0.8;1.2];

% Conditional survival probabilities (sure death at the end of age N_j)
Params.sj=[0.95*ones(1,N_j-1),0];
DiscountFactorParamNames={'beta','sj'};

% Age weights consistent with the survival probabilities
Params.mewj=ones(1,N_j);
for jj=2:N_j
    Params.mewj(jj)=Params.mewj(jj-1)*Params.sj(jj-1);
end
Params.mewj=Params.mewj/sum(Params.mewj);

% Everyone alive receives the lump-sum bequest transfer beq. The GE price beq
% is a length-N_i vector, which is how the toolkit knows it depends on ptype.
Params.beq=[0.03;0.05]; % initial guess

ReturnFn=@(aprime,a,z,r,w,fe,kappa_j,sigma,agej,Jr,pension,beq) ...
    ReturnFn_nod_z_noe_nosemiz_febeq(aprime,a,z,r,w,fe,kappa_j,sigma,agej,Jr,pension,beq);

FnsToEvaluate.bequestsleft=@(aprime,a,z,sj) (1-sj)*aprime;

% Bequests left equal bequests received, conditional on ptype
GeneralEqmEqns.BeqMarket=@(beq,bequestsleft) beq-bequestsleft;
GEPriceParamNames={'beq'};

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Transition: w falls from 1.1 to 1 (final equals the baseline Params.w=1)
w_initial=1.1;
ParamPathLocal.w=[linspace(w_initial,1,floor(T/2)),ones(1,T-floor(T/2))];

%% Stationary GE at the initial and final parameters: PType with GEptype
heteroagentoptions.GEptype={'BeqMarket'};
n_p=0;

% Initial stationary eqm (w=1.1)
Params.w=w_initial;
[p_eqm_init_PT,~]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
Params.beq=[p_eqm_init_PT.beq.(Names_i{1}); p_eqm_init_PT.beq.(Names_i{2})];
[~,Policy_init_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

% Check bequests left equal bequests received at the initial stationary eqm
AllStats_init_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(AgentDist_initial_PT,Policy_init_PT,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype bequests, initial stationary eqm, beq vs bequests left (type %s), should be close to zero: %2.8f \n',nm,abs(p_eqm_init_PT.beq.(nm)-AllStats_init_PT.bequestsleft.(nm).Mean))
end

% Final stationary eqm (w=1)
Params.w=1;
Params.beq=[0.03;0.05]; % reset the initial guess
[p_eqm_final_PT,~]=HeteroAgentStationaryEqm_Case1_FHorz_PType(n_d, n_a, n_z, N_j, Names_i, n_p, pi_z, d_grid, a_grid, z_grid, jequaloneDist, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, GEPriceParamNames, heteroagentoptions, simoptions, vfoptions);
Params.beq=[p_eqm_final_PT.beq.(Names_i{1}); p_eqm_final_PT.beq.(Names_i{2})];
[V_final_PT,Policy_final_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_final_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_final_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);

% Check bequests left equal bequests received at the final stationary eqm
AllStats_final_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_final_PT,Policy_final_PT,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype bequests, final stationary eqm, beq vs bequests left (type %s), should be close to zero: %2.8f \n',nm,abs(p_eqm_final_PT.beq.(nm)-AllStats_final_PT.bequestsleft.(nm).Mean))
end

% Sanity: the two types must have genuinely different beq (else the test could
% silently pass with the types accidentally identical)
fprintf('GEptype bequests, final stationary eqm, beq by type (these should differ, higher fe higher beq): %2.8f vs %2.8f \n',p_eqm_final_PT.beq.(Names_i{1}),p_eqm_final_PT.beq.(Names_i{2}))

%% Stationary GE: two independent one-type solves (no PType, scalar beq)
heteroagentoptions_solo=struct();
p_eqm_init_solo=zeros(N_i,1); p_eqm_final_solo=zeros(N_i,1);
V_final_solo=cell(N_i,1); AgentDist_initial_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.fe=Params.fe(ii);

    % Initial stationary eqm (w=1.1)
    Params_ii.w=w_initial;
    Params_ii.beq=0.04; % initial guess
    [p_eqm_ii,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_ii, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions_solo, simoptions, vfoptions);
    p_eqm_init_solo(ii)=p_eqm_ii.beq;
    Params_ii.beq=p_eqm_ii.beq;
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_solo{ii}=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);

    % Final stationary eqm (w=1)
    Params_ii.w=1;
    Params_ii.beq=0.04; % reset the initial guess
    [p_eqm_ii,~]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames, n_d, n_a, n_z, N_j, n_p, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_ii, DiscountFactorParamNames, [], [], [], GEPriceParamNames, heteroagentoptions_solo, simoptions, vfoptions);
    p_eqm_final_solo(ii)=p_eqm_ii.beq;
    Params_ii.beq=p_eqm_ii.beq;
    [V_final_solo{ii},~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
end

for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype bequests, initial stationary eqm, PType vs solo beq (type %s), should be close to zero: %2.8f \n',nm,abs(p_eqm_init_PT.beq.(nm)-p_eqm_init_solo(ii)))
    fprintf('GEptype bequests, final   stationary eqm, PType vs solo beq (type %s), should be close to zero: %2.8f \n',nm,abs(p_eqm_final_PT.beq.(nm)-p_eqm_final_solo(ii)))
end

%% Initial guess for the price path (shared by both fastOLG blocks)
PricePath0=struct();
PricePath0.beq.(Names_i{1})=linspace(p_eqm_init_PT.beq.(Names_i{1}),p_eqm_final_PT.beq.(Names_i{1}),T)';
PricePath0.beq.(Names_i{2})=linspace(p_eqm_init_PT.beq.(Names_i{2}),p_eqm_final_PT.beq.(Names_i{2}),T)';

%% GE transition path with GEptype, fastOLG=1
transpathoptions_fast.GEnewprice=3;
transpathoptions_fast.GEnewprice3.howtoupdate={'BeqMarket','beq',0,0.5}; % new beq = old beq - 0.5*(beq-bequestsleft)
transpathoptions_fast.GEptype={'BeqMarket'};
transpathoptions_fast.fastOLG=1;

PricePath_PT_fast=TransitionPath_Case1_FHorz_PType(PricePath0, ParamPathLocal, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions_fast, simoptions, vfoptions);

% Verify bequests left equal bequests received along the converged path
[~,PolicyPath_fast]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath_PT_fast, ParamPathLocal, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_fast, vfoptions);
AgentDistPath_fast=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath_PT_fast, ParamPathLocal, PolicyPath_fast, AgeWeightParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params, transpathoptions_fast, simoptions);
AggVarsPath_fast=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_fast, PolicyPath_fast, PricePath_PT_fast, ParamPathLocal, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_fast, simoptions);

for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype bequests (fastOLG=1), beq path vs bequests left path (type %s), should be close to zero: %2.8f \n',nm,max(abs(PricePath_PT_fast.beq.(nm)(:)-AggVarsPath_fast.bequestsleft.(nm).Mean(:))))
end

% Two independent one-type GE transitions (no PType, scalar beq)
transpathoptions_fast_solo.GEnewprice=3;
transpathoptions_fast_solo.GEnewprice3.howtoupdate={'BeqMarket','beq',0,0.5};
transpathoptions_fast_solo.fastOLG=1;

PricePath_solo_fast=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.fe=Params.fe(ii);
    Params_ii.beq=p_eqm_final_solo(ii);
    PricePath0_ii=struct();
    PricePath0_ii.beq=linspace(p_eqm_init_solo(ii),p_eqm_final_solo(ii),T)';
    PricePath_solo_fast{ii}=TransitionPath_Case1_FHorz(PricePath0_ii, ParamPathLocal, T, V_final_solo{ii}, AgentDist_initial_solo{ii}, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_ii, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions_fast_solo, simoptions, vfoptions);

    fprintf('GEptype bequests (fastOLG=1), PType vs solo beq path (type %s), should be close to zero: %2.8f \n',Names_i{ii},max(abs(PricePath_PT_fast.beq.(Names_i{ii})(:)-PricePath_solo_fast{ii}.beq(:))))
end

%% GE transition path with GEptype, fastOLG=0
transpathoptions_slow.GEnewprice=3;
transpathoptions_slow.GEnewprice3.howtoupdate={'BeqMarket','beq',0,0.5}; % new beq = old beq - 0.5*(beq-bequestsleft)
transpathoptions_slow.GEptype={'BeqMarket'};
transpathoptions_slow.fastOLG=0;

PricePath_PT_slow=TransitionPath_Case1_FHorz_PType(PricePath0, ParamPathLocal, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions_slow, simoptions, vfoptions);

% Verify bequests left equal bequests received along the converged path
[~,PolicyPath_slow]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath_PT_slow, ParamPathLocal, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_slow, vfoptions);
AgentDistPath_slow=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath_PT_slow, ParamPathLocal, PolicyPath_slow, AgeWeightParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params, transpathoptions_slow, simoptions);
AggVarsPath_slow=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_slow, PolicyPath_slow, PricePath_PT_slow, ParamPathLocal, Params, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_slow, simoptions);

for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype bequests (fastOLG=0), beq path vs bequests left path (type %s), should be close to zero: %2.8f \n',nm,max(abs(PricePath_PT_slow.beq.(nm)(:)-AggVarsPath_slow.bequestsleft.(nm).Mean(:))))
end

% Two independent one-type GE transitions (no PType, scalar beq)
transpathoptions_slow_solo.GEnewprice=3;
transpathoptions_slow_solo.GEnewprice3.howtoupdate={'BeqMarket','beq',0,0.5};
transpathoptions_slow_solo.fastOLG=0;

PricePath_solo_slow=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params;
    Params_ii.fe=Params.fe(ii);
    Params_ii.beq=p_eqm_final_solo(ii);
    PricePath0_ii=struct();
    PricePath0_ii.beq=linspace(p_eqm_init_solo(ii),p_eqm_final_solo(ii),T)';
    PricePath_solo_slow{ii}=TransitionPath_Case1_FHorz(PricePath0_ii, ParamPathLocal, T, V_final_solo{ii}, AgentDist_initial_solo{ii}, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_ii, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions_slow_solo, simoptions, vfoptions);

    fprintf('GEptype bequests (fastOLG=0), PType vs solo beq path (type %s), should be close to zero: %2.8f \n',Names_i{ii},max(abs(PricePath_PT_slow.beq.(Names_i{ii})(:)-PricePath_solo_slow{ii}.beq(:))))
end

%% fastOLG vs slowOLG
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('GEptype bequests, fastOLG vs slowOLG beq path (type %s), should be close to zero: %2.8f \n',nm,max(abs(PricePath_PT_fast.beq.(nm)(:)-PricePath_PT_slow.beq.(nm)(:))))
end

%% Graph the bequest paths
fig=figure(figure_c);
subplot(2,1,1); plot(1:1:T,PricePath_PT_fast.beq.(Names_i{1})(:), 1:1:T,PricePath_solo_fast{1}.beq(:), 1:1:T,PricePath_PT_slow.beq.(Names_i{1})(:))
title(['beq path, type ',Names_i{1}])
legend('PType fastOLG','solo fastOLG','PType slowOLG')
subplot(2,1,2); plot(1:1:T,PricePath_PT_fast.beq.(Names_i{2})(:), 1:1:T,PricePath_solo_fast{2}.beq(:), 1:1:T,PricePath_PT_slow.beq.(Names_i{2})(:))
title(['beq path, type ',Names_i{2}])
legend('PType fastOLG','solo fastOLG','PType slowOLG')

output=struct();

end
