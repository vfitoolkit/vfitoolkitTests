function output=CoreFHorzTPathPType_GETPath_1iter(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test 4: Run the GE transition path (TransitionPath_Case1_FHorz_PType) with
% transpathoptions.maxiter=2, so it returns the path after one update. With N_i=2
% identical types (and ptypeweights summing to one) the aggregates equal the
% no-PType aggregates, so the one-iteration price path update must equal the
% no-PType TransitionPath_Case1_FHorz one. Checked with fastOLG=1 and
% fastOLG=0. Only z is used (no d, no e, no semiz).

n_d=0;
d_grid=[];

% ReturnFn declares kappa_j (not kappa_j_pt), so both types are identical.
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

% Trivial period-T V and Policy (as in CoreFHorzTPathTests)
V_final=zeros([n_a,n_z,N_j],'gpuArray');
Policy_final=ones([1,n_a,n_z,N_j],'gpuArray');
names_PT={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_PT=struct(); Policy_final_PT=struct();
for ii=1:N_i
    V_final_PT.(names_PT{ii})=V_final;
    Policy_final_PT.(names_PT{ii})=Policy_final;
end

%% Create the initial agent distributions (from the first period of the policy path)
transpathoptions.fastOLG=1;
[~,PolicyPath]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,PolicyPath(:,:,:,:,1),n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
clear PolicyPath

[~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params, n_d, n_a, n_z, N_j, N_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
Policy_init_PT=struct();
for ii=1:N_i
    Policy_init_PT.(names_PT{ii})=PolicyPath_PT.(names_PT{ii})(:,:,:,:,1);
end
AgentDist_initial_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_init_PT,n_d,n_a,n_z,N_j,N_i,pi_z,Params,simoptions);
clear PolicyPath_PT

%% One iteration of the GE transition path, with fastOLG=1
% maxiter=2: since the shooting algorithm returns the path whose general eqm conditions it evaluated,
% maxiter=1 would return the starting path unchanged and compare it with itself. The general eqm
% condition depends on earnings, which depend on w on the ParamPath, so the update does too.
transpathoptions.maxiter=2;
transpathoptions.GEnewprice=3; % the shooting algorithm (GEnewprice has no default: 1 is quasi-Newton, 2 Anderson, 3 shooting)
transpathoptions.GEnewprice3.howtoupdate={'dummy','r',0,0.1};

GeneralEqmEqns.dummy=@(r,earnings) r-0.02-0.01*earnings;

transpathoptions.fastOLG=1;
PricePath_PT1=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
PricePath_NoPT1=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of GE TPath (fastOLG=1), PType vs no-PType, this should be zero: %.3e \n',max(abs(PricePath_PT1.r(:)-PricePath_NoPT1.r(:))))

%% And again with fastOLG=0
transpathoptions.fastOLG=0;
PricePath_PT0=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
PricePath_NoPT0=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

fprintf('One iter of GE TPath (fastOLG=0), PType vs no-PType, this should be zero: %.3e \n',max(abs(PricePath_PT0.r(:)-PricePath_NoPT0.r(:))))

fprintf('One iter of GE TPath, PType with/without fastOLG, this should be zero: %.3e \n',max(abs(PricePath_PT1.r(:)-PricePath_PT0.r(:))))

%% The same with Anderson acceleration (GEnewprice=2)
% warmup=0 and maxiter=3, so iterations 2 and 3 are Anderson steps (at maxiter=2 Anderson would still be
% the shooting algorithm here, as iteration 1 has no history to mix)
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
transpathoptionsA.fastOLG=1;
PricePathA_PT1=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
PricePathA_NoPT1=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptionsA.fastOLG=0;
PricePathA_PT0=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
% The shooting algorithm for the same number of iterations, to show the Anderson steps were taken
transpathoptionsS3=transpathoptions;
transpathoptionsS3.maxiter=3;
transpathoptionsS3.fastOLG=1;
PricePathS3_PT1=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsS3, simoptions, vfoptions);

fprintf('Anderson at maxiter=3 GE TPath (fastOLG=1), PType vs no-PType, this should be zero: %.3e \n',max(abs(PricePathA_PT1.r(:)-PricePathA_NoPT1.r(:))))
fprintf('Anderson at maxiter=3 GE TPath, PType with/without fastOLG, this should be zero: %.3e \n',max(abs(PricePathA_PT1.r(:)-PricePathA_PT0.r(:))))
fprintf('Anderson vs shooting at maxiter=3 GE TPath (fastOLG=1), PType, this should NOT be zero: %.3e \n',max(abs(PricePathA_PT1.r(:)-PricePathS3_PT1.r(:))))

output=struct();

end
