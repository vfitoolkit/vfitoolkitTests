function output=CoreInfHorzTPathPType_NivsNames(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,N_i,Names_i,transpathoptionsbaseline)
% TPath with N_i (numeric) vs Names_i (cell) should give identical output. The two types are genuinely
% different (sigma differs, given as a vector with one entry per ptype). Only z (no d, no e).
%   (1) the standalone commands: VPath, PolicyPath, AgentDistPath, AggVarsPath (by ptype and pooled)
%   (2) one general eqm solve (shooting, maxiter=2), common FnsToEvaluate
%   (3) the same with FnsToEvaluate.earnings given per ptype, as a structure keyed by the ptype names
%       (ptype001/ptype002 for N_i, the Names_i otherwise)
% V_final is the stationary value fn at the final prices and AgentDist_initial the stationary agent
% distribution at the period-1 prices, each type with its own sigma.

n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;

sigma_pt=[2;3];
Params_PT=Params;
Params_PT.sigma=sigma_pt; % one entry per ptype

names_A={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_A=struct(); Policy_final_A=struct(); AgentDist_initial_A=struct();
V_final_B=struct(); Policy_final_B=struct(); AgentDist_initial_B=struct();
for ii=1:N_i
    Params_ii=Params;
    Params_ii.sigma=sigma_pt(ii);
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_ii=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptions,Params_ii,[]);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_ii,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    V_final_A.(names_A{ii})=V_ii;
    Policy_final_A.(names_A{ii})=Policy_ii;
    AgentDist_initial_A.(names_A{ii})=AgentDist_ii;
    V_final_B.(Names_i{ii})=V_ii;
    Policy_final_B.(Names_i{ii})=Policy_ii;
    AgentDist_initial_B.(Names_i{ii})=AgentDist_ii;
end
AgentDist_initial_A.ptweights=Params.ptypeweights;
AgentDist_initial_B.ptweights=Params.ptypeweights;

%% (1) The standalone commands
[VPath_A,PolicyPath_A]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath, T, V_final_A, Policy_final_A, Params_PT, n_d,n_a,n_z, N_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
AgentDistPath_A=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_A, PricePath, ParamPath, PolicyPath_A, n_d,n_a,n_z, N_i, pi_z, T, Params_PT, transpathoptionsbaseline, simoptions);
AggVarsPath_A=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_A, PolicyPath_A, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, N_i, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptions);

[VPath_B,PolicyPath_B]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath, T, V_final_B, Policy_final_B, Params_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
AgentDistPath_B=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_B, PricePath, ParamPath, PolicyPath_B, n_d,n_a,n_z, Names_i, pi_z, T, Params_PT, transpathoptionsbaseline, simoptions);
AggVarsPath_B=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_B, PolicyPath_B, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, Names_i, d_grid, a_grid, z_grid, transpathoptionsbaseline, simoptions);

for ii=1:N_i
    nA=names_A{ii};
    nB=Names_i{ii};
    fprintf('N_i vs Names_i, VPath (type %d), this should be zero: %.3e \n',ii,max(abs(VPath_A.(nA)(:)-VPath_B.(nB)(:))))
    fprintf('N_i vs Names_i, PolicyPath (type %d), this should be zero: %.3e \n',ii,max(abs(PolicyPath_A.(nA)(:)-PolicyPath_B.(nB)(:))))
    fprintf('N_i vs Names_i, AgentDistPath (type %d), this should be zero: %.3e \n',ii,max(abs(AgentDistPath_A.(nA)(:)-AgentDistPath_B.(nB)(:))))
    fprintf('N_i vs Names_i, AggVarsPath assets (type %d), this should be zero: %.3e \n',ii,max(abs(AggVarsPath_A.assets.(nA).Mean(:)-AggVarsPath_B.assets.(nB).Mean(:))))
    fprintf('N_i vs Names_i, AggVarsPath earnings (type %d), this should be zero: %.3e \n',ii,max(abs(AggVarsPath_A.earnings.(nA).Mean(:)-AggVarsPath_B.earnings.(nB).Mean(:))))
end
fprintf('N_i vs Names_i, AggVarsPath assets pooled, this should be zero: %.3e \n',max(abs(AggVarsPath_A.assets.Mean(:)-AggVarsPath_B.assets.Mean(:))))
fprintf('N_i vs Names_i, AggVarsPath earnings pooled, this should be zero: %.3e \n',max(abs(AggVarsPath_A.earnings.Mean(:)-AggVarsPath_B.earnings.Mean(:))))
fprintf('N_i vs Names_i, for scale the two types have different VPath, this should NOT be zero: %.3e \n',max(abs(VPath_B.(Names_i{1})(:)-VPath_B.(Names_i{2})(:))))

%% (2) One general eqm solve, common FnsToEvaluate
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;

[PricePath_GA,GEcondnPath_GA]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_A, AgentDist_initial_A, n_d,n_a,n_z, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
[PricePath_GB,GEcondnPath_GB]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_B, AgentDist_initial_B, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
fprintf('N_i vs Names_i, general eqm, PricePath r, this should be zero: %.3e \n',max(abs(PricePath_GA.r(:)-PricePath_GB.r(:))))
fprintf('N_i vs Names_i, general eqm, GEcondnPath, this should be zero: %.3e \n',max(abs(GEcondnPath_GA.CapitalMarket(:)-GEcondnPath_GB.CapitalMarket(:))))
fprintf('N_i vs Names_i, general eqm, for scale the update moved r by %.3e (must not be zero, or the PricePath check above is empty) \n',max(abs(PricePath_GB.r(:)-PricePath.r(:))))

%% (3) FnsToEvaluate given per ptype (the same function for both, so the answer must not change)
FnsToEvaluate_A.assets=FnsToEvaluate.assets;
FnsToEvaluate_B.assets=FnsToEvaluate.assets;
for ii=1:N_i
    FnsToEvaluate_A.earnings.(names_A{ii})=FnsToEvaluate.earnings;
    FnsToEvaluate_B.earnings.(Names_i{ii})=FnsToEvaluate.earnings;
end
[PricePath_FA,GEcondnPath_FA]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_A, AgentDist_initial_A, n_d,n_a,n_z, N_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_A, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
fprintf('N_i, FnsToEvaluate per ptype vs common, general eqm, PricePath r, this should be zero: %.3e \n',max(abs(PricePath_FA.r(:)-PricePath_GA.r(:))))
fprintf('N_i, FnsToEvaluate per ptype vs common, general eqm, GEcondnPath, this should be zero: %.3e \n',max(abs(GEcondnPath_FA.CapitalMarket(:)-GEcondnPath_GA.CapitalMarket(:))))
[PricePath_FB,GEcondnPath_FB]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_B, AgentDist_initial_B, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate_B, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
fprintf('Names_i, FnsToEvaluate per ptype vs common, general eqm, PricePath r, this should be zero: %.3e \n',max(abs(PricePath_FB.r(:)-PricePath_GB.r(:))))
fprintf('Names_i, FnsToEvaluate per ptype vs common, general eqm, GEcondnPath, this should be zero: %.3e \n',max(abs(GEcondnPath_FB.CapitalMarket(:)-GEcondnPath_GB.CapitalMarket(:))))

output=struct();

end
