function output=CoreInfHorzTPathPType_NiTwoTypes_vsTwoSolves(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Two genuinely different types (sigma differs, given as a structure) against one solo solve per ptype
% with the commands without permanent types. Only z (no d, no e).
%   (A) common PricePath and ParamPath: the standalone commands, each ptype vs its solo solve, and the
%       pooled aggregates vs the ptweight-weighted solo ones
%   (B) the general eqm conditions the solver returns vs those evaluated by hand on the pooled aggregates
%       of (A). maxiter=1, so GEcondnPath is the conditions at the INPUT path, which is the path the
%       standalone commands were given. 'earnings' reads w (on the ParamPath), so this fails if the
%       aggregates inside the path do not see the period-t parameters.
%   (C) w on the ParamPath by ptype, as a structure: the standalone commands vs solo, and (B) again
%   (D) w on the ParamPath by ptype, as an N_i-by-T matrix: the same answers as (C)
%   (E) r on the PricePath by ptype, as a structure: the standalone commands vs solo

Names_i={'low','high'};
N_i=2;
n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
FnNames=fieldnames(FnsToEvaluate);
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);

ptweights=[0.6;0.4];
PTypeDistParamNames={'ptweights_test'};
Params_PT=Params;
Params_PT.ptweights_test=ptweights;
Params_PT.sigma=struct();
Params_PT.sigma.low=2;
Params_PT.sigma.high=3;

%% Per type: V_final (stationary at the final prices) and AgentDist_initial (stationary at the period-1 prices)
V_final_PT=struct(); Policy_final_PT=struct(); AgentDist_initial_PT=struct();
Params_solo=cell(N_i,1);
for ii=1:N_i
    Params_solo{ii}=Params;
    Params_solo{ii}.sigma=Params_PT.sigma.(Names_i{ii});
    Params_ii=Params_solo{ii};
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(Names_i{ii})=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptions,Params_ii,[]);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(Names_i{ii}),Policy_final_PT.(Names_i{ii})]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
end
AgentDist_initial_PT.ptweights=ptweights;

transpathoptions_sa=transpathoptionsbaseline;
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.maxiter=1; % GEcondnPath is the conditions at the path handed in
transpathoptions.verbose=0;

% The paths by ptype used in (C)-(E): typehigh has a 5% higher w, and an interest rate 0.005 higher
wpath_PT={ParamPath.w, 1.05*ParamPath.w};
rpath_PT={PricePath.r, PricePath.r+0.005};

%% (A) Common PricePath and ParamPath
[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_PT=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT, PricePath, ParamPath, PolicyPath_PT, n_d,n_a,n_z, Names_i, pi_z, T, Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, Names_i, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);

pooled=struct();
for ff=1:length(FnNames)
    pooled.(FnNames{ff})=zeros(1,T);
end
for ii=1:N_i
    nm=Names_i{ii};
    [VPath_ii,PolicyPath_ii]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_PT.(nm), Policy_final_PT.(nm), Params_solo{ii}, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_ii=AgentDistOnTransPath_InfHorz(AgentDist_initial_PT.(nm), PricePath, ParamPath, PolicyPath_ii, n_d,n_a,n_z, pi_z, T, Params_solo{ii}, transpathoptions_sa, simoptions);
    AggVarsPath_ii=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath_ii, PolicyPath_ii, PricePath, ParamPath, Params_solo{ii}, T, n_d, n_a, n_z, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    fprintf('Two types vs two solo solves, common paths, VPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(VPath_PT.(nm)(:))-gather(VPath_ii(:)))))
    fprintf('Two types vs two solo solves, common paths, PolicyPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(PolicyPath_PT.(nm)(:))-gather(PolicyPath_ii(:)))))
    fprintf('Two types vs two solo solves, common paths, AgentDistPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(AgentDistPath_PT.(nm)(:))-gather(AgentDistPath_ii(:)))))
    for ff=1:length(FnNames)
        fprintf('Two types vs two solo solves, common paths, AggVarsPath %s (type %s), this should be zero: %.3e \n',FnNames{ff},nm,max(abs(gather(AggVarsPath_PT.(FnNames{ff}).(nm).Mean(:))-gather(AggVarsPath_ii.(FnNames{ff}).Mean(:)))))
        pooled.(FnNames{ff})=pooled.(FnNames{ff})+ptweights(ii)*gather(reshape(AggVarsPath_ii.(FnNames{ff}).Mean,1,T));
    end
end
for ff=1:length(FnNames)
    fprintf('Two types vs two solo solves, common paths, AggVarsPath %s pooled vs ptweight-weighted solo, this should be zero: %.3e \n',FnNames{ff},max(abs(gather(reshape(AggVarsPath_PT.(FnNames{ff}).Mean,1,T))-pooled.(FnNames{ff}))))
end
fprintf('Two types vs two solo solves, for scale the two types have different AgentDistPath, this should NOT be zero: %.3e \n',max(abs(gather(AgentDistPath_PT.low(:))-gather(AgentDistPath_PT.high(:)))))

%% (B) In-path general eqm conditions vs those evaluated by hand on the standalone aggregates
[~,GEcondnPath_PT]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
GEcondn_byhand=zeros(1,T-1);
for tt=1:T-1
    GEcondn_byhand(tt)=GeneralEqmEqns.CapitalMarket(PricePath.r(tt),pooled.assets(tt),pooled.earnings(tt));
end
GEPT=gather(GEcondnPath_PT.CapitalMarket(:))';
fprintf('Two types, common paths, GE in-path vs standalone, GEcondnPath, this should be zero: %.3e \n',max(abs(GEPT(1:T-1)-GEcondn_byhand)))
% for scale: the same conditions with earnings at the INITIAL w, which is what in-path aggregates that
% missed the ParamPath would give
GEcondn_w1=zeros(1,T-1);
for tt=1:T-1
    GEcondn_w1(tt)=GeneralEqmEqns.CapitalMarket(PricePath.r(tt),pooled.assets(tt),pooled.earnings(tt)*Params.w/ParamPath.w(tt));
end
fprintf('Two types, common paths, for scale earnings at the initial w would move GEcondnPath by %.3e (must not be zero, or the check above cannot see a stale w) \n',max(abs(GEcondn_w1-GEcondn_byhand)))

%% (C) w by ptype, as a structure
ParamPath_C=struct();
ParamPath_C.w=struct();
for ii=1:N_i
    ParamPath_C.w.(Names_i{ii})=wpath_PT{ii};
end
pooledC=struct();
for ff=1:length(FnNames)
    pooledC.(FnNames{ff})=zeros(1,T);
end
[VPath_C,PolicyPath_C]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath_C, T, V_final_PT, Policy_final_PT, Params_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_C=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT, PricePath, ParamPath_C, PolicyPath_C, n_d,n_a,n_z, Names_i, pi_z, T, Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_C=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_C, PolicyPath_C, PricePath, ParamPath_C, Params_PT, T, n_d, n_a, n_z, Names_i, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    ParamPath_ii.w=wpath_PT{ii};
    [VPath_ii,PolicyPath_ii]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath_ii, T, V_final_PT.(nm), Policy_final_PT.(nm), Params_solo{ii}, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_ii=AgentDistOnTransPath_InfHorz(AgentDist_initial_PT.(nm), PricePath, ParamPath_ii, PolicyPath_ii, n_d,n_a,n_z, pi_z, T, Params_solo{ii}, transpathoptions_sa, simoptions);
    AggVarsPath_ii=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath_ii, PolicyPath_ii, PricePath, ParamPath_ii, Params_solo{ii}, T, n_d, n_a, n_z, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    fprintf('Two types vs two solo solves, w by ptype (structure), VPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(VPath_C.(nm)(:))-gather(VPath_ii(:)))))
    fprintf('Two types vs two solo solves, w by ptype (structure), AgentDistPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(AgentDistPath_C.(nm)(:))-gather(AgentDistPath_ii(:)))))
    fprintf('Two types vs two solo solves, w by ptype (structure), AggVarsPath earnings (type %s), this should be zero: %.3e \n',nm,max(abs(gather(AggVarsPath_C.earnings.(nm).Mean(:))-gather(AggVarsPath_ii.earnings.Mean(:)))))
    for ff=1:length(FnNames)
        pooledC.(FnNames{ff})=pooledC.(FnNames{ff})+ptweights(ii)*gather(reshape(AggVarsPath_ii.(FnNames{ff}).Mean,1,T));
    end
end
fprintf('Two types vs two solo solves, w by ptype (structure), for scale typehigh VPath with its own w vs the common w, this should NOT be zero: %.3e \n',max(abs(gather(VPath_C.high(:))-gather(VPath_PT.high(:)))))
[~,GEcondnPath_C]=TransitionPath_InfHorz_PType(PricePath, ParamPath_C, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
GEcondn_byhandC=zeros(1,T-1);
for tt=1:T-1
    GEcondn_byhandC(tt)=GeneralEqmEqns.CapitalMarket(PricePath.r(tt),pooledC.assets(tt),pooledC.earnings(tt));
end
GEPTC=gather(GEcondnPath_C.CapitalMarket(:))';
fprintf('Two types, w by ptype (structure), GE in-path vs solo standalone, GEcondnPath, this should be zero: %.3e \n',max(abs(GEPTC(1:T-1)-GEcondn_byhandC)))

%% (D) w by ptype, as an N_i-by-T matrix: must give the same as the structure in (C)
ParamPath_D=struct();
ParamPath_D.w=[wpath_PT{1}; wpath_PT{2}];
[VPath_D,PolicyPath_D]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath_D, T, V_final_PT, Policy_final_PT, Params_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_D=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT, PricePath, ParamPath_D, PolicyPath_D, n_d,n_a,n_z, Names_i, pi_z, T, Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_D=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_D, PolicyPath_D, PricePath, ParamPath_D, Params_PT, T, n_d, n_a, n_z, Names_i, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Two types, w by ptype, N_i-by-T matrix vs structure, VPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(VPath_D.(nm)(:))-gather(VPath_C.(nm)(:)))))
    fprintf('Two types, w by ptype, N_i-by-T matrix vs structure, AgentDistPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(AgentDistPath_D.(nm)(:))-gather(AgentDistPath_C.(nm)(:)))))
end
fprintf('Two types, w by ptype, N_i-by-T matrix vs structure, AggVarsPath earnings pooled, this should be zero: %.3e \n',max(abs(gather(AggVarsPath_D.earnings.Mean(:))-gather(AggVarsPath_C.earnings.Mean(:)))))
[~,GEcondnPath_D]=TransitionPath_InfHorz_PType(PricePath, ParamPath_D, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
GEPTD=gather(GEcondnPath_D.CapitalMarket(:))';
fprintf('Two types, w by ptype, N_i-by-T matrix vs solo standalone, general eqm, GEcondnPath, this should be zero: %.3e \n',max(abs(GEPTD(1:T-1)-GEcondn_byhandC)))

%% (E) r by ptype, as a structure, in the standalone commands
PricePath_E=struct();
PricePath_E.r=struct();
for ii=1:N_i
    PricePath_E.r.(Names_i{ii})=rpath_PT{ii};
end
[VPath_E,PolicyPath_E]=ValueFnOnTransPath_InfHorz_PType(PricePath_E, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_E=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT, PricePath_E, ParamPath, PolicyPath_E, n_d,n_a,n_z, Names_i, pi_z, T, Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_E=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_E, PolicyPath_E, PricePath_E, ParamPath, Params_PT, T, n_d, n_a, n_z, Names_i, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    PricePath_ii.r=rpath_PT{ii};
    [VPath_ii,PolicyPath_ii]=ValueFnOnTransPath_InfHorz(PricePath_ii, ParamPath, T, V_final_PT.(nm), Policy_final_PT.(nm), Params_solo{ii}, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_ii=AgentDistOnTransPath_InfHorz(AgentDist_initial_PT.(nm), PricePath_ii, ParamPath, PolicyPath_ii, n_d,n_a,n_z, pi_z, T, Params_solo{ii}, transpathoptions_sa, simoptions);
    AggVarsPath_ii=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath_ii, PolicyPath_ii, PricePath_ii, ParamPath, Params_solo{ii}, T, n_d, n_a, n_z, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    fprintf('Two types vs two solo solves, r by ptype (structure), VPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(VPath_E.(nm)(:))-gather(VPath_ii(:)))))
    fprintf('Two types vs two solo solves, r by ptype (structure), AgentDistPath (type %s), this should be zero: %.3e \n',nm,max(abs(gather(AgentDistPath_E.(nm)(:))-gather(AgentDistPath_ii(:)))))
    fprintf('Two types vs two solo solves, r by ptype (structure), AggVarsPath assets (type %s), this should be zero: %.3e \n',nm,max(abs(gather(AggVarsPath_E.assets.(nm).Mean(:))-gather(AggVarsPath_ii.assets.Mean(:)))))
end

output=struct();

end
