function output=CoreInfHorzTPathPType_ShocksByType(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,Names_i)
% Exogenous shocks that differ by ptype on a transition path, one ptype's varying over the path:
%   zvar: n_z points, the z_grid scaled by a factor that falls from 1.1 to 1 over the first half of the
%         path, given as a [n_z,T] grid (pi_z constant)
%   zfix: 3 points, its own constant z_grid and pi_z
% The exogenous shock setup decides per ptype whether the z path is trivial (constant) and writes the
% answer into transpathoptions, so ptypes that share one transpathoptions would all get whichever ptype
% was set up last. Run with Names_i in both orders, so that each ptype is the last one once.
% Checked against one solo solve per ptype:
%   (A) the standalone commands at a common price path: VPath, AgentDistPath, AggVarsPath
%   (B) the general eqm path, with the capital market condition by ptype (GEptype) and an interest
%       rate per ptype, shooting for maxiter=2
% V_final is the stationary value fn at the final prices and z_grid, AgentDist_initial the stationary
% agent distribution at the period-1 prices and z_grid.

N_i=2;
n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();
ordername=strjoin(Names_i,',');

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);

%% The shocks of each ptype
n_z_PT=struct(); z_grid_PT=struct(); pi_z_PT=struct();
zmult=[linspace(1.1,1,floor(T/2)),ones(1,T-floor(T/2))];
n_z_PT.zvar=n_z;
z_grid_PT.zvar=z_grid*zmult; % [n_z,T]
pi_z_PT.zvar=pi_z;
[z_grid_fix,pi_z_fix]=discretizeAR1_FarmerToda(0,0.8,0.05,3);
n_z_PT.zfix=3;
z_grid_PT.zfix=exp(z_grid_fix);
pi_z_PT.zfix=pi_z_fix;
% the period-1 and the final grid, for the stationary problems at the two ends of the path
z_grid_1=struct('zvar',z_grid_PT.zvar(:,1),'zfix',z_grid_PT.zfix);
z_grid_T=struct('zvar',z_grid_PT.zvar(:,T),'zfix',z_grid_PT.zfix);

%% Per type: V_final and AgentDist_initial
V_final_PT=struct(); Policy_final_PT=struct(); AgentDist_initial_PT=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z_PT.(nm),d_grid,a_grid,z_grid_1.(nm),pi_z_PT.(nm),ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(nm)=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z_PT.(nm),pi_z_PT.(nm),simoptions,Params_ii,[]);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(nm),Policy_final_PT.(nm)]=ValueFnIter_InfHorz(n_d,n_a,n_z_PT.(nm),d_grid,a_grid,z_grid_T.(nm),pi_z_PT.(nm),ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
end
AgentDist_initial_PT.ptweights=[0.5;0.5];
PTypeDistParamNames={'ptweights_test'};
Params_PT=Params;
Params_PT.ptweights_test=[0.5;0.5];

transpathoptions_sa=struct();

%% (A) The standalone commands
[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d,n_a,n_z_PT, Names_i, d_grid,a_grid,z_grid_PT, pi_z_PT, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_PT=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT, PricePath, ParamPath, PolicyPath_PT, n_d,n_a,n_z_PT, Names_i, pi_z_PT, T, Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z_PT, Names_i, d_grid, a_grid, z_grid_PT, transpathoptions_sa, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    [VPath_ii,PolicyPath_ii]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final_PT.(nm), Policy_final_PT.(nm), Params, n_d,n_a,n_z_PT.(nm), d_grid,a_grid,z_grid_PT.(nm), pi_z_PT.(nm), DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_ii=AgentDistOnTransPath_InfHorz(AgentDist_initial_PT.(nm), PricePath, ParamPath, PolicyPath_ii, n_d,n_a,n_z_PT.(nm), pi_z_PT.(nm), T, Params, transpathoptions_sa, simoptions);
    AggVarsPath_ii=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath_ii, PolicyPath_ii, PricePath, ParamPath, Params, T, n_d, n_a, n_z_PT.(nm), d_grid, a_grid,z_grid_PT.(nm), transpathoptions_sa, simoptions);
    fprintf('Shocks by type (Names_i order %s), standalone vs solo, VPath (type %s), this should be zero: %.3e \n',ordername,nm,max(abs(gather(VPath_PT.(nm)(:))-gather(VPath_ii(:)))))
    fprintf('Shocks by type (Names_i order %s), standalone vs solo, AgentDistPath (type %s), this should be zero: %.3e \n',ordername,nm,max(abs(gather(AgentDistPath_PT.(nm)(:))-gather(AgentDistPath_ii(:)))))
    fprintf('Shocks by type (Names_i order %s), standalone vs solo, AggVarsPath earnings (type %s), this should be zero: %.3e \n',ordername,nm,max(abs(gather(AggVarsPath_PT.earnings.(nm).Mean(:))-gather(AggVarsPath_ii.earnings.Mean(:)))))
end

%% (B) The general eqm path, by ptype
PricePath_PT=struct();
PricePath_PT.r=struct();
for ii=1:N_i
    PricePath_PT.r.(Names_i{ii})=PricePath.r;
end
transpathoptions=struct();
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.GEptype={'CapitalMarket'};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;
[PricePath_GE,GEcondnPath_GE]=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z_PT, Names_i, d_grid,a_grid,z_grid_PT, pi_z_PT, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
transpathoptions_solo=rmfield(transpathoptions,'GEptype');
for ii=1:N_i
    nm=Names_i{ii};
    [PricePath_s,GEcondnPath_s]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), n_d,n_a,n_z_PT.(nm), d_grid,a_grid,z_grid_PT.(nm), pi_z_PT.(nm), ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_solo, simoptions, vfoptions, []);
    GEs=gather(GEcondnPath_s.CapitalMarket(:));
    GEp=gather(GEcondnPath_GE.CapitalMarket.(nm)(:));
    fprintf('Shocks by type (Names_i order %s), GE PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',ordername,nm,max(abs(gather(PricePath_GE.r.(nm)(:))-gather(PricePath_s.r(:)))))
    fprintf('Shocks by type (Names_i order %s), GE PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',ordername,nm,max(abs(GEp(1:T-1)-GEs(1:T-1))))
end

% for scale: the zvar solo path with its z_grid held at the final grid, which is what a ptype that was
% given the 'trivial z path' of the other ptype would see
[~,GEcondnPath_c]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_PT.zvar, AgentDist_initial_PT.zvar, n_d,n_a,n_z, d_grid,a_grid,z_grid_T.zvar, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions_solo, simoptions, vfoptions, []);
GEc=gather(GEcondnPath_c.CapitalMarket(:));
GEzvar=gather(GEcondnPath_GE.CapitalMarket.zvar(:));
fprintf('Shocks by type (Names_i order %s), for scale zvar with a constant z_grid would move GEcondnPath by %.3e (must not be zero) \n',ordername,max(abs(GEc(1:T-1)-GEzvar(1:T-1))))

output=struct();

end
