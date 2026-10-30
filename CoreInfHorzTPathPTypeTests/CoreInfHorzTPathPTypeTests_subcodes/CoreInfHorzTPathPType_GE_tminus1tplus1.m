function output=CoreInfHorzTPathPType_GE_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% _tminus1 and _tplus1 on the InfHorz general eqm transition path with permanent types. The InfHorz
% counterpart of CoreFHorzTPathPType_GE_tminus1tplus1. Two types (sigma differs, given as a structure),
% each with its own interest rate and a general eqm condition conditional on ptype
% (transpathoptions.GEptype), so the types do not interact and the PType solve must equal two solo solves
% without permanent types. The general eqm condition uses r_tplus1 and the aggregate assets_tminus1, and a
% FnsToEvaluate uses r_tminus1, r_tplus1 and w_tminus1 (so each ptype sees them on its own, inside its own
% path solve). The initial values depend on ptype, given as a structure (r) and as a vector of length N_i
% (assets); w is common. The interest rate by ptype is a structure with its fields in the opposite order
% to Names_i. The aggregate assets_tminus1 is the mean conditional on the ptype.
%   (1) shooting for maxiter=2, PType vs two solo solves
%   (2) Anderson with memory=0 vs shooting at maxiter=3 (Anderson is not separable by ptype, so not vs
%       solo), and with the default memory it must differ
%   (3) the standalone commands at the input price path, EvalFnOnTransPath_AggVars_InfHorz_PType with
%       FnsToEvaluate that use the lags and leads, each ptype vs the command without permanent types on
%       its own paths with its own initial values
% Only z (no d, no e).

Names_i={'low','high'};
N_i=2;
n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
FnsToEvaluate.interest=@(aprime,a,z,r_tminus1,r_tplus1,w_tminus1) a*r_tminus1+aprime*r_tplus1+w_tminus1*z;

GeneralEqmEqns.CapitalMarket=@(r,r_tplus1,assets_tminus1,interest,earnings) r-0.02-0.005*(r_tplus1-r)-0.01*assets_tminus1/(1+earnings)+0.001*interest;

Params_PT=Params;
Params_PT.sigma=struct();
Params_PT.sigma.low=2;
Params_PT.sigma.high=3;
PTypeDistParamNames={'ptweights_GEtest'};
Params_PT.ptweights_GEtest=[0.7;0.3];

initialvalues_r=[0.04,0.045];
initialvalues_assets=[0.5,0.6];

%% Per type: the initial agent distribution (stationary at the period-1 prices) and V_final (stationary at the final prices)
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
AgentDist_initial_PT.ptweights=[0.7;0.3];

%% (1) One PType solve
% A price that depends on ptype, as a structure, deliberately with its fields in the opposite order to
% Names_i: the solver must read them by name
PricePath_PT=struct();
PricePath_PT.r=struct('high',PricePath.r,'low',PricePath.r);

transpathoptions=struct();
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.GEptype={'CapitalMarket'};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;
transpathoptions.initialvalues.r.low=initialvalues_r(1);
transpathoptions.initialvalues.r.high=initialvalues_r(2);
transpathoptions.initialvalues.w=1.2;
transpathoptions.initialvalues.assets=initialvalues_assets'; % one per ptype

[PricePath_GE,GEcondnPath_GE]=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

%% Two solo solves, without permanent types
transpathoptions_solo=cell(N_i,1);
for ii=1:N_i
    nm=Names_i{ii};
    transpathoptions_solo{ii}=struct();
    transpathoptions_solo{ii}.GEnewprice=3;
    transpathoptions_solo{ii}.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
    transpathoptions_solo{ii}.maxiter=2;
    transpathoptions_solo{ii}.verbose=0;
    transpathoptions_solo{ii}.initialvalues.r=initialvalues_r(ii);
    transpathoptions_solo{ii}.initialvalues.w=1.2;
    transpathoptions_solo{ii}.initialvalues.assets=initialvalues_assets(ii);
    [PricePath_s,GEcondnPath_s]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, transpathoptions_solo{ii}, simoptions, vfoptions, []);
    GEs=gather(GEcondnPath_s.CapitalMarket(:));
    GEp=gather(GEcondnPath_GE.CapitalMarket.(nm)(:));
    fprintf('GE _tminus1/_tplus1 by ptype, PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',nm,max(abs(gather(PricePath_GE.r.(nm)(:))-gather(PricePath_s.r(:)))))
    fprintf('GE _tminus1/_tplus1 by ptype, PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',nm,max(abs(GEp(1:T-1)-GEs(1:T-1))))
end
% (the initial values differ by ptype, so a solver that gave every ptype the same one would fail the above)

%% (2) Anderson acceleration (GEnewprice=2)
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
transpathoptionsA.anderson.memory=0;
PricePathA0=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptionsA=rmfield(transpathoptionsA,'anderson');
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
PricePathA5=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
transpathoptionsS3=transpathoptions;
transpathoptionsS3.maxiter=3;
PricePathS3=TransitionPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, PTypeDistParamNames, transpathoptionsS3, simoptions, vfoptions);
fprintf('GE _tminus1/_tplus1 by ptype, Anderson with memory=0 vs shooting, both at maxiter=3, both types, this should be zero: %.3e \n',max([max(abs(gather(PricePathA0.r.low(:))-gather(PricePathS3.r.low(:)))),max(abs(gather(PricePathA0.r.high(:))-gather(PricePathS3.r.high(:))))]))
fprintf('GE _tminus1/_tplus1 by ptype, Anderson with the default memory vs shooting, both at maxiter=3, both types, this should NOT be zero: %.3e \n',max([max(abs(gather(PricePathA5.r.low(:))-gather(PricePathS3.r.low(:)))),max(abs(gather(PricePathA5.r.high(:))-gather(PricePathS3.r.high(:))))]))

%% (3) The standalone commands with lags and leads, initial values by ptype
FnsToEvaluateL.assets=FnsToEvaluate.assets;
FnsToEvaluateL.interest=FnsToEvaluate.interest;
FnsToEvaluateL.y=@(aprime,a,z,assets_tminus1) a*assets_tminus1;
FnNamesL=fieldnames(FnsToEvaluateL);
transpathoptions_sa=struct();
transpathoptions_sa.initialvalues=transpathoptions.initialvalues;
[~,PolicyPath_PT]=ValueFnOnTransPath_InfHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d,n_a,n_z, Names_i, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_PT=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT, PricePath_PT, ParamPath, PolicyPath_PT, n_d,n_a,n_z, Names_i, pi_z, T, Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluateL, AgentDistPath_PT, PolicyPath_PT, PricePath_PT, ParamPath, Params_PT, T, n_d, n_a, n_z, Names_i, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);
for ii=1:N_i
    nm=Names_i{ii};
    transpathoptions_ii=struct();
    transpathoptions_ii.initialvalues.r=initialvalues_r(ii);
    transpathoptions_ii.initialvalues.w=1.2;
    transpathoptions_ii.initialvalues.assets=initialvalues_assets(ii);
    AggVarsPath_ii=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluateL, AgentDistPath_PT.(nm), PolicyPath_PT.(nm), PricePath, ParamPath, Params_solo{ii}, T, n_d, n_a, n_z, d_grid, a_grid,z_grid, transpathoptions_ii, simoptions);
    for ff=1:length(FnNamesL)
        fprintf('_tminus1/_tplus1 by ptype, EvalFnOnTransPath_AggVars_PType vs without permanent types with its own initial values, %s (type %s), this should be zero: %.3e \n',FnNamesL{ff},nm,max(abs(gather(AggVarsPath_PT.(FnNamesL{ff}).(nm).Mean(:))-gather(AggVarsPath_ii.(FnNamesL{ff}).Mean(:)))))
    end
end

output=struct();

end
