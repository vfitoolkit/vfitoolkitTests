function output=CoreFHorzTPathPType_GE_tminus1tplus1(T,PricePath,ParamPath,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,fastOLG)
% _tminus1 and _tplus1 on the general eqm transition path with permanent types. Two types (kappa_j
% differs, given as a structure), each with its own interest rate and a general eqm condition
% conditional on ptype (transpathoptions.GEptype), so the types do not interact and the PType solve
% must equal two solo solves without permanent types. The general eqm condition uses r_tplus1 and the
% aggregate assets_tminus1, and a FnsToEvaluate uses r_tminus1, r_tplus1 and w_tminus1 (so each ptype
% sees them on its own, inside its own path solve). The initial values depend on ptype, given as a
% structure (r) and as a vector of length N_i (assets); w is common. The interest rate by ptype is a
% structure with its fields in the opposite order to Names_i.
% Shooting for maxiter=2, with transpathoptions.fastOLG=fastOLG. Only z (no d, no e, no semiz).

Names_i={'low','high'};
N_i=2;

n_d=0;
d_grid=[];
vfoptions=struct();
simoptions=struct();

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate.interest=@(aprime,a,z,r_tminus1,r_tplus1,w_tminus1) a*r_tminus1+aprime*r_tplus1+w_tminus1*z;

GeneralEqmEqns.CapitalMarket=@(r,r_tplus1,assets_tminus1,interest,earnings) r-0.02-0.005*(r_tplus1-r)-0.01*assets_tminus1/(1+earnings)+0.001*interest;

jequaloneDist=zeros([n_a,n_z],'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

Params_PT=Params;
Params_PT.kappa_j=struct();
Params_PT.kappa_j.low=Params.kappa_j;
Params_PT.kappa_j.high=1.4*Params.kappa_j;
PTypeDistParamNames={'ptweights_GEtest'};
Params_PT.ptweights_GEtest=[0.7;0.3];

initialvalues_r=[0.04,0.045];
initialvalues_assets=[0.5,0.6];

%% Per type: the initial agent distribution (stationary at the period-1 prices) and V_final (stationary at the final prices)
V_final_PT=struct(); AgentDist_initial_PT=struct();
Params_solo=cell(N_i,1);
for ii=1:N_i
    Params_solo{ii}=Params;
    Params_solo{ii}.kappa_j=Params_PT.kappa_j.(Names_i{ii});
    Params_ii=Params_solo{ii};
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(Names_i{ii})=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(Names_i{ii}),~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
end
AgentDist_initial_PT.ptweights=[0.7;0.3];

%% One PType solve
% A price that depends on ptype, as a structure, deliberately with its fields in the opposite order to
% Names_i: the solver must read them by name
PricePath_PT=struct();
PricePath_PT.r=struct('high',PricePath.r,'low',PricePath.r);

transpathoptions=struct();
transpathoptions.fastOLG=fastOLG;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.GEptype={'CapitalMarket'};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;
transpathoptions.initialvalues.r.low=initialvalues_r(1);
transpathoptions.initialvalues.r.high=initialvalues_r(2);
transpathoptions.initialvalues.w=1.2;
transpathoptions.initialvalues.assets=initialvalues_assets'; % one per ptype

[PricePath_GE,GEcondnPath_GE]=TransitionPath_Case1_FHorz_PType(PricePath_PT, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

%% Two solo solves, without permanent types
for ii=1:N_i
    nm=Names_i{ii};
    transpathoptions_solo=struct();
    transpathoptions_solo.fastOLG=fastOLG;
    transpathoptions_solo.GEnewprice=3;
    transpathoptions_solo.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
    transpathoptions_solo.maxiter=2;
    transpathoptions_solo.verbose=0;
    transpathoptions_solo.initialvalues.r=initialvalues_r(ii);
    transpathoptions_solo.initialvalues.w=1.2;
    transpathoptions_solo.initialvalues.assets=initialvalues_assets(ii);
    [PricePath_solo,GEcondnPath_solo]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final_PT.(nm), AgentDist_initial_PT.(nm), jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_solo{ii}, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions_solo, simoptions, vfoptions);

    fprintf('GE _tminus1/_tplus1 by ptype (fastOLG=%i), PType vs solo, PricePath r (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(PricePath_GE.r.(nm)(:)-PricePath_solo.r(:))))
    fprintf('GE _tminus1/_tplus1 by ptype (fastOLG=%i), PType vs solo, GEcondnPath (type %s), this should be zero: %.3e \n',fastOLG,nm,max(abs(GEcondnPath_GE.CapitalMarket.(nm)(:)-GEcondnPath_solo.CapitalMarket(:))))
end

output=struct();

end
