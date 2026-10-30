function output=CoreInfHorzTPathPType_NiAllIdentical_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,transpathoptionsbaseline,combo)
% Permanent types that are all identical must give exactly what the commands without permanent types
% give, for one of the 8 combinations of d, z and e (combo, e.g. 'nod_z_e').
%   (1) the standalone commands, N_i=2 identical types: VPath, PolicyPath, AgentDistPath of each ptype,
%       and the pooled AggVarsPath (ptweights 0.6 and 0.4 sum to one)
%   (2) the general eqm path, shooting for maxiter=2 (returns the price path after one update together
%       with the general eqm conditions evaluated on it), N_i=1 and N_i=2 identical types
% The FnsToEvaluate 'earnings' reads w, which is on the ParamPath, and the general eqm condition reads
% earnings, so the aggregates computed inside the path must see the period-t parameters, not the
% initial ones. V_final is the stationary value fn at the final prices, AgentDist_initial the stationary
% agent distribution at the period-1 prices.

vfoptions=struct();
simoptions=struct();
if contains(combo,'_z_')
    % n_z, z_grid and pi_z as given
else
    n_z=0;
    z_grid=[];
    pi_z=[];
end
if contains(combo,'_e')
    vfoptions.n_e=vfoptionsbaseline.n_e;
    vfoptions.e_grid=vfoptionsbaseline.e_grid;
    vfoptions.pi_e=vfoptionsbaseline.pi_e;
    simoptions.n_e=simoptionsbaseline.n_e;
    simoptions.e_grid=simoptionsbaseline.e_grid;
    simoptions.pi_e=simoptionsbaseline.pi_e;
end
usesemiz=contains(combo,'_semiz');
if usesemiz
    vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
    vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
    vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
    vfoptions.l_dsemiz=vfoptionsbaseline.l_dsemiz;
    simoptions.n_semiz=simoptionsbaseline.n_semiz;
    simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
    simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
    simoptions.l_dsemiz=simoptionsbaseline.l_dsemiz;
    simoptions.d_grid=d_grid;
end
if strcmp(combo(1:3),'nod') && ~usesemiz % without d1, a semiz model still has its d2 (n_d and d_grid as given)
    n_d=0;
    d_grid=[];
end

switch combo
    case 'nod_noz_noe'
        ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.earnings=@(aprime,a,w) w;
    case 'd_noz_noe'
        ReturnFn=@(d,aprime,a,r,w,sigma,eta,varphi) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a) a;
        FnsToEvaluate.earnings=@(d,aprime,a,w) w*d;
    case 'nod_z_noe'
        ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w) w*z;
    case 'd_z_noe'
        ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*z*d;
    case 'nod_noz_e'
        ReturnFn=@(aprime,a,e,r,w,sigma) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,e) a;
        FnsToEvaluate.earnings=@(aprime,a,e,w) w*e;
    case 'd_noz_e'
        ReturnFn=@(d,aprime,a,e,r,w,sigma,eta,varphi) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,e) a;
        FnsToEvaluate.earnings=@(d,aprime,a,e,w) w*d*e;
    case 'nod_z_e'
        ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(aprime,a,z,e,w) w*z*e;
    case 'd_z_e'
        ReturnFn=@(d,aprime,a,z,e,r,w,sigma,eta,varphi) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,e,w) w*z*e*d;
    case 'nod1_noz_noe_semiz'
        ReturnFn=@(d2,aprime,a,semiz,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,sigma,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,w) w*semiz;
    case 'd1_noz_noe_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,w) w*d1*semiz;
    case 'nod1_z_noe_semiz'
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,sigma,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w) w*z*semiz;
    case 'd1_z_noe_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,z,w) w*d1*z*semiz;
    case 'nod1_noz_e_semiz'
        ReturnFn=@(d2,aprime,a,semiz,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,sigma,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,e) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,e,w) w*e*semiz;
    case 'd1_noz_e_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz,e) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,e,w) w*d1*e*semiz;
    case 'nod1_z_e_semiz'
        ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,sigma,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z,e) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,e,w) w*z*e*semiz;
    case 'd1_z_e_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_e_semiz(d1,d2,aprime,a,semiz,z,e,r,w,sigma,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz,z,e) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,z,e,w) w*d1*z*e*semiz;
    otherwise
        error('CoreInfHorzTPathPType_NiAllIdentical_vsNoPType: unknown combo %s',combo)
end

%% Initial agent distribution (stationary, at the period-1 prices) and V_final (stationary, at the final prices)
Params_init=Params;
Params_init.r=PricePath.r(1);
Params_init.w=ParamPath.w(1);
[~,Policy_init]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_init,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_init,n_d,n_a,n_z,pi_z,simoptions,Params_init,[]);
Params_final=Params;
Params_final.r=PricePath.r(end);
Params_final.w=ParamPath.w(end);
[V_final,Policy_final]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_final,DiscountFactorParamNames,[],vfoptions);

names_PT={'ptype001','ptype002'}; % auto-names used when N_i is numeric
V_final_PT1.ptype001=V_final;
Policy_final_PT1.ptype001=Policy_final;
AgentDist_initial_PT1.ptype001=AgentDist_initial;
AgentDist_initial_PT1.ptweights=1;
V_final_PT2=struct(); Policy_final_PT2=struct(); AgentDist_initial_PT2=struct();
for ii=1:2
    V_final_PT2.(names_PT{ii})=V_final;
    Policy_final_PT2.(names_PT{ii})=Policy_final;
    AgentDist_initial_PT2.(names_PT{ii})=AgentDist_initial;
end
AgentDist_initial_PT2.ptweights=[0.6;0.4];
PTypeDistParamNames={'ptweights_test'};
Params1=Params;
Params1.ptweights_test=1;
Params2=Params;
Params2.ptweights_test=[0.6;0.4];

%% (1) The standalone commands
transpathoptions_sa=transpathoptionsbaseline;
[VPath,PolicyPath]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPath, n_d,n_a,n_z, pi_z, T, Params, transpathoptions_sa, simoptions);
AggVarsPath=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPath, PolicyPath, PricePath, ParamPath, Params, T, n_d, n_a, n_z, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);

[VPath_PT,PolicyPath_PT]=ValueFnOnTransPath_InfHorz_PType(PricePath, ParamPath, T, V_final_PT2, Policy_final_PT2, Params2, n_d,n_a,n_z, 2, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_PT=AgentDistOnTransPath_InfHorz_PType(AgentDist_initial_PT2, PricePath, ParamPath, PolicyPath_PT, n_d,n_a,n_z, 2, pi_z, T, Params2, transpathoptions_sa, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_InfHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params2, T, n_d, n_a, n_z, 2, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);

for ii=1:2
    nm=names_PT{ii};
    fprintf('All identical vs no-PType (%s), VPath (type %d), this should be zero: %.3e \n',combo,ii,max(abs(gather(VPath_PT.(nm)(:))-gather(VPath(:)))))
    fprintf('All identical vs no-PType (%s), PolicyPath (type %d), this should be zero: %.3e \n',combo,ii,max(abs(gather(PolicyPath_PT.(nm)(:))-gather(PolicyPath(:)))))
    fprintf('All identical vs no-PType (%s), AgentDistPath (type %d), this should be zero: %.3e \n',combo,ii,max(abs(gather(AgentDistPath_PT.(nm)(:))-gather(AgentDistPath(:)))))
end
fprintf('All identical vs no-PType (%s), AggVarsPath assets pooled, this should be zero: %.3e \n',combo,max(abs(gather(AggVarsPath_PT.assets.Mean(:))-gather(AggVarsPath.assets.Mean(:)))))
fprintf('All identical vs no-PType (%s), AggVarsPath earnings pooled, this should be zero: %.3e \n',combo,max(abs(gather(AggVarsPath_PT.earnings.Mean(:))-gather(AggVarsPath.earnings.Mean(:)))))

%% (2) The general eqm path
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;

[PricePath_NoPT,GEcondnPath_NoPT]=TransitionPath_InfHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);
[PricePath_PT1,GEcondnPath_PT1]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_PT1, AgentDist_initial_PT1, n_d,n_a,n_z, 1, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params1, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);
[PricePath_PT2,GEcondnPath_PT2]=TransitionPath_InfHorz_PType(PricePath, ParamPath, T, V_final_PT2, AgentDist_initial_PT2, n_d,n_a,n_z, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

% Compared over periods 1 to T-1, the ones the path evaluates
GEnoPT=gather(GEcondnPath_NoPT.CapitalMarket(:));
GEPT1=gather(GEcondnPath_PT1.CapitalMarket(:));
GEPT2=gather(GEcondnPath_PT2.CapitalMarket(:));
fprintf('GE vs no-PType (%s), shooting, N_i=1, PricePath r, this should be zero: %.3e \n',combo,max(abs(gather(PricePath_PT1.r(:))-gather(PricePath_NoPT.r(:)))))
fprintf('GE vs no-PType (%s), shooting, N_i=1, GEcondnPath, this should be zero: %.3e \n',combo,max(abs(GEPT1(1:T-1)-GEnoPT(1:T-1))))
fprintf('GE vs no-PType (%s), shooting, N_i=2 identical types, PricePath r, this should be zero: %.3e \n',combo,max(abs(gather(PricePath_PT2.r(:))-gather(PricePath_NoPT.r(:)))))
fprintf('GE vs no-PType (%s), shooting, N_i=2 identical types, GEcondnPath, this should be zero: %.3e \n',combo,max(abs(GEPT2(1:T-1)-GEnoPT(1:T-1))))
fprintf('GE vs no-PType (%s), for scale the update moved r by %.3e (must not be zero, or the PricePath checks above are empty) \n',combo,max(abs(gather(PricePath_NoPT.r(:))-PricePath.r(:))))

output=struct();

end
