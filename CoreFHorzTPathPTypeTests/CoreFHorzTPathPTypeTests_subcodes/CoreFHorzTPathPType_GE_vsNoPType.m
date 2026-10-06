function output=CoreFHorzTPathPType_GE_vsNoPType(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline,combo,fastOLG)
% The general eqm transition path with permanent types (TransitionPath_Case1_FHorz_PType) against the
% one without (TransitionPath_Case1_FHorz), for one of the 16 combinations of d (d1), z, e and semiz,
% with transpathoptions.fastOLG=fastOLG. combo is the name of the model, as in the CoreFHorzTPath
% subcodes (e.g. 'nod_z_e_nosemiz', 'd1_noz_e_semiz'); n_d and d_grid are the ones that model uses.
%
% Each solve runs the shooting algorithm for maxiter=2, so it returns the price path after one update
% together with the general eqm conditions evaluated on it. The model has a general eqm condition in
% the interest rate that depends on assets and on earnings, and earnings depend on w, which is on the
% ParamPath, so the in-path aggregates must see the period-t parameters.
%   (1) N_i=1 is the model without permanent types: PricePath and GEcondnPath must be identical.
%   (2) N_i=2 identical types, ptweights 0.6 and 0.4: the same, up to the rounding in 0.6*x+0.4*x.
%   (3) N_i=2 genuinely different types (named, kappa_j differs by type, given as a structure): the
%       GEcondnPath returned by the solver must equal the general eqm condition evaluated by hand on
%       the aggregates of the standalone commands (ValueFnOnTransPath_Case1_FHorz_PType,
%       AgentDistOnTransPath_Case1_FHorz_PType, EvalFnOnTransPath_AggVars_Case1_FHorz_PType). This one
%       runs maxiter=1 and evaluates the standalone commands at the INPUT price path: GEcondnPath is
%       always the conditions on the last path evaluated, but when maxiter runs out the returned price
%       path has had one more update (as in CoreFHorzTPath_demographics and _tminus1tplus1).
% The initial agent distribution is the stationary one at the period-1 prices, and V_final the
% stationary value fn at the final prices (which the path ends at).

vfoptions=struct();
simoptions=struct();

usesemiz=contains(combo,'_semiz');
if contains(combo,'_z_')
    % n_z, z_grid and pi_z as given
else
    n_z=0;
    z_grid=[];
    pi_z=[];
end
if contains(combo,'_e_')
    vfoptions.n_e=vfoptionsbaseline.n_e;
    vfoptions.e_grid=vfoptionsbaseline.e_grid;
    vfoptions.pi_e=vfoptionsbaseline.pi_e;
    simoptions.n_e=simoptionsbaseline.n_e;
    simoptions.e_grid=simoptionsbaseline.e_grid;
    simoptions.pi_e=simoptionsbaseline.pi_e;
end
if usesemiz
    vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
    vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
    vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
    simoptions.n_semiz=simoptionsbaseline.n_semiz;
    simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
    simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
    simoptions.d_grid=d_grid;
end
if strcmp(combo(1:3),'nod')
    if ~usesemiz
        n_d=0;
        d_grid=[];
    end
end

%% The model pieces of this combination
switch combo
    case 'nod_noz_noe_nosemiz'
        ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.earnings=@(aprime,a,w,kappa_j) w*kappa_j;
        jequaloneDist=zeros(n_a,1,'gpuArray');
        jequaloneDist(1)=1;
    case 'd_noz_noe_nosemiz'
        ReturnFn=@(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluate.assets=@(d,aprime,a) a;
        FnsToEvaluate.earnings=@(d,aprime,a,w,kappa_j) w*kappa_j*d;
        jequaloneDist=zeros(n_a,1,'gpuArray');
        jequaloneDist(1)=1;
    case 'nod_z_noe_nosemiz'
        ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
        jequaloneDist=zeros([n_a,n_z],'gpuArray');
        jequaloneDist(1,ceil(n_z/2))=1;
    case 'd_z_noe_nosemiz'
        ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w,kappa_j) w*kappa_j*z*d;
        jequaloneDist=zeros([n_a,n_z],'gpuArray');
        jequaloneDist(1,ceil(n_z/2))=1;
    case 'nod_noz_e_nosemiz'
        ReturnFn=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate.assets=@(aprime,a,e) a;
        FnsToEvaluate.earnings=@(aprime,a,e,w,kappa_j) w*kappa_j*e;
        jequaloneDist=zeros([n_a,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_e/2))=1;
    case 'd_noz_e_nosemiz'
        ReturnFn=@(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_noz_e_nosemiz(d,aprime,a,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluate.assets=@(d,aprime,a,e) a;
        FnsToEvaluate.earnings=@(d,aprime,a,e,w,kappa_j) w*kappa_j*d*e;
        jequaloneDist=zeros([n_a,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_e/2))=1;
    case 'nod_z_e_nosemiz'
        ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
        FnsToEvaluate.assets=@(aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
        jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
    case 'd_z_e_nosemiz'
        ReturnFn=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
        FnsToEvaluate.assets=@(d,aprime,a,z,e) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,e,w,kappa_j) w*kappa_j*z*e*d;
        jequaloneDist=zeros([n_a,n_z,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
    case 'nod1_noz_noe_semiz'
        ReturnFn=@(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,w,kappa_j) w*kappa_j*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2))=1;
    case 'd1_noz_noe_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_noe_semiz(d1,d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,w,kappa_j) w*kappa_j*d1*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2))=1;
    case 'nod1_z_noe_semiz'
        ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;
    case 'd1_z_noe_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_noe_semiz(d1,d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz,z) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*d1*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;
    case 'nod1_noz_e_semiz'
        ReturnFn=@(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,e) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,e,w,kappa_j) w*kappa_j*e*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(vfoptions.n_e/2))=1;
    case 'd1_noz_e_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_noz_e_semiz(d1,d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz,e) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,e,w,kappa_j) w*kappa_j*e*d1*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(vfoptions.n_e/2))=1;
    case 'nod1_z_e_semiz'
        ReturnFn=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d2,aprime,a,semiz,z,e) a;
        FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,e,w,kappa_j) w*kappa_j*z*e*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
    case 'd1_z_e_semiz'
        ReturnFn=@(d1,d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost) ReturnFn_d1_z_e_semiz(d1,d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,eta,varphi,uempbenefit,searcheffortcost);
        FnsToEvaluate.assets=@(d1,d2,aprime,a,semiz,z,e) a;
        FnsToEvaluate.earnings=@(d1,d2,aprime,a,semiz,z,e,w,kappa_j) w*kappa_j*z*e*d1*semiz;
        jequaloneDist=zeros([n_a,vfoptions.n_semiz,n_z,vfoptions.n_e],'gpuArray');
        jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2),ceil(vfoptions.n_e/2))=1;
    otherwise
        error('CoreFHorzTPathPType_GE_vsNoPType: unknown combo %s',combo)
end

%% General eqm: one condition, in the interest rate
GeneralEqmEqns.CapitalMarket=@(r,assets,earnings) r-0.02-0.01*assets/(1+earnings);

transpathoptions=struct();
transpathoptions.fastOLG=fastOLG;
transpathoptions.GEnewprice=3; % the shooting algorithm
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.5};
transpathoptions.maxiter=2; % returns the price path after one update, and the general eqm conditions evaluated on it
transpathoptions.verbose=0;

%% Initial agent distribution (stationary, at the period-1 prices) and V_final (stationary, at the final prices)
Params_init=Params;
Params_init.r=PricePath.r(1);
Params_init.w=ParamPath.w(1);
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_init,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params_init,simoptions);
Params_final=Params;
Params_final.r=PricePath.r(end);
Params_final.w=ParamPath.w(end);
[V_final,~]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_final,DiscountFactorParamNames,[],vfoptions);

%% Without permanent types
[PricePath_NoPT,GEcondnPath_NoPT]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);

%% (1) N_i=1
PTypeDistParamNames={'ptweights_GEtest'};
Params.ptweights_GEtest=1;
V_final_PT1.ptype001=V_final;
AgentDist_initial_PT1.ptype001=AgentDist_initial;
AgentDist_initial_PT1.ptweights=1;
[PricePath_PT1,GEcondnPath_PT1]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT1, AgentDist_initial_PT1, jequaloneDist, n_d, n_a, n_z, N_j, 1, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

fprintf('GE vs no-PType (%s, fastOLG=%i), N_i=1, PricePath r, this should be zero: %.3e \n',combo,fastOLG,max(abs(PricePath_PT1.r(:)-PricePath_NoPT.r(:))))
fprintf('GE vs no-PType (%s, fastOLG=%i), N_i=1, GEcondnPath, this should be zero: %.3e \n',combo,fastOLG,max(abs(GEcondnPath_PT1.CapitalMarket(:)-GEcondnPath_NoPT.CapitalMarket(:))))

%% (2) N_i=2 identical types
Params.ptweights_GEtest=[0.6;0.4];
V_final_PT2.ptype001=V_final;
V_final_PT2.ptype002=V_final;
AgentDist_initial_PT2.ptype001=AgentDist_initial;
AgentDist_initial_PT2.ptype002=AgentDist_initial;
AgentDist_initial_PT2.ptweights=[0.6;0.4];
[PricePath_PT2,GEcondnPath_PT2]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT2, AgentDist_initial_PT2, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions, simoptions, vfoptions);

fprintf('GE vs no-PType (%s, fastOLG=%i), N_i=2 identical types, PricePath r, this should be zero: %.3e \n',combo,fastOLG,max(abs(PricePath_PT2.r(:)-PricePath_NoPT.r(:))))
fprintf('GE vs no-PType (%s, fastOLG=%i), N_i=2 identical types, GEcondnPath, this should be zero: %.3e \n',combo,fastOLG,max(abs(GEcondnPath_PT2.CapitalMarket(:)-GEcondnPath_NoPT.CapitalMarket(:))))

%% (3) N_i=2 different types: in-path general eqm conditions vs the standalone commands
Names_i={'typeA','typeB'};
Params_PT=Params;
Params_PT.kappa_j=struct();
Params_PT.kappa_j.typeA=Params.kappa_j;
Params_PT.kappa_j.typeB=1.3*Params.kappa_j;
V_final_PT=struct(); Policy_final_PT=struct(); AgentDist_initial_PT=struct();
for ii=1:2
    Params_ii=Params;
    Params_ii.kappa_j=Params_PT.kappa_j.(Names_i{ii});
    Params_ii.r=PricePath.r(1);
    Params_ii.w=ParamPath.w(1);
    [~,Policy_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial_PT.(Names_i{ii})=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    Params_ii.r=PricePath.r(end);
    Params_ii.w=ParamPath.w(end);
    [V_final_PT.(Names_i{ii}),Policy_final_PT.(Names_i{ii})]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
end
AgentDist_initial_PT.ptweights=[0.6;0.4];

transpathoptions3=transpathoptions;
transpathoptions3.maxiter=1; % GEcondnPath is the conditions at the path handed in
[~,GEcondnPath_PT]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, AgentDist_initial_PT, jequaloneDist, n_d, n_a, n_z, N_j, Names_i, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params_PT, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptions3, simoptions, vfoptions);

% The standalone commands at the same (input) price path
transpathoptions_sa=struct();
transpathoptions_sa.fastOLG=fastOLG;
[~,PolicyPath_PT]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT, Policy_final_PT, Params_PT, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);

GEcondn_byhand=zeros(1,T-1);
for tt=1:T-1
    GEcondn_byhand(tt)=GeneralEqmEqns.CapitalMarket(PricePath.r(tt),AggVarsPath_PT.assets.Mean(tt),AggVarsPath_PT.earnings.Mean(tt));
end
fprintf('GE in-path vs standalone (%s, fastOLG=%i), two different types, GEcondnPath, this should be zero: %.3e \n',combo,fastOLG,max(abs(gather(GEcondnPath_PT.CapitalMarket(1:T-1))-gather(GEcondn_byhand))))

output=struct();

end
