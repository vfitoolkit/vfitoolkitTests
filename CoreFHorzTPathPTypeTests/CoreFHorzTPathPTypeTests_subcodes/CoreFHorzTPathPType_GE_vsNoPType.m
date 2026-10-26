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
[V_final,Policy_final]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_final,DiscountFactorParamNames,[],vfoptions);

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

%% (1A)/(2A) The same with Anderson acceleration (GEnewprice=2)
% Anderson only mixes once it has a history: iteration 1 is always a plain shooting step, and by
% default (anderson.warmup=2) so is iteration 2, so at the maxiter=2 used above Anderson would be the
% shooting algorithm and these checks would test nothing new. warmup=0 and maxiter=3 make iterations 2
% and 3 Anderson steps. The safeguard is off (the transition path default, set here anyway) so that those
% steps are always taken, as a rejected step is a plain shooting step; it is tested in
% with_CoreFHorzTPathPTypeAlgoTests.
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=transpathoptions.GEnewprice3.howtoupdate;
transpathoptionsA.maxiter=3;
transpathoptionsA.anderson.warmup=0;
transpathoptionsA.anderson.safeguard=0;
[PricePathA_NoPT,GEcondnPathA_NoPT]=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions);
ParamsA1=Params;
ParamsA1.ptweights_GEtest=1;
[PricePathA_PT1,GEcondnPathA_PT1]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT1, AgentDist_initial_PT1, jequaloneDist, n_d, n_a, n_z, N_j, 1, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsA1, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
[PricePathA_PT2,GEcondnPathA_PT2]=TransitionPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT2, AgentDist_initial_PT2, jequaloneDist, n_d, n_a, n_z, N_j, 2, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, PTypeDistParamNames, transpathoptionsA, simoptions, vfoptions);
% The shooting algorithm for the same number of iterations, to show the Anderson steps were taken
transpathoptionsS3=transpathoptions;
transpathoptionsS3.maxiter=3;
PricePathS3_NoPT=TransitionPath_Case1_FHorz(PricePath, ParamPath, T, V_final, AgentDist_initial, jequaloneDist, n_d, n_a, n_z, N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS3, simoptions, vfoptions);

fprintf('GE vs no-PType (%s, fastOLG=%i), Anderson at maxiter=3, N_i=1, PricePath r, this should be zero: %.3e \n',combo,fastOLG,max(abs(PricePathA_PT1.r(:)-PricePathA_NoPT.r(:))))
fprintf('GE vs no-PType (%s, fastOLG=%i), Anderson at maxiter=3, N_i=1, GEcondnPath, this should be zero: %.3e \n',combo,fastOLG,max(abs(GEcondnPathA_PT1.CapitalMarket(:)-GEcondnPathA_NoPT.CapitalMarket(:))))
fprintf('GE vs no-PType (%s, fastOLG=%i), Anderson at maxiter=3, N_i=2 identical types, PricePath r, this should be zero: %.3e \n',combo,fastOLG,max(abs(PricePathA_PT2.r(:)-PricePathA_NoPT.r(:))))
fprintf('GE vs no-PType (%s, fastOLG=%i), Anderson at maxiter=3, N_i=2 identical types, GEcondnPath, this should be zero: %.3e \n',combo,fastOLG,max(abs(GEcondnPathA_PT2.CapitalMarket(:)-GEcondnPathA_NoPT.CapitalMarket(:))))
fprintf('GE vs no-PType (%s, fastOLG=%i), Anderson vs shooting, both at maxiter=3 without permanent types, PricePath r, this should NOT be zero: %.3e \n',combo,fastOLG,max(abs(PricePathA_NoPT.r(:)-PricePathS3_NoPT.r(:))))

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
AgentDistPath_PT=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT, jequaloneDist, PricePath, ParamPath, PolicyPath_PT, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,Names_i,pi_z, T,Params_PT, transpathoptions_sa, simoptions);
AggVarsPath_PT=EvalFnOnTransPath_AggVars_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);

GEcondn_byhand=zeros(1,T-1);
for tt=1:T-1
    GEcondn_byhand(tt)=GeneralEqmEqns.CapitalMarket(PricePath.r(tt),AggVarsPath_PT.assets.Mean(tt),AggVarsPath_PT.earnings.Mean(tt));
end
fprintf('GE in-path vs standalone (%s, fastOLG=%i), two different types, GEcondnPath, this should be zero: %.3e \n',combo,fastOLG,max(abs(gather(GEcondnPath_PT.CapitalMarket(1:T-1))-gather(GEcondn_byhand))))

%% Statistics along the path: EvalFnOnTransPath_AllStats_Case1_FHorz_PType and LifeCycleProfiles_TransPath_FHorz_Case1_PType
% Only for fastOLG=1: the stats commands are handed the agent dist and policy paths, which are laid out
% the same whatever fastOLG they were computed with, so running them again for fastOLG=0 adds nothing.
%   (PA1)/(PL1) N_i=1, and N_i=2 identical types, vs the commands without permanent types on the no-PType
%               paths: the pooled stats, and those of each ptype, every stat; the life-cycle profiles, all
%               seven outputs, both cohort layouts, with NaN in exactly the same places
%   (PA2)/(PL2) two different types: each ptype vs the commands without permanent types on that ptype's
%               own paths (its own kappa_j); and everything (pooled and by ptype) vs the stationary PType
%               commands period by period. The pooled Gini must differ from that of a single ptype, or
%               the pooled stats were not pooled.
%   (PA3)       two different types: the pooled Mean vs EvalFnOnTransPath_AggVars_Case1_FHorz_PType
% The comparison walks every field of the reference output, to any depth (MoreInequality is a struct, the
% comments are text), and finds the same field in the output being checked: for a ptype, the ptype name
% follows the FnsToEvaluate name. A field the reference has and the output does not is counted with the
% NaN mismatches.
if fastOLG==1
    statsname=combo;
    OutputTypes={'Mean','Median','Variance','LorenzCurve','Gini','QuantileCutoffs','QuantileMeans'};
    fnnames=fieldnames(FnsToEvaluate);

    % The standalone paths: without permanent types, N_i=1, N_i=2 identical types (two different types are above)
    [~,PolicyPath_NoPT]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_final, Policy_final, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_NoPT=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PolicyPath_NoPT, AgeWeightParamNames,n_d,n_a,n_z,N_j,pi_z, T,Params, transpathoptions_sa, simoptions);
    Policy_final_PT1.ptype001=Policy_final;
    [~,PolicyPath_PT1]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT1, Policy_final_PT1, ParamsA1, n_d, n_a, n_z, N_j, 1, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_PT1=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT1, jequaloneDist, PricePath, ParamPath, PolicyPath_PT1, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,1,pi_z, T,ParamsA1, transpathoptions_sa, simoptions);
    Policy_final_PT2.ptype001=Policy_final;
    Policy_final_PT2.ptype002=Policy_final;
    [~,PolicyPath_PT2]=ValueFnOnTransPath_Case1_FHorz_PType(PricePath, ParamPath, T, V_final_PT2, Policy_final_PT2, Params, n_d, n_a, n_z, N_j, 2, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions_sa, vfoptions);
    AgentDistPath_PT2=AgentDistOnTransPath_Case1_FHorz_PType(AgentDist_initial_PT2, jequaloneDist, PricePath, ParamPath, PolicyPath_PT2, AgeWeightParamNames, PTypeDistParamNames,n_d,n_a,n_z,N_j,2,pi_z, T,Params, transpathoptions_sa, simoptions);

    % The commands without permanent types: on the no-PType paths, and on each of the two different types' own paths
    AllStats_NoPT=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPath_NoPT, PolicyPath_NoPT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    LCP_NoPT=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPath_NoPT, PolicyPath_NoPT, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);
    AllStats_solo=struct(); LCP_solo=struct();
    for ii=1:2
        Params_ii=Params;
        Params_ii.kappa_j=Params_PT.kappa_j.(Names_i{ii});
        AllStats_solo.(Names_i{ii})=EvalFnOnTransPath_AllStats_Case1_FHorz(FnsToEvaluate, AgentDistPath_PT.(Names_i{ii}), PolicyPath_PT.(Names_i{ii}), PricePath, ParamPath, Params_ii, T, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
        LCP_solo.(Names_i{ii})=LifeCycleProfiles_TransPath_FHorz_Case1(FnsToEvaluate, AgentDistPath_PT.(Names_i{ii}), PolicyPath_PT.(Names_i{ii}), PricePath, ParamPath, Params_ii, T, n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, transpathoptions_sa, simoptions);
    end

    % The commands with permanent types
    AllStats_PT1=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT1, PolicyPath_PT1, PricePath, ParamPath, ParamsA1, T, n_d, n_a, n_z, N_j, 1, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    LCP_PT1=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate, AgentDistPath_PT1, PolicyPath_PT1, PricePath, ParamPath, ParamsA1, T, n_d, n_a, n_z, N_j, 1, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    AllStats_PT2=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT2, PolicyPath_PT2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, 2, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    LCP_PT2=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate, AgentDistPath_PT2, PolicyPath_PT2, PricePath, ParamPath, Params, T, n_d, n_a, n_z, N_j, 2, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    AllStats_PT=EvalFnOnTransPath_AllStats_Case1_FHorz_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);
    LCP_PT=LifeCycleProfiles_TransPath_FHorz_Case1_PType(FnsToEvaluate, AgentDistPath_PT, PolicyPath_PT, PricePath, ParamPath, Params_PT, T, n_d, n_a, n_z, N_j, Names_i, d_grid, a_grid,z_grid, transpathoptions_sa, simoptions);

    % (PA1)/(PA2) AllStats: {label, reference (no-PType layout), output with permanent types, ptype name ('' for pooled)}
    cmpcases={'(PA1) N_i=1, pooled', AllStats_NoPT, AllStats_PT1, '';
              '(PA1) N_i=1, ptype001', AllStats_NoPT, AllStats_PT1, 'ptype001';
              '(PA1) N_i=2 identical types, pooled', AllStats_NoPT, AllStats_PT2, '';
              '(PA1) N_i=2 identical types, ptype001', AllStats_NoPT, AllStats_PT2, 'ptype001';
              '(PA1) N_i=2 identical types, ptype002', AllStats_NoPT, AllStats_PT2, 'ptype002';
              '(PA2) two different types, typeA vs its own paths', AllStats_solo.typeA, AllStats_PT, 'typeA';
              '(PA2) two different types, typeB vs its own paths', AllStats_solo.typeB, AllStats_PT, 'typeB'};
    for cc=1:size(cmpcases,1)
        REF=cmpcases{cc,2};
        OUT=cmpcases{cc,3};
        nm=cmpcases{cc,4};
        maxdiff=0; nanmismatch=0;
        queue={};
        REFnames=fieldnames(REF);
        for kk=1:length(REFnames)
            queue{end+1}=REFnames(kk);
        end
        while ~isempty(queue)
            fieldpath=queue{1};
            queue(1)=[];
            O=getfield(REF,fieldpath{:});
            if isstruct(O)
                Onames=fieldnames(O);
                for kk=1:length(Onames)
                    queue{end+1}=[fieldpath,Onames(kk)];
                end
            else
                % the same field in OUT (for a ptype, the ptype name follows the FnsToEvaluate name)
                if isempty(nm)
                    Ppath=fieldpath;
                else
                    Ppath=[fieldpath(1),{nm},fieldpath(2:end)];
                end
                P=OUT;
                found=1;
                for kk=1:length(Ppath)
                    if isstruct(P) && isfield(P,Ppath{kk})
                        P=P.(Ppath{kk});
                    else
                        found=0;
                        break
                    end
                end
                if found==0
                    nanmismatch=nanmismatch+1;
                elseif iscell(O) % text (the comments): must be the same text
                    nanmismatch=nanmismatch+(~isequal(P,O));
                else
                    maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                    nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
                end
            end
        end
        fprintf('PType stats, %s: %s, EvalFnOnTransPath_AllStats_PType vs without permanent types, every stat, this should be zero: %.3e \n',statsname,cmpcases{cc,1},maxdiff)
        fprintf('PType stats, %s: %s, EvalFnOnTransPath_AllStats_PType vs without permanent types, every stat, positions NaN (or text) in one but not the other, or missing, this should be zero: %i \n',statsname,cmpcases{cc,1},nanmismatch)
    end

    % (PL1)/(PL2) LifeCycleProfiles: the same cases
    cmpcases={'(PL1) N_i=1, pooled', LCP_NoPT, LCP_PT1, '';
              '(PL1) N_i=1, ptype001', LCP_NoPT, LCP_PT1, 'ptype001';
              '(PL1) N_i=2 identical types, pooled', LCP_NoPT, LCP_PT2, '';
              '(PL1) N_i=2 identical types, ptype001', LCP_NoPT, LCP_PT2, 'ptype001';
              '(PL1) N_i=2 identical types, ptype002', LCP_NoPT, LCP_PT2, 'ptype002';
              '(PL2) two different types, typeA vs its own paths', LCP_solo.typeA, LCP_PT, 'typeA';
              '(PL2) two different types, typeB vs its own paths', LCP_solo.typeB, LCP_PT, 'typeB'};
    for cc=1:size(cmpcases,1)
        REF=cmpcases{cc,2};
        OUT=cmpcases{cc,3};
        nm=cmpcases{cc,4};
        maxdiff=0; nanmismatch=0;
        for ff=1:length(fnnames)
            if isempty(nm)
                OUTff=OUT.(fnnames{ff});
            else
                OUTff=OUT.(fnnames{ff}).(nm);
            end
            for oo=1:length(OutputTypes)
                PI=gather(OUTff.initialages.(OutputTypes{oo}));
                PB=gather(OUTff.bornduringtranstion.(OutputTypes{oo}));
                OI=gather(REF.(fnnames{ff}).initialages.(OutputTypes{oo}));
                OB=gather(REF.(fnnames{ff}).bornduringtranstion.(OutputTypes{oo}));
                maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
            end
        end
        fprintf('PType stats, %s: %s, LifeCycleProfiles_TransPath_PType vs without permanent types, all seven outputs, initialages and bornduringtranstion, this should be zero: %.3e \n',statsname,cmpcases{cc,1},maxdiff)
        fprintf('PType stats, %s: %s, LifeCycleProfiles_TransPath_PType vs without permanent types, entries NaN in one but not the other (the alive pattern), this should be zero: %i \n',statsname,cmpcases{cc,1},nanmismatch)
    end

    % (PA2)/(PL2) two different types: everything (pooled and by ptype) vs the stationary PType commands period by period
    PricePathNamesS=fieldnames(PricePath);
    ParamPathNamesS=fieldnames(ParamPath);
    AScell=cell(1,T);
    LCPcell=cell(1,T);
    for tt=1:T
        Params_tt=Params_PT;
        for kk=1:length(PricePathNamesS)
            Params_tt.(PricePathNamesS{kk})=PricePath.(PricePathNamesS{kk})(tt);
        end
        for kk=1:length(ParamPathNamesS)
            Params_tt.(ParamPathNamesS{kk})=ParamPath.(ParamPathNamesS{kk})(tt);
        end
        StationaryDist_tt=struct(); Policy_tt=struct();
        for ii=1:2
            sliceAgentDist=repmat({':'},1,ndims(AgentDistPath_PT.(Names_i{ii}))-1);
            slicePolicy=repmat({':'},1,ndims(PolicyPath_PT.(Names_i{ii}))-1);
            StationaryDist_tt.(Names_i{ii})=AgentDistPath_PT.(Names_i{ii})(sliceAgentDist{:},tt);
            Policy_tt.(Names_i{ii})=PolicyPath_PT.(Names_i{ii})(slicePolicy{:},tt);
        end
        StationaryDist_tt.ptweights=AgentDist_initial_PT.ptweights;
        AScell{tt}=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_tt,Policy_tt,FnsToEvaluate,Params_tt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
        LCPcell{tt}=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_tt,Policy_tt,FnsToEvaluate,Params_tt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
    end
    maxdiff=0; nanmismatch=0;
    for tt=1:T
        AS=AScell{tt};
        queue={{}};
        while ~isempty(queue)
            fieldpath=queue{1};
            queue(1)=[];
            if isempty(fieldpath)
                O=AS;
            else
                O=getfield(AS,fieldpath{:});
            end
            if isstruct(O)
                Onames=fieldnames(O);
                for kk=1:length(Onames)
                    queue{end+1}=[fieldpath,Onames(kk)];
                end
            elseif ischar(O) || isstring(O) || iscell(O) % text (the comments): must be the same text
                P=getfield(AllStats_PT,fieldpath{:});
                nanmismatch=nanmismatch+(~isequal(P{tt},O));
            else
                O=gather(O);
                P=getfield(AllStats_PT,fieldpath{:});
                P=P(:,tt);
                maxdiff=max([maxdiff; abs(P(:)-O(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(P(:))~=isnan(O(:)));
            end
        end
    end
    fprintf('PType stats, %s: (PA2) two different types, EvalFnOnTransPath_AllStats_PType vs the stationary PType AllStats period by period, pooled and by ptype, every stat, this should be zero: %.3e \n',statsname,maxdiff)
    fprintf('PType stats, %s: (PA2) two different types, EvalFnOnTransPath_AllStats_PType vs the stationary PType AllStats period by period, pooled and by ptype, positions NaN (or text) in one but not the other, this should be zero: %i \n',statsname,nanmismatch)
    fprintf('PType stats, %s: (PA2) two different types, pooled assets Gini minus that of typeA, max over the path, this should NOT be zero: %.3e \n',statsname,max(abs(AllStats_PT.assets.Gini-AllStats_PT.assets.typeA.Gini)))
    % The life-cycle profiles: the cohort bookkeeping of the per-period stationary profiles, pooled and by ptype
    scalarstat=[1,1,1,0,1,0,0]; % Mean, Median, Variance, Gini are one number per age
    groupnames={'','typeA','typeB'};
    maxdiff=0; nanmismatch=0;
    for gg=1:length(groupnames)
        for ff=1:length(fnnames)
            if isempty(groupnames{gg})
                OUTff=LCP_PT.(fnnames{ff});
            else
                OUTff=LCP_PT.(fnnames{ff}).(groupnames{gg});
            end
            for oo=1:length(OutputTypes)
                PI=gather(OUTff.initialages.(OutputTypes{oo}));
                PB=gather(OUTff.bornduringtranstion.(OutputTypes{oo}));
                OI=nan(size(PI));
                OB=nan(size(PB));
                % Initial cohorts: aged jj in period 1, so aged jj+tt-1 in period tt, alive while that is at most N_j
                for tt=1:min(T,N_j)
                    if isempty(groupnames{gg})
                        temp=gather(LCPcell{tt}.(fnnames{ff}).(OutputTypes{oo}));
                    else
                        temp=gather(LCPcell{tt}.(fnnames{ff}).(groupnames{gg}).(OutputTypes{oo}));
                    end
                    for jj=1:N_j-tt+1
                        if scalarstat(oo)==1
                            OI(jj,jj+tt-1)=temp(jj+tt-1);
                        else
                            OI(:,jj,jj+tt-1)=temp(:,jj+tt-1);
                        end
                    end
                end
                % Cohorts born during the transition: aged 1 in period ii, so aged age in period ii+age-1
                for ii=1:T
                    for age=1:N_j
                        tt=ii+age-1;
                        if tt<=T
                            if isempty(groupnames{gg})
                                temp=gather(LCPcell{tt}.(fnnames{ff}).(OutputTypes{oo}));
                            else
                                temp=gather(LCPcell{tt}.(fnnames{ff}).(groupnames{gg}).(OutputTypes{oo}));
                            end
                            if scalarstat(oo)==1
                                OB(ii,age)=temp(age);
                            else
                                OB(:,ii,age)=temp(:,age);
                            end
                        end
                    end
                end
                maxdiff=max([maxdiff; abs(PI(:)-OI(:)); abs(PB(:)-OB(:))],[],'omitnan');
                nanmismatch=nanmismatch+sum(isnan(PI(:))~=isnan(OI(:)))+sum(isnan(PB(:))~=isnan(OB(:)));
            end
        end
    end
    fprintf('PType stats, %s: (PL2) two different types, LifeCycleProfiles_TransPath_PType vs the per-period stationary PType profiles, pooled and by ptype, all seven outputs, this should be zero: %.3e \n',statsname,maxdiff)
    fprintf('PType stats, %s: (PL2) two different types, LifeCycleProfiles_TransPath_PType vs the per-period stationary PType profiles, entries NaN in one but not the other (the alive pattern), this should be zero: %i \n',statsname,nanmismatch)

    % (PA3) the pooled Mean vs the AggVars command
    fprintf('PType stats, %s: (PA3) two different types, EvalFnOnTransPath_AllStats_PType pooled Mean vs EvalFnOnTransPath_AggVars_PType, this should be zero: %.3e \n',statsname,max([max(abs(AllStats_PT.assets.Mean-gather(AggVarsPath_PT.assets.Mean))),max(abs(AllStats_PT.earnings.Mean-gather(AggVarsPath_PT.earnings.Mean)))]))

    clear AScell LCPcell AllStats_NoPT LCP_NoPT AllStats_solo LCP_solo AllStats_PT1 LCP_PT1 AllStats_PT2 LCP_PT2 AllStats_PT LCP_PT
    clear PolicyPath_NoPT AgentDistPath_NoPT PolicyPath_PT1 AgentDistPath_PT1 PolicyPath_PT2 AgentDistPath_PT2
end

output=struct();

end
