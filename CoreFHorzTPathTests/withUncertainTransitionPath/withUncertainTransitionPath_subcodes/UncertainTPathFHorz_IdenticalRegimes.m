function output=UncertainTPathFHorz_IdenticalRegimes(T,PricePath,ParamPath,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline,vfoptionsbaseline)
% Identical regimes: two S with the same prices and parameters, on a full-support tree, so every history is the
% deterministic path. Every node's V, Policy, AgentDist and AggVars must equal the deterministic path's in its period, as
% identities. Probabilities are dyadic (0.5), so that the expectation over two identical children is exactly the child.
% Written 2026-11-05, before the first run.
% Cases: nod z (fastOLG=1 and 0), d z, nod noz, nod z e, nod z with the grid interpolation layer, nod z with
% divide-and-conquer.

pi_S=[0.5,0.5;0.5,0.5]; T1=3; T2=5;
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
fprintf('Uncertain FHorz identical regimes: %i histories, %i nodes \n',Tree.nHist,Tree.nNodes)
casenames={'nod z','nod z fastOLG=0','d z','nod noz','nod z e','nod z GI','nod z DC'};
for cc=1:length(casenames)
    vfoptions=struct();
    simoptions=struct();
    transpathoptions=transpathoptionsbaseline;
    n_d_c=0; d_grid_c=[];
    n_z_c=n_z; z_grid_c=z_grid; pi_z_c=pi_z;
    jequaloneDist=zeros([n_a,n_z]);
    jequaloneDist(1,ceil(n_z/2))=1;
    ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
    clear FnsToEvaluate
    FnsToEvaluate.assets=@(aprime,a,z) a;
    FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
    switch casenames{cc}
        case 'nod z fastOLG=0'
            transpathoptions.fastOLG=0;
        case 'd z'
            n_d_c=n_d; d_grid_c=d_grid;
            ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);
            clear FnsToEvaluate
            FnsToEvaluate.assets=@(d,aprime,a,z) a;
            FnsToEvaluate.hours=@(d,aprime,a,z) d;
        case 'nod noz'
            n_z_c=0; z_grid_c=[]; pi_z_c=[];
            jequaloneDist=zeros([n_a,1]);
            jequaloneDist(1)=1;
            ReturnFn=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
            clear FnsToEvaluate
            FnsToEvaluate.assets=@(aprime,a) a;
        case 'nod z e'
            vfoptions.n_e=vfoptionsbaseline.n_e; vfoptions.e_grid=vfoptionsbaseline.e_grid; vfoptions.pi_e=vfoptionsbaseline.pi_e;
            simoptions.n_e=vfoptions.n_e; simoptions.e_grid=vfoptions.e_grid; simoptions.pi_e=vfoptions.pi_e;
            jequaloneDist=zeros([n_a,n_z,vfoptions.n_e]);
            jequaloneDist(1,ceil(n_z/2),:)=vfoptions.pi_e;
            ReturnFn=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
            clear FnsToEvaluate
            FnsToEvaluate.assets=@(aprime,a,z,e) a;
        case 'nod z GI'
            vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
            simoptions.gridinterplayer=1; simoptions.ngridinterp=5;
        case 'nod z DC'
            vfoptions.divideandconquer=1;
    end
    AggNames=fieldnames(FnsToEvaluate);
    % The deterministic path, and its final eqm and initial dist
    Params_T=Params; Params_T.r=PricePath.r(T); Params_T.w=ParamPath.w(T);
    [V_fin,Policy_fin]=ValueFnIter_Case1_FHorz(n_d_c,n_a,n_z_c,N_j,d_grid_c,a_grid,z_grid_c,pi_z_c,ReturnFn,Params_T,DiscountFactorParamNames,[],vfoptions);
    Params_1=Params; Params_1.r=PricePath.r(1); Params_1.w=ParamPath.w(1);
    [~,Policy_1]=ValueFnIter_Case1_FHorz(n_d_c,n_a,n_z_c,N_j,d_grid_c,a_grid,z_grid_c,pi_z_c,ReturnFn,Params_1,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_1,n_d_c,n_a,n_z_c,N_j,pi_z_c,Params_1,simoptions);
    [VD,PD]=ValueFnOnTransPath_Case1_FHorz(PricePath, ParamPath, T, V_fin, Policy_fin, Params, n_d_c, n_a, n_z_c, N_j, d_grid_c, a_grid,z_grid_c, pi_z_c, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    DD=AgentDistOnTransPath_Case1_FHorz(AgentDist_initial, jequaloneDist, PricePath, ParamPath, PD, AgeWeightParamNames,n_d_c,n_a,n_z_c,N_j,pi_z_c, T,Params, transpathoptions, simoptions);
    AD=EvalFnOnTransPath_AggVars_Case1_FHorz(FnsToEvaluate, DD, PD, PricePath, ParamPath, Params, T, n_d_c, n_a, n_z_c, N_j, d_grid_c, a_grid,z_grid_c, transpathoptions, simoptions);
    % The tree, with S=1 and S=2 identical
    clear PricePathU ParamPathU
    PricePathU.r=[PricePath.r; PricePath.r];
    ParamPathU.w=[ParamPath.w; ParamPath.w];
    V_final=struct('S1',V_fin,'S2',V_fin);
    Policy_final=struct('S1',Policy_fin,'S2',Policy_fin);
    [VU,PU]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d_c,n_a,n_z_c,N_j, d_grid_c,a_grid,z_grid_c, pi_z_c, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    DU=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d_c,n_a,n_z_c,N_j,pi_z_c, Params, transpathoptions, simoptions);
    AU=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluate, DU, PU, PricePathU, ParamPathU, Params, T, T1, T2, pi_S, 1, n_d_c,n_a,n_z_c,N_j, d_grid_c,a_grid,z_grid_c, transpathoptions, simoptions);
    VU=gather(reshape(VU,[],Tree.nNodes)); VD=gather(reshape(VD,[],T)); VD=VD(:,Tree.nodedate);
    PU=gather(reshape(PU,[],Tree.nNodes)); PD=gather(reshape(PD,[],T)); PD=PD(:,Tree.nodedate);
    DU=gather(reshape(DU,[],Tree.nNodes)); DD=gather(reshape(DD,[],T)); DD=DD(:,Tree.nodedate);
    dd=abs(VU-VD); dd(VU==VD)=0; dd(isnan(dd))=Inf;
    fprintf('Uncertain FHorz identical regimes, %s: every node''s V is the deterministic path''s in its period, this should be zero: %.3e \n',casenames{cc},max(dd,[],'all'))
    fprintf('Uncertain FHorz identical regimes, %s: Policy, this should be zero: %.3e \n',casenames{cc},max(abs(PU-PD),[],'all'))
    fprintf('Uncertain FHorz identical regimes, %s: AgentDist, this should be zero: %.3e \n',casenames{cc},max(abs(DU-DD),[],'all'))
    aggerr=0;
    for ff=1:length(AggNames)
        temp=gather(AD.(AggNames{ff}).Mean(:))';
        aggerr=max(aggerr,max(abs(AU.(AggNames{ff}).Mean-temp(Tree.nodedate'))));
    end
    fprintf('Uncertain FHorz identical regimes, %s: AggVars, this should be zero: %.3e \n',casenames{cc},aggerr)
end

output=struct();

end
