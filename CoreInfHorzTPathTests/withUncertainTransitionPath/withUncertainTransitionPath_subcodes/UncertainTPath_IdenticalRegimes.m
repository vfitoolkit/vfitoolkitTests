function output=UncertainTPath_IdenticalRegimes(T,PricePath,ParamPath,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptions,vfoptionsbaseline,simoptionsbaseline)
% Two (or three) regimes with IDENTICAL parameters, and full support over a window of three periods, so the tree has
% 8 (or 27) histories which are all the same economy: every node must equal the deterministic path in its period.
% This is the exact oracle for the expectation at a branching node (its weights sum to one), for the forward pass
% (each history carries the whole distribution, it is not split by probability), and for the numbering, on a tree
% that really branches. It is run on each kind of model the deterministic commands handle, since the tree reaches
% all of them through the deterministic commands:
%   nod z, d z, nod noz, nod z e, nod z with the grid interpolation layer, d z with divide-and-conquer, nod z with N_S=3
% The probabilities are DYADIC (0.5; and 0.5, 0.25, 0.25) so the expectation over identical children is the child's
% value exactly in floating point: sum_i p_i*V=V holds exactly only then, and a ULP could flip a Policy tie.
% Finally, for scale, the same tree with different regimes must NOT equal the deterministic path.

T1=3; T2=5;

for vv=1:7
    vfoptions=struct();
    simoptions=struct();
    FnsToEvaluate=struct();
    n_dV=0; d_gridV=[]; n_zV=n_z; z_gridV=z_grid; pi_zV=pi_z;
    pi_S=[0.5,0.5;0.5,0.5];
    if vv==1
        lab='nod z';
        ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.nextassets=@(aprime,a,z) aprime;
    elseif vv==2
        lab='d z';
        n_dV=n_d; d_gridV=d_grid;
        ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*d*z;
    elseif vv==3
        lab='nod noz';
        n_zV=0; z_gridV=[]; pi_zV=[];
        ReturnFn=@(aprime,a,r,w,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.nextassets=@(aprime,a) aprime;
    elseif vv==4
        lab='nod z e';
        vfoptions.n_e=vfoptionsbaseline.n_e; vfoptions.e_grid=vfoptionsbaseline.e_grid; vfoptions.pi_e=vfoptionsbaseline.pi_e;
        simoptions.n_e=simoptionsbaseline.n_e; simoptions.e_grid=simoptionsbaseline.e_grid; simoptions.pi_e=simoptionsbaseline.pi_e;
        ReturnFn=@(aprime,a,z,e,r,w,sigma) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,z,e) a;
        FnsToEvaluate.eval=@(aprime,a,z,e) e;
    elseif vv==5
        lab='nod z GI';
        vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
        simoptions.gridinterplayer=1; simoptions.ngridinterp=5;
        ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.nextassets=@(aprime,a,z) aprime;
    elseif vv==6
        lab='d z DC';
        n_dV=n_d; d_gridV=d_grid;
        vfoptions.divideandconquer=1; vfoptions.level1n=11;
        ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*d*z;
    elseif vv==7
        lab='nod z, N_S=3';
        pi_S=[0.5,0.25,0.25;0.25,0.5,0.25;0.25,0.25,0.5];
        ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.nextassets=@(aprime,a,z) aprime;
    end
    N_S=size(pi_S,1);

    % The same economy in every S. r is given N_S-by-T (identical rows), w as 1-by-T, so both readings are used.
    PricePathU.r=repmat(PricePath.r,N_S,1);
    ParamPathU.w=ParamPath.w;

    ParamsF=Params; ParamsF.r=PricePath.r(T); ParamsF.w=ParamPath.w(T);
    [V_finalD,Policy_finalD]=ValueFnIter_InfHorz(n_dV,n_a,n_zV,d_gridV,a_grid,z_gridV,pi_zV,ReturnFn,ParamsF,DiscountFactorParamNames,[],vfoptions);
    V_final=struct(); Policy_final=struct();
    for SS=1:N_S
        V_final.(['S',num2str(SS)])=V_finalD;
        Policy_final.(['S',num2str(SS)])=Policy_finalD;
    end
    ParamsI=Params; ParamsI.r=PricePath.r(1); ParamsI.w=ParamPath.w(1);
    [~,Policy_initial]=ValueFnIter_InfHorz(n_dV,n_a,n_zV,d_gridV,a_grid,z_gridV,pi_zV,ReturnFn,ParamsI,DiscountFactorParamNames,[],vfoptions);
    AgentDist_initial=StationaryDist_InfHorz(Policy_initial,n_dV,n_a,n_zV,pi_zV,simoptions,ParamsI,[]);

    % Deterministic
    [VPathD,PolicyPathD]=ValueFnOnTransPath_InfHorz(PricePath, ParamPath, T, V_finalD, Policy_finalD, Params, n_dV,n_a,n_zV, d_gridV,a_grid,z_gridV, pi_zV, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    AgentDistPathD=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePath, ParamPath, PolicyPathD, n_dV,n_a,n_zV,pi_zV,T,Params,transpathoptions,simoptions);
    AggVarsPathD=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathD, PolicyPathD, PricePath, ParamPath, Params, T, n_dV, n_a, n_zV, d_gridV, a_grid, z_gridV, transpathoptions, simoptions);

    % On the tree
    [VPathU,PolicyPathU]=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_dV,n_a,n_zV, d_gridV,a_grid,z_gridV, pi_zV, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
    AgentDistPathU=AgentDistOnUncertainTransPath_InfHorz(AgentDist_initial, PricePathU, ParamPathU, PolicyPathU, T, T1, T2, pi_S, 1, n_dV,n_a,n_zV,pi_zV, Params, transpathoptions, simoptions);
    AggVarsTreeU=EvalFnOnUncertainTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathU, PolicyPathU, PricePathU, ParamPathU, Params, T, T1, T2, pi_S, 1, n_dV,n_a,n_zV, d_gridV,a_grid,z_gridV, transpathoptions, simoptions);

    Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
    nNodes=Tree.nNodes;
    fprintf('Identical regimes %s: the tree has N_S^3 histories, should be one: %i \n',lab,Tree.nHist==N_S^3)

    % Node k against period nodedate(k) of the deterministic path
    A=gather(reshape(VPathU,[],nNodes)); B=gather(reshape(VPathD,[],T)); B=B(:,Tree.nodedate);
    dd=abs(A-B); dd(A==B)=0; dd(isnan(dd))=Inf;
    fprintf('Identical regimes %s: V on every node vs the deterministic path in its period, this should be zero: %.3e \n',lab,max(dd,[],'all'))
    A=gather(reshape(PolicyPathU,[],nNodes)); B=gather(reshape(PolicyPathD,[],T)); B=B(:,Tree.nodedate);
    fprintf('Identical regimes %s: Policy on every node, this should be zero: %.3e \n',lab,max(abs(A-B),[],'all'))
    A=gather(reshape(AgentDistPathU,[],nNodes)); B=gather(reshape(AgentDistPathD,[],T)); B=B(:,Tree.nodedate);
    fprintf('Identical regimes %s: AgentDist on every node (each history carries the whole distribution), this should be zero: %.3e \n',lab,max(abs(A-B),[],'all'))
    aggerr=0;
    AggNames=fieldnames(FnsToEvaluate);
    for ff=1:length(AggNames)
        B=gather(AggVarsPathD.(AggNames{ff}).Mean(:))';
        aggerr=max(aggerr,max(abs(gather(AggVarsTreeU.(AggNames{ff}).Mean)-B(Tree.nodedate'))));
    end
    fprintf('Identical regimes %s: AggVars on every node, this should be zero: %.3e \n',lab,aggerr)

    if vv==1
        % For scale: regimes that differ (S=2 has a 10 percent lower wage) must NOT give the deterministic path
        ParamPathX.w=[ParamPath.w; 0.9*ParamPath.w];
        V_finalX=V_final;
        ParamsF2=ParamsF; ParamsF2.w=0.9*ParamPath.w(T);
        [V_finalX.S2,Policy_finalX2]=ValueFnIter_InfHorz(n_dV,n_a,n_zV,d_gridV,a_grid,z_gridV,pi_zV,ReturnFn,ParamsF2,DiscountFactorParamNames,[],vfoptions);
        Policy_finalX=Policy_final; Policy_finalX.S2=Policy_finalX2;
        VPathX=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathX, T, T1, T2, pi_S, 1, V_finalX, Policy_finalX, Params, n_dV,n_a,n_zV, d_gridV,a_grid,z_gridV, pi_zV, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
        A=gather(reshape(VPathX,[],nNodes)); B=gather(reshape(VPathD,[],T)); B=B(:,Tree.nodedate);
        dd=abs(A-B); dd(A==B)=0; dd(isnan(dd))=0;
        fprintf('Identical regimes %s: for scale, with the S=2 wage 10 percent lower, V vs the deterministic path, should NOT be zero: %.3e \n',lab,max(dd,[],'all'))
        % And in period 1, before any S is drawn, the expectation of the lower wage already lowers V
        fprintf('Identical regimes %s: for scale, V in period 1 is lower than the deterministic one when S=2 is possible, should be one: %i \n',lab,all(A(:,1)<=B(:,1)) && any(A(:,1)<B(:,1)))
    end
end

output=struct();

end
