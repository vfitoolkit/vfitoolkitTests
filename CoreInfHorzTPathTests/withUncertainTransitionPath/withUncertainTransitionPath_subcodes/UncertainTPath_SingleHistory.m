function output=UncertainTPath_SingleHistory(T,rS,wS,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptions)
% A tree with a single history is a deterministic transition path, so the uncertain commands must reproduce
% ValueFnOnTransPath_InfHorz, AgentDistOnTransPath_InfHorz and EvalFnOnTransPath_AggVars_InfHorz on the same path
% exactly. Its nodes are numbered 1,...,T, so the comparison is node t against period t.
%   (a) N_S=1, with d (and with the grid interpolation layer)
%   (b) N_S=2 and pi_S=[0,1;0,1] with T1=T2=9 (uncertainty faced in period 9): S switches to 2 in period 10 for sure. The deterministic path has the
%       S=1 values before period 10 and the S=2 values from 10 on, and the S=2 final eqm.
%   (c) the same with a 3-D pi_S: identity, then [0,1;0,1], then identity, over T1=3..T2=5, so the switch is in period 5
%   (d) deterministic stepping across THREE states: N_S=3, T1=5, T2=7, the transition from 5 to 6 takes S=1 to 2, from
%       6 to 7 is the identity, and from 7 to 8 takes S=2 to 3 (the rows a history is never in are identity rows). So
%       S=1 up to period 5, S=2 in periods 6 and 7, S=3 from 8 on; against the deterministic path on that S path, with
%       the S=3 final eqm.

vfoptions=struct();
simoptions=struct();
ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
FnsToEvaluate.assets=@(d,aprime,a,z) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*d*z;
FnsToEvaluate.nextassets=@(d,aprime,a,z) aprime;

% A third regime, for (d): lower r and wage again
rS=[rS; rS(1,:)-0.01];
wS=[wS; 0.8*wS(1,:)];

for cc=1:5
    if cc==1
        lab='(a) N_S=1';
        pi_S=1; T1=3; T2=5; Spath=ones(1,T);
        vfoptionsC=vfoptions; simoptionsC=simoptions;
    elseif cc==2
        lab='(a) N_S=1, GI';
        pi_S=1; T1=3; T2=5; Spath=ones(1,T);
        vfoptionsC=vfoptions; vfoptionsC.gridinterplayer=1; vfoptionsC.ngridinterp=5;
        simoptionsC=simoptions; simoptionsC.gridinterplayer=1; simoptionsC.ngridinterp=5;
    elseif cc==3
        lab='(b) sure switch, 2-D pi_S';
        pi_S=[0,1;0,1]; T1=9; T2=9; Spath=[ones(1,9),2*ones(1,T-9)];
        vfoptionsC=vfoptions; simoptionsC=simoptions;
    elseif cc==4
        lab='(c) sure switch, 3-D pi_S';
        pi_S=cat(3,eye(2),[0,1;0,1],eye(2)); T1=3; T2=5; Spath=[ones(1,4),2*ones(1,T-4)];
        vfoptionsC=vfoptions; simoptionsC=simoptions;
    elseif cc==5
        lab='(d) stepping 1 to 2 to 3, 3-D pi_S';
        pi_S=cat(3,[0,1,0;0,1,0;0,0,1],eye(3),[1,0,0;0,0,1;0,0,1]); T1=5; T2=7; Spath=[ones(1,5),2,2,3*ones(1,T-7)];
        vfoptionsC=vfoptions; simoptionsC=simoptions;
    end
    N_S=size(pi_S,1);

    % The prices and parameters on the tree: N_S-by-T (with N_S=1, the S=1 rows)
    PricePathU.r=rS(1:N_S,:);
    ParamPathU.w=wS(1:N_S,:);
    % The same as a deterministic path: the value of S_t in each period t
    PricePathD.r=rS(Spath+3*(0:1:T-1));
    ParamPathD.w=wS(Spath+3*(0:1:T-1));
    Sfinal=Spath(T);

    % The final eqm of each S (the final value of its path), and the initial agent dist (the stationary dist of S=1)
    V_final=struct(); Policy_final=struct();
    for SS=1:N_S
        ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
        [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptionsC);
    end
    ParamsS=Params; ParamsS.r=rS(1,1); ParamsS.w=wS(1,1);
    [~,Policy_initial]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptionsC);
    AgentDist_initial=StationaryDist_InfHorz(Policy_initial,n_d,n_a,n_z,pi_z,simoptionsC,ParamsS,[]);

    % Deterministic
    [VPathD,PolicyPathD]=ValueFnOnTransPath_InfHorz(PricePathD, ParamPathD, T, V_final.(['S',num2str(Sfinal)]), Policy_final.(['S',num2str(Sfinal)]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsC);
    AgentDistPathD=AgentDistOnTransPath_InfHorz(AgentDist_initial, PricePathD, ParamPathD, PolicyPathD, n_d,n_a,n_z,pi_z,T,Params,transpathoptions,simoptionsC);
    AggVarsPathD=EvalFnOnTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathD, PolicyPathD, PricePathD, ParamPathD, Params, T, n_d, n_a, n_z, d_grid, a_grid, z_grid, transpathoptions, simoptionsC);

    % On the tree
    [VPathU,PolicyPathU]=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsC);
    AgentDistPathU=AgentDistOnUncertainTransPath_InfHorz(AgentDist_initial, PricePathU, ParamPathU, PolicyPathU, T, T1, T2, pi_S, 1, n_d,n_a,n_z,pi_z, Params, transpathoptions, simoptionsC);
    AggVarsTreeU=EvalFnOnUncertainTransPath_AggVars_InfHorz(FnsToEvaluate, AgentDistPathU, PolicyPathU, PricePathU, ParamPathU, Params, T, T1, T2, pi_S, 1, n_d,n_a,n_z, d_grid,a_grid,z_grid, transpathoptions, simoptionsC);

    Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
    fprintf('Single history %s: the tree has one history and T nodes, should be one: %i \n',lab,Tree.nHist==1 && Tree.nNodes==T)
    fprintf('Single history %s: the output shapes are the deterministic ones, should be one: %i \n',lab,isequal(size(VPathU),size(VPathD)) && isequal(size(PolicyPathU),size(PolicyPathD)) && isequal(size(AgentDistPathU),size(AgentDistPathD)))

    A=gather(VPathU(:)); B=gather(VPathD(:));
    dd=abs(A-B); dd(A==B)=0; dd(isnan(dd))=Inf;
    fprintf('Single history %s: V vs the deterministic path, this should be zero: %.3e \n',lab,max(dd))
    fprintf('Single history %s: Policy vs the deterministic path, this should be zero: %.3e \n',lab,max(abs(gather(PolicyPathU(:))-gather(PolicyPathD(:)))))
    fprintf('Single history %s: AgentDist vs the deterministic path, this should be zero: %.3e \n',lab,max(abs(gather(AgentDistPathU(:))-gather(AgentDistPathD(:)))))
    aggerr=0;
    AggNames=fieldnames(FnsToEvaluate);
    for ff=1:length(AggNames)
        aggerr=max(aggerr,max(abs(gather(AggVarsTreeU.(AggNames{ff}).Mean(:))-gather(AggVarsPathD.(AggNames{ff}).Mean(:)))));
    end
    fprintf('Single history %s: AggVars vs the deterministic path, this should be zero: %.3e \n',lab,aggerr)
end

output=struct();

end
