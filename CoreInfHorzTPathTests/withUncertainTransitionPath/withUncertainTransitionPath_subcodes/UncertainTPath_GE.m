function output=UncertainTPath_GE(T,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% UncertainTransitionPath_InfHorz, the general eqm on the tree, by shooting (stage U2). Written 2026-11-01, before
% the first run. The model is the households of CoreInfHorzTPath_localsearch_shooting with a Cobb-Douglas firm:
% K is assets and N is efficiency labour z, r and w clear the capital and labour markets. The regimes differ in TFP,
% firmA, which enters only the general eqm conditions (S=1: 0.5, S=2: 0.49, think of a permanent productivity loss).
% The loss is 2 percent and not 10: the first run (2026-11-01) used 0.45, from whose stationary eqm the default
% stationary solve (fminsearch from the S=1 prices) stopped at a residual of 0.019, and a path whose period T eqm is
% that far off cannot meet the path tolerance.
%   (a) a SURE switch to S=2 in period 2 is a deterministic transition: the tree (one history) against
%       TransitionPath_InfHorz from the same initial guess, with the same final eqm
%   (a2) deterministic stepping across three states (1 to 2 to 3), the same against TransitionPath_InfHorz
%   (b) identical regimes with a temporary TFP boost, on a full-support tree: every history is the deterministic
%       path, so all nodes in a period have the same prices, and they are TransitionPath_InfHorz's
%   (c) an absorbing switch (S=2 for good, from any period 3 to 10, with probability 1/4 each period), as in
%       ABCP2024: converged; the general eqm conditions re-evaluated by hand at the returned prices are within
%       tolerance and are the ones reported; period T is the final eqm of each S; the possibility of S=2 moves
%       the prices before any switch (should NOT be zero)
%   (d) restarting from the solution (transpathoptions.PricePathTree0) converges at once to the same prices
%   (e) refusals, each by its message
% Shooting is a fixed point, so (a) and (b) are equal only up to the convergence tolerance (the two solvers stop at
% different iterations, and evaluate the conditions on CPU and on GPU); they are checked against 10 times
% toleranceGEprices.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluate.K=@(aprime,a,z) a;
FnsToEvaluate.N=@(aprime,a,z) z;
GeneralEqmEqns.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

%% The initial stationary eqm (S=1, firmA=0.5)
Params.r=0.04; Params.w=1;
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptions, vfoptions);
fprintf('Uncertain GE: initial stationary eqm r=%2.8f, w=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_init,n_d,n_a,n_z,pi_z,simoptions,Params,[]);
% The S=2 final eqm, solved here exactly as UncertainTransitionPath_InfHorz will solve it (same start, same options),
% so that the path tolerance can be set above the worse of the two stationary residuals: below it, the conditions near
% period T could never be met and shooting would run to maxiter
Params2=Params; Params2.firmA=0.49;
[~,GEcondns2]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, DiscountFactorParamNames, [], [], [], {'r','w'}, struct(), simoptions, vfoptions);
GEfloor=max([abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket),abs(GEcondns2.CapitalMarket),abs(GEcondns2.LabourMarket)]);
fprintf('Uncertain GE: the S=2 stationary eqm has conditions %.3e %.3e \n',GEcondns2.CapitalMarket,GEcondns2.LabourMarket)
fprintf('Uncertain GE: both stationary eqms are solved to 2e-3, so the path tolerance below is meaningful, should be one: %i \n',GEfloor<=2e-3)

transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=max(5*GEfloor,10^(-4));
tolP=transpathoptions.toleranceGEprices;
fprintf('Uncertain GE: toleranceGEprices %.3e, toleranceGEcondns %.3e \n',transpathoptions.toleranceGEprices,transpathoptions.toleranceGEcondns)

PricePathShaper.r=min((1:1:T)/30,1);
PricePathShaper.w=min((1:1:T)/30,1);

%% (a) A sure switch to S=2 in period 2, against TransitionPath_InfHorz
ParamPathA.firmA=[0.5*ones(1,T); 0.49*ones(1,T)];
pi_SA=[0,1;0,1]; T1=1; T2=1; % uncertainty faced in period 1, so the switch lands in period 2
[PricePathTreeA,SummaryA]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathA, T, T1, T2, pi_SA, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions);
fprintf('Uncertain GE (a): a sure switch is a single history, should be one: %i \n',SummaryA.Tree.nHist==1)
fprintf('Uncertain GE (a): converged before maxiter, should be one: %i \n',SummaryA.iterations<transpathoptions.maxiter)
% The deterministic path, from the same guess (period 1 is S=1, which has no final eqm here, so the guess is the
% initial price there) and with the same final eqm
pf2=SummaryA.finaleqm.S2.p_eqm;
PricePath0D.r=[Params.r, Params.r+(pf2.r-Params.r)*PricePathShaper.r(2:T)];
PricePath0D.w=[Params.w, Params.w+(pf2.w-Params.w)*PricePathShaper.w(2:T)];
ParamPathD.firmA=[0.5, 0.49*ones(1,T-1)];
ParamsD=Params;
PricePathD=TransitionPath_InfHorz(PricePath0D, ParamPathD, T, SummaryA.finaleqm.S2.V_final, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, ParamsD, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);
diffA=max([max(abs(PricePathTreeA.r-gather(PricePathD.r(:))'),[],'all'),max(abs(PricePathTreeA.w-gather(PricePathD.w(:))'),[],'all')]);
fprintf('Uncertain GE (a): the tree vs TransitionPath_InfHorz, prices within 10*toleranceGEprices, should be one: %i \n',diffA<=10*tolP)
fprintf('Uncertain GE (a): for information, that difference is %.3e, after %i tree iterations \n',diffA,SummaryA.iterations)

%% (a2) Deterministic stepping across three states, 1 to 2 to 3, against TransitionPath_InfHorz on that S path
% T1=5, T2=7: from 5 to 6 S=1 goes to 2, from 6 to 7 nothing changes, from 7 to 8 S=2 goes to 3 (the rows a history
% is never in are identity rows). So S=1 up to period 5, S=2 in 6 and 7, S=3 from 8 on. Only S=3 has a final eqm, so
% the guess is the initial price wherever S is 1 or 2.
ParamPathA2.firmA=[0.5*ones(1,T); 0.495*ones(1,T); 0.49*ones(1,T)];
pi_SA2=cat(3,[0,1,0;0,1,0;0,0,1],eye(3),[1,0,0;0,0,1;0,0,1]); T1=5; T2=7;
Spath=[ones(1,5),2,2,3*ones(1,T-7)];
[PricePathTreeA2,SummaryA2]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathA2, T, T1, T2, pi_SA2, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions);
fprintf('Uncertain GE (a2): stepping 1 to 2 to 3 is a single history, should be one: %i \n',SummaryA2.Tree.nHist==1)
fprintf('Uncertain GE (a2): and its S path is 1 to period 5, 2 in 6 and 7, 3 from 8, should be one: %i \n',isequal(SummaryA2.Tree.historyS,Spath))
fprintf('Uncertain GE (a2): converged before maxiter, should be one: %i \n',SummaryA2.iterations<transpathoptions.maxiter)
pf3=SummaryA2.finaleqm.S3.p_eqm;
PricePath0D2.r=Params.r+(Spath==3).*(pf3.r-Params.r).*PricePathShaper.r;
PricePath0D2.w=Params.w+(Spath==3).*(pf3.w-Params.w).*PricePathShaper.w;
ParamPathD2.firmA=ParamPathA2.firmA(Spath+3*(0:1:T-1));
PricePathD2=TransitionPath_InfHorz(PricePath0D2, ParamPathD2, T, SummaryA2.finaleqm.S3.V_final, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);
diffA2=max([max(abs(PricePathTreeA2.r-gather(PricePathD2.r(:))'),[],'all'),max(abs(PricePathTreeA2.w-gather(PricePathD2.w(:))'),[],'all')]);
fprintf('Uncertain GE (a2): the tree vs TransitionPath_InfHorz, prices within 10*toleranceGEprices, should be one: %i \n',diffA2<=10*tolP)
fprintf('Uncertain GE (a2): for information, that difference is %.3e, after %i tree iterations \n',diffA2,SummaryA2.iterations)

%% (b) Identical regimes with a temporary TFP boost, on a full-support tree
boost=0.5*(1+0.05*max(0,1-(1:1:T)/20));
ParamPathB.firmA=[boost; boost];
pi_SB=[0.5,0.5;0.5,0.5]; T1=3; T2=5;
[PricePathTreeB,SummaryB]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathB, T, T1, T2, pi_SB, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions);
TreeB=SummaryB.Tree;
spreadB=0;
for tt=1:T
    nodes_t=find(TreeB.nodedate==tt);
    spreadB=max([spreadB,max(PricePathTreeB.r(nodes_t))-min(PricePathTreeB.r(nodes_t)),max(PricePathTreeB.w(nodes_t))-min(PricePathTreeB.w(nodes_t))]);
end
fprintf('Uncertain GE (b): identical regimes, %i histories, every node in a period has the same prices, this should be zero: %.3e \n',TreeB.nHist,spreadB)
ParamPathBD.firmA=boost;
% the guess the tree used (both S have the same final eqm, which is the initial one up to the stationary solver)
PricePath0BD.r=Params.r+(SummaryB.finaleqm.S1.p_eqm.r-Params.r)*PricePathShaper.r;
PricePath0BD.w=Params.w+(SummaryB.finaleqm.S1.p_eqm.w-Params.w)*PricePathShaper.w;
PricePathBD=TransitionPath_InfHorz(PricePath0BD, ParamPathBD, T, SummaryB.finaleqm.S1.V_final, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, []);
rD=gather(PricePathBD.r(:))'; wD=gather(PricePathBD.w(:))';
diffB=max([max(abs(PricePathTreeB.r-rD(TreeB.nodedate')),[],'all'),max(abs(PricePathTreeB.w-wD(TreeB.nodedate')),[],'all')]);
fprintf('Uncertain GE (b): the tree vs TransitionPath_InfHorz in each period, within 10*toleranceGEprices, should be one: %i \n',diffB<=10*tolP)
fprintf('Uncertain GE (b): for information, that difference is %.3e; the boost moves r by %.3e from the initial eqm \n',diffB,max(abs(rD-Params.r)))

%% (c) An absorbing switch, as in ABCP2024
P=0.75;
pi_SC=[P,1-P;0,1]; T1=2; T2=9; % agents face it in periods 2 to 9, so the switch can land in any of 3 to 10
ParamPathC=ParamPathA;
[PricePathTreeC,SummaryC]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions);
TreeC=SummaryC.Tree;
fprintf('Uncertain GE (c): %i histories, converged before maxiter, should be one: %i \n',TreeC.nHist,SummaryC.iterations<transpathoptions.maxiter)
% The general eqm conditions re-evaluated by hand at the returned prices
[~,PolicyC]=ValueFnOnUncertainTransPath_InfHorz(PricePathTreeC, ParamPathC, T, T1, T2, pi_SC, 1, struct('S1',SummaryC.finaleqm.S1.V_final,'S2',SummaryC.finaleqm.S2.V_final), struct('S1',SummaryC.finaleqm.S1.Policy_final,'S2',SummaryC.finaleqm.S2.Policy_final), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DistC=AgentDistOnUncertainTransPath_InfHorz(AgentDist_initial, PricePathTreeC, ParamPathC, PolicyC, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,pi_z, Params, transpathoptions, simoptions);
AggC=EvalFnOnUncertainTransPath_AggVars_InfHorz(FnsToEvaluate, DistC, PolicyC, PricePathTreeC, ParamPathC, Params, T, T1, T2, pi_SC, 1, n_d,n_a,n_z, d_grid,a_grid,z_grid, transpathoptions, simoptions);
[~,~,ParamPathNodesC]=UncertainTransPathTree(T,T1,T2,pi_SC,1,[],ParamPathC);
GEhand=zeros(2,TreeC.nNodes);
for kk=1:TreeC.nNodes
    r=PricePathTreeC.r(kk); w=PricePathTreeC.w(kk); K=AggC.K.Mean(kk); N=AggC.N.Mean(kk); A=ParamPathNodesC.firmA(kk);
    GEhand(1,kk)=GeneralEqmEqns.CapitalMarket(r,K,N,Params.firmalpha,Params.firmdelta,A);
    GEhand(2,kk)=GeneralEqmEqns.LabourMarket(w,K,N,Params.firmalpha,A);
end
notlast=(TreeC.nodedate<T)';
fprintf('Uncertain GE (c): the general eqm conditions by hand at the returned prices, within toleranceGEcondns at every node before T, should be one: %i \n',max(sqrt(sum(GEhand(:,notlast).^2,1)))<=transpathoptions.toleranceGEcondns)
fprintf('Uncertain GE (c): and they are the ones reported, this should be zero: %.3e \n',max(abs(GEhand-[SummaryC.GEcondnTree.CapitalMarket;SummaryC.GEcondnTree.LabourMarket]),[],'all'))
lastnodes=find(TreeC.nodedate==T)';
heldT=0;
for kk=lastnodes
    pfk=SummaryC.finaleqm.(['S',num2str(TreeC.nodeS(kk))]).p_eqm;
    heldT=max([heldT,abs(PricePathTreeC.r(kk)-pfk.r),abs(PricePathTreeC.w(kk)-pfk.w)]);
end
fprintf('Uncertain GE (c): period T is the final eqm of each S, this should be zero: %.3e \n',heldT)
% Before any switch is possible, the deterministic S=1 path would be the initial eqm (a null reform); the possibility
% of the TFP loss is what moves path 0
hh0=find(all(TreeC.historyS==1,2));
nodes0=TreeC.histories(hh0,:);
fprintf('Uncertain GE (c): for scale, the possibility of S=2 moves the never-switch prices in period 1 away from the initial eqm, should NOT be zero: %.3e \n',max(abs(PricePathTreeC.r(nodes0(1))-Params.r),abs(PricePathTreeC.w(nodes0(1))-Params.w)))
fprintf('Uncertain GE (c): for information, %i iterations, price distance %.3e, general eqm distance %.3e \n',SummaryC.iterations,SummaryC.PricePathDist,SummaryC.GEcondnDist)

%% (d) Restarting from the solution
transpathoptionsD=transpathoptions;
transpathoptionsD.PricePathTree0=PricePathTreeC;
[PricePathTreeC2,SummaryC2]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptionsD, simoptions, vfoptions);
fprintf('Uncertain GE (d): restarted from the solution it converges in one iteration, should be one: %i \n',SummaryC2.iterations==1)
fprintf('Uncertain GE (d): to the same prices, this should be zero: %.3e \n',max([abs(PricePathTreeC2.r-PricePathTreeC.r),abs(PricePathTreeC2.w-PricePathTreeC.w)]))

%% (e) Refusals, each by its message
err=zeros(1,7);
tpE=transpathoptions; tpE=rmfield(tpE,'GEnewprice');
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, tpE, simoptions, vfoptions);
catch ME
    err(1)=contains(ME.message,'GEnewprice must be set');
end
tpE=transpathoptions; tpE.GEnewprice=1;
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, tpE, simoptions, vfoptions);
catch ME
    err(2)=contains(ME.message,'GEnewprice=1 (quasi-Newton) is not implemented');
end
tpE=transpathoptions; tpE.GEnewprice=2;
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, tpE, simoptions, vfoptions);
catch ME
    err(3)=contains(ME.message,'GEnewprice=2 (Anderson acceleration) is not yet implemented');
end
tpE=transpathoptions; tpE.GEnewprice3=rmfield(tpE.GEnewprice3,'howtoupdate');
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, tpE, simoptions, vfoptions);
catch ME
    err(4)=contains(ME.message,'requires transpathoptions.GEnewprice3.howtoupdate');
end
tpE=transpathoptions; tpE.constrainpositive={'w'};
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, tpE, simoptions, vfoptions);
catch ME
    err(5)=contains(ME.message,'price constraints');
end
GeneralEqmEqnsE=GeneralEqmEqns;
GeneralEqmEqnsE.CapitalMarket=@(r,K_tminus1,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K_tminus1^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqnsE, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions);
catch ME
    err(6)=contains(ME.message,'_tminus1 and _tplus1 are not yet implemented');
end
vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptionsE);
catch ME
    err(7)=contains(ME.message,'localsearch is not yet implemented');
end
fprintf('Uncertain GE (e): GEnewprice unset is refused by name, should be one: %i \n',err(1))
fprintf('Uncertain GE (e): GEnewprice=1 is refused by name, should be one: %i \n',err(2))
fprintf('Uncertain GE (e): GEnewprice=2 is refused by name (not yet), should be one: %i \n',err(3))
fprintf('Uncertain GE (e): GEnewprice=3 without howtoupdate is refused by name, should be one: %i \n',err(4))
fprintf('Uncertain GE (e): price constraints are refused by name (not yet), should be one: %i \n',err(5))
fprintf('Uncertain GE (e): a _tminus1 in the transition general eqm conditions is refused by name (not yet), should be one: %i \n',err(6))
fprintf('Uncertain GE (e): local search is refused by name (not yet), should be one: %i \n',err(7))

output=struct();

end
