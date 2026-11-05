function output=UncertainTPathFHorz_tplus1IntEqns(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% UncertainTransitionPath_FHorz with _tplus1, _tminus1 and intermediate eqns (stage F3). Written 2026-11-06, before the
% first run. The FHorz port of UncertainTPath_tplus1 (CoreInfHorzUncertainTPathTests doPart(9)) and of the intermediate
% eqns half of UncertainTPath_LSIntEqns (doPart(11)). The model is that of UncertainTPathFHorz_GE (the Cobb-Douglas
% firm OLG model, the regimes differing in TFP: S=1 0.5, S=2 0.52) with an asset price q that the households do not see:
%   FnsToEvaluate.D=r_tplus1*aprime, next period's interest on today's saving (a FnsToEvaluate that uses _tplus1)
%   AssetPrice: q-(D+q_tplus1)/(1+r_tplus1) on the path, NONLINEAR in the _tplus1 prices; q-K in the final eqms
% The rule (see UncertainTransitionPath_FHorz): at a branching node a condition that uses _tplus1 holds in
% expectation over the children, sum_c nodecondprob(c)*condition_c, with child c's prices and the aggregates given c
% (AggVarsGivenChildTree); at a node with one child it is the deterministic condition.
%   (a) a SURE switch (one history) against TransitionPath_Case1_FHorz on the same S path
%   (b) identical regimes on a full-support tree: every node in a period has the same prices, and they are
%       TransitionPath_Case1_FHorz's
%   (c) the absorbing tree, shooting: converged; at every node the condition rebuilt by hand from the children is the
%       one reported, and within tolerance before T; at a branching node AggVarsTree.D is the expectation of
%       AggVarsGivenChildTree.D; below a node with one child AggVarsGivenChildTree.D is that node's D; and
%       AggVarsGivenChildTree.D by hand from the node's dist and policy (every age) and the child's r
%   (d) one evaluation (maxiter=1) at a tree whose children at the first branching node are far apart. The
%       DISCRIMINATOR: for a condition linear in the _tplus1 inputs the expectation of the condition and the condition
%       at the expected prices agree; for AssetPrice they do not, and the reported one is the expectation. Also
%       _tminus1 of an aggregate, a price and a parameter in a condition, against the parent node by hand.
%   (e) Anderson on the absorbing tree: converged, and the conditions rebuilt by hand are the ones reported
%   (f) EvalFnOnUncertainTransPath_AggVars_FHorz with _tminus1 of a price, a parameter and an aggregate in
%       FnsToEvaluate: at every node, the parent node's value (transpathoptions.initialvalues at the period 1 node)
%   (g) refusals, each by its message
% INTERMEDIATE EQNS (transpathoptions.intermediateEqns), on the model without q:
%   (h) output Y as an intermediate eqn (CapitalMarket and LabourMarket written in Y) on the absorbing tree: the same
%       prices as without it to 10*toleranceGEprices, and the conditions rebuilt by hand (Y by hand) are the ones reported
%   (i) a SURE switch against TransitionPath_Case1_FHorz with the same intermediate eqns
%   (j) one evaluation (maxiter=1) with a chain of intermediate eqns that uses _tplus1, Rnext=1+r_tplus1 and
%       Dnext=1/Rnext, in a condition ql-Dnext: at a branching node it holds in expectation over the children (by
%       hand), which differs from the condition at the expected r (the DISCRIMINATOR); the chain also tests that an
%       intermediate eqn that uses another one that uses _tplus1 is evaluated per child
% Two shooting solves stop at different iterations, so they are held to 10 times toleranceGEprices, as in
% UncertainTPathFHorz_GE.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
heteroagentoptions=struct();
heteroagentoptions.verbose=0;

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluate0.K=@(aprime,a,z) a;
FnsToEvaluate0.N=@(aprime,a,z,kappa_j) kappa_j*z;
GeneralEqmEqns0.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns0.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

FnsToEvaluate=FnsToEvaluate0;
FnsToEvaluate.D=@(aprime,a,z,r_tplus1) r_tplus1*aprime;
% the final (stationary) eqms, and the path
GeneralEqmEqns=GeneralEqmEqns0;
GeneralEqmEqns.AssetPrice=@(q,K) q-K;
GeneralEqmEqns_Transition=GeneralEqmEqns0;
GeneralEqmEqns_Transition.AssetPrice=@(q,D,q_tplus1,r_tplus1) q-(D+q_tplus1)/(1+r_tplus1);

%% The initial stationary eqm (S=1, firmA=0.5), without q, then q=K there
[p_eqm,~,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate0,GeneralEqmEqns0,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptions,vfoptions);
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
AggVars0=EvalFnOnAgentDist_AggVars_FHorz_Case1(AgentDist_initial, Policy_init, FnsToEvaluate0, Params, [], n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, simoptions);
K0=AggVars0.K.Mean;
Params.q=K0;
% The stationary solves evaluate every FnsToEvaluate, and D needs an r_tplus1; the final eqms do not use D, and on the
% path r_tplus1 is set at every node
Params.r_tplus1=Params.r;
fprintf('Uncertain FHorz _tplus1: initial stationary eqm r=%2.8f, w=%2.8f, q=K=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,Params.q,GEcondns.CapitalMarket,GEcondns.LabourMarket)
% The S=2 final eqm, solved exactly as UncertainTransitionPath_FHorz will solve it, for information
Params2=Params; Params2.firmA=0.52;
[~,~,GEcondns2]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params2,DiscountFactorParamNames,[],[],[],{'r','w','q'},heteroagentoptions,simoptions,vfoptions);
GEfloor=max([abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket),abs(GEcondns2.CapitalMarket),abs(GEcondns2.LabourMarket),abs(GEcondns2.AssetPrice)]);
fprintf('Uncertain FHorz _tplus1: the S=2 stationary eqm has conditions %.3e %.3e %.3e \n',GEcondns2.CapitalMarket,GEcondns2.LabourMarket,GEcondns2.AssetPrice)
% NOT a gate, as in UncertainTPath_tplus1: the period T-1 prices are free and absorb the final eqm residual, so the path
% tolerance is fixed at 1e-3
fprintf('Uncertain FHorz _tplus1: for information, the stationary eqms are solved to %.3e \n',GEfloor)

% q gets factor 1, so its update is q=(D+q_tplus1)/(1+r_tplus1) at the old prices: a backward recursion, done in at
% most T iterations, while r and w take longer
howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1; 'AssetPrice','q',0,1};
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate=howtoupdate;
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=10^(-3);
tolP=transpathoptions.toleranceGEprices;
tolC=transpathoptions.toleranceGEcondns;
fprintf('Uncertain FHorz _tplus1: toleranceGEprices %.3e, toleranceGEcondns %.3e \n',tolP,tolC)

PricePathShaper.r=min((1:1:T)/20,1);
PricePathShaper.w=min((1:1:T)/20,1);
PricePathShaper.q=min((1:1:T)/20,1);

%% (a) A sure switch to S=2 in period 2, against TransitionPath_Case1_FHorz
ParamPathA.firmA=[0.5*ones(1,T); 0.52*ones(1,T)];
pi_SA=[0,1;0,1]; T1=1; T2=1;
[PricePathTreeA,SummaryA]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathA, T, T1, T2, pi_SA, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz _tplus1 (a): a sure switch is a single history, should be one: %i \n',SummaryA.Tree.nHist==1)
fprintf('Uncertain FHorz _tplus1 (a): converged before maxiter, should be one: %i \n',SummaryA.iterations<transpathoptions.maxiter)
pf2=SummaryA.finaleqm.S2.p_eqm;
PricePath0D=struct();
for pname={'r','w','q'}
    PricePath0D.(pname{1})=[Params.(pname{1}), Params.(pname{1})+(pf2.(pname{1})-Params.(pname{1}))*PricePathShaper.(pname{1})(2:T)];
end
ParamPathD.firmA=[0.5, 0.52*ones(1,T-1)];
PricePathD=TransitionPath_Case1_FHorz(PricePath0D, ParamPathD, T, SummaryA.finaleqm.S2.V_final, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);
diffA=max([max(abs(PricePathTreeA.r-gather(PricePathD.r(:))')),max(abs(PricePathTreeA.w-gather(PricePathD.w(:))')),max(abs(PricePathTreeA.q-gather(PricePathD.q(:))'))]);
fprintf('Uncertain FHorz _tplus1 (a): the tree vs TransitionPath_Case1_FHorz, prices (r, w and q) within 10*toleranceGEprices, should be one: %i \n',diffA<=10*tolP)
fprintf('Uncertain FHorz _tplus1 (a): for information, that difference is %.3e, after %i tree iterations; q moves by %.3e over the path \n',diffA,SummaryA.iterations,max(PricePathTreeA.q)-min(PricePathTreeA.q))

%% (b) Identical regimes with a temporary TFP boost, on a full-support tree (dyadic probabilities)
boost=0.5*(1+0.05*max(0,1-(1:1:T)/20));
ParamPathB.firmA=[boost; boost];
pi_SB=[0.5,0.5;0.5,0.5]; T1=3; T2=5;
[PricePathTreeB,SummaryB]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathB, T, T1, T2, pi_SB, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
TreeB=SummaryB.Tree;
spreadB=0;
for tt=1:T
    nodes_t=find(TreeB.nodedate==tt);
    for pname={'r','w','q'}
        spreadB=max(spreadB,max(PricePathTreeB.(pname{1})(nodes_t))-min(PricePathTreeB.(pname{1})(nodes_t)));
    end
end
fprintf('Uncertain FHorz _tplus1 (b): identical regimes, %i histories, every node in a period has the same prices, this should be zero: %.3e \n',TreeB.nHist,spreadB)
ParamPathBD.firmA=boost;
PricePath0BD=struct();
for pname={'r','w','q'}
    PricePath0BD.(pname{1})=Params.(pname{1})+(SummaryB.finaleqm.S1.p_eqm.(pname{1})-Params.(pname{1}))*PricePathShaper.(pname{1});
end
PricePathBD=TransitionPath_Case1_FHorz(PricePath0BD, ParamPathBD, T, SummaryB.finaleqm.S1.V_final, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);
diffB=0;
for pname={'r','w','q'}
    pD=gather(PricePathBD.(pname{1})(:))';
    diffB=max(diffB,max(abs(PricePathTreeB.(pname{1})-pD(TreeB.nodedate'))));
end
fprintf('Uncertain FHorz _tplus1 (b): the tree vs TransitionPath_Case1_FHorz in each period, within 10*toleranceGEprices, should be one: %i \n',diffB<=10*tolP)
fprintf('Uncertain FHorz _tplus1 (b): for information, that difference is %.3e \n',diffB)

%% (c) The absorbing tree, shooting
P=0.75;
pi_SC=[P,1-P;0,1]; T1=2; T2=5;
ParamPathC=ParamPathA;
[PricePathTreeC,SummaryC]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
TreeC=SummaryC.Tree;
nNodesC=TreeC.nNodes;
notlastC=(TreeC.nodedate<T)';
branchingC=find(TreeC.nodenchild>1)';
fprintf('Uncertain FHorz _tplus1 (c): %i histories, %i branching nodes, converged before maxiter, should be one: %i \n',TreeC.nHist,length(branchingC),SummaryC.iterations<transpathoptions.maxiter)
% The conditions rebuilt by hand from the children
AggC=SummaryC.AggVarsTree;
DgivenC=SummaryC.AggVarsGivenChildTree.D;
[~,~,ParamPathNodesC]=UncertainTransPathTree(T,T1,T2,pi_SC,1,[],ParamPathC);
GEhand=zeros(3,nNodesC);
for kk=1:nNodesC
    r=PricePathTreeC.r(kk); w=PricePathTreeC.w(kk); q=PricePathTreeC.q(kk); K=AggC.K.Mean(kk); N=AggC.N.Mean(kk); A=ParamPathNodesC.firmA(kk);
    GEhand(1,kk)=GeneralEqmEqns0.CapitalMarket(r,K,N,Params.firmalpha,Params.firmdelta,A);
    GEhand(2,kk)=GeneralEqmEqns0.LabourMarket(w,K,N,Params.firmalpha,A);
    if TreeC.nodenchild(kk)==0
        GEhand(3,kk)=GeneralEqmEqns_Transition.AssetPrice(q,AggC.D.Mean(kk),q,r);
    else
        children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
        wc=TreeC.nodecondprob(children);
        if TreeC.nodenchild(kk)==1
            wc=1;
        end
        for cc=1:length(children)
            GEhand(3,kk)=GEhand(3,kk)+wc(cc)*GeneralEqmEqns_Transition.AssetPrice(q,DgivenC(children(cc)),PricePathTreeC.q(children(cc)),PricePathTreeC.r(children(cc)));
        end
    end
end
GErep=[SummaryC.GEcondnTree.CapitalMarket;SummaryC.GEcondnTree.LabourMarket;SummaryC.GEcondnTree.AssetPrice];
fprintf('Uncertain FHorz _tplus1 (c): the conditions rebuilt by hand (AssetPrice as the expectation over the children) are the ones reported, this should be zero: %.3e \n',max(abs(GEhand-GErep),[],'all'))
fprintf('Uncertain FHorz _tplus1 (c): and are within toleranceGEcondns at every node before T, should be one: %i \n',max(sqrt(sum(GEhand(:,notlastC).^2,1)))<=tolC)
% AggVarsTree.D at a branching node is the expectation of AggVarsGivenChildTree.D
expD=0;
for kk=branchingC
    children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
    expD=max(expD,abs(AggC.D.Mean(kk)-sum(reshape(TreeC.nodecondprob(children),1,[]).*DgivenC(children))));
end
fprintf('Uncertain FHorz _tplus1 (c): at every branching node AggVarsTree.D is the expectation of AggVarsGivenChildTree.D over the children, this should be zero: %.3e \n',expD)
singlechildC=find(TreeC.nodenchild==1)';
fprintf('Uncertain FHorz _tplus1 (c): below a node with one child, AggVarsGivenChildTree.D is that node''s AggVarsTree.D, this should be zero: %.3e \n',max(abs(DgivenC(TreeC.nodefirstchild(singlechildC))-AggC.D.Mean(singlechildC))))
fprintf('Uncertain FHorz _tplus1 (c): the period 1 node has no parent, so AggVarsGivenChildTree.D is NaN there, should be one: %i \n',isnan(DgivenC(1)))
% AggVarsGivenChildTree.D by hand: the child's r times the node's aggregate saving, over every age (the agent dist
% already carries the age weights)
DistC=reshape(gather(SummaryC.AgentDistPathTree),[N_a*N_z*N_j,nNodesC]);
PolC=reshape(gather(SummaryC.PolicyPathTree),[N_a*N_z*N_j,nNodesC]);
handD=0;
for kk=find(TreeC.nodenchild>0)'
    Aprime_k=sum(DistC(:,kk).*a_grid(PolC(:,kk)));
    children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
    for cc=children
        handD=max(handD,abs(DgivenC(cc)-PricePathTreeC.r(cc)*Aprime_k));
    end
end
fprintf('Uncertain FHorz _tplus1 (c): AggVarsGivenChildTree.D by hand from the node''s dist and policy and the child''s r, this should be zero: %.3e \n',handD)
kb=branchingC(1);
childrenkb=TreeC.nodefirstchild(kb)+(0:TreeC.nodenchild(kb)-1);
fprintf('Uncertain FHorz _tplus1 (c): for scale, at the first branching node the two children''s r differ by, should NOT be zero: %.3e \n',abs(diff(PricePathTreeC.r(childrenkb))))
fprintf('Uncertain FHorz _tplus1 (c): for information, %i iterations, price distance %.3e, general eqm distance %.3e \n',SummaryC.iterations,SummaryC.PricePathDist,SummaryC.GEcondnDist)

%% (d) One evaluation at a tree whose children at the first branching node are far apart
% AssetPrice (nonlinear in the _tplus1 inputs), AssetLinear (linear in them, price ql), and Lag (_tminus1 of an
% aggregate, a price and a parameter, price Kl). maxiter=1 evaluates the tree at PricePathTree0 and stops before any
% update. Only the evaluation is checked here, so the accuracy of the final eqms does not matter.
GeneralEqmEqnsD=GeneralEqmEqns;
GeneralEqmEqnsD.AssetLinear=@(ql,q,r,K) ql-0.96*(r*K+q);
GeneralEqmEqnsD.Lag=@(Kl,K,r,firmA) Kl-(K+r*firmA);
GeneralEqmEqns_TransitionD=GeneralEqmEqns_Transition;
GeneralEqmEqns_TransitionD.AssetLinear=@(ql,D,q_tplus1) ql-0.96*(D+q_tplus1);
GeneralEqmEqns_TransitionD.Lag=@(Kl,K_tminus1,r_tminus1,firmA_tminus1) Kl-(K_tminus1+r_tminus1*firmA_tminus1);
ParamsD=Params;
ParamsD.ql=0.96*(Params.r*K0+K0);
ParamsD.Kl=K0+Params.r*Params.firmA;
PricePathShaperD=PricePathShaper;
PricePathShaperD.ql=PricePathShaper.q;
PricePathShaperD.Kl=PricePathShaper.q;
transpathoptionsD=transpathoptions;
transpathoptionsD.GEnewprice3.howtoupdate=[howtoupdate; {'AssetLinear','ql',0,1; 'Lag','Kl',0,1}];
transpathoptionsD.maxiter=1;
transpathoptionsD.initialvalues.K=K0;
transpathoptionsD.initialvalues.r=Params.r;
transpathoptionsD.initialvalues.firmA=Params.firmA;
PricePathTree0=PricePathTreeC;
PricePathTree0.ql=0.96*(PricePathTreeC.r.*AggC.K.Mean+PricePathTreeC.q);
PricePathTree0.Kl=AggC.K.Mean+0.01;
kS2=childrenkb(TreeC.nodeS(childrenkb)==2);
PricePathTree0.r(kS2)=PricePathTree0.r(kS2)+0.05;
PricePathTree0.q(kS2)=1.5*PricePathTree0.q(kS2);
transpathoptionsD.PricePathTree0=PricePathTree0;
[PricePathTreeD,SummaryD]=UncertainTransitionPath_FHorz(PricePathShaperD, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsD, GeneralEqmEqns_TransitionD, ParamsD, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsD, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz _tplus1 (d): maxiter=1 returns the prices it evaluated, this should be zero: %.3e \n',max([max(abs(PricePathTreeD.r(notlastC)-PricePathTree0.r(notlastC))),max(abs(PricePathTreeD.q(notlastC)-PricePathTree0.q(notlastC)))]))
AggD=SummaryD.AggVarsTree;
DgivenD=SummaryD.AggVarsGivenChildTree.D;
handAsset=zeros(1,nNodesC); handLinear=zeros(1,nNodesC); handLag=zeros(1,nNodesC);
for kk=1:nNodesC
    q=PricePathTreeD.q(kk); ql=PricePathTreeD.ql(kk);
    if TreeC.nodenchild(kk)==0
        children=kk; wc=1; Dc=AggD.D.Mean(kk);
    else
        children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
        if TreeC.nodenchild(kk)==1
            wc=1;
        else
            wc=TreeC.nodecondprob(children);
        end
        Dc=DgivenD(children);
    end
    for cc=1:length(children)
        handAsset(kk)=handAsset(kk)+wc(cc)*GeneralEqmEqns_TransitionD.AssetPrice(q,Dc(cc),PricePathTreeD.q(children(cc)),PricePathTreeD.r(children(cc)));
        handLinear(kk)=handLinear(kk)+wc(cc)*GeneralEqmEqns_TransitionD.AssetLinear(ql,Dc(cc),PricePathTreeD.q(children(cc)));
    end
    parentnode=TreeC.nodeparent(kk);
    if parentnode>0
        handLag(kk)=GeneralEqmEqns_TransitionD.Lag(PricePathTreeD.Kl(kk),AggD.K.Mean(parentnode),PricePathTreeD.r(parentnode),ParamPathNodesC.firmA(parentnode));
    else
        handLag(kk)=GeneralEqmEqns_TransitionD.Lag(PricePathTreeD.Kl(kk),K0,Params.r,Params.firmA);
    end
end
fprintf('Uncertain FHorz _tplus1 (d): AssetPrice, the expectation over the children by hand is the one reported, this should be zero: %.3e \n',max(abs(handAsset-SummaryD.GEcondnTree.AssetPrice)))
fprintf('Uncertain FHorz _tplus1 (d): AssetLinear, the same, this should be zero: %.3e \n',max(abs(handLinear-SummaryD.GEcondnTree.AssetLinear)))
fprintf('Uncertain FHorz _tplus1 (d): Lag, _tminus1 of K, r and firmA are the parent node''s (initialvalues at the period 1 node), this should be zero: %.3e \n',max(abs(handLag-SummaryD.GEcondnTree.Lag)))
% The discriminator, at the first branching node: the condition at the expected prices (and expected D)
wkb=reshape(TreeC.nodecondprob(childrenkb),1,[]);
EqD=sum(wkb.*PricePathTreeD.q(childrenkb));
ErD=sum(wkb.*PricePathTreeD.r(childrenkb));
EDD=sum(wkb.*DgivenD(childrenkb));
plugAsset=GeneralEqmEqns_TransitionD.AssetPrice(PricePathTreeD.q(kb),EDD,EqD,ErD);
plugLinear=GeneralEqmEqns_TransitionD.AssetLinear(PricePathTreeD.ql(kb),EDD,EqD);
fprintf('Uncertain FHorz _tplus1 (d): AssetLinear is linear in the _tplus1 inputs, so at the branching node the condition at the expected prices is the expected condition, this should be zero: %.3e \n',abs(plugLinear-SummaryD.GEcondnTree.AssetLinear(kb)))
fprintf('Uncertain FHorz _tplus1 (d): AssetPrice is not, so the condition at the expected prices differs from the one reported by more than 1e-4, should be one: %i \n',abs(plugAsset-SummaryD.GEcondnTree.AssetPrice(kb))>1e-4)
fprintf('Uncertain FHorz _tplus1 (d): for information, that difference is %.3e (the reported condition is %.3e) \n',abs(plugAsset-SummaryD.GEcondnTree.AssetPrice(kb)),SummaryD.GEcondnTree.AssetPrice(kb))

%% (e) Anderson on the absorbing tree
transpathoptionsA=transpathoptions;
transpathoptionsA=rmfield(transpathoptionsA,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=howtoupdate;
transpathoptionsA.maxiter=200;
[PricePathTreeE,SummaryE]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz _tplus1 (e): Anderson converged, should be one: %i \n',SummaryE.GEcondnDist<=tolC)
fprintf('Uncertain FHorz _tplus1 (e): and took more than one iteration (the initial guess is not already within tolerance), should be one: %i \n',SummaryE.iterations>1)
AggE=SummaryE.AggVarsTree;
DgivenE=SummaryE.AggVarsGivenChildTree.D;
handE=zeros(1,nNodesC);
for kk=1:nNodesC
    q=PricePathTreeE.q(kk);
    if TreeC.nodenchild(kk)==0
        handE(kk)=GeneralEqmEqns_Transition.AssetPrice(q,AggE.D.Mean(kk),q,PricePathTreeE.r(kk));
    else
        children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
        wc=TreeC.nodecondprob(children);
        if TreeC.nodenchild(kk)==1
            wc=1;
        end
        for cc=1:length(children)
            handE(kk)=handE(kk)+wc(cc)*GeneralEqmEqns_Transition.AssetPrice(q,DgivenE(children(cc)),PricePathTreeE.q(children(cc)),PricePathTreeE.r(children(cc)));
        end
    end
end
fprintf('Uncertain FHorz _tplus1 (e): AssetPrice rebuilt by hand from the children is the one reported, this should be zero: %.3e \n',max(abs(handE-SummaryE.GEcondnTree.AssetPrice)))
diffE=max([max(abs(PricePathTreeE.r-PricePathTreeC.r)),max(abs(PricePathTreeE.w-PricePathTreeC.w)),max(abs(PricePathTreeE.q-PricePathTreeC.q))]);
fprintf('Uncertain FHorz _tplus1 (e): for information, Anderson vs shooting prices differ by %.3e, after %i Anderson iterations \n',diffE,SummaryE.iterations)

%% (f) _tminus1 in FnsToEvaluate, EvalFnOnUncertainTransPath_AggVars_FHorz on the tree of (c)
FnsToEvaluateF.K=@(aprime,a,z) a;
FnsToEvaluateF.Dlag=@(aprime,a,z,r_tminus1) r_tminus1*a;
FnsToEvaluateF.Alag=@(aprime,a,z,firmA_tminus1) firmA_tminus1*a;
FnsToEvaluateF.Klag=@(aprime,a,z,K_tminus1) K_tminus1;
transpathoptionsF=transpathoptionsbaseline;
transpathoptionsF.initialvalues.r=Params.r;
transpathoptionsF.initialvalues.firmA=Params.firmA;
transpathoptionsF.initialvalues.K=K0;
AggF=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluateF, SummaryC.AgentDistPathTree, SummaryC.PolicyPathTree, PricePathTreeC, ParamPathC, Params, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptionsF, simoptions);
rpar=zeros(1,nNodesC); Apar=zeros(1,nNodesC); Kpar=zeros(1,nNodesC);
for kk=1:nNodesC
    parentnode=TreeC.nodeparent(kk);
    if parentnode>0
        rpar(kk)=PricePathTreeC.r(parentnode); Apar(kk)=ParamPathNodesC.firmA(parentnode); Kpar(kk)=AggF.K.Mean(parentnode);
    else
        rpar(kk)=Params.r; Apar(kk)=Params.firmA; Kpar(kk)=K0;
    end
end
fprintf('Uncertain FHorz _tplus1 (f): r_tminus1*a aggregates to the parent node''s r times K, this should be zero: %.3e \n',max(abs(AggF.Dlag.Mean-rpar.*AggF.K.Mean)))
fprintf('Uncertain FHorz _tplus1 (f): firmA_tminus1*a aggregates to the parent node''s firmA times K, this should be zero: %.3e \n',max(abs(AggF.Alag.Mean-Apar.*AggF.K.Mean)))
fprintf('Uncertain FHorz _tplus1 (f): K_tminus1 aggregates to the parent node''s K, this should be zero: %.3e \n',max(abs(AggF.Klag.Mean-Kpar)))
fprintf('Uncertain FHorz _tplus1 (f): and the K of (f) is the K of (c), this should be zero: %.3e \n',max(abs(AggF.K.Mean-AggC.K.Mean)))
% the oracle has teeth: S=2 raises firmA at the switch, so a child's firmA_tminus1 differs from its own firmA
switchnodes=find(TreeC.nodeS'==2 & [0,TreeC.nodeS(TreeC.nodeparent(2:end))']==1);
fprintf('Uncertain FHorz _tplus1 (f): for scale, at the switch nodes firmA_tminus1 differs from firmA by, should NOT be zero: %.3e \n',max(abs(Apar(switchnodes)-ParamPathNodesC.firmA(switchnodes))))

%% (g) Refusals, each by its message
GEcall=@(tpo,GEtrans) UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GEtrans, Params, DiscountFactorParamNames, AgeWeightParamNames, tpo, simoptions, vfoptions, heteroagentoptions);
err=zeros(1,5);
tpE=transpathoptions; tpE.stockvars={'q'};
try; GEcall(tpE,GeneralEqmEqns_Transition); catch ME; err(1)=contains(ME.message,'transpathoptions.stockvars is not yet implemented'); end
GeneralEqmEqns_TransitionE=GeneralEqmEqns_Transition;
GeneralEqmEqns_TransitionE.CapitalMarket=@(r,K_tminus1,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K_tminus1^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
try; GEcall(transpathoptions,GeneralEqmEqns_TransitionE); catch ME; err(2)=contains(ME.message,'a _tminus1 input is used, so transpathoptions.initialvalues must be given'); end
GeneralEqmEqns_TransitionE=GeneralEqmEqns_Transition;
GeneralEqmEqns_TransitionE.AssetPrice=@(q,D_tminus1,q_tplus1,r_tplus1) q-(D_tminus1+q_tplus1)/(1+r_tplus1);
tpE=transpathoptions; tpE.initialvalues.D=Params.r*K0;
try; GEcall(tpE,GeneralEqmEqns_TransitionE); catch ME; err(3)=contains(ME.message,'D_tminus1 is used, but FnsToEvaluate.D uses _tplus1'); end
FnsToEvaluateE=FnsToEvaluate;
FnsToEvaluateE.Dlag=@(aprime,a,z,D_tminus1) D_tminus1;
tpE=transpathoptionsbaseline; tpE.initialvalues.D=Params.r*K0;
try; EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluateE, SummaryC.AgentDistPathTree, SummaryC.PolicyPathTree, PricePathTreeC, ParamPathC, Params, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, tpE, simoptions); catch ME; err(4)=contains(ME.message,'D_tminus1 is used, but FnsToEvaluate.D uses _tplus1'); end
FnsToEvaluateE=FnsToEvaluate0;
FnsToEvaluateE.Anext=@(aprime,a,z,firmA_tplus1) firmA_tplus1*a;
try; EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluateE, SummaryC.AgentDistPathTree, SummaryC.PolicyPathTree, PricePathTreeC, ParamPathC, Params, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptionsbaseline, simoptions); catch ME; err(5)=contains(ME.message,'_tplus1 input that is NOT a price'); end
fprintf('Uncertain FHorz _tplus1 (g): transpathoptions.stockvars is refused by name (not yet), should be one: %i \n',err(1))
fprintf('Uncertain FHorz _tplus1 (g): a _tminus1 in the conditions without transpathoptions.initialvalues is refused by name, should be one: %i \n',err(2))
fprintf('Uncertain FHorz _tplus1 (g): D_tminus1 in the conditions, with D a FnsToEvaluate that uses _tplus1, is refused by name, should be one: %i \n',err(3))
fprintf('Uncertain FHorz _tplus1 (g): the same in FnsToEvaluate, by EvalFnOnUncertainTransPath_AggVars_FHorz, should be one: %i \n',err(4))
fprintf('Uncertain FHorz _tplus1 (g): _tplus1 of a parameter is refused by name, should be one: %i \n',err(5))

%% (h) Output Y as an intermediate eqn, on the model without q
howtoupdate0={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptions0=transpathoptions;
transpathoptions0.GEnewprice3.howtoupdate=howtoupdate0;
transpathoptions0.toleranceGEcondns=2*10^(-3); % the two-price stationary solves stop well inside this (UncertainTPathFHorz_GE)
PricePathShaper0.r=PricePathShaper.r;
PricePathShaper0.w=PricePathShaper.w;
Params0=rmfield(Params,{'q','r_tplus1'});
PricePathTree0H=UncertainTransitionPath_FHorz(PricePathShaper0, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate0, GeneralEqmEqns0, GeneralEqmEqns0, Params0, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions0, simoptions, vfoptions, heteroagentoptions);
GeneralEqmEqnsI.CapitalMarket=@(r,Y,K,firmalpha,firmdelta) r-(firmalpha*Y/K-firmdelta);
GeneralEqmEqnsI.LabourMarket=@(w,Y,N,firmalpha) w-((1-firmalpha)*Y/N);
transpathoptionsI=transpathoptions0;
transpathoptionsI.intermediateEqns.Y=@(K,N,firmA,firmalpha) firmA*(K^firmalpha)*(N^(1-firmalpha));
heteroagentoptionsI=heteroagentoptions;
heteroagentoptionsI.intermediateEqns.Y=transpathoptionsI.intermediateEqns.Y;
[PricePathTreeH,SummaryH]=UncertainTransitionPath_FHorz(PricePathShaper0, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate0, GeneralEqmEqnsI, GeneralEqmEqnsI, Params0, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions, heteroagentoptionsI);
fprintf('Uncertain FHorz intermediate (h): converged before maxiter, should be one: %i \n',SummaryH.iterations<transpathoptionsI.maxiter)
diffH=max([max(abs(PricePathTreeH.r-PricePathTree0H.r)),max(abs(PricePathTreeH.w-PricePathTree0H.w))]);
fprintf('Uncertain FHorz intermediate (h): the same model written with Y, the same prices within 10*toleranceGEprices, should be one: %i \n',diffH<=10*tolP)
fprintf('Uncertain FHorz intermediate (h): for information, that difference is %.3e, after %i iterations \n',diffH,SummaryH.iterations)
AggH=SummaryH.AggVarsTree;
GEhandH=zeros(2,nNodesC);
for kk=1:nNodesC
    K=AggH.K.Mean(kk); N=AggH.N.Mean(kk);
    Y=transpathoptionsI.intermediateEqns.Y(K,N,ParamPathNodesC.firmA(kk),Params.firmalpha);
    GEhandH(1,kk)=GeneralEqmEqnsI.CapitalMarket(PricePathTreeH.r(kk),Y,K,Params.firmalpha,Params.firmdelta);
    GEhandH(2,kk)=GeneralEqmEqnsI.LabourMarket(PricePathTreeH.w(kk),Y,N,Params.firmalpha);
end
GErepH=[SummaryH.GEcondnTree.CapitalMarket;SummaryH.GEcondnTree.LabourMarket];
fprintf('Uncertain FHorz intermediate (h): the conditions rebuilt by hand (Y by hand) are the ones reported, this should be zero: %.3e \n',max(abs(GEhandH-GErepH),[],'all'))

%% (i) A sure switch to S=2 in period 2, against TransitionPath_Case1_FHorz with the same intermediate eqns
[PricePathTreeI,SummaryI]=UncertainTransitionPath_FHorz(PricePathShaper0, ParamPathC, T, 1, 1, pi_SA, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate0, GeneralEqmEqnsI, GeneralEqmEqnsI, Params0, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions, heteroagentoptionsI);
fprintf('Uncertain FHorz intermediate (i): a sure switch is a single history, converged before maxiter, should be one: %i \n',SummaryI.Tree.nHist==1 && SummaryI.iterations<transpathoptionsI.maxiter)
pf2=SummaryI.finaleqm.S2.p_eqm;
clear PricePath0D
PricePath0D.r=[Params.r, Params.r+(pf2.r-Params.r)*PricePathShaper.r(2:T)];
PricePath0D.w=[Params.w, Params.w+(pf2.w-Params.w)*PricePathShaper.w(2:T)];
PricePathDI=TransitionPath_Case1_FHorz(PricePath0D, ParamPathD, T, SummaryI.finaleqm.S2.V_final, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate0, GeneralEqmEqnsI, Params0, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsI, simoptions, vfoptions);
diffI=max([max(abs(PricePathTreeI.r-gather(PricePathDI.r(:))')),max(abs(PricePathTreeI.w-gather(PricePathDI.w(:))'))]);
fprintf('Uncertain FHorz intermediate (i): the tree vs TransitionPath_Case1_FHorz with intermediateEqns, prices within 10*toleranceGEprices, should be one: %i \n',diffI<=10*tolP)
fprintf('Uncertain FHorz intermediate (i): for information, that difference is %.3e, after %i tree iterations \n',diffI,SummaryI.iterations)

%% (j) One evaluation with a chain of intermediate eqns that uses _tplus1
GeneralEqmEqnsJ=GeneralEqmEqnsI;
GeneralEqmEqnsJ.Discount=@(ql,r) ql-1/(1+r); % the final eqms
GeneralEqmEqns_TransitionJ=GeneralEqmEqnsI;
GeneralEqmEqns_TransitionJ.Discount=@(ql,Dnext) ql-Dnext;
transpathoptionsJ=transpathoptionsI;
transpathoptionsJ.intermediateEqns.Rnext=@(r_tplus1) 1+r_tplus1;
transpathoptionsJ.intermediateEqns.Dnext=@(Rnext) 1/Rnext;
transpathoptionsJ.GEnewprice3.howtoupdate=[howtoupdate0; {'Discount','ql',0,1}];
transpathoptionsJ.maxiter=1;
ParamsJ=Params0;
ParamsJ.ql=1/(1+Params.r);
PricePathShaperJ=PricePathShaper0;
PricePathShaperJ.ql=PricePathShaper0.r;
PricePathTree0J=PricePathTreeH;
PricePathTree0J.ql=1./(1+PricePathTreeH.r);
PricePathTree0J.r(kS2)=PricePathTree0J.r(kS2)+0.05; % the children of the first branching node far apart
transpathoptionsJ.PricePathTree0=PricePathTree0J;
[PricePathTreeJ,SummaryJ]=UncertainTransitionPath_FHorz(PricePathShaperJ, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate0, GeneralEqmEqnsJ, GeneralEqmEqns_TransitionJ, ParamsJ, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsJ, simoptions, vfoptions, heteroagentoptionsI);
fprintf('Uncertain FHorz intermediate (j): maxiter=1 returns the prices it evaluated, this should be zero: %.3e \n',max(abs(PricePathTreeJ.r(notlastC)-PricePathTree0J.r(notlastC))))
handJ=zeros(1,nNodesC);
for kk=1:nNodesC
    ql=PricePathTreeJ.ql(kk);
    if TreeC.nodenchild(kk)==0
        children=kk; wc=1;
    else
        children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
        if TreeC.nodenchild(kk)==1
            wc=1;
        else
            wc=TreeC.nodecondprob(children);
        end
    end
    for cc=1:length(children)
        Rnext=transpathoptionsJ.intermediateEqns.Rnext(PricePathTreeJ.r(children(cc)));
        Dnext=transpathoptionsJ.intermediateEqns.Dnext(Rnext);
        handJ(kk)=handJ(kk)+wc(cc)*GeneralEqmEqns_TransitionJ.Discount(ql,Dnext);
    end
end
fprintf('Uncertain FHorz intermediate (j): Discount, the expectation over the children of the chain Rnext, Dnext by hand is the one reported, this should be zero: %.3e \n',max(abs(handJ-SummaryJ.GEcondnTree.Discount)))
plugJ=GeneralEqmEqns_TransitionJ.Discount(PricePathTreeJ.ql(kb),1/(1+sum(wkb.*PricePathTreeJ.r(childrenkb))));
fprintf('Uncertain FHorz intermediate (j): Dnext is nonlinear in r_tplus1, so at the branching node the condition at the expected r differs from the one reported by more than 1e-6, should be one: %i \n',abs(plugJ-SummaryJ.GEcondnTree.Discount(kb))>1e-6)
fprintf('Uncertain FHorz intermediate (j): for information, that difference is %.3e \n',abs(plugJ-SummaryJ.GEcondnTree.Discount(kb)))
AggJ=SummaryJ.AggVarsTree;
GEhandJ=zeros(1,nNodesC);
for kk=1:nNodesC
    K=AggJ.K.Mean(kk); N=AggJ.N.Mean(kk);
    Y=transpathoptionsJ.intermediateEqns.Y(K,N,ParamPathNodesC.firmA(kk),Params.firmalpha);
    GEhandJ(kk)=GeneralEqmEqnsI.CapitalMarket(PricePathTreeJ.r(kk),Y,K,Params.firmalpha,Params.firmdelta);
end
fprintf('Uncertain FHorz intermediate (j): and CapitalMarket, with Y at the node, by hand is the one reported, this should be zero: %.3e \n',max(abs(GEhandJ-SummaryJ.GEcondnTree.CapitalMarket)))

output=struct();

end
