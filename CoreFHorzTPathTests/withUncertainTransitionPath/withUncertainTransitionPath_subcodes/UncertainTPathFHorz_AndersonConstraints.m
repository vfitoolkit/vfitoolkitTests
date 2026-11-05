function output=UncertainTPathFHorz_AndersonConstraints(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% UncertainTransitionPath_FHorz with Anderson acceleration (GEnewprice=2) and with price constraints (stage F3).
% Written 2026-11-06, before the first run. The FHorz port of UncertainTPath_AndersonConstraints
% (CoreInfHorzUncertainTPathTests doPart(8)). The model is that of UncertainTPathFHorz_GE (the Cobb-Douglas firm OLG
% model, the regimes differing in TFP: S=1 0.5, S=2 0.52).
%   (a) Anderson on a sure switch (a single history) against TransitionPath_Case1_FHorz with Anderson, from the same guess
%   (b) Anderson on the absorbing tree (switch to S=2 possible in periods 3 to 6): converged, and the general eqm
%       conditions re-evaluated by hand at the returned prices are within tolerance and are the ones reported; for
%       information, its prices against the shooting solve of the same tree
%   (c) constraints that never bind (r and w positive) change nothing, shooting and Anderson: identities
%   (d) a constraint that binds, shooting, one update by hand: the update of maxiter=2 is the fraction-to-the-boundary
%       rule applied to the unconstrained proposal at the guess (maxiter=1 gives the guess and its conditions)
%   (e) refusals of the constraint setup, each by its message

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluate.K=@(aprime,a,z) a;
FnsToEvaluate.N=@(aprime,a,z,kappa_j) kappa_j*z;
GeneralEqmEqns.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
heteroagentoptions=struct();
heteroagentoptions.verbose=0;

%% The initial stationary eqm, and the S=2 one (for the tolerance, as in UncertainTPathFHorz_GE)
[p_eqm,~,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptions,vfoptions);
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
Params2=Params; Params2.firmA=0.52;
[~,~,GEcondns2]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params2,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptions,vfoptions);
GEfloor=max([abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket),abs(GEcondns2.CapitalMarket),abs(GEcondns2.LabourMarket)]);
fprintf('Uncertain FHorz Anderson/constraints: both stationary eqms are solved to 2e-3, so the path tolerance below is meaningful, should be one: %i \n',GEfloor<=2e-3)

howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptionsS=transpathoptionsbaseline;
transpathoptionsS.GEnewprice=3;
transpathoptionsS.GEnewprice3.howtoupdate=howtoupdate;
transpathoptionsS.maxiter=1000;
transpathoptionsS.verbose=0;
transpathoptionsS.toleranceGEprices=10^(-4);
transpathoptionsS.toleranceGEcondns=max(5*GEfloor,10^(-4));
transpathoptionsA=transpathoptionsbaseline;
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=howtoupdate;
transpathoptionsA.maxiter=200;
transpathoptionsA.verbose=0;
transpathoptionsA.toleranceGEprices=transpathoptionsS.toleranceGEprices;
transpathoptionsA.toleranceGEcondns=transpathoptionsS.toleranceGEcondns;
tolC=transpathoptionsS.toleranceGEcondns;
fprintf('Uncertain FHorz Anderson/constraints: toleranceGEcondns %.3e \n',tolC)

PricePathShaper.r=min((1:1:T)/20,1);
PricePathShaper.w=min((1:1:T)/20,1);
ParamPathC.firmA=[0.5*ones(1,T); 0.52*ones(1,T)];

%% (a) Anderson on a sure switch, against TransitionPath_Case1_FHorz with Anderson
pi_SA=[0,1;0,1]; T1=1; T2=1; % the switch lands in period 2
[PricePathTreeA,SummaryA]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SA, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz Anderson (a): converged, should be one: %i \n',SummaryA.GEcondnDist<=tolC)
pf2=SummaryA.finaleqm.S2.p_eqm;
PricePath0D.r=[Params.r, Params.r+(pf2.r-Params.r)*PricePathShaper.r(2:T)];
PricePath0D.w=[Params.w, Params.w+(pf2.w-Params.w)*PricePathShaper.w(2:T)];
ParamPathD.firmA=[0.5, 0.52*ones(1,T-1)];
PricePathD=TransitionPath_Case1_FHorz(PricePath0D, ParamPathD, T, SummaryA.finaleqm.S2.V_final, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions);
diffA=max([max(abs(PricePathTreeA.r-gather(PricePathD.r(:))')),max(abs(PricePathTreeA.w-gather(PricePathD.w(:))'))]);
fprintf('Uncertain FHorz Anderson (a): the tree vs TransitionPath_Case1_FHorz with Anderson, prices within 10*toleranceGEprices, should be one: %i \n',diffA<=10*transpathoptionsA.toleranceGEprices)
fprintf('Uncertain FHorz Anderson (a): for information, that difference is %.3e, after %i Anderson iterations on the tree \n',diffA,SummaryA.iterations)

%% (b) Anderson on the absorbing tree
P=0.75;
pi_SC=[P,1-P;0,1]; T1=2; T2=5;
[PricePathTreeB,SummaryB]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptions, heteroagentoptions);
TreeB=SummaryB.Tree;
fprintf('Uncertain FHorz Anderson (b): %i histories, converged, should be one: %i \n',TreeB.nHist,SummaryB.GEcondnDist<=tolC)
V_finalB=struct('S1',SummaryB.finaleqm.S1.V_final,'S2',SummaryB.finaleqm.S2.V_final);
Policy_finalB=struct('S1',SummaryB.finaleqm.S1.Policy_final,'S2',SummaryB.finaleqm.S2.Policy_final);
[~,PolicyB]=ValueFnOnUncertainTransPath_FHorz(PricePathTreeB, ParamPathC, T, T1, T2, pi_SC, 1, V_finalB, Policy_finalB, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsA, vfoptions);
DistB=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathTreeB, ParamPathC, PolicyB, AgeWeightParamNames, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptionsA, simoptions);
AggB=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluate, DistB, PolicyB, PricePathTreeB, ParamPathC, Params, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptionsA, simoptions);
[~,~,ParamPathNodesB]=UncertainTransPathTree(T,T1,T2,pi_SC,1,[],ParamPathC);
GEhand=zeros(2,TreeB.nNodes);
for kk=1:TreeB.nNodes
    r=PricePathTreeB.r(kk); w=PricePathTreeB.w(kk); K=AggB.K.Mean(kk); N=AggB.N.Mean(kk); A=ParamPathNodesB.firmA(kk);
    GEhand(1,kk)=GeneralEqmEqns.CapitalMarket(r,K,N,Params.firmalpha,Params.firmdelta,A);
    GEhand(2,kk)=GeneralEqmEqns.LabourMarket(w,K,N,Params.firmalpha,A);
end
notlast=(TreeB.nodedate<T)';
fprintf('Uncertain FHorz Anderson (b): the general eqm conditions by hand at the returned prices, within toleranceGEcondns at every node before T, should be one: %i \n',max(sqrt(sum(GEhand(:,notlast).^2,1)))<=tolC)
fprintf('Uncertain FHorz Anderson (b): and they are the ones reported, this should be zero: %.3e \n',max(abs(GEhand-[SummaryB.GEcondnTree.CapitalMarket;SummaryB.GEcondnTree.LabourMarket]),[],'all'))
[PricePathTreeBS,SummaryBS]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsS, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz Anderson (b): for information, Anderson %i iterations against shooting %i; the two solutions differ by at most %.3e in r and %.3e in w \n',SummaryB.iterations,SummaryBS.iterations,max(abs(PricePathTreeB.r-PricePathTreeBS.r)),max(abs(PricePathTreeB.w-PricePathTreeBS.w)))

%% (c) Constraints that never bind change nothing
tpSc=transpathoptionsS; tpSc.constrainpositive={'r','w'};
PricePathTreeSc=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, tpSc, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz constraints (c): shooting, r and w constrained positive (never binds) vs unconstrained, this should be zero: %.3e \n',max([abs(PricePathTreeSc.r-PricePathTreeBS.r),abs(PricePathTreeSc.w-PricePathTreeBS.w)]))
tpAc=transpathoptionsA; tpAc.constrainpositive={'r','w'};
PricePathTreeAc=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, tpAc, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz constraints (c): Anderson, the same, this should be zero: %.3e \n',max([abs(PricePathTreeAc.r-PricePathTreeB.r),abs(PricePathTreeAc.w-PricePathTreeB.w)]))

%% (d) A constraint that binds, shooting, one update by hand
% maxiter=1 returns the guess and the general eqm conditions there; maxiter=2 returns the guess after one update.
% The bounds on r are the range of the guess itself (widened by 1e-7, so the guess is strictly inside), so every node
% whose proposed r leaves that range is cut back; the count is printed and must not be zero, or the check is empty.
tp1=transpathoptionsS; tp1.maxiter=1;
[PricePathTree1,Summary1]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, tp1, simoptions, vfoptions, heteroagentoptions);
lbr=min(PricePathTree1.r)-1e-7; ubr=max(PricePathTree1.r)+1e-7;
tp2=transpathoptionsS; tp2.maxiter=2; tp2.constrainAtoB={'r'}; tp2.constrainAtoBlimits.r=[lbr,ubr]; tp2.constraintstepfraction=0.5;
PricePathTree2=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, tp2, simoptions, vfoptions, heteroagentoptions);
notlast1=(Summary1.Tree.nodedate<T)';
% The shooting proposal at the guess (add=0, factor=0.1: price minus 0.1 times its condition), then the rule
rprop=PricePathTree1.r-0.1*Summary1.GEcondnTree.CapitalMarket;
wprop=PricePathTree1.w-0.1*Summary1.GEcondnTree.LabourMarket;
rprop(~notlast1)=PricePathTree1.r(~notlast1); wprop(~notlast1)=PricePathTree1.w(~notlast1); % period T is never updated
below=(rprop<=lbr); above=(rprop>=ubr);
rhand=rprop;
rhand(below)=PricePathTree1.r(below)+0.5*(lbr-PricePathTree1.r(below));
rhand(above)=PricePathTree1.r(above)+0.5*(ubr-PricePathTree1.r(above));
fprintf('Uncertain FHorz constraints (d): the bound on r cuts back %i of the %i proposed node prices, should NOT be zero: %i \n',sum(below|above),sum(notlast1),sum(below|above))
fprintf('Uncertain FHorz constraints (d): one constrained update vs the rule by hand, r, this should be zero: %.3e \n',max(abs(PricePathTree2.r-rhand)))
fprintf('Uncertain FHorz constraints (d): and w (unconstrained) is the plain proposal, this should be zero: %.3e \n',max(abs(PricePathTree2.w-wprop)))
fprintf('Uncertain FHorz constraints (d): the constrained r stays strictly inside its bounds, should be one: %i \n',all(PricePathTree2.r>lbr & PricePathTree2.r<ubr))

%% (e) Refusals of the constraint setup, each by its message
GEcall=@(tpo) UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, tpo, simoptions, vfoptions, heteroagentoptions);
err=zeros(1,4);
tpE=transpathoptionsS; tpE.constrainpositive={'notaprice'};
try; GEcall(tpE); catch ME; err(1)=contains(ME.message,'which is not a price in PricePathShaper'); end
tpE=transpathoptionsS; tpE.constrainAtoB={'r'};
try; GEcall(tpE); catch ME; err(2)=contains(ME.message,'so you must set transpathoptions.constrainAtoBlimits.r'); end
tpE=transpathoptionsS; tpE.constrainpositive={'r'}; tpE.constrain0to1={'r'};
try; GEcall(tpE); catch ME; err(3)=contains(ME.message,'appears in more than one of'); end
tpE=transpathoptionsS; tpE.constrainpositive={'r'}; tpE.constraintstepfraction=1;
try; GEcall(tpE); catch ME; err(4)=contains(ME.message,'must be strictly between 0 and 1'); end
fprintf('Uncertain FHorz constraints (e): a constrained name that is not a price is refused by name, should be one: %i \n',err(1))
fprintf('Uncertain FHorz constraints (e): constrainAtoB without constrainAtoBlimits is refused by name, should be one: %i \n',err(2))
fprintf('Uncertain FHorz constraints (e): a price in two constraints is refused by name, should be one: %i \n',err(3))
fprintf('Uncertain FHorz constraints (e): a step fraction of one is refused by name, should be one: %i \n',err(4))

output=struct();

end
