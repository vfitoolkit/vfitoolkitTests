function output=UncertainTPathFHorz_GE(T,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% UncertainTransitionPath_FHorz, the general eqm on the tree, by shooting (stage F2). Written 2026-11-06, before the
% first run. The FHorz port of UncertainTPath_GE (CoreInfHorzUncertainTPathTests doPart(7)). The model is the
% Cobb-Douglas firm OLG model of MRFHorzTPath_byhand (nod z): K is assets and N is efficiency labour kappa_j*z, r and w
% clear the capital and labour markets. The regimes differ in TFP, firmA, which enters only the general eqm conditions
% (S=1: 0.5, S=2: 0.52, the long-run level the multiple reveal bank moves to).
%   (a) a SURE switch to S=2 in period 2 is a deterministic transition: the tree (one history) against
%       TransitionPath_Case1_FHorz from the same initial guess, with the same final eqm
%   (b) identical regimes with a temporary TFP boost, on a full-support tree: every history is the deterministic
%       path, so all nodes in a period have the same prices, and they are TransitionPath_Case1_FHorz's
%   (c) an absorbing switch (S=2 for good, from any period 3 to 6, with probability 1/4 each period): converged; the
%       general eqm conditions re-evaluated by hand at the returned prices are within tolerance and are the ones
%       reported; period T is the final eqm of each S; the possibility of S=2 moves the prices before any switch
%       (should NOT be zero)
%   (d) restarting from the solution (transpathoptions.PricePathTree0) converges at once to the same prices
%   (e) refusals, each by its message: the options, what remains refused of the stage F3 features, and the features
%       of stage F4 not yet built
% Shooting is a fixed point, so (a) and (b) are equal only up to the convergence tolerance (the two solvers stop at
% different iterations); they are checked against 10 times toleranceGEprices.

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

%% The initial stationary eqm (S=1, firmA=0.5)
[p_eqm,~,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptions,vfoptions);
fprintf('Uncertain FHorz GE: initial stationary eqm r=%2.8f, w=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);
% The S=2 final eqm, solved here exactly as UncertainTransitionPath_FHorz will solve it (same start, same options), so
% that the path tolerance can be set above the worse of the two stationary residuals
Params2=Params; Params2.firmA=0.52;
[~,~,GEcondns2]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params2,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptions,vfoptions);
GEfloor=max([abs(GEcondns.CapitalMarket),abs(GEcondns.LabourMarket),abs(GEcondns2.CapitalMarket),abs(GEcondns2.LabourMarket)]);
fprintf('Uncertain FHorz GE: the S=2 stationary eqm has conditions %.3e %.3e \n',GEcondns2.CapitalMarket,GEcondns2.LabourMarket)
fprintf('Uncertain FHorz GE: both stationary eqms are solved to 2e-3, so the path tolerance below is meaningful, should be one: %i \n',GEfloor<=2e-3)

transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=max(5*GEfloor,10^(-4));
tolP=transpathoptions.toleranceGEprices;
fprintf('Uncertain FHorz GE: toleranceGEprices %.3e, toleranceGEcondns %.3e \n',transpathoptions.toleranceGEprices,transpathoptions.toleranceGEcondns)

PricePathShaper.r=min((1:1:T)/20,1);
PricePathShaper.w=min((1:1:T)/20,1);

%% (a) A sure switch to S=2 in period 2, against TransitionPath_Case1_FHorz
ParamPathA.firmA=[0.5*ones(1,T); 0.52*ones(1,T)];
pi_SA=[0,1;0,1]; T1=1; T2=1; % uncertainty faced in period 1, so the switch lands in period 2
[PricePathTreeA,SummaryA]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathA, T, T1, T2, pi_SA, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz GE (a): a sure switch is a single history, should be one: %i \n',SummaryA.Tree.nHist==1)
fprintf('Uncertain FHorz GE (a): converged before maxiter, should be one: %i \n',SummaryA.iterations<transpathoptions.maxiter)
% The deterministic path, from the same guess (period 1 is S=1, which has no final eqm here, so the guess is the
% initial price there) and with the same final eqm
pf2=SummaryA.finaleqm.S2.p_eqm;
PricePath0D.r=[Params.r, Params.r+(pf2.r-Params.r)*PricePathShaper.r(2:T)];
PricePath0D.w=[Params.w, Params.w+(pf2.w-Params.w)*PricePathShaper.w(2:T)];
ParamPathD.firmA=[0.5, 0.52*ones(1,T-1)];
PricePathD=TransitionPath_Case1_FHorz(PricePath0D, ParamPathD, T, SummaryA.finaleqm.S2.V_final, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);
diffA=max([max(abs(PricePathTreeA.r-gather(PricePathD.r(:))'),[],'all'),max(abs(PricePathTreeA.w-gather(PricePathD.w(:))'),[],'all')]);
fprintf('Uncertain FHorz GE (a): the tree vs TransitionPath_Case1_FHorz, prices within 10*toleranceGEprices, should be one: %i \n',diffA<=10*tolP)
fprintf('Uncertain FHorz GE (a): for information, that difference is %.3e, after %i tree iterations \n',diffA,SummaryA.iterations)
fprintf('Uncertain FHorz GE (a): for scale, the TFP rise moves r on the path by, should NOT be zero: %.3e \n',max(abs(PricePathTreeA.r-Params.r)))

%% (b) Identical regimes with a temporary TFP boost, on a full-support tree
boost=0.5*(1+0.05*max(0,1-(1:1:T)/20));
ParamPathB.firmA=[boost; boost];
pi_SB=[0.5,0.5;0.5,0.5]; T1=3; T2=5;
[PricePathTreeB,SummaryB]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathB, T, T1, T2, pi_SB, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
TreeB=SummaryB.Tree;
spreadB=0;
for tt=1:T
    nodes_t=find(TreeB.nodedate==tt);
    spreadB=max([spreadB,max(PricePathTreeB.r(nodes_t))-min(PricePathTreeB.r(nodes_t)),max(PricePathTreeB.w(nodes_t))-min(PricePathTreeB.w(nodes_t))]);
end
fprintf('Uncertain FHorz GE (b): identical regimes, %i histories, every node in a period has the same prices, this should be zero: %.3e \n',TreeB.nHist,spreadB)
ParamPathBD.firmA=boost;
% the guess the tree used (both S have the same final eqm, which is the initial one up to the stationary solver)
PricePath0BD.r=Params.r+(SummaryB.finaleqm.S1.p_eqm.r-Params.r)*PricePathShaper.r;
PricePath0BD.w=Params.w+(SummaryB.finaleqm.S1.p_eqm.w-Params.w)*PricePathShaper.w;
PricePathBD=TransitionPath_Case1_FHorz(PricePath0BD, ParamPathBD, T, SummaryB.finaleqm.S1.V_final, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions);
rD=gather(PricePathBD.r(:))'; wD=gather(PricePathBD.w(:))';
diffB=max([max(abs(PricePathTreeB.r-rD(TreeB.nodedate')),[],'all'),max(abs(PricePathTreeB.w-wD(TreeB.nodedate')),[],'all')]);
fprintf('Uncertain FHorz GE (b): the tree vs TransitionPath_Case1_FHorz in each period, within 10*toleranceGEprices, should be one: %i \n',diffB<=10*tolP)
fprintf('Uncertain FHorz GE (b): for information, that difference is %.3e; the boost moves r by %.3e from the initial eqm \n',diffB,max(abs(rD-Params.r)))

%% (c) An absorbing switch
P=0.75;
pi_SC=[P,1-P;0,1]; T1=2; T2=5; % agents face it in periods 2 to 5, so the switch can land in any of 3 to 6
ParamPathC=ParamPathA;
[PricePathTreeC,SummaryC]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
TreeC=SummaryC.Tree;
fprintf('Uncertain FHorz GE (c): %i histories, converged before maxiter, should be one: %i \n',TreeC.nHist,SummaryC.iterations<transpathoptions.maxiter)
% The general eqm conditions re-evaluated by hand at the returned prices
V_finalC=struct('S1',SummaryC.finaleqm.S1.V_final,'S2',SummaryC.finaleqm.S2.V_final);
Policy_finalC=struct('S1',SummaryC.finaleqm.S1.Policy_final,'S2',SummaryC.finaleqm.S2.Policy_final);
[~,PolicyC]=ValueFnOnUncertainTransPath_FHorz(PricePathTreeC, ParamPathC, T, T1, T2, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
DistC=AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathTreeC, ParamPathC, PolicyC, AgeWeightParamNames, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions);
AggC=EvalFnOnUncertainTransPath_AggVars_FHorz(FnsToEvaluate, DistC, PolicyC, PricePathTreeC, ParamPathC, Params, T, T1, T2, pi_SC, 1, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, transpathoptions, simoptions);
[~,~,ParamPathNodesC]=UncertainTransPathTree(T,T1,T2,pi_SC,1,[],ParamPathC);
GEhand=zeros(2,TreeC.nNodes);
for kk=1:TreeC.nNodes
    r=PricePathTreeC.r(kk); w=PricePathTreeC.w(kk); K=AggC.K.Mean(kk); N=AggC.N.Mean(kk); A=ParamPathNodesC.firmA(kk);
    GEhand(1,kk)=GeneralEqmEqns.CapitalMarket(r,K,N,Params.firmalpha,Params.firmdelta,A);
    GEhand(2,kk)=GeneralEqmEqns.LabourMarket(w,K,N,Params.firmalpha,A);
end
notlast=(TreeC.nodedate<T)';
fprintf('Uncertain FHorz GE (c): the general eqm conditions by hand at the returned prices, within toleranceGEcondns at every node before T, should be one: %i \n',max(sqrt(sum(GEhand(:,notlast).^2,1)))<=transpathoptions.toleranceGEcondns)
fprintf('Uncertain FHorz GE (c): and they are the ones reported, this should be zero: %.3e \n',max(abs(GEhand-[SummaryC.GEcondnTree.CapitalMarket;SummaryC.GEcondnTree.LabourMarket]),[],'all'))
lastnodes=find(TreeC.nodedate==T)';
heldT=0;
for kk=lastnodes
    pfk=SummaryC.finaleqm.(['S',num2str(TreeC.nodeS(kk))]).p_eqm;
    heldT=max([heldT,abs(PricePathTreeC.r(kk)-pfk.r),abs(PricePathTreeC.w(kk)-pfk.w)]);
end
fprintf('Uncertain FHorz GE (c): period T is the final eqm of each S, this should be zero: %.3e \n',heldT)
% The deterministic S=1 path would be the initial eqm (a null reform); the possibility of the TFP rise is what moves
% path 0. Not in period 1: K there is the initial dist's and N is exogenous, so r and w in period 1 are the initial
% eqm's whatever is expected (the first run, 2026-11-06, measured 1.4e-5 there, below toleranceGEprices). Over the
% whole never-switch history savings respond, and the move must beat the solver's own noise. Each condition is a price
% minus a marginal product, so a converged price is within about toleranceGEcondns of clearing its market: that is the
% noise (toleranceGEprices bounds only the last step; the second run measured 8.0e-4, 3.3 times toleranceGEcondns)
hh0=find(all(TreeC.historyS==1,2));
nodes0=TreeC.histories(hh0,1:T-1);
move0=max(max(abs(PricePathTreeC.r(nodes0)-Params.r)),max(abs(PricePathTreeC.w(nodes0)-Params.w)));
fprintf('Uncertain FHorz GE (c): for scale, the possibility of S=2 moves the never-switch prices away from the initial eqm, should NOT be zero: %.3e \n',move0)
fprintf('Uncertain FHorz GE (c): and by more than 2*toleranceGEcondns, so it is not solver noise, should be one: %i \n',move0>2*transpathoptions.toleranceGEcondns)
fprintf('Uncertain FHorz GE (c): for information, %i iterations, price distance %.3e, general eqm distance %.3e \n',SummaryC.iterations,SummaryC.PricePathDist,SummaryC.GEcondnDist)

%% (d) Restarting from the solution
transpathoptionsD=transpathoptions;
transpathoptionsD.PricePathTree0=PricePathTreeC;
[PricePathTreeC2,SummaryC2]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsD, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain FHorz GE (d): restarted from the solution it converges in one iteration, should be one: %i \n',SummaryC2.iterations==1)
fprintf('Uncertain FHorz GE (d): to the same prices, this should be zero: %.3e \n',max([abs(PricePathTreeC2.r-PricePathTreeC.r),abs(PricePathTreeC2.w-PricePathTreeC.w)]))

%% (e) Refusals, each by its message (all but the constraint one are refused before the final eqms are solved)
GEcall=@(tpo,vfo,GEtrans,FnsE,DFn) UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsE, GeneralEqmEqns, GEtrans, Params, DFn, AgeWeightParamNames, tpo, simoptions, vfo, heteroagentoptions);
nE=12;
err=zeros(1,nE);
msg=cell(1,nE);
tpE=transpathoptions; tpE=rmfield(tpE,'GEnewprice');
try; GEcall(tpE,vfoptions,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(1)=contains(ME.message,'GEnewprice must be set'); end
msg{1}='GEnewprice unset';
tpE=transpathoptions; tpE.GEnewprice=1;
try; GEcall(tpE,vfoptions,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(2)=contains(ME.message,'GEnewprice=1 (quasi-Newton) is not implemented'); end
msg{2}='GEnewprice=1 (quasi-Newton)';
% Anderson, price constraints, _tminus1, _tplus1, intermediateEqns and discount factors read at the child USED TO BE
% REFUSED here (stage F2); they are implemented in F3 (parts 9-11), so these now check the refusals that remain for them
tpE=transpathoptions; tpE.GEnewprice=2;
try; GEcall(tpE,vfoptions,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(3)=contains(ME.message,'requires transpathoptions.GEnewprice2.howtoupdate'); end
msg{3}='GEnewprice=2 without GEnewprice2.howtoupdate';
tpE=transpathoptions; tpE.GEnewprice3=rmfield(tpE.GEnewprice3,'howtoupdate');
try; GEcall(tpE,vfoptions,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(4)=contains(ME.message,'requires transpathoptions.GEnewprice3.howtoupdate'); end
msg{4}='GEnewprice=3 without howtoupdate';
tpE=transpathoptions; tpE.constrainAtoB={'w'}; tpE.constrainAtoBlimits.w=[0,0.1];
try; GEcall(tpE,vfoptions,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(5)=contains(ME.message,'is not strictly inside its constraint'); end
msg{5}='price constraints whose bounds exclude the initial guess';
GeneralEqmEqnsE=GeneralEqmEqns;
GeneralEqmEqnsE.CapitalMarket=@(r,K_tminus1,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K_tminus1^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
try; GEcall(transpathoptions,vfoptions,GeneralEqmEqnsE,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(6)=contains(ME.message,'a _tminus1 input is used, so transpathoptions.initialvalues must be given'); end
msg{6}='a _tminus1 without transpathoptions.initialvalues';
FnsE=FnsToEvaluate; FnsE.D=@(aprime,a,z,r_tplus1) r_tplus1*aprime;
GeneralEqmEqnsE=GeneralEqmEqns; GeneralEqmEqnsE.CapitalMarket=@(r,K,N,D_tminus1,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta)+0*D_tminus1;
tpE=transpathoptions; tpE.initialvalues.D=0;
try; GEcall(tpE,vfoptions,GeneralEqmEqnsE,FnsE,DiscountFactorParamNames); catch ME; err(7)=contains(ME.message,'D_tminus1 is used, but FnsToEvaluate.D uses _tplus1'); end
msg{7}='D_tminus1 of a FnsToEvaluate D that uses _tplus1';
try; GEcall(transpathoptions,vfoptions,GeneralEqmEqns,FnsToEvaluate,{'beta','r','r_tplus1'}); catch ME; err(8)=contains(ME.message,'would be applied twice, at the node and at the child'); end
msg{8}='r and r_tplus1 together in DiscountFactorParamNames';
try; GEcall(transpathoptions,vfoptions,GeneralEqmEqns,FnsToEvaluate,{'beta','sigma_tplus1'}); catch ME; err(9)=contains(ME.message,'is not a price in PricePathShaper or a parameter in ParamPath'); end
msg{9}='sigma_tplus1 (sigma is not on the path) in DiscountFactorParamNames';
vfoE=vfoptions; vfoE.localsearch=1; vfoE.nlocalsearch=1;
try; GEcall(transpathoptions,vfoE,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(10)=contains(ME.message,'vfoptions.localsearch is not yet implemented'); end
msg{10}='local search (stage F4)';
tpE=transpathoptions; tpE.fastOLG=2;
try; GEcall(tpE,vfoptions,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(11)=contains(ME.message,'transpathoptions.fastOLG=2 is not yet implemented'); end
msg{11}='fastOLG=2 (stage F4)';
vfoE=vfoptions; vfoE.exoticpreferences='QuasiHyperbolic';
try; GEcall(transpathoptions,vfoE,GeneralEqmEqns,FnsToEvaluate,DiscountFactorParamNames); catch ME; err(12)=contains(ME.message,'vfoptions.exoticpreferences is not yet implemented for uncertain transition paths'); end
msg{12}='exoticpreferences';
for ee=1:nE
    fprintf('Uncertain FHorz GE (e): UncertainTransitionPath_FHorz refuses %s, by name, should be one: %i \n',msg{ee},err(ee))
end

output=struct();

end
