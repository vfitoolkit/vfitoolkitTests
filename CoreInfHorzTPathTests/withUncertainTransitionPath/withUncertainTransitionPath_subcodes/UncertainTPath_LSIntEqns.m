function output=UncertainTPath_LSIntEqns(T,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% Local search and intermediate eqns on the tree, UncertainTransitionPath_InfHorz. Written 2026-11-05, before the first
% run. The model is that of UncertainTPath_GE (households with a Cobb-Douglas firm, the regimes differing in TFP:
% S=1 0.5, S=2 0.49), on its absorbing tree (S=2 for good, possible from period 3 to 10, probability 1/4 each period).
% LOCAL SEARCH (vfoptions.localsearch=1): restricted sweeps around a reference per node, certified by a standard sweep.
%   (a) shooting with local search against shooting without: both converged, the same prices to 10*toleranceGEprices;
%       restricted sweeps were taken and certified (the counts); and the summary is of a STANDARD sweep, so its V and
%       Policy are those of ValueFnOnUncertainTransPath_InfHorz with no local search at the returned prices, exactly
%   (b) Anderson with local search: the same
%   (c) ValueFnOnUncertainTransPath_InfHorz on its own: localsearch=1 is ignored with a warning (the same V as 0,
%       exactly); localsearch=2 forces a stay-put window, on a single history exactly ValueFnOnTransPath_InfHorz with
%       localsearch=2, and different from the unrestricted V (the restriction has teeth)
% INTERMEDIATE EQNS (transpathoptions.intermediateEqns):
%   (d) the same model with output Y as an intermediate eqn (CapitalMarket and LabourMarket written in Y): the same
%       prices as (a) to 10*toleranceGEprices, and the conditions rebuilt by hand (Y by hand) are the ones reported
%   (e) a SURE switch (one history) against TransitionPath_InfHorz with the same intermediate eqns
%   (f) one evaluation (maxiter=1) with a chain of intermediate eqns that uses _tplus1, Rnext=1+r_tplus1 and
%       Dnext=1/Rnext, in a condition ql-Dnext: at a branching node it holds in expectation over the children
%       (by hand), which differs from the condition at the expected r (the DISCRIMINATOR); the chain also tests that
%       an intermediate eqn that uses another one that uses _tplus1 is evaluated per child
%   (g) refusals, each by its message
% Shooting stops at a tolerance, so two solves of the same model are held to 10 times toleranceGEprices, as in
% UncertainTPath_GE.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;

FnsToEvaluate.K=@(aprime,a,z) a;
FnsToEvaluate.N=@(aprime,a,z) z;
GeneralEqmEqns.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

%% The initial stationary eqm (S=1, firmA=0.5)
heteroagentoptions=struct();
Params.r=0.04; Params.w=1;
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a, n_z, 0, pi_z, d_grid, a_grid, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptions, simoptions, vfoptions);
fprintf('Uncertain LS/intermediate: initial stationary eqm r=%2.8f, w=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_init,n_d,n_a,n_z,pi_z,simoptions,Params,[]);

howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate=howtoupdate;
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=2*10^(-3); % the two-price stationary solves stop at about 1e-3 (UncertainTPath_GE)
tolP=transpathoptions.toleranceGEprices;
tolC=transpathoptions.toleranceGEcondns;
fprintf('Uncertain LS/intermediate: toleranceGEprices %.3e, toleranceGEcondns %.3e \n',tolP,tolC)
PricePathShaper.r=min((1:1:T)/30,1);
PricePathShaper.w=min((1:1:T)/30,1);

ParamPathC.firmA=[0.5*ones(1,T); 0.49*ones(1,T)];
pi_SC=[0.75,0.25;0,1]; T1=2; T2=9;

vfoptionsLS=vfoptions;
vfoptionsLS.localsearch=1;
vfoptionsLS.nlocalsearch=3;

%% (a) Shooting, without and with local search
[PricePathTreeA0,SummaryA0]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
[PricePathTreeA,SummaryA]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptionsLS, heteroagentoptions);
TreeC=SummaryA.Tree;
nNodesC=TreeC.nNodes;
fprintf('Uncertain LS (a): %i histories; without local search converged before maxiter, should be one: %i \n',TreeC.nHist,SummaryA0.iterations<transpathoptions.maxiter)
fprintf('Uncertain LS (a): with local search converged before maxiter, should be one: %i \n',SummaryA.iterations<transpathoptions.maxiter)
fprintf('Uncertain LS (a): restricted sweeps were taken, should be one: %i \n',SummaryA.localsearch.nrestricted>0)
fprintf('Uncertain LS (a): and certified by at least one standard verification sweep, should be one: %i \n',SummaryA.localsearch.nverify>=1)
fprintf('Uncertain LS (a): the last sweep, which the summary is of, was a standard one, should be one: %i \n',SummaryA.localsearch.lastsweeprestricted==0)
diffA=max([max(abs(PricePathTreeA.r-PricePathTreeA0.r)),max(abs(PricePathTreeA.w-PricePathTreeA0.w))]);
fprintf('Uncertain LS (a): with vs without local search, prices within 10*toleranceGEprices, should be one: %i \n',diffA<=10*tolP)
fprintf('Uncertain LS (a): for information, that difference is %.3e; %i iterations with local search (%i restricted, %i verification), %i without \n',diffA,SummaryA.iterations,SummaryA.localsearch.nrestricted,SummaryA.localsearch.nverify,SummaryA0.iterations)
% The summary is of a standard sweep: V and Policy at the returned prices with no local search
V_finalA=struct(); Policy_finalA=struct();
for SS=1:2
    V_finalA.(['S',num2str(SS)])=SummaryA.finaleqm.(['S',num2str(SS)]).V_final;
    Policy_finalA.(['S',num2str(SS)])=SummaryA.finaleqm.(['S',num2str(SS)]).Policy_final;
end
[VA,PA]=ValueFnOnUncertainTransPath_InfHorz(PricePathTreeA, ParamPathC, T, T1, T2, pi_SC, 1, V_finalA, Policy_finalA, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VA=gather(reshape(VA,[N_a*N_z,nNodesC])); VS=gather(reshape(SummaryA.VPathTree,[N_a*N_z,nNodesC]));
dd=abs(VA-VS); dd(VA==VS)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain LS (a): the summary V is the unrestricted V at the returned prices, this should be zero: %.3e \n',max(dd,[],'all'))
fprintf('Uncertain LS (a): and so is the summary Policy, this should be zero: %.3e \n',max(abs(gather(PA(:))-gather(SummaryA.PolicyPathTree(:)))))

%% (b) Anderson with local search
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=howtoupdate;
transpathoptionsA.maxiter=200;
[PricePathTreeB,SummaryB]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptionsA, simoptions, vfoptionsLS, heteroagentoptions);
fprintf('Uncertain LS (b): Anderson with local search converged, should be one: %i \n',SummaryB.GEcondnDist<=tolC)
fprintf('Uncertain LS (b): and took more than one iteration, should be one: %i \n',SummaryB.iterations>1)
fprintf('Uncertain LS (b): and was certified by at least one unrestricted sweep, should be one: %i \n',SummaryB.localsearch.nverify>=1)
[VB,PB]=ValueFnOnUncertainTransPath_InfHorz(PricePathTreeB, ParamPathC, T, T1, T2, pi_SC, 1, V_finalA, Policy_finalA, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VB=gather(reshape(VB,[N_a*N_z,nNodesC])); VSB=gather(reshape(SummaryB.VPathTree,[N_a*N_z,nNodesC]));
dd=abs(VB-VSB); dd(VB==VSB)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain LS (b): the summary V is the unrestricted V at the returned prices, this should be zero: %.3e \n',max(dd,[],'all'))
fprintf('Uncertain LS (b): and so is the summary Policy, this should be zero: %.3e \n',max(abs(gather(PB(:))-gather(SummaryB.PolicyPathTree(:)))))
diffB=max([max(abs(PricePathTreeB.r-PricePathTreeA.r)),max(abs(PricePathTreeB.w-PricePathTreeA.w))]);
fprintf('Uncertain LS (b): for information, Anderson vs shooting (both with local search) prices differ by %.3e, after %i Anderson iterations \n',diffB,SummaryB.iterations)

%% (c) ValueFnOnUncertainTransPath_InfHorz on its own
vfoptionsC=vfoptions; vfoptionsC.localsearch=1; vfoptionsC.nlocalsearch=1;
lastwarn('');
VC1=ValueFnOnUncertainTransPath_InfHorz(PricePathTreeA, ParamPathC, T, T1, T2, pi_SC, 1, V_finalA, Policy_finalA, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsC);
fprintf('Uncertain LS (c): localsearch=1 warns that it is ignored, should be one: %i \n',contains(lastwarn,'localsearch=1 is IGNORED here'))
VC1=gather(reshape(VC1,[N_a*N_z,nNodesC]));
dd=abs(VC1-VA); dd(VC1==VA)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain LS (c): and gives the unrestricted V, this should be zero: %.3e \n',max(dd,[],'all'))
% localsearch=2 on a single history (S stays 1), against ValueFnOnTransPath_InfHorz with localsearch=2
hh0=find(all(TreeC.historyS==1,2));
nodes0=TreeC.histories(hh0,:);
clear PricePathH ParamPathH PricePathD ParamPathD
PricePathH.r=PricePathTreeA.r(nodes0);
PricePathH.w=PricePathTreeA.w(nodes0);
ParamPathH.firmA=0.5*ones(1,T);
vfoptionsC.localsearch=2;
VC2=ValueFnOnUncertainTransPath_InfHorz(PricePathH, ParamPathH, T, 1, 1, eye(2), 1, V_finalA, Policy_finalA, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsC);
VD2=ValueFnOnTransPath_InfHorz(PricePathH, ParamPathH, T, V_finalA.S1, Policy_finalA.S1, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsC);
VC2=gather(reshape(VC2,[N_a*N_z,T])); VD2=gather(reshape(VD2,[N_a*N_z,T]));
dd=abs(VC2-VD2); dd(VC2==VD2)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain LS (c): localsearch=2 on a single history is ValueFnOnTransPath_InfHorz with localsearch=2, this should be zero: %.3e \n',max(dd,[],'all'))
dd=abs(VC2-VA(:,nodes0)); dd(VC2==VA(:,nodes0))=0; dd(isnan(dd))=0;
fprintf('Uncertain LS (c): for scale, the forced stay-put window (nlocalsearch=1) against the unrestricted V on that history, should NOT be zero: %.3e \n',max(dd,[],'all'))

%% (d) Output Y as an intermediate eqn
GeneralEqmEqnsI.CapitalMarket=@(r,Y,K,firmalpha,firmdelta) r-(firmalpha*Y/K-firmdelta);
GeneralEqmEqnsI.LabourMarket=@(w,Y,N,firmalpha) w-((1-firmalpha)*Y/N);
transpathoptionsI=transpathoptions;
transpathoptionsI.intermediateEqns.Y=@(K,N,firmA,firmalpha) firmA*(K^firmalpha)*(N^(1-firmalpha));
heteroagentoptionsI=heteroagentoptions;
heteroagentoptionsI.intermediateEqns.Y=transpathoptionsI.intermediateEqns.Y;
[PricePathTreeD,SummaryD]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI, GeneralEqmEqnsI, Params, DiscountFactorParamNames, transpathoptionsI, simoptions, vfoptions, heteroagentoptionsI);
fprintf('Uncertain intermediate (d): converged before maxiter, should be one: %i \n',SummaryD.iterations<transpathoptions.maxiter)
diffD=max([max(abs(PricePathTreeD.r-PricePathTreeA0.r)),max(abs(PricePathTreeD.w-PricePathTreeA0.w))]);
fprintf('Uncertain intermediate (d): the same model as (a) written with Y, the same prices within 10*toleranceGEprices, should be one: %i \n',diffD<=10*tolP)
fprintf('Uncertain intermediate (d): for information, that difference is %.3e, after %i iterations \n',diffD,SummaryD.iterations)
[~,~,ParamNodesC]=UncertainTransPathTree(T,T1,T2,pi_SC,1,[],ParamPathC);
AggD=SummaryD.AggVarsTree;
GEhandD=zeros(2,nNodesC);
for kk=1:nNodesC
    K=AggD.K.Mean(kk); N=AggD.N.Mean(kk);
    Y=transpathoptionsI.intermediateEqns.Y(K,N,ParamNodesC.firmA(kk),Params.firmalpha);
    GEhandD(1,kk)=GeneralEqmEqnsI.CapitalMarket(PricePathTreeD.r(kk),Y,K,Params.firmalpha,Params.firmdelta);
    GEhandD(2,kk)=GeneralEqmEqnsI.LabourMarket(PricePathTreeD.w(kk),Y,N,Params.firmalpha);
end
GErepD=[SummaryD.GEcondnTree.CapitalMarket;SummaryD.GEcondnTree.LabourMarket];
fprintf('Uncertain intermediate (d): the conditions rebuilt by hand (Y by hand) are the ones reported, this should be zero: %.3e \n',max(abs(GEhandD-GErepD),[],'all'))

%% (e) A sure switch to S=2 in period 2, against TransitionPath_InfHorz with the same intermediate eqns
pi_SE=[0,1;0,1];
[PricePathTreeE,SummaryE]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, 1, 1, pi_SE, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI, GeneralEqmEqnsI, Params, DiscountFactorParamNames, transpathoptionsI, simoptions, vfoptions, heteroagentoptionsI);
fprintf('Uncertain intermediate (e): a sure switch is a single history, converged before maxiter, should be one: %i \n',SummaryE.Tree.nHist==1 && SummaryE.iterations<transpathoptions.maxiter)
pf2=SummaryE.finaleqm.S2.p_eqm;
clear PricePath0D ParamPathD
PricePath0D.r=[Params.r, Params.r+(pf2.r-Params.r)*PricePathShaper.r(2:T)];
PricePath0D.w=[Params.w, Params.w+(pf2.w-Params.w)*PricePathShaper.w(2:T)];
ParamPathD.firmA=[0.5, 0.49*ones(1,T-1)];
PricePathDE=TransitionPath_InfHorz(PricePath0D, ParamPathD, T, SummaryE.finaleqm.S2.V_final, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsI, Params, DiscountFactorParamNames, transpathoptionsI, simoptions, vfoptions, []);
diffE=max([max(abs(PricePathTreeE.r-gather(PricePathDE.r(:))')),max(abs(PricePathTreeE.w-gather(PricePathDE.w(:))'))]);
fprintf('Uncertain intermediate (e): the tree vs TransitionPath_InfHorz with intermediateEqns, prices within 10*toleranceGEprices, should be one: %i \n',diffE<=10*tolP)
fprintf('Uncertain intermediate (e): for information, that difference is %.3e, after %i tree iterations \n',diffE,SummaryE.iterations)

%% (f) One evaluation with a chain of intermediate eqns that uses _tplus1
GeneralEqmEqnsF=GeneralEqmEqnsI;
GeneralEqmEqnsF.Discount=@(ql,r) ql-1/(1+r); % the final eqms
GeneralEqmEqns_TransitionF=GeneralEqmEqnsI;
GeneralEqmEqns_TransitionF.Discount=@(ql,Dnext) ql-Dnext;
transpathoptionsF=transpathoptionsI;
transpathoptionsF.intermediateEqns.Rnext=@(r_tplus1) 1+r_tplus1;
transpathoptionsF.intermediateEqns.Dnext=@(Rnext) 1/Rnext;
transpathoptionsF.GEnewprice3.howtoupdate=[howtoupdate; {'Discount','ql',0,1}];
transpathoptionsF.maxiter=1;
ParamsF=Params;
ParamsF.ql=1/(1+Params.r);
PricePathShaperF=PricePathShaper;
PricePathShaperF.ql=PricePathShaper.r;
PricePathTree0=PricePathTreeD;
PricePathTree0.ql=1./(1+PricePathTreeD.r);
kb=find(TreeC.nodenchild>1,1);
childrenkb=TreeC.nodefirstchild(kb)+(0:TreeC.nodenchild(kb)-1);
kS2=childrenkb(TreeC.nodeS(childrenkb)==2);
PricePathTree0.r(kS2)=PricePathTree0.r(kS2)+0.05; % the children of the first branching node far apart
transpathoptionsF.PricePathTree0=PricePathTree0;
[PricePathTreeF,SummaryF]=UncertainTransitionPath_InfHorz(PricePathShaperF, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqnsF, GeneralEqmEqns_TransitionF, ParamsF, DiscountFactorParamNames, transpathoptionsF, simoptions, vfoptions, heteroagentoptionsI);
notlastC=(TreeC.nodedate<T)';
fprintf('Uncertain intermediate (f): maxiter=1 returns the prices it evaluated, this should be zero: %.3e \n',max(abs(PricePathTreeF.r(notlastC)-PricePathTree0.r(notlastC))))
handF=zeros(1,nNodesC);
for kk=1:nNodesC
    ql=PricePathTreeF.ql(kk);
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
        Rnext=transpathoptionsF.intermediateEqns.Rnext(PricePathTreeF.r(children(cc)));
        Dnext=transpathoptionsF.intermediateEqns.Dnext(Rnext);
        handF(kk)=handF(kk)+wc(cc)*GeneralEqmEqns_TransitionF.Discount(ql,Dnext);
    end
end
fprintf('Uncertain intermediate (f): Discount, the expectation over the children of the chain Rnext, Dnext by hand is the one reported, this should be zero: %.3e \n',max(abs(handF-SummaryF.GEcondnTree.Discount)))
wkb=reshape(TreeC.nodecondprob(childrenkb),1,[]);
plugF=GeneralEqmEqns_TransitionF.Discount(PricePathTreeF.ql(kb),1/(1+sum(wkb.*PricePathTreeF.r(childrenkb))));
fprintf('Uncertain intermediate (f): Dnext is nonlinear in r_tplus1, so at the branching node the condition at the expected r differs from the one reported by more than 1e-6, should be one: %i \n',abs(plugF-SummaryF.GEcondnTree.Discount(kb))>1e-6)
fprintf('Uncertain intermediate (f): for information, that difference is %.3e \n',abs(plugF-SummaryF.GEcondnTree.Discount(kb)))
AggF=SummaryF.AggVarsTree;
GEhandF=zeros(1,nNodesC);
for kk=1:nNodesC
    K=AggF.K.Mean(kk); N=AggF.N.Mean(kk);
    Y=transpathoptionsF.intermediateEqns.Y(K,N,ParamNodesC.firmA(kk),Params.firmalpha);
    GEhandF(kk)=GeneralEqmEqnsI.CapitalMarket(PricePathTreeF.r(kk),Y,K,Params.firmalpha,Params.firmdelta);
end
fprintf('Uncertain intermediate (f): and CapitalMarket, with Y at the node, by hand is the one reported, this should be zero: %.3e \n',max(abs(GEhandF-SummaryF.GEcondnTree.CapitalMarket)))

%% (g) Refusals, each by its message
err=zeros(1,3);
vfoptionsE=vfoptions; vfoptionsE.localsearch=1;
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptionsE, heteroagentoptions);
catch ME
    err(1)=contains(ME.message,'When using vfoptions.localsearch=1 you must set vfoptions.nlocalsearch');
end
vfoptionsE.nlocalsearch=N_a;
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptionsE, heteroagentoptions);
catch ME
    err(2)=contains(ME.message,'vfoptions.nlocalsearch is too large for the grid');
end
vfoptionsE.nlocalsearch=1; vfoptionsE.nlocalsearchup=0;
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, Params, DiscountFactorParamNames, transpathoptions, simoptions, vfoptionsE, heteroagentoptions);
catch ME
    err(3)=contains(ME.message,'nlocalsearchup and vfoptions.nlocalsearchdown must both be at least 1');
end
fprintf('Uncertain LS (g): localsearch=1 without nlocalsearch is refused by name, should be one: %i \n',err(1))
fprintf('Uncertain LS (g): an nlocalsearch too large for the grid is refused by name, should be one: %i \n',err(2))
fprintf('Uncertain LS (g): nlocalsearchup=0 is refused by name, should be one: %i \n',err(3))

output=struct();

end
