function output=UncertainTPathFHorz_SDF(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% A stochastic discount factor on the FHorz tree: name_tplus1 in DiscountFactorParamNames (stage F3). Written
% 2026-11-06, before the first run. The FHorz port of UncertainTPath_SDF (CoreInfHorzUncertainTPathTests doPart(10)).
% The rule (see ValueFnOnUncertainTransPath_FHorz): 'name' is read at the node, 'name_tplus1' at the next node; at a
% branching node it goes inside the expectation, next period value fn sum_c nodecondprob(c)*m_c*V_c, at every age.
% The households of the parent bank; mS is a price and betaS a parameter, both depending on S, and the children of a
% branching node have different values of each.
%   (b) a single history, given prices: m_tplus1 (a price) and betaS_tplus1 (a parameter) against
%       ValueFnOnTransPath_Case1_FHorz with the discount factor shifted back one period by hand, as an identity
%   (c) the absorbing tree, given prices, m_tplus1: every node before T is one step from its children, rebuilt by hand
%       with a 2-period ValueFnOnTransPath_Case1_FHorz; at a branching node from sum_c p_c m_c V_c, at a node with one
%       child c with the discount beta*m_c. As an identity, every age.
%   (d) the DISCRIMINATOR, at the first branching node: folding m out of the expectation, (sum p m)(sum p V), and
%       leaving it out, sum p V, both give a different V
%   (e) betaS with no _tplus1 is still read at the node: at a branching node V is one step from sum_c p_c V_c with
%       the node's own betaS, as an identity
%   (f) general eqm, households who discount with beta*m_tplus1, m a price with condition m-C_tminus1/C on the path
%       and m-1 in the final eqms: converged; the final eqm is the stationary eqm solved by hand with {'beta','m'};
%       the conditions rebuilt by hand are the ones reported and within tolerance; Anderson the same
%   (g) refusals, each by its message
% (f) uses the grid interpolation layer and damps m to factor 0.2 from the start: UncertainTPath_SDF needed both
% (C_tminus1/C turns the jumps of C, a step function of the prices on a pure discretization, into jumps of m, and a
% full step on m overshoots through the households' discount).

vfoptions=struct();
simoptions=struct();
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);
NV=N_a*N_z*N_j;
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;

tvec=1:1:T;
mS=[1+0.02*sin(tvec/3); 1-0.03*cos(tvec/4)]; % a price that depends on S
betaS=[0.95*ones(1,T); 0.93+0.005*sin(tvec/5)]; % a parameter that depends on S
PricePathU.r=rS;
PricePathU.m=mS;
ParamPathU.w=wS;
ParamPathU.betaS=betaS;
Params.m=1;
Params.betaS=Params.beta;

%% (b) A single history, given prices: S=1 in period 1, S=2 from period 2
pi_SB=[0,1;0,1]; T1=1; T2=1;
Spath=[1,2*ones(1,T-1)];
idx=Spath+2*(tvec-1); % linear index of (S_t,t) in an N_S-by-T array
clear PricePathD ParamPathD
PricePathD.r=rS(idx);
ParamPathD.w=wS(idx);
% m_tplus1
DFnames={'beta','m_tplus1'};
ParamsS=Params; ParamsS.r=rS(2,T); ParamsS.w=wS(2,T); ParamsS.m=mS(2,T);
[V_finalB.S2,Policy_finalB.S2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'beta','m'},[],vfoptions);
VtreeB=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_SB, 1, V_finalB, Policy_finalB, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DFnames, ReturnFn, transpathoptionsbaseline, vfoptions);
mD=mS(idx);
ParamPathD.mlead=[mD(2:T),1]; % period t discounts with m_{t+1}; the period T entry is never used
VD=ValueFnOnTransPath_Case1_FHorz(PricePathD, ParamPathD, T, V_finalB.S2, Policy_finalB.S2, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VtreeB=gather(reshape(VtreeB,[NV,T])); VD=gather(reshape(VD,[NV,T]));
dd=abs(VtreeB-VD); dd(VtreeB==VD)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz SDF (b): a single history with m_tplus1 is ValueFnOnTransPath_Case1_FHorz with m shifted back by hand, V, this should be zero: %.3e \n',max(dd,[],'all'))
ParamPathD.mlead=mD; % the same, read at the node: the shift has teeth
VDnode=ValueFnOnTransPath_Case1_FHorz(PricePathD, ParamPathD, T, V_finalB.S2, Policy_finalB.S2, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VDnode=gather(reshape(VDnode,[NV,T]));
dd=abs(VtreeB-VDnode); dd(VtreeB==VDnode)=0; dd(isnan(dd))=0;
fprintf('Uncertain FHorz SDF (b): for scale, against m read at the node instead, V differs by, should NOT be zero: %.3e \n',max(dd,[],'all'))
% betaS_tplus1, a parameter
ParamsS=Params; ParamsS.r=rS(2,T); ParamsS.w=wS(2,T); ParamsS.betaS=betaS(2,T);
[V_finalB2.S2,Policy_finalB2.S2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'betaS'},[],vfoptions);
VtreeB2=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_SB, 1, V_finalB2, Policy_finalB2, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'betaS_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
betaD=betaS(idx);
ParamPathD.betalead=[betaD(2:T),1];
VD2=ValueFnOnTransPath_Case1_FHorz(PricePathD, ParamPathD, T, V_finalB2.S2, Policy_finalB2.S2, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'betalead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VtreeB2=gather(reshape(VtreeB2,[NV,T])); VD2=gather(reshape(VD2,[NV,T]));
dd=abs(VtreeB2-VD2); dd(VtreeB2==VD2)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz SDF (b): a single history with betaS_tplus1 (a parameter that depends on S) is ValueFnOnTransPath_Case1_FHorz with betaS shifted back by hand, V, this should be zero: %.3e \n',max(dd,[],'all'))

%% (c) The absorbing tree, given prices, m_tplus1: every node one step from its children
P=0.75;
pi_SC=[P,1-P;0,1]; T1=2; T2=5;
V_finalC=struct(); Policy_finalC=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T); ParamsS.m=mS(SS,T);
    [V_finalC.(['S',num2str(SS)]),Policy_finalC.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'beta','m'},[],vfoptions);
end
Policyshape=size(Policy_finalC.S1);
NP=prod(Policyshape);
[VPathC,PolicyPathC]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'beta','m_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
[TreeC,PriceNodesC,ParamNodesC]=UncertainTransPathTree(T,T1,T2,pi_SC,1,PricePathU,ParamPathU);
nNodesC=TreeC.nNodes;
VC=gather(reshape(VPathC,[NV,nNodesC]));
PC=gather(reshape(PolicyPathC,[NP,nNodesC]));
errbranch=0; errsingle=0; nbranch=0; nsingle=0;
for kk=find(TreeC.nodedate<T)'
    children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
    clear PricePathH ParamPathH
    PricePathH.r=PriceNodesC.r([kk,kk]);
    ParamPathH.w=ParamNodesC.w([kk,kk]);
    if TreeC.nodenchild(kk)>1
        % written as the command writes it, so that it is an identity
        Vnext=sum(VC(:,children).*(reshape(TreeC.nodecondprob(children),1,[]).*PriceNodesC.m(children)),2);
        ParamPathH.mlead=[1,1];
    else
        Vnext=VC(:,children);
        ParamPathH.mlead=PriceNodesC.m([children,children]);
    end
    Vh=ValueFnOnTransPath_Case1_FHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vnext,[n_a,n_z,N_j])), reshape(PC(:,children(1)),Policyshape), Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
    Vh=gather(reshape(Vh,[NV,2]));
    dd=abs(Vh(:,1)-VC(:,kk)); dd(Vh(:,1)==VC(:,kk))=0; dd(isnan(dd))=Inf;
    if TreeC.nodenchild(kk)>1
        errbranch=max(errbranch,max(dd)); nbranch=nbranch+1;
    else
        errsingle=max(errsingle,max(dd)); nsingle=nsingle+1;
    end
end
fprintf('Uncertain FHorz SDF (c): %i histories, %i branching nodes, %i nodes with one child \n',TreeC.nHist,nbranch,nsingle)
fprintf('Uncertain FHorz SDF (c): at every branching node V is one step from sum_c p_c m_c V_c (every age), this should be zero: %.3e \n',errbranch)
fprintf('Uncertain FHorz SDF (c): at every node with one child c, V is one step from V_c with the discount beta*m_c, this should be zero: %.3e \n',errsingle)

%% (d) The discriminator, at the first branching node
kb=find(TreeC.nodenchild>1,1);
childrenkb=TreeC.nodefirstchild(kb)+(0:TreeC.nodenchild(kb)-1);
pkb=reshape(TreeC.nodecondprob(childrenkb),1,[]);
mkb=PriceNodesC.m(childrenkb);
fprintf('Uncertain FHorz SDF (d): for scale, the two children of the first branching node have m differing by, should NOT be zero: %.3e \n',abs(diff(mkb)))
clear PricePathH ParamPathH
PricePathH.r=PriceNodesC.r([kb,kb]);
ParamPathH.w=ParamNodesC.w([kb,kb]);
ParamPathH.mlead=[1,1];
Vfold=sum(pkb.*mkb)*sum(VC(:,childrenkb).*pkb,2);
Vh=ValueFnOnTransPath_Case1_FHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vfold,[n_a,n_z,N_j])), reshape(PC(:,childrenkb(1)),Policyshape), Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
Vh=gather(reshape(Vh,[NV,2]));
dd=abs(Vh(:,1)-VC(:,kb)); dd(isnan(dd))=0;
fprintf('Uncertain FHorz SDF (d): folding m out of the expectation, (sum p m)(sum p V), gives a V that differs by more than 1e-8, should be one: %i \n',max(dd)>1e-8)
fprintf('Uncertain FHorz SDF (d): for information, that difference is %.3e \n',max(dd))
Vplain=sum(VC(:,childrenkb).*pkb,2);
Vh=ValueFnOnTransPath_Case1_FHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vplain,[n_a,n_z,N_j])), reshape(PC(:,childrenkb(1)),Policyshape), Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
Vh=gather(reshape(Vh,[NV,2]));
dd=abs(Vh(:,1)-VC(:,kb)); dd(isnan(dd))=0;
fprintf('Uncertain FHorz SDF (d): leaving m out, sum p V, gives a V that differs by more than 1e-8, should be one: %i \n',max(dd)>1e-8)

%% (e) betaS with no _tplus1 is read at the node
V_finalE=struct(); Policy_finalE=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T); ParamsS.betaS=betaS(SS,T);
    [V_finalE.(['S',num2str(SS)]),Policy_finalE.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'betaS'},[],vfoptions);
end
[VPathE,PolicyPathE]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_SC, 1, V_finalE, Policy_finalE, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'betaS'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VE=gather(reshape(VPathE,[NV,nNodesC]));
PE=gather(reshape(PolicyPathE,[NP,nNodesC]));
errE=0; scaleE=0;
for kk=find(TreeC.nodenchild>1)'
    children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
    clear PricePathH ParamPathH
    PricePathH.r=PriceNodesC.r([kk,kk]);
    ParamPathH.w=ParamNodesC.w([kk,kk]);
    ParamPathH.betaS=ParamNodesC.betaS([kk,kk]);
    Vnext=sum(VE(:,children).*reshape(TreeC.nodecondprob(children),1,[]),2);
    Vh=ValueFnOnTransPath_Case1_FHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vnext,[n_a,n_z,N_j])), reshape(PE(:,children(1)),Policyshape), Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, {'betaS'}, ReturnFn, transpathoptionsbaseline, vfoptions);
    Vh=gather(reshape(Vh,[NV,2]));
    dd=abs(Vh(:,1)-VE(:,kk)); dd(Vh(:,1)==VE(:,kk))=0; dd(isnan(dd))=Inf;
    errE=max(errE,max(dd));
    scaleE=max(scaleE,max(abs(ParamNodesC.betaS(children)-ParamNodesC.betaS(kk))));
end
fprintf('Uncertain FHorz SDF (e): betaS (no _tplus1) at a branching node is the node''s own, V one step from sum_c p_c V_c, this should be zero: %.3e \n',errE)
fprintf('Uncertain FHorz SDF (e): for scale, at the branching nodes a child''s betaS differs from the node''s by, should NOT be zero: %.3e \n',scaleE)

%% (f) General eqm: households who discount with beta*m_tplus1, m=C_tminus1/C on the path
vfoptionsF=vfoptions; vfoptionsF.gridinterplayer=1; vfoptionsF.ngridinterp=10;
simoptionsF=simoptions; simoptionsF.gridinterplayer=1; simoptionsF.ngridinterp=10;
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
FnsToEvaluate0.K=@(aprime,a,z) a;
FnsToEvaluate0.N=@(aprime,a,z,kappa_j) kappa_j*z;
GeneralEqmEqns0.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns0.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
FnsToEvaluate=FnsToEvaluate0;
% consumption; in the parent bank's setup Jr>N_j, so every age is on the working branch of the return fn
FnsToEvaluate.C=@(aprime,a,z,r,w,kappa_j) (1+r)*a+w*kappa_j*z-aprime;
GeneralEqmEqns=GeneralEqmEqns0;
GeneralEqmEqns.SDF=@(m) m-1;
GeneralEqmEqns_Transition=GeneralEqmEqns0;
GeneralEqmEqns_Transition.SDF=@(m,C,C_tminus1) m-C_tminus1/C;
DFnamesF={'beta','m_tplus1'};
heteroagentoptions=struct();
heteroagentoptions.verbose=0;

% The initial stationary eqm (S=1, firmA=0.5), m=1 so the discount is beta
ParamsF=rmfield(Params,'betaS');
ParamsF.m=1;
[p_eqm,~,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate0,GeneralEqmEqns0,ParamsF,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptionsF,vfoptionsF);
ParamsF.r=p_eqm.r;
ParamsF.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsF,DiscountFactorParamNames,[],vfoptionsF);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,ParamsF,simoptionsF);
AggVars0=EvalFnOnAgentDist_AggVars_FHorz_Case1(AgentDist_initial, Policy_init, FnsToEvaluate, ParamsF, [], n_d, n_a, n_z, N_j, d_grid, a_grid, z_grid, simoptionsF);
C0=AggVars0.C.Mean;
fprintf('Uncertain FHorz SDF (f): initial stationary eqm r=%2.8f, w=%2.8f, C=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,C0,GEcondns.CapitalMarket,GEcondns.LabourMarket)

howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1; 'SDF','m',0,0.2};
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate=howtoupdate;
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=10^(-3); % as UncertainTPath_SDF: the three-price stationary solves can stop at a few 1e-3
transpathoptions.initialvalues.C=C0;
tolC=transpathoptions.toleranceGEcondns;
PricePathShaper.r=min((1:1:T)/20,1);
PricePathShaper.w=min((1:1:T)/20,1);
PricePathShaper.m=min((1:1:T)/20,1);

ParamPathF.firmA=[0.5*ones(1,T); 0.52*ones(1,T)];
[PricePathTreeF,SummaryF]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, ParamsF, DFnamesF, AgeWeightParamNames, transpathoptions, simoptionsF, vfoptionsF, heteroagentoptions);
TreeF=SummaryF.Tree;
nNodesF=TreeF.nNodes;
notlastF=(TreeF.nodedate<T)';
fprintf('Uncertain FHorz SDF (f): %i histories, converged before maxiter, should be one: %i \n',TreeF.nHist,SummaryF.iterations<transpathoptions.maxiter)
% The final eqm is the stationary eqm with the _tplus1 removed, solved by hand exactly as the command solves it
Params2=ParamsF; Params2.firmA=0.52;
p_eqm2=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,Params2,{'beta','m'},[],[],[],{'r','w','m'},heteroagentoptions,simoptionsF,vfoptionsF);
fprintf('Uncertain FHorz SDF (f): the S=2 final eqm is the stationary eqm solved by hand with {''beta'',''m''}, this should be zero: %.3e \n',max(abs([SummaryF.finaleqm.S2.p_eqm.r-p_eqm2.r,SummaryF.finaleqm.S2.p_eqm.w-p_eqm2.w,SummaryF.finaleqm.S2.p_eqm.m-p_eqm2.m])))
% The conditions rebuilt by hand, SDF with the parent node's C (initialvalues at the period 1 node)
[~,~,ParamNodesF]=UncertainTransPathTree(T,T1,T2,pi_SC,1,[],ParamPathF);
AggF=SummaryF.AggVarsTree;
GEhandF=zeros(3,nNodesF);
for kk=1:nNodesF
    r=PricePathTreeF.r(kk); w=PricePathTreeF.w(kk); m=PricePathTreeF.m(kk); K=AggF.K.Mean(kk); N=AggF.N.Mean(kk); A=ParamNodesF.firmA(kk);
    GEhandF(1,kk)=GeneralEqmEqns0.CapitalMarket(r,K,N,Params.firmalpha,Params.firmdelta,A);
    GEhandF(2,kk)=GeneralEqmEqns0.LabourMarket(w,K,N,Params.firmalpha,A);
    parentnode=TreeF.nodeparent(kk);
    if parentnode>0
        GEhandF(3,kk)=GeneralEqmEqns_Transition.SDF(m,AggF.C.Mean(kk),AggF.C.Mean(parentnode));
    else
        GEhandF(3,kk)=GeneralEqmEqns_Transition.SDF(m,AggF.C.Mean(kk),C0);
    end
end
GErepF=[SummaryF.GEcondnTree.CapitalMarket;SummaryF.GEcondnTree.LabourMarket;SummaryF.GEcondnTree.SDF];
fprintf('Uncertain FHorz SDF (f): the conditions rebuilt by hand (SDF with the parent node''s C) are the ones reported, this should be zero: %.3e \n',max(abs(GEhandF-GErepF),[],'all'))
fprintf('Uncertain FHorz SDF (f): and are within toleranceGEcondns at every node before T, should be one: %i \n',max(sqrt(sum(GEhandF(:,notlastF).^2,1)))<=tolC)
GEnamesF={'CapitalMarket','LabourMarket','SDF'};
for gg=1:3
    [temp,kmax]=max(abs(GEhandF(gg,:)).*notlastF);
    fprintf('Uncertain FHorz SDF (f): for information, the largest %s condition before T is %.3e, at a node in period %i with S=%i \n',GEnamesF{gg},temp,TreeF.nodedate(kmax),TreeF.nodeS(kmax))
end
fprintf('Uncertain FHorz SDF (f): for information, the S=2 final eqm conditions are %.3e %.3e %.3e (CapitalMarket, LabourMarket, SDF) \n',SummaryF.finaleqm.S2.GeneralEqmCondnValues.CapitalMarket,SummaryF.finaleqm.S2.GeneralEqmCondnValues.LabourMarket,SummaryF.finaleqm.S2.GeneralEqmCondnValues.SDF)
kbF=find(TreeF.nodenchild>1,1);
childrenkbF=TreeF.nodefirstchild(kbF)+(0:TreeF.nodenchild(kbF)-1);
fprintf('Uncertain FHorz SDF (f): for scale, m moves away from one by, should NOT be zero: %.3e \n',max(abs(PricePathTreeF.m-1)))
fprintf('Uncertain FHorz SDF (f): for scale, the children of the first branching node have m differing by, should NOT be zero: %.3e \n',abs(diff(PricePathTreeF.m(childrenkbF))))
fprintf('Uncertain FHorz SDF (f): for information, %i iterations, price distance %.3e, general eqm distance %.3e \n',SummaryF.iterations,SummaryF.PricePathDist,SummaryF.GEcondnDist)
% Anderson
transpathoptionsA=transpathoptions;
transpathoptionsA=rmfield(transpathoptionsA,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=howtoupdate;
transpathoptionsA.maxiter=200;
[PricePathTreeFA,SummaryFA]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, ParamsF, DFnamesF, AgeWeightParamNames, transpathoptionsA, simoptionsF, vfoptionsF, heteroagentoptions);
fprintf('Uncertain FHorz SDF (f): Anderson converged, should be one: %i \n',SummaryFA.GEcondnDist<=tolC)
fprintf('Uncertain FHorz SDF (f): and took more than one iteration, should be one: %i \n',SummaryFA.iterations>1)
AggFA=SummaryFA.AggVarsTree;
handFA=zeros(1,nNodesF);
for kk=1:nNodesF
    parentnode=TreeF.nodeparent(kk);
    if parentnode>0
        handFA(kk)=GeneralEqmEqns_Transition.SDF(PricePathTreeFA.m(kk),AggFA.C.Mean(kk),AggFA.C.Mean(parentnode));
    else
        handFA(kk)=GeneralEqmEqns_Transition.SDF(PricePathTreeFA.m(kk),AggFA.C.Mean(kk),C0);
    end
end
fprintf('Uncertain FHorz SDF (f): Anderson, SDF rebuilt by hand is the one reported, this should be zero: %.3e \n',max(abs(handFA-SummaryFA.GEcondnTree.SDF)))
diffFA=max([max(abs(PricePathTreeFA.r-PricePathTreeF.r)),max(abs(PricePathTreeFA.w-PricePathTreeF.w)),max(abs(PricePathTreeFA.m-PricePathTreeF.m))]);
fprintf('Uncertain FHorz SDF (f): for information, Anderson vs shooting prices differ by %.3e, after %i Anderson iterations \n',diffFA,SummaryFA.iterations)

%% (g) Refusals, each by its message
VFcall=@(PP,DFn) ValueFnOnUncertainTransPath_FHorz(PricePathU, PP, T, 2, 5, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DFn, ReturnFn, transpathoptionsbaseline, vfoptions);
GEcall=@(PP,DFn) UncertainTransitionPath_FHorz(PricePathShaper, PP, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, ParamsF, DFn, AgeWeightParamNames, transpathoptions, simoptionsF, vfoptionsF, heteroagentoptions);
err=zeros(1,7);
try; VFcall(ParamPathU,{'beta','m','m_tplus1'}); catch ME; err(1)=contains(ME.message,'would be applied twice, at the node and at the child'); end
try; VFcall(ParamPathU,{'beta','sigma_tplus1'}); catch ME; err(2)=contains(ME.message,'is not in PricePath or ParamPath (a discount factor read at the child'); end
ParamPathE=ParamPathU; ParamPathE.DFtplus1_UncertainTPath=ones(1,T);
try; VFcall(ParamPathE,{'beta','m_tplus1'}); catch ME; err(3)=contains(ME.message,'DFtplus1_UncertainTPath is a name used internally'); end
try; GEcall(ParamPathF,{'beta','m','m_tplus1'}); catch ME; err(4)=contains(ME.message,'would be applied twice, at the node and at the child'); end
try; GEcall(ParamPathF,{'beta','sigma_tplus1'}); catch ME; err(5)=contains(ME.message,'is not a price in PricePathShaper or a parameter in ParamPath'); end
% a discount factor that depends on age is fine read at the node (it is the FHorz commands' usual form), but not at
% the child
ParamPathE=ParamPathF; ParamPathE.betaj=repmat(0.99*ones(N_j,1),1,T);
try; GEcall(ParamPathE,{'beta','m_tplus1','betaj_tplus1'}); catch ME; err(6)=contains(ME.message,'betaj is not a scalar; a discount factor read at the child must be a scalar'); end
ParamPathE=ParamPathU; ParamPathE.betaj=repmat(0.99*ones(N_j,1),1,T);
try; VFcall(ParamPathE,{'beta','betaj_tplus1'}); catch ME; err(7)=contains(ME.message,'betaj is not a scalar; a discount factor read at the child must be a scalar'); end
fprintf('Uncertain FHorz SDF (g): ValueFnOnUncertainTransPath_FHorz refuses m and m_tplus1 together, by name, should be one: %i \n',err(1))
fprintf('Uncertain FHorz SDF (g): ValueFnOnUncertainTransPath_FHorz refuses sigma_tplus1 (sigma is not on the path), by name, should be one: %i \n',err(2))
fprintf('Uncertain FHorz SDF (g): ValueFnOnUncertainTransPath_FHorz refuses a parameter with the internal name, by name, should be one: %i \n',err(3))
fprintf('Uncertain FHorz SDF (g): UncertainTransitionPath_FHorz refuses m and m_tplus1 together, by name, should be one: %i \n',err(4))
fprintf('Uncertain FHorz SDF (g): UncertainTransitionPath_FHorz refuses sigma_tplus1, by name, should be one: %i \n',err(5))
fprintf('Uncertain FHorz SDF (g): UncertainTransitionPath_FHorz refuses an age-dependent discount factor read at the child, by name, should be one: %i \n',err(6))
fprintf('Uncertain FHorz SDF (g): and so does ValueFnOnUncertainTransPath_FHorz, by name, should be one: %i \n',err(7))

output=struct();

end
