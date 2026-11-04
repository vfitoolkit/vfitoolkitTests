function output=UncertainTPath_SDF(T,rS,wS,n_a,n_z,a_grid,z_grid,pi_z,n_a_GE,a_grid_GE,Params,DiscountFactorParamNames,transpathoptionsbaseline)
% A stochastic discount factor on the tree: name_tplus1 in DiscountFactorParamNames. Written 2026-11-04, before the
% first run. The rule (see ValueFnOnUncertainTransPath_InfHorz): 'name' is read at the node, 'name_tplus1' at the next
% node; at a branching node it goes inside the expectation, next period value fn sum_c nodecondprob(c)*m_c*V_c.
% The households of the parent bank; mS is a price and betaS a parameter, both depending on S, and the children of a
% branching node have different values of each.
%   (a) the commands with no _tplus1 in DiscountFactorParamNames are parts 1-9 (nothing to do here)
%   (b) a single history, given prices: m_tplus1 (a price) and betaS_tplus1 (a parameter) against
%       ValueFnOnTransPath_InfHorz with the discount factor shifted back one period by hand, as an identity
%   (c) the absorbing tree (as ABCP2024), given prices, m_tplus1: every node before T is one step from its children,
%       rebuilt by hand with a 2-period ValueFnOnTransPath_InfHorz; at a branching node from sum_c p_c m_c V_c, at a
%       node with one child c with the discount beta*m_c. As an identity.
%   (d) the DISCRIMINATOR, at the first branching node: folding m out of the expectation, (sum p m)(sum p V), and
%       leaving it out, sum p V, both give a different V
%   (e) betaS with no _tplus1 is still read at the node: at a branching node V is one step from sum_c p_c V_c with
%       the node's own betaS, as an identity
%   (f) general eqm, households who discount with beta*m_tplus1, m a price with condition m-C_tminus1/C on the path
%       and m-1 in the final eqms: converged; the final eqm is the stationary eqm solved by hand with {'beta','m'};
%       the conditions rebuilt by hand are the ones reported and within tolerance; Anderson the same
%   (g) refusals, each by its message

vfoptions=struct();
simoptions=struct();
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);

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
[V_finalB.S2,Policy_finalB.S2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'beta','m'},[],vfoptions);
VtreeB=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_SB, 1, V_finalB, Policy_finalB, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DFnames, ReturnFn, transpathoptionsbaseline, vfoptions);
mD=mS(idx);
ParamPathD.mlead=[mD(2:T),1]; % period t discounts with m_{t+1}; the period T entry is never used
VD=ValueFnOnTransPath_InfHorz(PricePathD, ParamPathD, T, V_finalB.S2, Policy_finalB.S2, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VtreeB=gather(reshape(VtreeB,[N_a*N_z,T])); VD=gather(reshape(VD,[N_a*N_z,T]));
dd=abs(VtreeB-VD); dd(VtreeB==VD)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain SDF (b): a single history with m_tplus1 is ValueFnOnTransPath_InfHorz with m shifted back by hand, V, this should be zero: %.3e \n',max(dd,[],'all'))
ParamPathD.mlead=mD; % the same, read at the node: the shift has teeth
VDnode=ValueFnOnTransPath_InfHorz(PricePathD, ParamPathD, T, V_finalB.S2, Policy_finalB.S2, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VDnode=gather(reshape(VDnode,[N_a*N_z,T]));
dd=abs(VtreeB-VDnode); dd(VtreeB==VDnode)=0; dd(isnan(dd))=0;
fprintf('Uncertain SDF (b): for scale, against m read at the node instead, V differs by, should NOT be zero: %.3e \n',max(dd,[],'all'))
% betaS_tplus1, a parameter
ParamsS=Params; ParamsS.r=rS(2,T); ParamsS.w=wS(2,T); ParamsS.betaS=betaS(2,T);
[V_finalB2.S2,Policy_finalB2.S2]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'betaS'},[],vfoptions);
VtreeB2=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_SB, 1, V_finalB2, Policy_finalB2, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'betaS_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
betaD=betaS(idx);
ParamPathD.betalead=[betaD(2:T),1];
VD2=ValueFnOnTransPath_InfHorz(PricePathD, ParamPathD, T, V_finalB2.S2, Policy_finalB2.S2, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'betalead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VtreeB2=gather(reshape(VtreeB2,[N_a*N_z,T])); VD2=gather(reshape(VD2,[N_a*N_z,T]));
dd=abs(VtreeB2-VD2); dd(VtreeB2==VD2)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain SDF (b): a single history with betaS_tplus1 (a parameter that depends on S) is ValueFnOnTransPath_InfHorz with betaS shifted back by hand, V, this should be zero: %.3e \n',max(dd,[],'all'))

%% (c) The absorbing tree, given prices, m_tplus1: every node one step from its children
P=0.75;
pi_SC=[P,1-P;0,1]; T1=2; T2=5;
V_finalC=struct(); Policy_finalC=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T); ParamsS.m=mS(SS,T);
    [V_finalC.(['S',num2str(SS)]),Policy_finalC.(['S',num2str(SS)])]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'beta','m'},[],vfoptions);
end
[VPathC,PolicyPathC]=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','m_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
[TreeC,PriceNodesC,ParamNodesC]=UncertainTransPathTree(T,T1,T2,pi_SC,1,PricePathU,ParamPathU);
nNodesC=TreeC.nNodes;
VC=gather(reshape(VPathC,[N_a*N_z,nNodesC]));
PC=gather(reshape(PolicyPathC,[1,N_a*N_z,nNodesC]));
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
    Vh=ValueFnOnTransPath_InfHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vnext,[N_a,N_z])), reshape(PC(:,:,children(1)),[1,N_a,N_z]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
    Vh=gather(reshape(Vh,[N_a*N_z,2]));
    dd=abs(Vh(:,1)-VC(:,kk)); dd(Vh(:,1)==VC(:,kk))=0; dd(isnan(dd))=Inf;
    if TreeC.nodenchild(kk)>1
        errbranch=max(errbranch,max(dd)); nbranch=nbranch+1;
    else
        errsingle=max(errsingle,max(dd)); nsingle=nsingle+1;
    end
end
fprintf('Uncertain SDF (c): %i histories, %i branching nodes, %i nodes with one child \n',TreeC.nHist,nbranch,nsingle)
fprintf('Uncertain SDF (c): at every branching node V is one step from sum_c p_c m_c V_c, this should be zero: %.3e \n',errbranch)
fprintf('Uncertain SDF (c): at every node with one child c, V is one step from V_c with the discount beta*m_c, this should be zero: %.3e \n',errsingle)

%% (d) The discriminator, at the first branching node
kb=find(TreeC.nodenchild>1,1);
childrenkb=TreeC.nodefirstchild(kb)+(0:TreeC.nodenchild(kb)-1);
pkb=reshape(TreeC.nodecondprob(childrenkb),1,[]);
mkb=PriceNodesC.m(childrenkb);
fprintf('Uncertain SDF (d): for scale, the two children of the first branching node have m differing by, should NOT be zero: %.3e \n',abs(diff(mkb)))
clear PricePathH ParamPathH
PricePathH.r=PriceNodesC.r([kb,kb]);
ParamPathH.w=ParamNodesC.w([kb,kb]);
ParamPathH.mlead=[1,1];
Vfold=sum(pkb.*mkb)*sum(VC(:,childrenkb).*pkb,2);
Vh=ValueFnOnTransPath_InfHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vfold,[N_a,N_z])), reshape(PC(:,:,childrenkb(1)),[1,N_a,N_z]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
Vh=gather(reshape(Vh,[N_a*N_z,2]));
dd=abs(Vh(:,1)-VC(:,kb)); dd(isnan(dd))=0;
fprintf('Uncertain SDF (d): folding m out of the expectation, (sum p m)(sum p V), gives a V that differs by more than 1e-8, should be one: %i \n',max(dd)>1e-8)
fprintf('Uncertain SDF (d): for information, that difference is %.3e \n',max(dd))
Vplain=sum(VC(:,childrenkb).*pkb,2);
Vh=ValueFnOnTransPath_InfHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vplain,[N_a,N_z])), reshape(PC(:,:,childrenkb(1)),[1,N_a,N_z]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','mlead'}, ReturnFn, transpathoptionsbaseline, vfoptions);
Vh=gather(reshape(Vh,[N_a*N_z,2]));
dd=abs(Vh(:,1)-VC(:,kb)); dd(isnan(dd))=0;
fprintf('Uncertain SDF (d): leaving m out, sum p V, gives a V that differs by more than 1e-8, should be one: %i \n',max(dd)>1e-8)

%% (e) betaS with no _tplus1 is read at the node
V_finalE=struct(); Policy_finalE=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T); ParamsS.betaS=betaS(SS,T);
    [V_finalE.(['S',num2str(SS)]),Policy_finalE.(['S',num2str(SS)])]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,{'betaS'},[],vfoptions);
end
[VPathE,PolicyPathE]=ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_SC, 1, V_finalE, Policy_finalE, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'betaS'}, ReturnFn, transpathoptionsbaseline, vfoptions);
VE=gather(reshape(VPathE,[N_a*N_z,nNodesC]));
PE=gather(reshape(PolicyPathE,[1,N_a*N_z,nNodesC]));
errE=0; scaleE=0;
for kk=find(TreeC.nodenchild>1)'
    children=TreeC.nodefirstchild(kk)+(0:TreeC.nodenchild(kk)-1);
    clear PricePathH ParamPathH
    PricePathH.r=PriceNodesC.r([kk,kk]);
    ParamPathH.w=ParamNodesC.w([kk,kk]);
    ParamPathH.betaS=ParamNodesC.betaS([kk,kk]);
    Vnext=sum(VE(:,children).*reshape(TreeC.nodecondprob(children),1,[]),2);
    Vh=ValueFnOnTransPath_InfHorz(PricePathH, ParamPathH, 2, gpuArray(reshape(Vnext,[N_a,N_z])), reshape(PE(:,:,children(1)),[1,N_a,N_z]), Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'betaS'}, ReturnFn, transpathoptionsbaseline, vfoptions);
    Vh=gather(reshape(Vh,[N_a*N_z,2]));
    dd=abs(Vh(:,1)-VE(:,kk)); dd(Vh(:,1)==VE(:,kk))=0; dd(isnan(dd))=Inf;
    errE=max(errE,max(dd));
    scaleE=max(scaleE,max(abs(ParamNodesC.betaS(children)-ParamNodesC.betaS(kk))));
end
fprintf('Uncertain SDF (e): betaS (no _tplus1) at a branching node is the node''s own, V one step from sum_c p_c V_c, this should be zero: %.3e \n',errE)
fprintf('Uncertain SDF (e): for scale, at the branching nodes a child''s betaS differs from the node''s by, should NOT be zero: %.3e \n',scaleE)

%% (f) General eqm: households who discount with beta*m_tplus1, m=C_tminus1/C on the path
% With the grid interpolation layer, from here on (the initial eqm, the tree, the final eqms, and (g)). The second run
% (2026-11-05, pure discretization, m damped to 0.5) still stalled: the SDF condition at 2.0e-3 in period 36 (S=2), m
% moving by half that each iteration, and the S=2 final eqm at 3.2e-3. On the n_a=201 grid C is a step function of the
% prices, and the ratio C_tminus1/C turns its jumps into jumps of m.
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=10;
simoptions.gridinterplayer=1; simoptions.ngridinterp=10;
Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
FnsToEvaluate0.K=@(aprime,a,z) a;
FnsToEvaluate0.N=@(aprime,a,z) z;
GeneralEqmEqns0.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns0.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));
FnsToEvaluate=FnsToEvaluate0;
FnsToEvaluate.C=@(aprime,a,z,r,w) (1+r)*a+w*z-aprime;
GeneralEqmEqns=GeneralEqmEqns0;
GeneralEqmEqns.SDF=@(m) m-1;
GeneralEqmEqns_Transition=GeneralEqmEqns0;
GeneralEqmEqns_Transition.SDF=@(m,C,C_tminus1) m-C_tminus1/C;
DFnamesF={'beta','m_tplus1'};

% The initial stationary eqm (S=1, firmA=0.5), m=1 so the discount is beta
heteroagentoptions=struct();
Params.r=0.04; Params.w=1; Params.m=1;
[p_eqm,GEcondns]=HeteroAgentStationaryEqm_InfHorz(n_d, n_a_GE, n_z, 0, pi_z, d_grid, a_grid_GE, z_grid, ReturnFn, FnsToEvaluate0, GeneralEqmEqns0, Params, DiscountFactorParamNames, [], [], [], {'r','w'}, heteroagentoptions, simoptions, vfoptions);
Params.r=p_eqm.r;
Params.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_InfHorz(n_d,n_a_GE,n_z,d_grid,a_grid_GE,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_InfHorz(Policy_init,n_d,n_a_GE,n_z,pi_z,simoptions,Params,[]);
AggVars0=EvalFnOnAgentDist_AggVars_InfHorz(AgentDist_initial, Policy_init, FnsToEvaluate, Params, [], n_d, n_a_GE, n_z, d_grid, a_grid_GE, z_grid, simoptions);
C0=AggVars0.C.Mean;
fprintf('Uncertain SDF (f): initial stationary eqm r=%2.8f, w=%2.8f, C=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,C0,GEcondns.CapitalMarket,GEcondns.LabourMarket)

% m is damped. The first run (2026-11-05) used factor 1, m=C_tminus1/C at the old prices, and shooting stalled at
% 4.05e-3 for 1000 iterations with the price move equal to the SDF residual (m chasing its target). m feeds straight back
% through the households' discount, beta*m_tplus1: a higher m_{t+1} raises saving in t and lowers C_t, which moves both
% m_t and m_{t+1}, so a full step overshoots. The prints after the solve say which condition is left, where.
% Third run (2026-11-05, GI, factor 0.5): shooting still stalled with the SDF condition at 1.3e-3 (period 33, S=2) and m
% moving by half of it each iteration, while Anderson on the same shooting map converged in 28 iterations: the fixed
% point is there, and 0.5 is still too long a step for shooting. So 0.2.
howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1; 'SDF','m',0,0.2};
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate=howtoupdate;
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=10^(-3); % as UncertainTPath_tplus1: the three-price stationary solves stop at a few 1e-3
transpathoptions.initialvalues.C=C0;
tolC=transpathoptions.toleranceGEcondns;
PricePathShaper.r=min((1:1:T)/30,1);
PricePathShaper.w=min((1:1:T)/30,1);
PricePathShaper.m=min((1:1:T)/30,1);

ParamPathF.firmA=[0.5*ones(1,T); 0.49*ones(1,T)];
T1=2; T2=9;
[PricePathTreeF,SummaryF]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a_GE,n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DFnamesF, transpathoptions, simoptions, vfoptions, heteroagentoptions);
TreeF=SummaryF.Tree;
nNodesF=TreeF.nNodes;
notlastF=(TreeF.nodedate<T)';
fprintf('Uncertain SDF (f): %i histories, converged before maxiter, should be one: %i \n',TreeF.nHist,SummaryF.iterations<transpathoptions.maxiter)
% The final eqm is the stationary eqm with the _tplus1 removed, solved by hand exactly as the command solves it
Params2=Params; Params2.firmA=0.49;
p_eqm2=HeteroAgentStationaryEqm_InfHorz(n_d, n_a_GE, n_z, 0, pi_z, d_grid, a_grid_GE, z_grid, ReturnFn, FnsToEvaluate, GeneralEqmEqns, Params2, {'beta','m'}, [], [], [], {'r','w','m'}, heteroagentoptions, simoptions, vfoptions);
fprintf('Uncertain SDF (f): the S=2 final eqm is the stationary eqm solved by hand with {''beta'',''m''}, this should be zero: %.3e \n',max(abs([SummaryF.finaleqm.S2.p_eqm.r-p_eqm2.r,SummaryF.finaleqm.S2.p_eqm.w-p_eqm2.w,SummaryF.finaleqm.S2.p_eqm.m-p_eqm2.m])))
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
fprintf('Uncertain SDF (f): the conditions rebuilt by hand (SDF with the parent node''s C) are the ones reported, this should be zero: %.3e \n',max(abs(GEhandF-GErepF),[],'all'))
fprintf('Uncertain SDF (f): and are within toleranceGEcondns at every node before T, should be one: %i \n',max(sqrt(sum(GEhandF(:,notlastF).^2,1)))<=tolC)
GEnamesF={'CapitalMarket','LabourMarket','SDF'};
for gg=1:3
    [temp,kmax]=max(abs(GEhandF(gg,:)).*notlastF);
    fprintf('Uncertain SDF (f): for information, the largest %s condition before T is %.3e, at a node in period %i with S=%i \n',GEnamesF{gg},temp,TreeF.nodedate(kmax),TreeF.nodeS(kmax))
end
fprintf('Uncertain SDF (f): for information, the S=2 final eqm conditions are %.3e %.3e %.3e (CapitalMarket, LabourMarket, SDF) \n',SummaryF.finaleqm.S2.GeneralEqmCondnValues.CapitalMarket,SummaryF.finaleqm.S2.GeneralEqmCondnValues.LabourMarket,SummaryF.finaleqm.S2.GeneralEqmCondnValues.SDF)
kbF=find(TreeF.nodenchild>1,1);
childrenkbF=TreeF.nodefirstchild(kbF)+(0:TreeF.nodenchild(kbF)-1);
fprintf('Uncertain SDF (f): for scale, m moves away from one by, should NOT be zero: %.3e \n',max(abs(PricePathTreeF.m-1)))
fprintf('Uncertain SDF (f): for scale, the children of the first branching node have m differing by, should NOT be zero: %.3e \n',abs(diff(PricePathTreeF.m(childrenkbF))))
fprintf('Uncertain SDF (f): for information, %i iterations, price distance %.3e, general eqm distance %.3e \n',SummaryF.iterations,SummaryF.PricePathDist,SummaryF.GEcondnDist)
% Anderson
transpathoptionsA=transpathoptions;
transpathoptionsA=rmfield(transpathoptionsA,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=howtoupdate;
transpathoptionsA.maxiter=200;
[PricePathTreeFA,SummaryFA]=UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a_GE,n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DFnamesF, transpathoptionsA, simoptions, vfoptions, heteroagentoptions);
fprintf('Uncertain SDF (f): Anderson converged, should be one: %i \n',SummaryFA.GEcondnDist<=tolC)
fprintf('Uncertain SDF (f): and took more than one iteration, should be one: %i \n',SummaryFA.iterations>1)
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
fprintf('Uncertain SDF (f): Anderson, SDF rebuilt by hand is the one reported, this should be zero: %.3e \n',max(abs(handFA-SummaryFA.GEcondnTree.SDF)))
diffFA=max([max(abs(PricePathTreeFA.r-PricePathTreeF.r)),max(abs(PricePathTreeFA.w-PricePathTreeF.w)),max(abs(PricePathTreeFA.m-PricePathTreeF.m))]);
fprintf('Uncertain SDF (f): for information, Anderson vs shooting prices differ by %.3e, after %i Anderson iterations \n',diffFA,SummaryFA.iterations)

%% (g) Refusals, each by its message
err=zeros(1,8);
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, 2, 5, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','m','m_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
catch ME
    err(1)=contains(ME.message,'would be applied twice, at the node and at the child');
end
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, 2, 5, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','sigma_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
catch ME
    err(2)=contains(ME.message,'is not in PricePath or ParamPath (a discount factor read at the child');
end
vfoptionsE=vfoptions; vfoptionsE.exoticpreferences='EpsteinZin';
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, 2, 5, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','m_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptionsE);
catch ME
    err(3)=contains(ME.message,'is not implemented with vfoptions.exoticpreferences');
end
ParamPathE=ParamPathU; ParamPathE.DFtplus1_UncertainTPath=ones(1,T);
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathE, T, 2, 5, pi_SC, 1, V_finalC, Policy_finalC, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, {'beta','m_tplus1'}, ReturnFn, transpathoptionsbaseline, vfoptions);
catch ME
    err(4)=contains(ME.message,'DFtplus1_UncertainTPath is a name used internally');
end
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a_GE,n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, {'beta','m','m_tplus1'}, transpathoptions, simoptions, vfoptions, heteroagentoptions);
catch ME
    err(5)=contains(ME.message,'would be applied twice, at the node and at the child');
end
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a_GE,n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, {'beta','sigma_tplus1'}, transpathoptions, simoptions, vfoptions, heteroagentoptions);
catch ME
    err(6)=contains(ME.message,'is not a price in PricePathShaper or a parameter in ParamPath');
end
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a_GE,n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, Params, DFnamesF, transpathoptions, simoptions, vfoptionsE, heteroagentoptions);
catch ME
    err(7)=contains(ME.message,'is not implemented with vfoptions.exoticpreferences');
end
ParamsE=Params; ParamsE.beta=[0.95,0.96];
try
    UncertainTransitionPath_InfHorz(PricePathShaper, ParamPathF, T, T1, T2, pi_SC, 1, AgentDist_initial, n_d,n_a_GE,n_z, d_grid,a_grid_GE,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns_Transition, ParamsE, {'beta'}, transpathoptions, simoptions, vfoptions, heteroagentoptions);
catch ME
    err(8)=contains(ME.message,'in DiscountFactorParamNames, is not a scalar');
end
fprintf('Uncertain SDF (g): ValueFnOnUncertainTransPath_InfHorz refuses m and m_tplus1 together, by name, should be one: %i \n',err(1))
fprintf('Uncertain SDF (g): ValueFnOnUncertainTransPath_InfHorz refuses sigma_tplus1 (sigma is not on the path), by name, should be one: %i \n',err(2))
fprintf('Uncertain SDF (g): ValueFnOnUncertainTransPath_InfHorz refuses m_tplus1 with exoticpreferences, by name, should be one: %i \n',err(3))
fprintf('Uncertain SDF (g): ValueFnOnUncertainTransPath_InfHorz refuses a parameter with the internal name, by name, should be one: %i \n',err(4))
fprintf('Uncertain SDF (g): UncertainTransitionPath_InfHorz refuses m and m_tplus1 together, by name, should be one: %i \n',err(5))
fprintf('Uncertain SDF (g): UncertainTransitionPath_InfHorz refuses sigma_tplus1, by name, should be one: %i \n',err(6))
fprintf('Uncertain SDF (g): UncertainTransitionPath_InfHorz refuses m_tplus1 with exoticpreferences, by name, should be one: %i \n',err(7))
fprintf('Uncertain SDF (g): UncertainTransitionPath_InfHorz refuses a vector discount factor in Parameters, by name, should be one: %i \n',err(8))

output=struct();

end
