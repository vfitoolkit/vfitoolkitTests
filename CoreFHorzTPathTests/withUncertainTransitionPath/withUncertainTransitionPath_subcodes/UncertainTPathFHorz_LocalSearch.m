function output=UncertainTPathFHorz_LocalSearch(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptionsbaseline)
% Local search on the FHorz tree, and transpathoptions.fastOLG=2 (stage F4). Written 2026-11-06, before the first run.
% The FHorz port of the local search half of UncertainTPath_LSIntEqns (CoreInfHorzUncertainTPathTests doPart(11)).
% The model of UncertainTPathFHorz_GE (the Cobb-Douglas firm OLG model, the regimes differing in TFP: S=1 0.5,
% S=2 0.52), on its absorbing tree (S=2 for good, possible from period 3 to 6, probability 1/4 each period).
% LOCAL SEARCH (vfoptions.localsearch=1): restricted sweeps around a reference per node, certified by a standard sweep.
%   (a) shooting with local search against shooting without: both converged, the same prices to 10*toleranceGEprices;
%       restricted sweeps were taken and certified (the counts); and the summary is of a STANDARD sweep, so its V and
%       Policy are those of ValueFnOnUncertainTransPath_FHorz with no local search at the returned prices, exactly
%   (b) Anderson with local search: the same
%   (c) transpathoptions.fastOLG=2 with local search, shooting: the same, against (a)
%   (d) ValueFnOnUncertainTransPath_FHorz on its own, given prices: localsearch=1 is ignored with a warning (the same V
%       as 0, exactly); localsearch=2 forces a stay-put window, on a single history exactly
%       ValueFnOnTransPath_Case1_FHorz with localsearch=2, and different from the unrestricted V (the restriction has teeth)
%   (e) the INTERNAL reference option of ValueFnOnTransPath_Case1_FHorz, given prices, for fastOLG 0, 1 and 2: a seed
%       sweep builds a reference; a restricted sweep at the same prices centred on it is the unrestricted answer
%       exactly, moves no reference, and hits no window edge. The DISCRIMINATOR for fastOLG=2 (seeds solved slowOLG,
%       put in the fastOLG layout of the restricted sweeps): the slowOLG seed of fastOLG=0, not converted, moves it.
%   (f) refusals, each by its message
% Shooting stops at a tolerance, so two solves of the same model are held to 10 times toleranceGEprices, as in
% UncertainTPathFHorz_GE.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
N_a=prod(n_a); N_z=prod(n_z);
NV=N_a*N_z*N_j;
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist=zeros([n_a,n_z],'gpuArray'); jequaloneDist(1,ceil(n_z/2))=1;
heteroagentoptions=struct();
heteroagentoptions.verbose=0;

Params.firmalpha=0.36;
Params.firmdelta=0.05;
Params.firmA=0.5;
FnsToEvaluate.K=@(aprime,a,z) a;
FnsToEvaluate.N=@(aprime,a,z,kappa_j) kappa_j*z;
GeneralEqmEqns.CapitalMarket=@(r,K,N,firmalpha,firmdelta,firmA) r-(firmalpha*firmA*(K^(firmalpha-1))*(N^(1-firmalpha))-firmdelta);
GeneralEqmEqns.LabourMarket=@(w,K,N,firmalpha,firmA) w-((1-firmalpha)*firmA*(K^firmalpha)*(N^(-firmalpha)));

%% The initial stationary eqm (S=1, firmA=0.5)
ParamsGE=Params;
[p_eqm,~,GEcondns]=HeteroAgentStationaryEqm_Case1_FHorz(jequaloneDist,AgeWeightParamNames,n_d,n_a,n_z,N_j,0,pi_z,d_grid,a_grid,z_grid,ReturnFn,FnsToEvaluate,GeneralEqmEqns,ParamsGE,DiscountFactorParamNames,[],[],[],{'r','w'},heteroagentoptions,simoptions,vfoptions);
fprintf('Uncertain FHorz LS: initial stationary eqm r=%2.8f, w=%2.8f, conditions %.3e %.3e \n',p_eqm.r,p_eqm.w,GEcondns.CapitalMarket,GEcondns.LabourMarket)
ParamsGE.r=p_eqm.r;
ParamsGE.w=p_eqm.w;
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsGE,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,ParamsGE,simoptions);

howtoupdate={'CapitalMarket','r',0,0.1; 'LabourMarket','w',0,0.1};
transpathoptions=transpathoptionsbaseline;
transpathoptions.GEnewprice=3;
transpathoptions.GEnewprice3.howtoupdate=howtoupdate;
transpathoptions.maxiter=1000;
transpathoptions.verbose=0;
transpathoptions.toleranceGEprices=10^(-4);
transpathoptions.toleranceGEcondns=2*10^(-3); % the two-price stationary solves stop well inside this (UncertainTPathFHorz_GE)
tolP=transpathoptions.toleranceGEprices;
tolC=transpathoptions.toleranceGEcondns;
fprintf('Uncertain FHorz LS: toleranceGEprices %.3e, toleranceGEcondns %.3e \n',tolP,tolC)
PricePathShaper.r=min((1:1:T)/20,1);
PricePathShaper.w=min((1:1:T)/20,1);

ParamPathC.firmA=[0.5*ones(1,T); 0.52*ones(1,T)];
pi_SC=[0.75,0.25;0,1]; T1=2; T2=5;

vfoptionsLS=vfoptions;
vfoptionsLS.localsearch=1;
vfoptionsLS.nlocalsearch=3;

%% (a) Shooting, without and with local search
[PricePathTreeA0,SummaryA0]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptions, heteroagentoptions);
[PricePathTreeA,SummaryA]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions, simoptions, vfoptionsLS, heteroagentoptions);
TreeC=SummaryA.Tree;
nNodesC=TreeC.nNodes;
fprintf('Uncertain FHorz LS (a): %i histories, both converged before maxiter, should be one: %i \n',TreeC.nHist,SummaryA0.iterations<transpathoptions.maxiter && SummaryA.iterations<transpathoptions.maxiter)
diffA=max([max(abs(PricePathTreeA.r-PricePathTreeA0.r)),max(abs(PricePathTreeA.w-PricePathTreeA0.w))]);
fprintf('Uncertain FHorz LS (a): with and without local search, the same prices within 10*toleranceGEprices, should be one: %i \n',diffA<=10*tolP)
fprintf('Uncertain FHorz LS (a): for information, that difference is %.3e; %i iterations without, %i with \n',diffA,SummaryA0.iterations,SummaryA.iterations)
fprintf('Uncertain FHorz LS (a): restricted sweeps were taken, should be one: %i \n',SummaryA.localsearch.nrestricted>0)
fprintf('Uncertain FHorz LS (a): and certified by at least one standard sweep, should be one: %i \n',SummaryA.localsearch.nverify>=1)
fprintf('Uncertain FHorz LS (a): the last sweep was a standard one, should be one: %i \n',SummaryA.localsearch.lastsweeprestricted==0)
V_finalC=struct('S1',SummaryA.finaleqm.S1.V_final,'S2',SummaryA.finaleqm.S2.V_final);
Policy_finalC=struct('S1',SummaryA.finaleqm.S1.Policy_final,'S2',SummaryA.finaleqm.S2.Policy_final);
[VhandA,PhandA]=ValueFnOnUncertainTransPath_FHorz(PricePathTreeA, ParamPathC, T, T1, T2, pi_SC, 1, V_finalC, Policy_finalC, ParamsGE, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VhandA=gather(reshape(VhandA,[NV,nNodesC])); VA=gather(reshape(SummaryA.VPathTree,[NV,nNodesC]));
dd=abs(VhandA-VA); dd(VhandA==VA)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz LS (a): the summary V is the unrestricted V at the returned prices, this should be zero: %.3e \n',max(dd,[],'all'))
fprintf('Uncertain FHorz LS (a): and the summary Policy, this should be zero: %.3e \n',max(abs(gather(PhandA(:))-gather(SummaryA.PolicyPathTree(:)))))

%% (b) Anderson with local search
transpathoptionsA=rmfield(transpathoptions,'GEnewprice3');
transpathoptionsA.GEnewprice=2;
transpathoptionsA.GEnewprice2.howtoupdate=howtoupdate;
transpathoptionsA.maxiter=200;
[PricePathTreeB,SummaryB]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptionsA, simoptions, vfoptionsLS, heteroagentoptions);
fprintf('Uncertain FHorz LS (b): Anderson with local search converged, should be one: %i \n',SummaryB.GEcondnDist<=tolC)
fprintf('Uncertain FHorz LS (b): and was certified by at least one unrestricted sweep, should be one: %i \n',SummaryB.localsearch.nverify>=1)
[VhandB,PhandB]=ValueFnOnUncertainTransPath_FHorz(PricePathTreeB, ParamPathC, T, T1, T2, pi_SC, 1, V_finalC, Policy_finalC, ParamsGE, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VhandB=gather(reshape(VhandB,[NV,nNodesC])); VB=gather(reshape(SummaryB.VPathTree,[NV,nNodesC]));
dd=abs(VhandB-VB); dd(VhandB==VB)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz LS (b): the summary V is the unrestricted V at the returned prices, this should be zero: %.3e \n',max(dd,[],'all'))
fprintf('Uncertain FHorz LS (b): and the summary Policy, this should be zero: %.3e \n',max(abs(gather(PhandB(:))-gather(SummaryB.PolicyPathTree(:)))))
fprintf('Uncertain FHorz LS (b): for information, Anderson with local search vs shooting without differ by %.3e, after %i Anderson iterations \n',max([max(abs(PricePathTreeB.r-PricePathTreeA0.r)),max(abs(PricePathTreeB.w-PricePathTreeA0.w))]),SummaryB.iterations)

%% (c) transpathoptions.fastOLG=2 with local search, shooting
transpathoptions2=transpathoptions;
transpathoptions2.fastOLG=2;
[PricePathTreeC2,SummaryC2]=UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, transpathoptions2, simoptions, vfoptionsLS, heteroagentoptions);
fprintf('Uncertain FHorz LS (c): fastOLG=2 converged before maxiter, with restricted sweeps certified, should be one: %i \n',SummaryC2.iterations<transpathoptions.maxiter && SummaryC2.localsearch.nrestricted>0 && SummaryC2.localsearch.nverify>=1)
diffC2=max([max(abs(PricePathTreeC2.r-PricePathTreeA.r)),max(abs(PricePathTreeC2.w-PricePathTreeA.w))]);
fprintf('Uncertain FHorz LS (c): fastOLG=2 vs fastOLG=1, both with local search, the same prices within 10*toleranceGEprices, should be one: %i \n',diffC2<=10*tolP)
fprintf('Uncertain FHorz LS (c): for information, that difference is %.3e, after %i iterations \n',diffC2,SummaryC2.iterations)
VhandC2=ValueFnOnUncertainTransPath_FHorz(PricePathTreeC2, ParamPathC, T, T1, T2, pi_SC, 1, V_finalC, Policy_finalC, ParamsGE, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VhandC2=gather(reshape(VhandC2,[NV,nNodesC])); VC2=gather(reshape(SummaryC2.VPathTree,[NV,nNodesC]));
dd=abs(VhandC2-VC2); dd(VhandC2==VC2)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz LS (c): the summary V is the unrestricted V (fastOLG=1) at the returned prices, this should be zero: %.3e \n',max(dd,[],'all'))

%% Given prices, for (d) and (e): the parent bank's regimes
V_finalU=struct(); Policy_finalU=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_finalU.(['S',num2str(SS)]),Policy_finalU.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end
PricePathU.r=rS;
ParamPathU.w=wS;

%% (d) ValueFnOnUncertainTransPath_FHorz on its own
VU0=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_SC, 1, V_finalU, Policy_finalU, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
vfo1=vfoptions; vfo1.localsearch=1; vfo1.nlocalsearch=1;
lastwarn('');
warning('off','all'); % lastwarn still records it; keeps the diary readable
VU1=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_SC, 1, V_finalU, Policy_finalU, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfo1);
warning('on','all');
warnmsg=lastwarn;
fprintf('Uncertain FHorz LS (d): localsearch=1 is ignored with a warning, should be one: %i \n',contains(warnmsg,'vfoptions.localsearch=1 is IGNORED here'))
fprintf('Uncertain FHorz LS (d): and the V is the unrestricted one, this should be zero: %.3e \n',max(abs(gather(VU1(:))-gather(VU0(:)))))
% localsearch=2, a stay-put window of one point either side, on a single history
pi_SA=[0,1;0,1];
vfo2=vfoptions; vfo2.localsearch=2; vfo2.nlocalsearch=1;
VUA2=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, pi_SA, 1, V_finalU, Policy_finalU, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfo2);
VUA0=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, 1, 1, pi_SA, 1, V_finalU, Policy_finalU, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoptions);
clear PricePathD ParamPathD
PricePathD.r=[rS(1,1),rS(2,2:T)];
ParamPathD.w=[wS(1,1),wS(2,2:T)];
VD2=ValueFnOnTransPath_Case1_FHorz(PricePathD, ParamPathD, T, V_finalU.S2, Policy_finalU.S2, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfo2);
VUA2=gather(reshape(VUA2,[NV,T])); VUA0=gather(reshape(VUA0,[NV,T])); VD2=gather(reshape(VD2,[NV,T]));
dd=abs(VUA2-VD2); dd(VUA2==VD2)=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz LS (d): localsearch=2 on a single history is ValueFnOnTransPath_Case1_FHorz with localsearch=2, this should be zero: %.3e \n',max(dd,[],'all'))
dd=abs(VUA2-VUA0); dd(VUA2==VUA0)=0; dd(isnan(dd))=0;
fprintf('Uncertain FHorz LS (d): for scale, the stay-put window differs from the unrestricted V by, should NOT be zero: %.3e \n',max(dd,[],'all'))

%% (e) The internal reference option of ValueFnOnTransPath_Case1_FHorz
PricePathE.r=rS(1,:);
ParamPathE.w=wS(1,:);
fastOLGcases=[0,1,2];
movesE=zeros(1,3); errVE=zeros(1,3); errPE=zeros(1,3); edgesE=zeros(1,3); seedrefs=cell(1,3);
for cc=1:3
    tpE=transpathoptionsbaseline; tpE.fastOLG=fastOLGcases(cc);
    vfoS=vfoptions; vfoS.localsearch=0; vfoS.localsearchseed=1;
    [VE0,PE0,refE]=ValueFnOnTransPath_Case1_FHorz(PricePathE, ParamPathE, T, V_finalU.S1, Policy_finalU.S1, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, tpE, vfoS);
    seedrefs{cc}=refE;
    vfoR=vfoptions; vfoR.localsearch=2; vfoR.nlocalsearch=3; vfoR.aprimeReferencePath=refE;
    [VE1,PE1,refE1,lsE]=ValueFnOnTransPath_Case1_FHorz(PricePathE, ParamPathE, T, V_finalU.S1, Policy_finalU.S1, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, tpE, vfoR);
    VE0=gather(VE0(:)); VE1=gather(VE1(:));
    dd=abs(VE1-VE0); dd(VE1==VE0)=0; dd(isnan(dd))=Inf;
    errVE(cc)=max(dd);
    errPE(cc)=max(abs(gather(PE1(:))-gather(PE0(:))));
    movesE(cc)=max(abs(double(gather(refE1(:)))-double(gather(refE(:)))));
    edgesE(cc)=lsE.nedge;
end
for cc=1:3
    fprintf('Uncertain FHorz LS (e): fastOLG=%i, a restricted sweep centred on the seed at the same prices is the unrestricted V, this should be zero: %.3e \n',fastOLGcases(cc),errVE(cc))
    fprintf('Uncertain FHorz LS (e): fastOLG=%i, and Policy, this should be zero: %.3e \n',fastOLGcases(cc),errPE(cc))
    fprintf('Uncertain FHorz LS (e): fastOLG=%i, it moves no reference and hits no window edge, should be one: %i \n',fastOLGcases(cc),movesE(cc)==0 && edgesE(cc)==0)
end
% The discriminator: under fastOLG=2 the restricted sweep reads the fastOLG layout, so the unconverted slowOLG seed of
% fastOLG=0 is the wrong reference there
tpE=transpathoptionsbaseline; tpE.fastOLG=2;
vfoR=vfoptions; vfoR.localsearch=2; vfoR.nlocalsearch=3; vfoR.aprimeReferencePath=seedrefs{1};
[~,~,refEbad]=ValueFnOnTransPath_Case1_FHorz(PricePathE, ParamPathE, T, V_finalU.S1, Policy_finalU.S1, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, tpE, vfoR);
fprintf('Uncertain FHorz LS (e): for scale, under fastOLG=2 the slowOLG seed, not put in the fastOLG layout, moves the reference by, should NOT be zero: %i \n',max(abs(double(gather(refEbad(:)))-double(gather(seedrefs{1}(:))))))
fprintf('Uncertain FHorz LS (e): and the fastOLG=2 seed is the fastOLG=1 seed, the slowOLG one put in the fastOLG layout, should be one: %i \n',isequal(gather(seedrefs{3}),gather(seedrefs{2})))

%% (f) Refusals, each by its message
GEcall=@(tpo,vfo) UncertainTransitionPath_FHorz(PricePathShaper, ParamPathC, T, T1, T2, pi_SC, 1, AgentDist_initial, jequaloneDist, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, ReturnFn, FnsToEvaluate, GeneralEqmEqns, GeneralEqmEqns, ParamsGE, DiscountFactorParamNames, AgeWeightParamNames, tpo, simoptions, vfo, heteroagentoptions);
err=zeros(1,4);
vfoE=vfoptions; vfoE.localsearch=2; vfoE.nlocalsearch=1;
try; GEcall(transpathoptions,vfoE); catch ME; err(1)=contains(ME.message,'the value 2 is only meaningful to ValueFnOnTransPath_Case1_FHorz and ValueFnOnUncertainTransPath_FHorz'); end
vfoE=vfoptions; vfoE.localsearch=1; vfoE.nlocalsearch=N_a;
try; GEcall(transpathoptions,vfoE); catch ME; err(2)=contains(ME.message,'vfoptions.nlocalsearch is too large for the grid'); end
vfoE=vfoptions; vfoE.localsearch=1; vfoE.nlocalsearch=1; vfoE.nlocalsearchup=0;
try; GEcall(transpathoptions,vfoE); catch ME; err(3)=contains(ME.message,'nlocalsearchup and vfoptions.nlocalsearchdown must both be at least 1'); end
vfoE=vfoptions; vfoE.aprimeReferencePath=seedrefs{2};
try; ValueFnOnTransPath_Case1_FHorz(PricePathE, ParamPathE, T, V_finalU.S1, Policy_finalU.S1, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptionsbaseline, vfoE); catch ME; err(4)=contains(ME.message,'is only used by a restricted sweep, vfoptions.localsearch=2'); end
fprintf('Uncertain FHorz LS (f): UncertainTransitionPath_FHorz refuses localsearch=2, by name, should be one: %i \n',err(1))
fprintf('Uncertain FHorz LS (f): an nlocalsearch too large for the grid is refused by name, should be one: %i \n',err(2))
fprintf('Uncertain FHorz LS (f): nlocalsearchup=0 is refused by name, should be one: %i \n',err(3))
fprintf('Uncertain FHorz LS (f): ValueFnOnTransPath_Case1_FHorz refuses the internal reference without a restricted sweep, by name, should be one: %i \n',err(4))

output=struct();

end
