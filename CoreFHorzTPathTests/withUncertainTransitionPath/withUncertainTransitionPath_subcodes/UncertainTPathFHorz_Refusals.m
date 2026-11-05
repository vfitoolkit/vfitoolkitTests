function output=UncertainTPathFHorz_Refusals(T,rS,wS,n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,transpathoptions,vfoptionsbaseline)
% The refusals of the FHorz uncertain transition path commands (stage F1), each checked by its own message, as a set.
% And the warning on a square entry (N_j=T) of a deterministic FHorz path: given by the deterministic commands, and NOT
% given when the tree itself makes a segment N_j periods long (it builds that entry square on purpose), with the warning
% state restored afterwards. Written 2026-11-05, before the first run.

vfoptions=struct();
simoptions=struct();
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist=zeros([n_a,n_z]); jequaloneDist(1,ceil(n_z/2))=1;
pi_S=[0.75,0.25;0,1]; T1=2; T2=4;
PricePathU.r=rS;
ParamPathU.w=wS;
V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end
Params1=Params; Params1.r=rS(1,1); Params1.w=wS(1,1);
[~,Policy_init]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params1,DiscountFactorParamNames,[],vfoptions);
AgentDist_initial=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_init,n_d,n_a,n_z,N_j,pi_z,Params1,simoptions);
[~,PU]=ValueFnOnUncertainTransPath_FHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);

err=zeros(1,10);
msg=cell(1,10);
VFcall=@(vfo,piz,DFn,PP) ValueFnOnUncertainTransPath_FHorz(PricePathU, PP, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, piz, DFn, ReturnFn, transpathoptions, vfo);
vfoE=vfoptions; vfoE.n_semiz=2;
try; VFcall(vfoE,pi_z,DiscountFactorParamNames,ParamPathU); catch ME; err(1)=contains(ME.message,'semi-exogenous states (vfoptions.n_semiz) are not yet implemented'); end
msg{1}='ValueFnOnUncertainTransPath_FHorz refuses a semi-exogenous state';
vfoE=vfoptions; vfoE.experienceasset=1;
try; VFcall(vfoE,pi_z,DiscountFactorParamNames,ParamPathU); catch ME; err(2)=contains(ME.message,'experience assets are not yet implemented'); end
msg{2}='ValueFnOnUncertainTransPath_FHorz refuses an experience asset';
vfoE=vfoptions; vfoE.exoticpreferences='QuasiHyperbolic';
try; VFcall(vfoE,pi_z,DiscountFactorParamNames,ParamPathU); catch ME; err(3)=contains(ME.message,'vfoptions.exoticpreferences is not yet implemented'); end
msg{3}='ValueFnOnUncertainTransPath_FHorz refuses exoticpreferences';
try; VFcall(vfoptions,repmat(pi_z,1,1,N_j,T),DiscountFactorParamNames,ParamPathU); catch ME; err(4)=contains(ME.message,'z_grid and pi_z that vary over the path are not implemented'); end
msg{4}='ValueFnOnUncertainTransPath_FHorz refuses a pi_z that varies over the path (4-D)';
vfoE=vfoptions; vfoE.n_e=vfoptionsbaseline.n_e; vfoE.e_grid=vfoptionsbaseline.e_grid; vfoE.pi_e=repmat(vfoptionsbaseline.pi_e,1,N_j,T);
try; VFcall(vfoE,pi_z,DiscountFactorParamNames,ParamPathU); catch ME; err(5)=contains(ME.message,'an e (vfoptions.pi_e) that varies over the path is not implemented'); end
msg{5}='ValueFnOnUncertainTransPath_FHorz refuses a pi_e that varies over the path';
vfoE=vfoptions; vfoE.localsearch=3;
try; VFcall(vfoE,pi_z,DiscountFactorParamNames,ParamPathU); catch ME; err(6)=contains(ME.message,'vfoptions.localsearch must be 0, 1 (ignored here, with a warning) or 2 (forced)'); end
msg{6}='ValueFnOnUncertainTransPath_FHorz refuses localsearch=3';
ParamPathE=ParamPathU; ParamPathE.bvec=repmat(0.99*ones(N_j,1),1,T);
try; VFcall(vfoptions,pi_z,{'beta','bvec_tplus1'},ParamPathE); catch ME; err(7)=contains(ME.message,'is not a scalar; a discount factor read at the child must be a scalar'); end
msg{7}='ValueFnOnUncertainTransPath_FHorz refuses a discount factor read at the child that is a vector';
jeqE.S1=jequaloneDist;
try; AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jeqE, PricePathU, ParamPathU, PU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions); catch ME; err(8)=contains(ME.message,'jequaloneDist is given per S, but has no field S2'); end
msg{8}='AgentDistOnUncertainTransPath_FHorz refuses a per-S jequaloneDist missing an S';
try; AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, zeros(n_a,n_z,3), PricePathU, ParamPathU, PU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoptions); catch ME; err(9)=contains(ME.message,'jequaloneDist must be [n_a,n_z] (the same every period), [n_a,n_z,T]'); end
msg{9}='AgentDistOnUncertainTransPath_FHorz refuses a jequaloneDist of the wrong size';
simoE=simoptions; simoE.n_semiz=2;
try; AgentDistOnUncertainTransPath_FHorz(AgentDist_initial, jequaloneDist, PricePathU, ParamPathU, PU, AgeWeightParamNames, T, T1, T2, pi_S, 1, n_d,n_a,n_z,N_j,pi_z, Params, transpathoptions, simoE); catch ME; err(10)=contains(ME.message,'semi-exogenous states (simoptions.n_semiz) are not yet implemented'); end
msg{10}='AgentDistOnUncertainTransPath_FHorz refuses a semi-exogenous state';
for ee=1:10
    fprintf('Uncertain FHorz refusals: %s, by name, should be one: %i \n',msg{ee},err(ee))
end

%% The square-entry warning
% A deterministic FHorz path with T=N_j and an N_j-by-T parameter: warned, and read as N_j-by-T
TS=N_j;
clear PricePathS ParamPathS
PricePathS.r=rS(1,1:TS);
ParamPathS.w=wS(1,1:TS);
ParamPathS.kappa_j=repmat(Params.kappa_j(:),1,TS).*(1+0.01*(1:TS)); % N_j-by-T, square, varies over t
% what is printed is checked (evalc), not lastwarn: MATLAB updates lastwarn even for a warning that is turned off
txtS=evalc('[VS,~]=ValueFnOnTransPath_Case1_FHorz(PricePathS, ParamPathS, TS, V_final.S1, Policy_final.S1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);');
fprintf('Uncertain FHorz square entry: a deterministic path with N_j=T and an N_j-by-T parameter warns how it is read, should be one: %i \n',contains(txtS,'so it is square; it is being read as'))
% read as N_j-by-T: the same as the transpose given as T-by-N_j when T~=N_j cannot be compared, so compare with the
% first period computed directly with that period's kappa_j (period TS-1 is one step from V_final)
ParamsC=Params; ParamsC.r=rS(1,TS-1); ParamsC.w=wS(1,TS-1); ParamsC.kappa_j=ParamPathS.kappa_j(:,TS-1)';
clear PricePathC ParamPathC
PricePathC.r=[ParamsC.r,ParamsC.r];
ParamPathC.w=[ParamsC.w,ParamsC.w];
ParamPathC.kappa_j=[ParamsC.kappa_j',ParamsC.kappa_j'];
VC=ValueFnOnTransPath_Case1_FHorz(PricePathC, ParamPathC, 2, V_final.S1, Policy_final.S1, Params, n_d, n_a, n_z, N_j, d_grid, a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
VS=gather(reshape(VS,[],TS)); VC=gather(reshape(VC,[],2));
dd=abs(VS(:,TS-1)-VC(:,1)); dd(VS(:,TS-1)==VC(:,1))=0; dd(isnan(dd))=Inf;
fprintf('Uncertain FHorz square entry: and it is read as N_j-by-T (period T-1 is the one step with that column of kappa_j), this should be zero: %.3e \n',max(dd))
% The tree making a segment N_j periods long itself: T=20, a branching node in period 5, so the children's segments
% are periods 6,...,20, 15=N_j periods; the user's kappa_j is N_j-by-20, not square
T20=20;
clear PricePath20 ParamPath20
PricePath20.r=rS(:,1:T20);
ParamPath20.w=wS(:,1:T20);
ParamPath20.kappa_j=repmat(Params.kappa_j(:),1,T20);
pi_S20=[0.5,0.5;0,1];
txt20=evalc('VFcall20=ValueFnOnUncertainTransPath_FHorz(PricePath20, ParamPath20, T20, 5, 5, pi_S20, 1, V_final, Policy_final, Params, n_d,n_a,n_z,N_j, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);');
Tree20=UncertainTransPathTree(T20,5,5,pi_S20,1);
seglengths=cellfun(@length,Tree20.segnodes);
fprintf('Uncertain FHorz square entry: for scale, the tree has a segment N_j periods long, should be one: %i \n',any(seglengths==N_j))
fprintf('Uncertain FHorz square entry: the tree does not warn about the square entry it builds itself, should be one: %i \n',~contains(txt20,'so it is square; it is being read as'))
wstate=warning('query','VFIToolkit:SquarePathEntry');
fprintf('Uncertain FHorz square entry: and the warning is on again afterwards, should be one: %i \n',strcmp(wstate.state,'on'))

output=struct();

end
