function output=UncertainTPath_Refusals(T,rS,wS,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,transpathoptions)
% The refusals of the uncertain transition path commands, each checked by its own message (a refusal whose message
% changes, or that stops firing, shows up here as a zero rather than passing in silence).

vfoptions=struct();
ReturnFn=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
n_d=0; d_grid=[];
pi_S=[0.75,0.25;0,1];
T1=3; T2=6;
PricePathU.r=rS;
ParamPathU.w=wS;
V_final=struct(); Policy_final=struct();
for SS=1:2
    ParamsS=Params; ParamsS.r=rS(SS,T); ParamsS.w=wS(SS,T);
    [V_final.(['S',num2str(SS)]),Policy_final.(['S',num2str(SS)])]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,ParamsS,DiscountFactorParamNames,[],vfoptions);
end

err=zeros(1,14);
% The tree
try
    UncertainTransPathTree(T,T1,T2,[0.75,0.2;0,1],1);
catch ME
    err(1)=contains(ME.message,'every row of pi_S must sum to one');
end
try
    UncertainTransPathTree(T,T1,T-1,pi_S,1); % the last draw would land in period T
catch ME
    err(2)=contains(ME.message,'1<=T1<=T2<=T-2');
end
try
    UncertainTransPathTree(T,0,T2,pi_S,1);
catch ME
    err(3)=contains(ME.message,'1<=T1<=T2<=T-2');
end
try
    UncertainTransPathTree(T,3,13,[0.5,0.5;0.5,0.5],1); % 2^11=2048 histories
catch ME
    err(4)=contains(ME.message,'more than transpathoptions.maxnhistories');
end
try
    UncertainTransPathTree(T,T1,T2,pi_S,3);
catch ME
    err(5)=contains(ME.message,'S_initial must be one of');
end
try
    UncertainTransPathTree(T,T1,T2,repmat(pi_S,1,1,2),1);
catch ME
    err(6)=contains(ME.message,'N_S-by-N_S-by-(T2-T1+1)');
end
ParamPathE.w=wS(:,1:T-1);
try
    UncertainTransPathTree(T,T1,T2,pi_S,1,PricePathU,ParamPathE);
catch ME
    err(7)=contains(ME.message,'ParamPath.w must be N_S-by-T');
end
PricePathE.r=[rS;rS];
try
    UncertainTransPathTree(T,T1,T2,pi_S,1,PricePathE,ParamPathU);
catch ME
    err(8)=contains(ME.message,'PricePath.r must be 1-by-nNodes');
end
% raising maxnhistories lets the big tree through
transpathoptionsH=transpathoptions; transpathoptionsH.maxnhistories=4096;
TreeH=UncertainTransPathTree(T,3,13,[0.5,0.5;0.5,0.5],1,[],[],transpathoptionsH);
err(9)=(TreeH.nHist==2048);
% The commands
vfoptionsE=vfoptions; vfoptionsE.localsearch=1; vfoptionsE.nlocalsearch=1;
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptionsE);
catch ME
    err(10)=contains(ME.message,'localsearch is not implemented for uncertain transition paths');
end
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_final, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, repmat(pi_z,1,1,T), DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
catch ME
    err(11)=contains(ME.message,'z_grid and pi_z that vary over the path are not implemented');
end
V_finalE=rmfield(V_final,'S2');
try
    ValueFnOnUncertainTransPath_InfHorz(PricePathU, ParamPathU, T, T1, T2, pi_S, 1, V_finalE, Policy_final, Params, n_d,n_a,n_z, d_grid,a_grid,z_grid, pi_z, DiscountFactorParamNames, ReturnFn, transpathoptions, vfoptions);
catch ME
    err(12)=contains(ME.message,'need a field S2');
end
FnsToEvaluateE.lagged=@(aprime,a,z,r_tminus1) r_tminus1*a;
Tree=UncertainTransPathTree(T,T1,T2,pi_S,1);
try
    EvalFnOnUncertainTransPath_AggVars_InfHorz(FnsToEvaluateE, zeros(n_a,n_z,Tree.nNodes), ones(1,n_a,n_z,Tree.nNodes), PricePathU, ParamPathU, Params, T, T1, T2, pi_S, 1, n_d,n_a,n_z, d_grid,a_grid,z_grid, transpathoptions, struct());
catch ME
    err(13)=contains(ME.message,'_tminus1 and _tplus1 inputs are not implemented');
end
try
    AgentDistOnUncertainTransPath_InfHorz(ones(n_a,n_z)/(n_a*n_z), PricePathU, ParamPathU, ones(1,n_a,n_z,T), T, T1, T2, pi_S, 1, n_d,n_a,n_z,pi_z, Params, transpathoptions, struct());
catch ME
    err(14)=contains(ME.message,'last dimension of PolicyPathTree must be');
end

fprintf('Refusals: pi_S rows that do not sum to one, should be one: %i \n',err(1))
fprintf('Refusals: T2=T-1 (the last draw would land in the final period), should be one: %i \n',err(2))
fprintf('Refusals: T1 below one, should be one: %i \n',err(3))
fprintf('Refusals: more histories than transpathoptions.maxnhistories, should be one: %i \n',err(4))
fprintf('Refusals: S_initial out of range, should be one: %i \n',err(5))
fprintf('Refusals: a 3-D pi_S whose third dimension is not T2-T1+1, should be one: %i \n',err(6))
fprintf('Refusals: a ParamPath entry neither N_S-by-T nor 1-by-T, should be one: %i \n',err(7))
fprintf('Refusals: a PricePath entry of none of the three forms, should be one: %i \n',err(8))
fprintf('Refusals: raising transpathoptions.maxnhistories lets the 2048-history tree through, should be one: %i \n',err(9))
fprintf('Refusals: vfoptions.localsearch, should be one: %i \n',err(10))
fprintf('Refusals: a pi_z that varies over the path, should be one: %i \n',err(11))
fprintf('Refusals: a V_final without a field for an S that can be reached in period T, should be one: %i \n',err(12))
fprintf('Refusals: a _tminus1 input to FnsToEvaluate, should be one: %i \n',err(13))
fprintf('Refusals: a PolicyPathTree with T periods instead of the nodes, should be one: %i \n',err(14))

output=struct();

end
