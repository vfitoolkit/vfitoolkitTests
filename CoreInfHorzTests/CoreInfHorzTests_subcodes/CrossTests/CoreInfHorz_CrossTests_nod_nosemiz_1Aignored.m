function output=CoreInfHorz_CrossTests_nod_nosemiz_1Aignored(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% CROSS-TEST: two endogenous states where the FIRST one does nothing, vs just one endogenous state.
% Without d.
%
% This is the mirror image of CoreInfHorz_CrossTests_nod_nosemiz_2Aignored: here it is a1 that the
% return function ignores, and a2 that is the actual asset. So it is again economically the same
% problem as the one endogenous state model (which is the fig 3 model), and V and the a2prime
% policy must be identical (with a1 just replicating them).
%
% Note: no grid interpolation layer here. GI interpolates the FIRST endogenous state, which in this
% test is the ignored dummy, so a GI version would not be testing anything meaningful.
%
% Note: a1prime is NOT compared. Since a1 is ignored, every a1prime is equally optimal, so which
% one gets chosen is just an arbitrary tie-break and carries no information.
% Note: policies are compared as VALUES (via PolicyInd2Val) rather than as indexes, as the
% index encodings differ between the one and two endogenous state cases.

DF=DiscountFactorParamNames;
N_z=prod(n_z);

n_d=0; d_grid=[];

% One endogenous state (this is the fig 3 model)
ReturnFn_1A=@(aprime,a,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,sigma);
FnsToEvaluate_1A.assets=@(aprime,a,z) a;

% Two endogenous states, the first of which is ignored by the return fn (a2 is the actual asset)
n_a1=3;
a1_grid=[1;2;3]; % values are irrelevant, a1 is ignored
n_a_2A=[n_a1,n_a];
a_grid_2A=[a1_grid; a_grid];
ReturnFn_2A=@(a1prime,a2prime,a1,a2,z,r,w,sigma) ReturnFn_nod_z_noe_nosemiz_1Aignored(a1prime,a2prime,a1,a2,z,r,w,sigma);
FnsToEvaluate_2A.assets=@(a1prime,a2prime,a1,a2,z) a2;

fprintf('\n================ CROSS-TEST: 2A with ignored FIRST asset vs 1A (nod) ================\n');

vfoptions=struct(); simoptions=struct();

[V1A,Policy1A]=ValueFnIter_InfHorz(n_d,n_a,   n_z,d_grid,a_grid,   z_grid,pi_z,ReturnFn_1A,Params,DF,[],vfoptions);
[V2A,Policy2A]=ValueFnIter_InfHorz(n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DF,[],vfoptions);

% V: the 2A value fn is just the 1A value fn, repeated across the (ignored) a1 dimension
V1Arep=repmat(reshape(V1A,[1,n_a,N_z]),n_a1,1,1);
fprintf('no GI : V, 2A(ignored a1) vs 1A, should be zero:            %.3e \n',max(abs(reshape(V2A,[n_a1,n_a,N_z])-V1Arep),[],'all'));

% Policy: compare the a2prime VALUES (a1prime is an arbitrary tie-break, so is not compared)
PolicyVals1A=PolicyInd2Val_InfHorz(Policy1A,n_d,n_a,   n_z,d_grid,a_grid,   vfoptions); % [1,n_a,n_z]
PolicyVals2A=PolicyInd2Val_InfHorz(Policy2A,n_d,n_a_2A,n_z,d_grid,a_grid_2A,vfoptions); % [2,n_a1,n_a2,n_z]
a2prime1A=repmat(reshape(PolicyVals1A(1,:,:),[1,n_a,N_z]),n_a1,1,1);
a2prime2A=reshape(PolicyVals2A(2,:,:,:),[n_a1,n_a,N_z]);
fprintf('no GI : a2prime policy values, 2A vs 1A, should be zero:    %.3e \n',max(abs(a2prime2A-a2prime1A),[],'all'));

% ValueFnFromPolicy on the 2A model should reproduce the 2A value fn
V2AfromPolicy=ValueFnFromPolicy_InfHorz(Policy2A,n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DF,vfoptions);
fprintf('no GI : ValueFnFromPolicy on 2A, should be zero:            %.3e \n',max(abs(V2AfromPolicy(:)-V2A(:))));

% StationaryDist: the a2-marginal of the 2A dist equals the 1A dist
% (a2 evolves independently of a1, so this holds whichever a1prime the tie-break picked)
StationaryDist1A=StationaryDist_InfHorz(Policy1A,n_d,n_a,   n_z,pi_z,simoptions,Params,[]);
StationaryDist2A=StationaryDist_InfHorz(Policy2A,n_d,n_a_2A,n_z,pi_z,simoptions,Params,[]);
StationaryDist2A_a2=reshape(sum(StationaryDist2A,1),[n_a,N_z]);
fprintf('no GI : StationaryDist a2-marginal, 2A vs 1A, ~zero:        %.3e \n',max(abs(StationaryDist2A_a2-reshape(StationaryDist1A,[n_a,N_z])),[],'all'));

% AllStats on assets should agree
AllStats1A=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1A,Policy1A,FnsToEvaluate_1A,Params,[],n_d,n_a,   n_z,d_grid,a_grid,   z_grid,simoptions);
AllStats2A=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist2A,Policy2A,FnsToEvaluate_2A,Params,[],n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,simoptions);
fprintf('no GI : AllStats assets Mean, 2A vs 1A, should be zero:     %.3e \n',abs(AllStats2A.assets.Mean-AllStats1A.assets.Mean));
fprintf('no GI : AllStats assets Gini, 2A vs 1A, should be zero:     %.3e \n',abs(AllStats2A.assets.Gini-AllStats1A.assets.Gini));

output=struct();

end
