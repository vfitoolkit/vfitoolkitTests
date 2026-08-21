function output=CoreInfHorz_CrossTests_d_nosemiz_2Aignored(n_d,n_a,n_a_big,n_z,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline)
% CROSS-TEST: two endogenous states where the SECOND one does nothing, vs just one endogenous state.
% With d.
%
% The 2A model adds a second endogenous state a2 that never appears in the return function, so
% it is economically the same problem as the one endogenous state model (which is the fig 4 model).
% Hence V, the d policy, and the a1prime policy must be identical (with a2 just replicating them).
% This checks the two-endogenous-state code paths against the (already tested) one endogenous
% state code paths. Done both without and with the grid interpolation layer, so it also checks
% the GI2A code path against the GI (one asset) code path.
%
% Note: a2prime is NOT compared. Since a2 is ignored, every a2prime is equally optimal, so which
% one gets chosen is just an arbitrary tie-break and carries no information.
% Note: policies are compared as VALUES (via PolicyInd2Val) rather than as indexes, as the
% index encodings differ between the one and two endogenous state cases (and under GI).

DF=DiscountFactorParamNames;
N_z=prod(n_z);

% One endogenous state (this is the fig 4 model)
ReturnFn_1A=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);
FnsToEvaluate_1A.assets=@(d,aprime,a,z) a;

% Two endogenous states, the second of which is ignored by the return fn
n_a2=3;
a2_grid=[1;2;3]; % values are irrelevant, a2 is ignored
n_a_2A=[n_a,n_a2];
a_grid_2A=[a_grid; a2_grid];
ReturnFn_2A=@(d,aprime,a2prime,a,a2,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz_2Aignored(d,aprime,a2prime,a,a2,z,r,w,sigma,eta,varphi);
FnsToEvaluate_2A.assets=@(d,aprime,a2prime,a,a2,z) a;

fprintf('\n================ CROSS-TEST: 2A with ignored second asset vs 1A (with d) ================\n');

for gi=0:1
    vfoptions=struct(); simoptions=struct();
    vfoptions.verbose_advice=0;
    if gi==0
        gistr='no GI ';
    else
        gistr='with GI';
        vfoptions.gridinterplayer=1; vfoptions.ngridinterp=5;
        simoptions.gridinterplayer=vfoptions.gridinterplayer; simoptions.ngridinterp=vfoptions.ngridinterp;
    end

    [V1A,Policy1A]=ValueFnIter_InfHorz(n_d,n_a,   n_z,d_grid,a_grid,   z_grid,pi_z,ReturnFn_1A,Params,DF,[],vfoptions);
    [V2A,Policy2A]=ValueFnIter_InfHorz(n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DF,[],vfoptions);

    % V: the 2A value fn is just the 1A value fn, repeated across the a2 dimension
    V1Arep=repmat(reshape(V1A,[n_a,1,N_z]),1,n_a2,1);
    fprintf('%s: V, 2A(ignored a2) vs 1A, should be zero:            %2.8f \n',gistr,max(abs(reshape(V2A,[n_a,n_a2,N_z])-V1Arep),[],'all'));

    % Policy: compare the d and aprime VALUES (a2prime is an arbitrary tie-break, so is not compared)
    PolicyVals1A=PolicyInd2Val_InfHorz(Policy1A,n_d,n_a,   n_z,d_grid,a_grid,   vfoptions); % [2,n_a,n_z]      (d,aprime)
    PolicyVals2A=PolicyInd2Val_InfHorz(Policy2A,n_d,n_a_2A,n_z,d_grid,a_grid_2A,vfoptions); % [3,n_a1,n_a2,n_z] (d,a1prime,a2prime)
    d1A=repmat(reshape(PolicyVals1A(1,:,:),[n_a,1,N_z]),1,n_a2,1);
    d2A=reshape(PolicyVals2A(1,:,:,:),[n_a,n_a2,N_z]);
    fprintf('%s: d policy values, 2A vs 1A, should be zero:          %2.8f \n',gistr,max(abs(d2A-d1A),[],'all'));
    aprime1A=repmat(reshape(PolicyVals1A(2,:,:),[n_a,1,N_z]),1,n_a2,1);
    aprime2A=reshape(PolicyVals2A(2,:,:,:),[n_a,n_a2,N_z]);
    fprintf('%s: aprime policy values, 2A vs 1A, should be zero:     %2.8f \n',gistr,max(abs(aprime2A-aprime1A),[],'all'));

    % ValueFnFromPolicy on the 2A model should reproduce the 2A value fn
    V2AfromPolicy=ValueFnFromPolicy_InfHorz(Policy2A,n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,pi_z,ReturnFn_2A,Params,DF,vfoptions);
    fprintf('%s: ValueFnFromPolicy on 2A, should be zero:            %2.8f \n',gistr,max(abs(V2AfromPolicy(:)-V2A(:))));

    % StationaryDist: the a1-marginal of the 2A dist equals the 1A dist
    % (a1 evolves independently of a2, so this holds whichever a2prime the tie-break picked)
    StationaryDist1A=StationaryDist_InfHorz(Policy1A,n_d,n_a,   n_z,pi_z,simoptions,Params,[]);
    StationaryDist2A=StationaryDist_InfHorz(Policy2A,n_d,n_a_2A,n_z,pi_z,simoptions,Params,[]);
    StationaryDist2A_a1=reshape(sum(StationaryDist2A,2),[n_a,N_z]);
    fprintf('%s: StationaryDist a1-marginal, 2A vs 1A, ~zero:        %2.8f \n',gistr,max(abs(StationaryDist2A_a1-reshape(StationaryDist1A,[n_a,N_z])),[],'all'));

    % AllStats on assets should agree
    AllStats1A=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1A,Policy1A,FnsToEvaluate_1A,Params,[],n_d,n_a,   n_z,d_grid,a_grid,   z_grid,simoptions);
    AllStats2A=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist2A,Policy2A,FnsToEvaluate_2A,Params,[],n_d,n_a_2A,n_z,d_grid,a_grid_2A,z_grid,simoptions);
    fprintf('%s: AllStats assets Mean, 2A vs 1A, should be zero:     %2.8f \n',gistr,abs(AllStats2A.assets.Mean-AllStats1A.assets.Mean));
    fprintf('%s: AllStats assets Gini, 2A vs 1A, should be zero:     %2.8f \n',gistr,abs(AllStats2A.assets.Gini-AllStats1A.assets.Gini));
end

output=struct();

end
