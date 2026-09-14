function output=CoreInfHorzInheritAsset_d1_noa1_z_noe(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% The one tier the toolkit supports: d1, no a1, z, no e.
%
% There is no divide-and-conquer here (not usable for InfHorz), no grid interpolation layer
% (ValueFnIter_InfHorz_InheritAsset has no gridinterplayer branch), and no Howards (the Howards
% block in ValueFnIter_InfHorz_InheritAsset_noa1_raw is commented out). That leaves lowmemory as
% the only solver option to sweep, so most of the work here is done by independent
% recomputation rather than by comparing two toolkit paths against each other.

n_d1=n_d(1);
n_d2=n_d(2);
N_d1=prod(n_d1);
N_d2=prod(n_d2);
N_a=prod(n_a);
N_z=prod(n_z);

d1_grid=d_grid(1:N_d1);
d2_grid=d_grid(N_d1+1:end);

ReturnFn=@(d1,d2,a,z,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(d1,d2,a,z,r,w,sigma,eta,varphi);

% Setup some FnsToEvaluate. Note there is no aprime input for an inheritance asset, so the
% FnsToEvaluate take (d1,d2,a,z).
FnsToEvaluate.assets=@(d1,d2,a,z) a;
FnsToEvaluate.earnings=@(d1,d2,a,z,w) w*z*d1;
FnsToEvaluate.bequest=@(d1,d2,a,z) d2;

%% Baseline VFI
vfoptions1=vfoptionsbaseline;
simoptions1=simoptionsbaseline;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
V1=gather(V1); Policy1=gather(Policy1);

%% Census of V
% max(abs(A-B)) IGNORES NaN, so every 'this should be zero' check below is silently vacuous at
% any entry where a NaN sits. The census is the only instrument that sees the NaN class at all,
% which is why it comes first and why it is printed rather than differenced.
fprintf('V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V1(:))),sum(V1(:)==-Inf),sum(V1(:)==Inf),sum(isnan(V1(:))),numel(V1));
fprintf('NaN in V, this should be zero: %i \n',sum(isnan(V1(:))))

%% Is the zero-weight guard actually reached?
% The interpolation weight a2primeProbs is exactly 1 wherever a2prime falls off the bottom of
% a2_grid and exactly 0 wherever it falls off the top. Those are the only weights that can put a
% zero against a -Inf node, so if there are none of them the NaN check above proves nothing.
z_gridvals=CreateGridvals(n_z,gpuArray(z_grid),1);
aprimeFnParamsVec=CreateVectorFromParams(Params,{'inheritrisk'});
[a2primeIndex,a2primeProbs]=CreateInheritanceAssetFnMatrix(vfoptions1.aprimeFn, n_d2, n_a, n_z, n_z, gpuArray(d2_grid), gpuArray(a_grid), z_gridvals, z_gridvals, aprimeFnParamsVec);
a2primeIndex=gather(a2primeIndex); a2primeProbs=gather(a2primeProbs);
fprintf('a2primeProbs: %i exactly 0, %i exactly 1 (of %i). Both should be non-zero, or the zero-weight guard is never reached \n', ...
    sum(a2primeProbs(:)==0),sum(a2primeProbs(:)==1),numel(a2primeProbs));

%% Recompute V from Policy, independently of the solver
% ValueFnFromPolicy_InfHorz does not support inheritance assets (it has no aprimeFn handling at
% all), so the bank does the recomputation itself. This is the strongest check in the file: it
% reads BOTH channels of Policy, so it catches a wrong d1 as well as a wrong d2, and it rebuilds
% the continuation from CreateInheritanceAssetFnMatrix rather than trusting the solver's own.
ReturnFnParamsVec=CreateVectorFromParams(Params,{'r','w','sigma','eta','varphi'});
d_gridvals=CreateGridvals(n_d,gpuArray(d_grid),1);
ReturnMatrix=CreateReturnFnMatrix_Case2_Disc(ReturnFn, n_d, n_a, n_z, d_gridvals, gpuArray(a_grid), z_gridvals, ReturnFnParamsVec);
ReturnMatrix=gather(reshape(ReturnMatrix,[N_d1,N_d2,N_a,N_z]));

% The continuation, built the same way the raw builds it (including the zero-weight guard)
addindexforzprime=N_a*(0:1:N_z-1);
Vlower=reshape(V1(a2primeIndex+addindexforzprime),[N_d2,N_z,N_z]); % (d2,zprime,z)
Vupper=reshape(V1(a2primeIndex+1+addindexforzprime),[N_d2,N_z,N_z]);
skipinterp=(Vlower==Vupper);
probs2=a2primeProbs;
probs2(skipinterp)=0;
EV=probs2.*Vlower+(1-probs2).*Vupper;
EV(probs2==0)=Vupper(probs2==0); % a zero weight against an infinite node gives 0*(-Inf)=NaN
EV(probs2==1)=Vlower(probs2==1);
EV=EV.*shiftdim(pi_z',-1);
EV(isnan(EV))=0; % remove nan created where value fn is -Inf but probability is zero
EV=reshape(sum(EV,2),[N_d2,N_z]); % (d2,z)

d1p=reshape(Policy1(1,:,:),[N_a,N_z]);
d2p=reshape(Policy1(2,:,:),[N_a,N_z]);
avec=repmat((1:N_a)',1,N_z);
zvec=repmat(1:N_z,N_a,1);
Fatpolicy=reshape(ReturnMatrix(d1p+(d2p-1)*N_d1+(avec-1)*N_d1*N_d2+(zvec-1)*N_d1*N_d2*N_a),[N_a,N_z]);
EVatpolicy=reshape(EV(d2p+(zvec-1)*N_d2),[N_a,N_z]);
V1recomputed=Fatpolicy+Params.beta*EVatpolicy;

dV=abs(V1(:)-V1recomputed(:));
fprintf('V recomputed from Policy: %i finite differences, %i NaN differences (of %i) \n', ...
    sum(isfinite(dV)),sum(isnan(dV)),numel(dV));
fprintf('V recomputed from Policy, this should be zero: %.3e \n',max(dV))

%% Is the reported d1 actually optimal given the reported d2?
% The refine step solves d1*(d2,a,z) first, then maximises over d2, then recovers d1 by linear
% indexing into a [N_d2,N_a,N_z] array. A stride mistake in that recovery returns a valid but
% wrong d1, which the value check above catches only because it reads d1 too. This check says
% directly what went wrong, which is worth having separately when it fires.
d1best=zeros(N_a,N_z);
for z_c=1:N_z
    for a_c=1:N_a
        [~,d1best(a_c,z_c)]=max(ReturnMatrix(:,d2p(a_c,z_c),a_c,z_c));
    end
end
fprintf('reported d1 matches argmax over d1 at the reported d2, this should be zero: %.3e \n',max(abs(d1p(:)-d1best(:))))

%% V0-independence
% The value fn iteration is a contraction, so the converged V cannot depend on the starting
% guess. It did once: a2primeProbs was zeroed IN PLACE on an array built outside the while loop,
% so the zeroed entries accumulated across iterations and the answer depended on the path taken.
% Fixed in e0aff864; this is the regression test for it.
vfoptions1_V0=vfoptions1;
vfoptions1_V0.V0=-10*ones(N_a,N_z);
[V1_V0,Policy1_V0]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_V0);
fprintf('different V0 gives the same V, this should be zero: %.3e \n',max(abs(V1(:)-gather(V1_V0(:)))))
fprintf('different V0 gives the same Policy, this should be zero: %.3e \n',max(abs(Policy1(:)-gather(Policy1_V0(:)))))

%% lowmemory: lowmemory=1 should give the same V and Policy as lowmemory=0
vfoptions1_lm=vfoptions1;
vfoptions1_lm.lowmemory=1;
[V1B,Policy1B]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_lm);
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(V1(:)-gather(V1B(:)))))
fprintf('lowmemory=1, this should be zero: %.3e \n',max(abs(Policy1(:)-gather(Policy1B(:)))))

%% PolicyInd2Val
PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);
PolicyVals1=gather(PolicyVals1);
fprintf('PolicyInd2Val ran, d1 values in [%g,%g] (grid is [%g,%g]) \n', ...
    min(min(PolicyVals1(1,:,:))),max(max(PolicyVals1(1,:,:))),min(d1_grid),max(d1_grid));
fprintf('PolicyInd2Val ran, d2 values in [%g,%g] (grid is [%g,%g]) \n', ...
    min(min(PolicyVals1(2,:,:))),max(max(PolicyVals1(2,:,:))),min(d2_grid),max(d2_grid));

%% Stationary distribution and statistics
StationaryDist1=StationaryDist_InfHorz(Policy1,n_d,n_a,n_z,pi_z,simoptions1,Params,[]);
fprintf('StationaryDist sums to one, this should be zero: %.3e \n',abs(sum(gather(StationaryDist1(:)))-1))
fprintf('StationaryDist has no negative mass, this should be zero: %i \n',sum(gather(StationaryDist1(:))<0))

%% Do the two aprime builders agree?
% There are two, and they deliberately use TRANSPOSED layouts: CreateInheritanceAssetFnMatrix
% (used by the solver) returns [N_d2,N_zprime,N_z], while CreateaprimePolicyInheritanceAsset
% (used by the panel) returns [N_a,N_z,N_zprime] and says so in a comment reading 'Following are
% different to how they are in CreateInheritanceAssetFnMatrix()'. Nothing else in the toolkit
% cross-checks them, and a dimension mix-up between the two is silent: it gives a wrong but
% entirely plausible asset transition. Evaluating the policy builder at the solved Policy has to
% reproduce the solver builder at the chosen d2.
[a2primeIndexP,a2primeProbsP]=CreateaprimePolicyInheritanceAsset(gpuArray(Policy1),vfoptions1.aprimeFn,length(n_d),n_d,0,n_a,n_z,n_z,gpuArray(d_grid),gpuArray(a_grid),z_gridvals,z_gridvals,aprimeFnParamsVec);
a2primeIndexP=gather(a2primeIndexP); a2primeProbsP=gather(a2primeProbsP);
dIdx=zeros(N_a,N_z); dPrb=zeros(N_a,N_z);
for z_c=1:N_z
    for zp_c=1:N_z
        dIdx(:,z_c)=max(dIdx(:,z_c),abs(a2primeIndexP(:,z_c,zp_c)-a2primeIndex(d2p(:,z_c),zp_c,z_c)));
        dPrb(:,z_c)=max(dPrb(:,z_c),abs(a2primeProbsP(:,z_c,zp_c)-a2primeProbs(d2p(:,z_c),zp_c,z_c)));
    end
end
fprintf('the two aprime builders agree at the solved Policy, these should be zero: index %.3e, prob %.3e \n',max(dIdx(:)),max(dPrb(:)))

%% Is the stationary distribution actually invariant under the inheritance-asset transition?
% This is the distribution-side counterpart of 'V recomputed from Policy' above, and the bank
% needs it: the cross-test below is DE-RISKED (its aprimeFn ignores zprime, which is what makes
% the match exact), so nothing else here exercises the zprime-dependent asset transition, and
% that is the whole point of an inheritance asset.
%
% Deliberately written as a dumb triple loop over (a,z,zprime) rather than vectorised. The
% toolkit builds this transition with a vectorised scatter; an independent check that shares that
% construction can share its bug, and 101*5*5 iterations cost nothing.
StationaryDist1g=gather(StationaryDist1);
pi_zg=gather(pi_z); % keep the loop off the GPU; scalar indexing a gpuArray 2525 times is pointless
Distnext=zeros(N_a,N_z);
for z_c=1:N_z
    for zp_c=1:N_z
        for a_c=1:N_a
            mass=StationaryDist1g(a_c,z_c)*pi_zg(z_c,zp_c);
            ilower=a2primeIndex(d2p(a_c,z_c),zp_c,z_c);
            plower=a2primeProbs(d2p(a_c,z_c),zp_c,z_c);
            Distnext(ilower,zp_c)=Distnext(ilower,zp_c)+mass*plower;
            Distnext(ilower+1,zp_c)=Distnext(ilower+1,zp_c)+mass*(1-plower);
        end
    end
end
fprintf('one-period-ahead dist sums to one, this should be zero: %.3e \n',abs(sum(Distnext(:))-1))
fprintf('StationaryDist is invariant under the bank''s own transition, this should be zero: %.3e \n', ...
    max(abs(Distnext(:)-StationaryDist1g(:))))
fprintf('  (mean of a under StationaryDist vs under its one-period-ahead: %.6f %.6f) \n', ...
    sum(sum(StationaryDist1g,2).*a_grid),sum(sum(Distnext,2).*a_grid))

AllStats1=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions1);
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions1);
fprintf('AllStats and AggVars means should agree, these should be zero: %.3e %.3e %.3e \n', ...
    abs(AllStats1.assets.Mean-AggVars1.assets.Mean), ...
    abs(AllStats1.earnings.Mean-AggVars1.earnings.Mean), ...
    abs(AllStats1.bequest.Mean-AggVars1.bequest.Mean))

ValuesOnGrid1=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions1); %#ok<NASGU>
fprintf('ValuesOnGrid ran \n')

%% SimPanelValues
% Monte Carlo, so these are roughly-equal checks and must never be compared as an exact zero
% (the parfor workers do not take the client's rng, so a panel is not reproducible run to run).
simoptionsPanel=simoptions1;
simoptionsPanel.numbersims=10^4;
simoptionsPanel.simperiods=100;
simoptionsPanel.burnin=0;
SimPanel1=SimPanelValues_InfHorz(StationaryDist1,Policy1,FnsToEvaluate,[],Params,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptionsPanel);
spA=SimPanel1.assets(:); spE=SimPanel1.earnings(:);
fprintf('SimPanelValues mean should roughly match AllStats mean (Monte Carlo) \n')
[AllStats1.assets.Mean, mean(spA)]
[AllStats1.earnings.Mean, mean(spE)]
fprintf('SimPanelValues std should roughly match AllStats std \n')
[AllStats1.assets.StdDeviation, std(spA,1)]

%% Panel by period: where does a mismatch enter?
% SimPanel1.assets is [simperiods,numbersims]. Period 1 is a straight draw from InitialDist with
% no transition applied, so it must match the stationary dist up to sampling error whatever the
% transition does. Period 2 is the first period that has been through exactly one transition, and
% the last period is where any per-period bias has fully compounded. Splitting them says whether
% a panel-vs-dist gap comes from the draw or from the transition, which pooling over all periods
% cannot tell you.
fprintf('panel assets by period (mean; AllStats mean is %.6f) \n',AllStats1.assets.Mean)
fprintf('  period 1 (pure draw, no transition): %.6f \n',mean(SimPanel1.assets(1,:)))
fprintf('  period 2 (one transition):           %.6f \n',mean(SimPanel1.assets(2,:)))
fprintf('  period %i (fully compounded):        %.6f \n',simoptionsPanel.simperiods,mean(SimPanel1.assets(end,:)))
fprintf('panel assets by period (std; AllStats std is %.6f) \n',AllStats1.assets.StdDeviation)
fprintf('  period 1 (pure draw, no transition): %.6f \n',std(SimPanel1.assets(1,:),1))
fprintf('  period 2 (one transition):           %.6f \n',std(SimPanel1.assets(2,:),1))
fprintf('  period %i (fully compounded):        %.6f \n',simoptionsPanel.simperiods,std(SimPanel1.assets(end,:),1))
fprintf('panel earnings period 1 vs AllStats (control: a pure function of z, so this isolates the asset transition) \n')
fprintf('  %.6f %.6f \n',AllStats1.earnings.Mean,mean(SimPanel1.earnings(1,:)))

%% Plot the stationary dist over the inheritance asset
fig=figure(figure_c); %#ok<NASGU>
plot(a_grid,cumsum(sum(gather(StationaryDist1),2)))
title('CDF of the inheritance asset'); xlabel('a2')

%%
output=struct(); % Not currently used for anything. Maybe will do so later.

end
