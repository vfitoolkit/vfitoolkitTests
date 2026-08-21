function output=CoreInfHorz_nod_z_noe_nosemiz_with2A(Params,DiscountFactorParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% Two-endogenous-state InfHorz test (fig 11), mirroring the fig-3 (nod, z) battery.
% Model: simplified Kitao (2008), no taxes, partial equilibrium (fixed r,w).
%   Endogenous states: a (assets), e (occupation: 0=worker, 1=entrepreneur)
%   Exogenous markov z: eta (labor productivity), theta (entrepreneurial ability)
% This exercises the InfHorz two-endogenous-state (GI2A) code paths.
%
% Note: like the published Kitao replication, an entrepreneur with (near) zero assets
% has no feasible positive consumption, so V takes -Inf at those states. The toolkit's
% iterated VFI and GI handle this (agents never optimally transition there).

DF=DiscountFactorParamNames;

%% Grids and model parameters (self-contained)
n_asset=101;
n_asset_big=501; % for the with/without GI moment check
n_e=2;           % occupation: worker/entrepreneur
n_eta=3;         % labor productivity
n_theta=2;       % entrepreneurial ability

n_a=[n_asset,n_e];
n_a_big=[n_asset_big,n_e];
n_z=[n_eta,n_theta];

% Assets: cubic-spaced, puts points near zero
assetmaxfactor=60;
asset_grid=assetmaxfactor*(linspace(0,1,n_asset).^3)';
asset_grid_big=assetmaxfactor*(linspace(0,1,n_asset_big).^3)';
e_grid=[0;1]; % 0=worker, 1=entrepreneur
a_grid=[asset_grid; e_grid];
a_grid_big=[asset_grid_big; e_grid];

% Exogenous shocks
[eta_grid,pi_eta]=discretizeAR1_FarmerToda(0,0.9,0.2,n_eta);
eta_grid=exp(eta_grid);
theta_grid=[0; 1.5];               % entrepreneurial ability (0 => no production)
pi_theta=[0.9,0.1; 0.2,0.8];       % persistent
z_grid=[eta_grid; theta_grid];
pi_z=kron(pi_theta,pi_eta);        % eta varies fastest (matches n_z=[n_eta,n_theta])

n_d=0; d_grid=[];

% Kitao parameters (no taxes)
Params.r=0.04; Params.w=1.4;
Params.alpha=0.36; Params.delta=0.06; Params.upsilon=0.88;
Params.upsilon1=Params.alpha*Params.upsilon;
Params.upsilon2=(1-Params.alpha)*Params.upsilon;
Params.leverage=0.5; % max borrowing leverage (was 'd' in the Kitao code)
Params.phi=0.05;     % extra borrowing cost
% Params.sigma, Params.beta come from the shared setup

ReturnFn=@(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi) ...
    ReturnFn_nod_z_noe_nosemiz_with2A(aprime,eprime,a,e,eta,theta,r,w,sigma,delta,upsilon1,upsilon2,leverage,phi);

% FnsToEvaluate (functions of (aprime,eprime,a,e,eta,theta))
FnsToEvaluate.assets=@(aprime,eprime,a,e,eta,theta) a;
FnsToEvaluate.entrepreneur=@(aprime,eprime,a,e,eta,theta) e; % fraction who are entrepreneurs

%% Baseline VFI
vfoptions=struct(); simoptions=struct();
vfoptions.verbose_advice=0; % 2-endo GI would otherwise sound the postGI advice
vfoptions1=vfoptions; simoptions1=simoptions;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1);

PolicyVals1=PolicyInd2Val_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,vfoptions1);

V1fromPolicy=ValueFnFromPolicy_InfHorz(Policy1,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions1);
fprintf('with2A: ValueFnFromPolicy, this should be zero: %2.8f \n',max(abs(V1fromPolicy(:)-V1(:))))

%% Grid-interpolation layer (GI2A: interpolates the asset dimension)
vfoptions3=vfoptions;
vfoptions3.gridinterplayer=1;
vfoptions3.ngridinterp=5;
simoptions3=simoptions;
simoptions3.gridinterplayer=vfoptions3.gridinterplayer;
simoptions3.ngridinterp=vfoptions3.ngridinterp;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3);

PolicyVals3=PolicyInd2Val_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,vfoptions3);

V3fromPolicy=ValueFnFromPolicy_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions3);
fprintf('with2A: ValueFnFromPolicy with grid interp, this should be zero: %2.8f \n',max(abs(V3fromPolicy(:)-V3(:))))

%% TEMPORARY DIAGNOSTIC: the line above is the only check in this bank that does not come out at ~1e-8.
% It is ~1e-5, and only for this model (the with-d version of it, and the 2Aignored cross-tests which
% use the same raws, are all fine), so it looks model specific rather than a broken code path. These
% prints are to find out whether it is one boundary state or something systematic. Delete once known.
dV=abs(V3fromPolicy(:)-V3(:));
sz=size(V3);
fprintf('  diag: %d of %d states differ by more than 1e-7 \n',sum(dV>1e-7),numel(dV));
fprintf('  diag: %d states of V3 are non-finite \n',sum(~isfinite(V3(:))));
[~,iworst]=max(dV);
% V3 is [n_a(1), n_a(2), n_z...] and n_z can itself be multi-dimensional, so decode over every
% dimension. Calling ind2sub with fewer outputs than there are dimensions silently folds the
% trailing ones into the last output, which is what made the previous run print 'z=6 of 3'.
sub=cell(1,numel(sz));
[sub{:}]=ind2sub(sz,iworst);
fprintf('  diag: worst at subscripts [%s] within size [%s] \n',num2str(cell2mat(sub)),num2str(sz));
fprintf('  diag: there V3=%g, V3fromPolicy=%g \n',V3(iworst),V3fromPolicy(iworst));
fprintf('  diag: Policy3 has %d channels, size [%s] \n',size(Policy3,1),num2str(size(Policy3)));
bad=(dV>1e-7);
if size(Policy3,1)>=4 && any(bad)
    P3=reshape(Policy3,4,[]); % [4, N_a*N_z], same linear ordering as V3(:)
    fprintf('  diag: Policy there = [L1=%g, a2prime=%g, L2=%g, L2flag=%g] \n',P3(1,iworst),P3(2,iworst),P3(3,iworst),P3(4,iworst));
    % Are the bad states all at an a1 edge, all at one L2, or all flagged?
    Policy3L1=reshape(Policy3(1,:,:,:),[],1);
    Policy3L2=reshape(Policy3(3,:,:,:),[],1);
    Policy3flag=reshape(Policy3(4,:,:,:),[],1);
    fprintf('  diag: among the differing states, L1 ranges %d to %d (a1 grid is 1 to %d) \n',min(Policy3L1(bad)),max(Policy3L1(bad)),sz(1));
    fprintf('  diag: among the differing states, L2 ranges %d to %d (L2 runs 1 to %d) \n',min(Policy3L2(bad)),max(Policy3L2(bad)),vfoptions3.ngridinterp+2);
    fprintf('  diag: among the differing states, how many have L2flag~=2: %d \n',sum(Policy3flag(bad)~=2));
    fprintf('  diag: total states in the whole model with L2flag~=2: %d \n',sum(Policy3flag~=2));
end
% Print the actual numbers at the one differing state, rather than reasoning about it from outside.
% The policy there is L1=3, a2prime=1, so the two interpolation points are joint a indices 3 and 4.
Vlin=reshape(V3,[],prod(sz(3:end))); % [N_a, N_z]
NZtot=prod(sz(3:end));
P3lin=reshape(Policy3,4,[]);
il=P3lin(1,iworst)+sz(1)*(P3lin(2,iworst)-1); % L1 + n_a(1)*(a2prime-1)
iu=il+1;
izworst=ceil(iworst/prod(sz(1:2))); % which z the worst state is in
fprintf('  diag: interpolation points are joint a indices %d and %d, worst state is in z %d of %d \n',il,iu,izworst,NZtot);
fprintf('  diag: V at lower point across z: %s \n',num2str(Vlin(il,:),'%14.7g'));
fprintf('  diag: V at upper point across z: %s \n',num2str(Vlin(iu,:),'%14.7g'));
fprintf('  diag: non-finite among those two rows: %d \n',sum(~isfinite([Vlin(il,:),Vlin(iu,:)])));
fprintf('  diag: the %d non-finite V3 states are at linear indices: %s \n',sum(~isfinite(V3(:))),num2str(find(~isfinite(V3(:)))'));
pizfull=reshape(pi_z,NZtot,[]);
fprintf('  diag: pi_z row for that z: %s \n',num2str(pizfull(izworst,:),'%14.7g'));
% If one z-prime contribution were being lost, the error would be beta*pi_z(z,zprime)*V(lower,zprime).
% If none of these is near the observed discrepancy then no dropped transition explains it.
fprintf('  diag: beta*pi_z*V(lower) per zprime: %s \n',num2str(Params.beta*pizfull(izworst,:).*Vlin(il,:),'%14.7g'));
fprintf('  diag: observed discrepancy for comparison: %2.10f \n',max(abs(V3fromPolicy(:)-V3(:))));

% FofPolicy inside ValueFnFromPolicy is built by running PolicyInd2Val_InfHorz on the policy and then
% evaluating the ReturnFn on those values, so the aprime VALUES are the remaining suspect. This state
% has L2=1, meaning the chosen aprime sits exactly on coarse grid point L1, so PolicyInd2Val has to
% return exactly asset_grid(L1) with no interpolation at all. That is checkable against an exact number.
PV=PolicyInd2Val_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,vfoptions3);
PVlin=reshape(PV,size(PV,1),[]);
L1w=P3lin(1,iworst); a2w=P3lin(2,iworst);
fprintf('  diag: PolicyInd2Val aprime at worst state: [%s] \n',num2str(PVlin(:,iworst)','%18.12g'));
fprintf('  diag: asset_grid(L1)=%0.12g, asset_grid(L1+1)=%0.12g, e_grid(a2prime)=%0.12g \n',asset_grid(L1w),asset_grid(L1w+1),e_grid(a2w));
fprintf('  diag: a1prime minus asset_grid(L1) = %g   (L2=1, so this should be exactly 0) \n',PVlin(1,iworst)-asset_grid(L1w));
fprintf('  diag: a2prime minus e_grid(a2prime) = %g \n',PVlin(2,iworst)-e_grid(a2w));
Fimp=V3(iworst)-Params.beta*sum(pizfull(izworst,:).*Vlin(il,:));
fprintf('  diag: F implied by V3 at that state = %0.12g \n',Fimp);

% The aprime values are exact, so on the face of it ReturnFn gets identical inputs and FofPolicy
% should match. Reproduce FofPolicy here exactly as ValueFnFromPolicy_InfHorz_GI builds it, and
% compare against the F that the value fn's own V implies. If they differ then it is the return
% evaluation after all (and given aprime is exact, that would point at the state values rather than
% the policy values). If they agree then the discrepancy is inside the iteration instead.
PVperm=permute(reshape(PV,[size(PV,1),prod(n_a),prod(n_z)]),[2,3,1]);
agv=CreateGridvals(n_a,a_grid,1);
[zgv,~,~]=ExogShockSetup_InfHorz(n_z,z_grid,pi_z,Params,vfoptions3,3);
RFPC=CreateCellFromParams(Params,ReturnFnParamNamesFn(ReturnFn,n_d,n_a,n_z,0,vfoptions3,Params));
FofP=EvalFnOnAgentDist_Grid(ReturnFn,RFPC,PVperm,2,n_a,n_z,agv,zgv);
fprintf('  diag: FofPolicy at worst state = %0.12g \n',FofP(iworst));
fprintf('  diag: F implied by V3          = %0.12g \n',Fimp);
fprintf('  diag: difference               = %g \n',FofP(iworst)-Fimp);
% For context, how many states have FofPolicy differing from the F implied by their own V by >1e-7.
% (only valid where the policy sits exactly on a coarse point, so just report the worst state's.)
fprintf('  diag: a and z values ReturnFn saw there: a=[%s], z=[%s] \n',num2str(agv(mod(iworst-1,prod(n_a))+1,:),'%14.8g'),num2str(zgv(izworst,:),'%14.8g'));

% Is the residual just where the ValueFnFromPolicy fixed-point iteration stopped? That loop ends on
% max(abs(VKron-VKronold)), and where V is -Inf those differences are NaN, which max quietly skips,
% so the stopping test only sees a subset of the states. If tightening the tolerance collapses the
% difference then there is nothing wrong here; if it does not, the residual is structural.
vfoptions3tight=vfoptions3; vfoptions3tight.tolerance=1e-12;
V3fromPolicyTight=ValueFnFromPolicy_InfHorz(Policy3,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,vfoptions3tight);
fprintf('  diag: same check but with vfoptions.tolerance=1e-12: %2.10f \n',max(abs(V3fromPolicyTight(:)-V3(:))))
fprintf('  diag: (it was %2.10f at the default tolerance) \n',max(abs(V3fromPolicy(:)-V3(:))))
clear V3fromPolicyTight vfoptions3tight

clear V1 V3 Policy1 Policy3 PolicyVals1 PolicyVals3 V1fromPolicy V3fromPolicy

%% Big asset grid: stationary dist and moments should be essentially the same with/without grid interp
[~,Policy1b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1);
StationaryDist1=StationaryDist_InfHorz(Policy1b,n_d,n_a_big,n_z,pi_z,simoptions1,Params,[]);
AllStats1=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);

[~,Policy3b]=ValueFnIter_InfHorz(n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3);
StationaryDist3=StationaryDist_InfHorz(Policy3b,n_d,n_a_big,n_z,pi_z,simoptions3,Params,[]);
AllStats3=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);
AggVars3=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist3,Policy3b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions3);

fprintf('with2A: with/without grid interp, should get much the same moments (for big asset grid) \n')
fprintf('with2A: StationaryDist with/without grid interp, close to zero: %2.8f \n',max(abs(StationaryDist1(:)-StationaryDist3(:))))
[AllStats1.assets.Mean,AllStats3.assets.Mean]
[AllStats1.entrepreneur.Mean,AllStats3.entrepreneur.Mean]  % fraction entrepreneurs
[AggVars1.assets.Mean,AggVars3.assets.Mean]

%% Check the remaining InfHorz dist commands run for two endogenous states
AutoCorr1=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptions1);
CrossSectionCorr1=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ValuesOnGrid1=EvalFnOnAgentDist_ValuesOnGrid_InfHorz(Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
ProbDensityFns1=EvalFnOnAgentDist_ProbDensityFn_InfHorz(StationaryDist1,Policy1b,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,simoptions1);
fprintf('with2A: AutoCorrTransProbs, CrossSectionCovarCorr, ValuesOnGrid, ProbDensityFn all ran \n')

%% SimTimeSeries: a long simulated series should reproduce the AgentDist mean (Monte Carlo)
simoptionsTS=simoptions1;
simoptionsTS.simperiods=10^5;
simoptionsTS.burnin=10^3;
TimeSeries1=SimTimeSeriesValues_InfHorz(Policy1b,FnsToEvaluate,Params,n_d,n_a_big,n_z,d_grid,a_grid_big,z_grid,pi_z,simoptionsTS);
fprintf('with2A: SimTimeSeries mean should roughly match AllStats/AggVars mean (Monte Carlo) \n')
[AllStats1.assets.Mean, AggVars1.assets.Mean, mean(TimeSeries1.assets(:))]
[AllStats1.entrepreneur.Mean, AggVars1.entrepreneur.Mean, mean(TimeSeries1.entrepreneur(:))]

%% Plot: asset distribution and entrepreneur share, with/without grid interp
fig=figure(figure_c);
% marginal over assets: StationaryDist is [n_asset,n_e,n_eta,n_theta], so sum over dimensions 2,3,4
assetdist1=sum(sum(sum(StationaryDist1,4),3),2);
assetdist3=sum(sum(sum(StationaryDist3,4),3),2);
subplot(2,1,1); plot(asset_grid_big,cumsum(assetdist1), asset_grid_big,cumsum(assetdist3))
title('with2A: CDF of assets (without vs with grid interp)'); legend('1','3')
subplot(2,1,2); plot(1:2,[AllStats1.entrepreneur.Mean,AllStats3.entrepreneur.Mean],'o')
title('with2A: entrepreneur share')

%% Howards iteration
% Howards improvement iterations are just an accelerator for the value function iteration, so
% turning them off (vfoptions.howards=0, which is then pure value function iteration) must give
% the same V and Policy.
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1,Policy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1);
[V1noH,Policy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions1_noH);
fprintf('howards=0 (pure VFI), this should be zero: %2.8f \n',max(abs(V1(:)-V1noH(:))))
fprintf('howards=0 (pure VFI), this should be zero: %2.8f \n',max(abs(Policy1(:)-Policy1noH(:))))

% Same again, with the grid interpolation layer
vfoptions3_noH=vfoptions3;
vfoptions3_noH.howards=0;
[V3,Policy3]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3);
[V3noH,Policy3noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DF,[],vfoptions3_noH);
fprintf('howards=0 (pure VFI, with GI), this should be zero: %2.8f \n',max(abs(V3(:)-V3noH(:))))
fprintf('howards=0 (pure VFI, with GI), this should be zero: %2.8f \n',max(abs(Policy3(:)-Policy3noH(:))))

clear V1 V3 V1noH V3noH Policy1 Policy3 Policy1noH Policy3noH

%%
output=struct();

end
