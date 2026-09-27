function output=EEInfHorz_d(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% endogenousexit=1, with d, GPU. Reaches ValueFnIter_InfHorz_EndogExit_Par2_raw.
%
% Also carries the d-vs-nod cross-test: a model with n_d=1 and d_grid=1 has a utilisation cost
% of exactly zero, so its return function is identical to the no-d one, and the two code paths
% (EndogExit_Par2_raw and EndogExit_nod_Par2_raw) must agree to machine zero.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn_d=@(d,aprime,a,z,p,alpha,tau,cf,empcap,psi) ReturnFn_EE_d(d,aprime,a,z,p,alpha,tau,cf,empcap,psi);
ReturnFn_nod=@(aprime,a,z,p,alpha,tau,cf,empcap) ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap);

FnsToEvaluate.employment=@(d,aprime,a,z) aprime;
FnsToEvaluate.output=@(d,aprime,a,z,alpha) z*((d*aprime)^alpha);

%% Baseline, with d
vfoptions1=vfoptionsbaseline;
simoptions1=simoptionsbaseline;
[V1,Policy1,ExitPolicy1]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfoptions1);
V1=gather(V1); Policy1=gather(Policy1); ExitPolicy1=gather(ExitPolicy1);

fprintf('V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V1(:))),sum(V1(:)==-Inf),sum(V1(:)==Inf),sum(isnan(V1(:))),numel(V1));
fprintf('NaN in V, this should be zero: %i \n',sum(isnan(V1(:))))
fprintf('ExitPolicy is binary, this should be zero: %i \n',sum(~(ExitPolicy1(:)==0 | ExitPolicy1(:)==1)))

%% Howards on/off
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1noH,Policy1noH,ExitPolicy1noH]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('howards=0 (pure VFI), this should be zero: V %.3e \n',max(abs(V1(:)-gather(V1noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: Policy %.3e \n',max(abs(Policy1(:)-gather(Policy1noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: ExitPolicy %.3e \n',max(abs(ExitPolicy1(:)-gather(ExitPolicy1noH(:)))))

%% keeppolicyonexit
vfoptions1_keep=vfoptions1;
vfoptions1_keep.keeppolicyonexit=1;
[~,Policy1keep,ExitPolicy1keep]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfoptions1_keep);
Policy1keep=gather(Policy1keep); ExitPolicy1keep=gather(ExitPolicy1keep);
stay=(ExitPolicy1keep==0);
Pd=reshape(Policy1(1,:,:),[N_a,N_z]); Pa=reshape(Policy1(2,:,:),[N_a,N_z]);
Pdk=reshape(Policy1keep(1,:,:),[N_a,N_z]); Pak=reshape(Policy1keep(2,:,:),[N_a,N_z]);
fprintf('keeppolicyonexit agrees where the firm stays, these should be zero: d %.3e, aprime %.3e \n', ...
    max(abs(Pd(stay)-Pdk(stay))),max(abs(Pa(stay)-Pak(stay))))
fprintf('default Policy is zero where the firm exits, these should be zero: d %.3e, aprime %.3e \n', ...
    max(abs(Pd(~stay))),max(abs(Pa(~stay))))

%% CROSS-TEST: n_d=1 at d=1 must reproduce the no-d path
% Same Bellman, two different raws. At d=1 the utilisation cost psi*(d-1)^2 is exactly zero, so
% ReturnFn_EE_d(1,...) is ReturnFn_EE_nod(...) term for term.
fprintf('\n----- cross-test: trivial d (n_d=1, d=1) vs no d ----- \n')
[Vtriv,Policytriv,ExitPolicytriv]=ValueFnIter_InfHorz(1,n_a,n_z,1,a_grid,z_grid,pi_z,ReturnFn_d,Params,DiscountFactorParamNames,[],vfoptions1);
[Vnod,Policynod,ExitPolicynod]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn_nod,Params,DiscountFactorParamNames,[],vfoptions1);
Vtriv=gather(Vtriv); Vnod=gather(Vnod);
Policytriv=gather(Policytriv); Policynod=gather(Policynod);
fprintf('trivial-d V census: %i finite, %i -Inf, %i NaN \n',sum(isfinite(Vtriv(:))),sum(Vtriv(:)==-Inf),sum(isnan(Vtriv(:))));
fprintf('no-d      V census: %i finite, %i -Inf, %i NaN \n',sum(isfinite(Vnod(:))),sum(Vnod(:)==-Inf),sum(isnan(Vnod(:))));
fprintf('states where exactly one of the two is finite, this should be zero: %i \n', ...
    sum(xor(isfinite(Vtriv(:)),isfinite(Vnod(:)))))
fprintf('CrossTest (trivial d vs no d), this should be zero: V %.3e \n',max(abs(Vtriv(:)-Vnod(:))))
fprintf('CrossTest (trivial d vs no d), this should be zero: aprime %.3e \n', ...
    max(abs(reshape(Policytriv(2,:,:),[N_a*N_z,1])-reshape(Policynod,[N_a*N_z,1]))))
fprintf('CrossTest (trivial d vs no d), this should be zero: ExitPolicy %.3e \n', ...
    max(abs(gather(ExitPolicytriv(:))-gather(ExitPolicynod(:)))))

%% Stationary distribution and statistics
Params.zeta=1-ExitPolicy1;
StationaryDist1=StationaryDist_InfHorz(Policy1,n_d,n_a,n_z,pi_z,simoptions1,Params,EntryExitParamNames);
fprintf('\nStationaryDist.pdf sums to one, this should be zero: %.3e \n',abs(sum(gather(StationaryDist1.pdf(:)))-1))
fprintf('StationaryDist.pdf has no NaN, this should be zero: %i \n',sum(isnan(gather(StationaryDist1.pdf(:)))))
AggVars1=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist1,Policy1,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions1,EntryExitParamNames);
fprintf('AggVars employment %g, output %g \n',AggVars1.employment.Aggregate,AggVars1.output.Aggregate)
fprintf('AggVars have no NaN, this should be zero: %i \n',sum(isnan([AggVars1.employment.Aggregate,AggVars1.output.Aggregate])))

%% Plot the exit decision
fig=figure(figure_c); %#ok<NASGU>
surf(ExitPolicy1)
title('Exit decision, with d (1 indicates exit)'); xlabel('z'); ylabel('a')

%%
output=struct();

end
