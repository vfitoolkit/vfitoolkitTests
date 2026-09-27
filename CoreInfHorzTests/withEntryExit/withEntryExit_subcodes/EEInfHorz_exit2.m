function output=EEInfHorz_exit2(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% endogenousexit=2: a mixture of no exit, endogenous exit and exogenous exit. GPU only (the
% parallel=0 and parallel=1 branches of ValueFnIter_InfHorz_EndogExit2 are commented out, so
% asking for them leaves VKron undefined rather than erroring cleanly; that is probed in the
% guards cross-test).
%
% Reaches ValueFnIter_InfHorz_EndogExit2_nod_Par2_raw and ValueFnIter_InfHorz_EndogExit2_Par2_raw.
%
% NOTE THE DIFFERENT ReturnToExitFn SIGNATURE. endogenousexit=1 builds the return to exit with
% CreateReturnToExitFnMatrix_Case1_Disc_Par2 from (a,z); endogenousexit=2 builds it with
% CreateReturnFnMatrix_Disc from (d,aprime,a,z) and then maximises over (d,aprime) to produce
% PolicyWhenExit. ReturnToExitFn_EE_aprime ignores aprime (this subcode has no d) and returns exactly
% ReturnToExitFn_EE, which is what makes the cross-test below an exact comparison.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,p,alpha,tau,cf,empcap) ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap);
ReturnFn_d=@(d,aprime,a,z,p,alpha,tau,cf,empcap,psi) ReturnFn_EE_d(d,aprime,a,z,p,alpha,tau,cf,empcap,psi);

FnsToEvaluate.employment=@(aprime,a,z) aprime;

%% CROSS-TEST: pure endogenous exit under endogenousexit=2 must reproduce endogenousexit=1
% exitprobabilities is [1-sum(named), named...], so exitprob_endog=1 and exitprob_exog=0 gives
% the weights [0,1,0]. With continuationcost=0 the surviving leg is
%     ExitPolicy.*FWhenExit+(1-ExitPolicy).*Vtemp
% which is exactly what endogenousexit=1 computes.
%
% This is also the only place in the bank that puts EXACTLY ZERO SCALAR weights against the
% other two legs. A scalar zero times a -Inf value is NaN in the same way an array one is, so
% this check is the probe for that shape as well as a cross-validation.
fprintf('\n----- cross-test: endogenousexit=2 with weights [0,1,0] vs endogenousexit=1 ----- \n')
fprintf('exitprobabilities are [%g,%g,%g]; the first and third being exactly zero is the point \n', ...
    1-(Params.exitprob_endog+Params.exitprob_exog),Params.exitprob_endog,Params.exitprob_exog)

vf2=vfoptionsbaseline;
vf2.endogenousexit=2;
vf2.ReturnToExitFn=@(aprime,a,z,tau,minexit) ReturnToExitFn_EE_aprime(aprime,a,z,tau,minexit); % (d,)aprime,a,z for exit2, and this subcode has no d
vf2.exitprobabilities={'exitprob_endog','exitprob_exog'};
vf2.endogenousexitcontinuationcost={'continuationcost'};
[V2,Policy2,PolicyWhenExit2,ExitPolicy2]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf2);
V2=gather(V2); Policy2=gather(Policy2); PolicyWhenExit2=gather(PolicyWhenExit2); ExitPolicy2=gather(ExitPolicy2);

vf1=vfoptionsbaseline;
[V1,Policy1,ExitPolicy1]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf1);
V1=gather(V1); Policy1=gather(Policy1); ExitPolicy1=gather(ExitPolicy1);

fprintf('exit2 V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V2(:))),sum(V2(:)==-Inf),sum(V2(:)==Inf),sum(isnan(V2(:))),numel(V2));
fprintf('exit1 V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V1(:))),sum(V1(:)==-Inf),sum(V1(:)==Inf),sum(isnan(V1(:))),numel(V1));
fprintf('NaN in exit2 V, this should be zero: %i \n',sum(isnan(V2(:))))
fprintf('states where exactly one of the two is finite, this should be zero: %i \n',sum(xor(isfinite(V2(:)),isfinite(V1(:)))))
fprintf('CrossTest (exit2 weights [0,1,0] vs exit1), this should be zero: V %.3e \n',max(abs(V2(:)-V1(:))))
fprintf('CrossTest (exit2 weights [0,1,0] vs exit1), this should be zero: ExitPolicy %.3e \n',max(abs(ExitPolicy2(:)-ExitPolicy1(:))))
% Policy has to be compared only where the firm STAYS. The two paths differ on purpose at the
% exiting states: endogenousexit=1 zeroes Policy there (the raws call it a 'deliberate' sentinel
% so that downstream mistakes throw rather than pass silently), whereas endogenousexit=2 leaves
% the policy in place, because with end-of-period exit timing the firm still acts in the period
% it exits. So the second line below is EXPECTED to be non-zero, and is printed as information.
stay=(ExitPolicy1==0);
P1=reshape(Policy1,[prod(n_a),prod(n_z)]); P2=reshape(Policy2,[prod(n_a),prod(n_z)]);
fprintf('CrossTest (exit2 weights [0,1,0] vs exit1), this should be zero: Policy where the firm stays %.3e \n',max(abs(P2(stay)-P1(stay))))
fprintf('  exit2 keeps the policy where exit1 zeroes it, so this one is EXPECTED non-zero: %.3e \n',max(abs(P2(~stay)-P1(~stay))))
fprintf('PolicyWhenExit has no NaN, this should be zero: %i \n',sum(isnan(PolicyWhenExit2(:))))

%% A genuine three-way mixture
% Now give all three legs positive weight, so nothing is degenerate and the exogenous-exit leg
% is actually doing something.
fprintf('\n----- genuine three-way exit mixture ----- \n')
Params3=Params;
Params3.exitprob_endog=0.6;
Params3.exitprob_exog=0.1;
Params3.continuationcost=0.05;
fprintf('exitprobabilities are [%g,%g,%g] \n',1-(Params3.exitprob_endog+Params3.exitprob_exog),Params3.exitprob_endog,Params3.exitprob_exog)

[V3,Policy3,PolicyWhenExit3,ExitPolicy3]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params3,DiscountFactorParamNames,[],vf2);
V3=gather(V3); ExitPolicy3=gather(ExitPolicy3);
fprintf('V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V3(:))),sum(V3(:)==-Inf),sum(V3(:)==Inf),sum(isnan(V3(:))),numel(V3));
fprintf('NaN in V, this should be zero: %i \n',sum(isnan(V3(:))))
fprintf('ExitPolicy is binary, this should be zero: %i \n',sum(~(ExitPolicy3(:)==0 | ExitPolicy3(:)==1)))
% Compare only where both are finite. V3 legitimately carries -Inf now (the exit2 fix turned the
% old NaN into the correct -Inf), and max(abs(finite-(-Inf))) is just Inf, which says nothing.
bothfinite=isfinite(V3(:)) & isfinite(V1(:));
fprintf('the mixture actually changed the answer, this should be non-zero: %.3e \n',max(abs(V3(bothfinite)-V1(bothfinite))))

%% Howards on/off, under the three-way mixture
vf2_noH=vf2;
vf2_noH.howards=0;
[V3noH,~,~,ExitPolicy3noH]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params3,DiscountFactorParamNames,[],vf2_noH);
fprintf('howards=0 (pure VFI), this should be zero: V %.3e \n',max(abs(V3(:)-gather(V3noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: ExitPolicy %.3e \n',max(abs(ExitPolicy3(:)-gather(ExitPolicy3noH(:)))))

%% With d
fprintf('\n----- endogenousexit=2 with d ----- \n')
vf2d=vf2;
[V4,Policy4,PolicyWhenExit4,ExitPolicy4]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn_d,Params3,DiscountFactorParamNames,[],vf2d);
V4=gather(V4); ExitPolicy4=gather(ExitPolicy4);
fprintf('V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V4(:))),sum(V4(:)==-Inf),sum(V4(:)==Inf),sum(isnan(V4(:))),numel(V4));
fprintf('NaN in V, this should be zero: %i \n',sum(isnan(V4(:))))
fprintf('ExitPolicy is binary, this should be zero: %i \n',sum(~(ExitPolicy4(:)==0 | ExitPolicy4(:)==1)))

%% Stationary distribution under endogenousexit=2
Params3.zeta=1-ExitPolicy3;
simoptions2=simoptionsbaseline;
simoptions2.endogenousexit=2;
simoptions2.exitprobabilities={'exitprob_endog','exitprob_exog'}; % StationaryDist needs these too, by name, exactly as vfoptions does
StationaryDist2=StationaryDist_InfHorz(Policy3,0,n_a,n_z,pi_z,simoptions2,Params3,EntryExitParamNames);
fprintf('\nStationaryDist.pdf sums to one, this should be zero: %.3e \n',abs(sum(gather(StationaryDist2.pdf(:)))-1))
fprintf('StationaryDist.pdf has no NaN, this should be zero: %i \n',sum(isnan(gather(StationaryDist2.pdf(:)))))
fprintf('StationaryDist.mass is finite and positive, this should be one: %i \n',isfinite(gather(StationaryDist2.mass)) && gather(StationaryDist2.mass)>0)

%% Plot the exit decision under the three-way mixture
fig=figure(figure_c); %#ok<NASGU>
surf(ExitPolicy3)
title('Exit decision, endogenousexit=2 (1 indicates exit)'); xlabel('z'); ylabel('a')

%%
output=struct();

end
