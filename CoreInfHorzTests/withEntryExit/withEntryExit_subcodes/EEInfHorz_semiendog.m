function output=EEInfHorz_semiendog(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,EntryExitParamNames,vfoptionsbaseline,simoptionsbaseline,figure_c)
% endogenousexit=1 with a SEMI-ENDOGENOUS shock, no d, GPU.
% Reaches ValueFnIter_InfHorz_EndogExit_SemiEndog_nod_Par2_raw.
%
% There is no 'with d' semi-endogenous raw: the dispatcher only looks for vfoptions.SemiEndogShock
% inside its n_d(1)==0 branch, so a semi-endogenous shock with d silently takes the plain path.
% That is probed in the guards cross-test rather than here.
%
% vfoptions.SemiEndogShock is [N_a,N_z,N_zprime]: the transition matrix of z is allowed to depend
% on the current endogenous state.

N_a=prod(n_a);
N_z=prod(n_z);

ReturnFn=@(aprime,a,z,p,alpha,tau,cf,empcap) ReturnFn_EE_nod(aprime,a,z,p,alpha,tau,cf,empcap);

FnsToEvaluate.employment=@(aprime,a,z) aprime;

%% CROSS-TEST first: a FAKE semi-endogenous shock must reproduce the plain path
% Build pi_z_semiendog so that it does not actually depend on a: every a-slice is pi_z itself.
% Then the semi-endogenous raw and the plain raw are solving the identical problem, and any
% difference is a defect in the semi-endogenous code path rather than in the economics.
fprintf('\n----- cross-test: fake (a-independent) semi-endogenous shock vs plain ----- \n')
pi_z_fake=repmat(shiftdim(gather(pi_z),-1),[N_a,1,1]); % [N_a,N_z,N_zprime], same for every a

vf_fake=vfoptionsbaseline;
vf_fake.SemiEndogShock=pi_z_fake;
[Vf,Policyf,ExitPolicyf]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf_fake);
Vf=gather(Vf); Policyf=gather(Policyf); ExitPolicyf=gather(ExitPolicyf);

vf_plain=vfoptionsbaseline;
[Vp,Policyp,ExitPolicyp]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vf_plain);
Vp=gather(Vp); Policyp=gather(Policyp); ExitPolicyp=gather(ExitPolicyp);

fprintf('fake-semiendog V census: %i finite, %i -Inf, %i NaN \n',sum(isfinite(Vf(:))),sum(Vf(:)==-Inf),sum(isnan(Vf(:))));
fprintf('plain          V census: %i finite, %i -Inf, %i NaN \n',sum(isfinite(Vp(:))),sum(Vp(:)==-Inf),sum(isnan(Vp(:))));
fprintf('states where exactly one of the two is finite, this should be zero: %i \n',sum(xor(isfinite(Vf(:)),isfinite(Vp(:)))))
fprintf('CrossTest (fake semiendog vs plain), this should be zero: V %.3e \n',max(abs(Vf(:)-Vp(:))))
fprintf('CrossTest (fake semiendog vs plain), this should be zero: Policy %.3e \n',max(abs(Policyf(:)-Policyp(:))))
fprintf('CrossTest (fake semiendog vs plain), this should be zero: ExitPolicy %.3e \n',max(abs(ExitPolicyf(:)-ExitPolicyp(:))))

%% A genuinely semi-endogenous shock
% Bigger firms are more persistent: mix pi_z towards the identity by an amount that rises with a.
fprintf('\n----- genuinely semi-endogenous shock ----- \n')
lambda=reshape(linspace(0,0.5,N_a),[N_a,1,1]); % weight on 'stay where you are', rising in a
pi_z_semiendog=(1-lambda).*repmat(shiftdim(gather(pi_z),-1),[N_a,1,1]) + lambda.*repmat(shiftdim(eye(N_z),-1),[N_a,1,1]);
fprintf('SemiEndogShock rows sum to one, this should be zero: %.3e \n',max(abs(sum(pi_z_semiendog,3)-1),[],'all'))

vfoptions1=vfoptionsbaseline;
vfoptions1.SemiEndogShock=pi_z_semiendog;
simoptions1=simoptionsbaseline;
simoptions1.SemiEndogShock=pi_z_semiendog;
[V1,Policy1,ExitPolicy1]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1);
V1=gather(V1); Policy1=gather(Policy1); ExitPolicy1=gather(ExitPolicy1);

fprintf('V census: %i finite, %i -Inf, %i +Inf, %i NaN (of %i) \n', ...
    sum(isfinite(V1(:))),sum(V1(:)==-Inf),sum(V1(:)==Inf),sum(isnan(V1(:))),numel(V1));
fprintf('NaN in V, this should be zero: %i \n',sum(isnan(V1(:))))
fprintf('ExitPolicy is binary, this should be zero: %i \n',sum(~(ExitPolicy1(:)==0 | ExitPolicy1(:)==1)))
fprintf('semi-endogenous shock actually changed the answer, this should be non-zero: %.3e \n',max(abs(V1(:)-Vp(:))))

%% Howards on/off
vfoptions1_noH=vfoptions1;
vfoptions1_noH.howards=0;
[V1noH,Policy1noH,ExitPolicy1noH]=ValueFnIter_InfHorz(0,n_a,n_z,[],a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions1_noH);
fprintf('howards=0 (pure VFI), this should be zero: V %.3e \n',max(abs(V1(:)-gather(V1noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: Policy %.3e \n',max(abs(Policy1(:)-gather(Policy1noH(:)))))
fprintf('howards=0 (pure VFI), this should be zero: ExitPolicy %.3e \n',max(abs(ExitPolicy1(:)-gather(ExitPolicy1noH(:)))))

%% Plot the exit decision
fig=figure(figure_c); %#ok<NASGU>
surf(ExitPolicy1)
title('Exit decision, semi-endogenous shock (1 indicates exit)'); xlabel('z'); ylabel('a')

%%
output=struct();

end
