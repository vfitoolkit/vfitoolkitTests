function output=CoreInfHorzInheritAsset_CrossTests_vsplain(n_d,n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,figure_c)
% CROSS-TEST: a de-risked inheritance asset must reproduce an ordinary one-asset model.
%
% Set aprimeFn to return d2 itself, ignoring the shock transition, and set d2_grid equal to
% a_grid. Then a2prime = d2 exactly, so choosing d2 IS choosing aprime and the inheritance-asset
% model is literally the plain model with d2 renamed. Everything must agree to machine zero:
% V, both Policy channels, the stationary distribution and the aggregates.
%
% This is the only instrument in the bank that tests the inheritance-asset machinery against
% code that is independently exercised (the plain InfHorz path is covered by CoreInfHorzTests),
% so it is where most of the confidence comes from.
%
% It also reaches both branches of the zero-weight guard, for free. a2prime lands exactly on a
% grid point every time, so the interpolation weight is exactly 1 everywhere except at the top
% grid point, where the toolkit's off-the-top-of-grid handling makes it exactly 0. Those are
% precisely the weights that used to produce 0*(-Inf)=NaN.

n_d1=n_d(1);
N_d1=prod(n_d1);
N_a=prod(n_a);
N_z=prod(n_z);

% d2_grid IS a_grid, so the two models have the same choice set for next period's asset
n_d2_x=n_a;
d1_grid=linspace(0,1,N_d1)';
d_grid_x=[d1_grid; a_grid];
n_d_x=[n_d1,n_d2_x];

aprimeFn_x=@(d2,z,zprime,inheritrisk) d2+0*inheritrisk*(z-zprime); % ignores the shock transition

ReturnFn_x=@(d1,d2,a,z,r,w,sigma,eta,varphi) ReturnFn_d1d2_a_z(d1,d2,a,z,r,w,sigma,eta,varphi);
ReturnFn_plain=@(d1,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_plain_d1_aprime_a_z(d1,aprime,a,z,r,w,sigma,eta,varphi);

FnsToEvaluate_x.assets=@(d1,d2,a,z) a;
FnsToEvaluate_x.earnings=@(d1,d2,a,z,w) w*z*d1;
FnsToEvaluate_plain.assets=@(d1,aprime,a,z) a;
FnsToEvaluate_plain.earnings=@(d1,aprime,a,z,w) w*z*d1;

fprintf('\n===== cross-test: de-risked inheritance asset vs plain one-asset model ===== \n')

%% The inheritance-asset version
vf_x=struct();
vf_x.inheritanceasset=1;
vf_x.aprimeFn=aprimeFn_x;
[Vx,Policyx]=ValueFnIter_InfHorz(n_d_x,n_a,n_z,d_grid_x,a_grid,z_grid,pi_z,ReturnFn_x,Params,DiscountFactorParamNames,[],vf_x);
Vx=gather(Vx); Policyx=gather(Policyx);

so_x=struct();
so_x.inheritanceasset=1;
so_x.aprimeFn=aprimeFn_x;
so_x.a_grid=a_grid;
so_x.d_grid=d_grid_x;
so_x.z_grid=z_grid;
StationaryDistx=StationaryDist_InfHorz(Policyx,n_d_x,n_a,n_z,pi_z,so_x,Params,[]);
AggVarsx=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDistx,Policyx,FnsToEvaluate_x,Params,[],n_d_x,n_a,n_z,d_grid_x,a_grid,z_grid,so_x);

%% The plain version
vf_p=struct();
[Vp,Policyp]=ValueFnIter_InfHorz(n_d1,n_a,n_z,d1_grid,a_grid,z_grid,pi_z,ReturnFn_plain,Params,DiscountFactorParamNames,[],vf_p);
Vp=gather(Vp); Policyp=gather(Policyp);

so_p=struct();
StationaryDistp=StationaryDist_InfHorz(Policyp,n_d1,n_a,n_z,pi_z,so_p,Params,[]);
AggVarsp=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDistp,Policyp,FnsToEvaluate_plain,Params,[],n_d1,n_a,n_z,d1_grid,a_grid,z_grid,so_p);

%% Compare
% Census first: if either V carries NaN the differences below are vacuous wherever it sits.
fprintf('inheritance V census: %i finite, %i -Inf, %i NaN (of %i) \n', ...
    sum(isfinite(Vx(:))),sum(Vx(:)==-Inf),sum(isnan(Vx(:))),numel(Vx));
fprintf('plain       V census: %i finite, %i -Inf, %i NaN (of %i) \n', ...
    sum(isfinite(Vp(:))),sum(Vp(:)==-Inf),sum(isnan(Vp(:))),numel(Vp));
fprintf('states where exactly one of the two is finite, this should be zero: %i \n', ...
    sum(xor(isfinite(Vx(:)),isfinite(Vp(:)))))

fprintf('CrossTest (de-risked inheritance vs plain), this should be zero: V %.3e \n',max(abs(Vx(:)-Vp(:))))
% Policy channel 1 is d1 in both models; channel 2 is d2 in one and aprime in the other, and
% they index the same grid by construction.
d1x=reshape(Policyx(1,:,:),[N_a,N_z]); d2x=reshape(Policyx(2,:,:),[N_a,N_z]);
d1p_=reshape(Policyp(1,:,:),[N_a,N_z]); apr=reshape(Policyp(2,:,:),[N_a,N_z]);
fprintf('CrossTest (de-risked inheritance vs plain), this should be zero: Policy d1 %.3e \n',max(abs(d1x(:)-d1p_(:))))
fprintf('CrossTest (de-risked inheritance vs plain), this should be zero: Policy d2 vs aprime %.3e \n',max(abs(d2x(:)-apr(:))))
fprintf('CrossTest (de-risked inheritance vs plain), this should be zero: AgentDist %.3e \n',max(abs(gather(StationaryDistx(:))-gather(StationaryDistp(:)))))
fprintf('CrossTest (de-risked inheritance vs plain), this should be zero: AggVars assets %.3e, earnings %.3e \n', ...
    abs(AggVarsx.assets.Mean-AggVarsp.assets.Mean),abs(AggVarsx.earnings.Mean-AggVarsp.earnings.Mean))

%% Plot both CDFs on top of each other
fig=figure(figure_c); %#ok<NASGU>
plot(a_grid,cumsum(sum(gather(StationaryDistx),2)),a_grid,cumsum(sum(gather(StationaryDistp),2)),'--')
title('CDF of assets: de-risked inheritance asset vs plain'); legend('inheritance','plain')

%%
output=struct();

end
