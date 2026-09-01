function output=TestFnsToEvaluate_d_z_noe_nosemiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Pure markov-z config -- exercises EvalFnOnAgentDist_AutoCorrTransProbs_FHorz in most
% detail (quantile TransitionProbs, deferred-error guards); e and semiz configs have their own Section J.
% Focus: smoke-test the 7 consumers on this third config, then exercise AutoCorrTransProbs_FHorz in detail.

fprintf('\n========== TestFnsToEvaluate_d_z_noe_nosemiz ==========\n')

vfoptions=struct();
simoptions=struct();

Params.J1=floor(N_j/2);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

ReturnFn=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

%% FnsToEvaluate (markov z, no e)
FnsToEvaluate.assets=@(d,aprime,a,z) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,w,kappa_j) w*kappa_j*z*d;
FnsToEvaluate.consumption=@(d,aprime,a,z,r,w,kappa_j) (1+r)*a + w*kappa_j*z*d - aprime;
FnsToEvaluate.one=@(d,aprime,a,z) 1;
FnsToEvaluate.Jnumbers=@(d,aprime,a,z,agej) agej;
FnsToEvaluate.retired=@(d,aprime,a,z,agej,Jr) (agej>=Jr);
FnsToEvaluate.atJ1=@(d,aprime,a,z,agej,J1) (agej==J1);
FnsToEvaluate.zvar=@(d,aprime,a,z) z;
FnNames=fieldnames(FnsToEvaluate);

% Counter incremented for every test that exceeds its tolerance (reported at end of subcode)
fail_count=0;
TOL_EXACT=1e-10;
TOL_TP=1e-8;     % TransitionProbs row-sum precision

%% Solve VFI + StationaryDist
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

%% Smoke-test the other 6 consumers (full cross-validation already done in the other two subcodes)
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
LifeCycle=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
SimPanelValues_FHorz_Case1(jequaloneDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf('-- All 6 of the other consumers ran without error.\n')

%% ===== AutoCorrTransProbs_FHorz: main call =====
fprintf('\n-- AutoCorrTransProbs_FHorz --\n')

simoptions_acp=simoptions;
simoptions_acp.transprobs={'assets','one','zvar'};      % request TransitionProbs for these
simoptions_acp.transprobquantiles=[];                    % default: per-age unique values
ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_acp);
fprintf('AutoCorrTransProbs_FHorz ran without error.\n')

% Shape checks
fprintf('Shape: Mean length=%d (expected %d)\n',length(ACP.assets.Mean),N_j); fail_count=fail_count+(length(ACP.assets.Mean)~=N_j);
fprintf('Shape: StdDeviation length=%d (expected %d)\n',length(ACP.assets.StdDeviation),N_j); fail_count=fail_count+(length(ACP.assets.StdDeviation)~=N_j);
fprintf('Shape: AutoCovariance length=%d (expected %d)\n',length(ACP.assets.AutoCovariance),N_j-1); fail_count=fail_count+(length(ACP.assets.AutoCovariance)~=N_j-1);
fprintf('Shape: AutoCorrelation length=%d (expected %d)\n',length(ACP.assets.AutoCorrelation),N_j-1); fail_count=fail_count+(length(ACP.assets.AutoCorrelation)~=N_j-1);

%% Consumer cross-checks
% AutoCorrTransProbs.Mean and .StdDeviation should equal LifeCycleProfiles.Mean and .StdDeviation
for ff=1:length(FnNames)
    fn=FnNames{ff};
    err_m=max(abs(gather(ACP.(fn).Mean)-gather(LifeCycle.(fn).Mean)));
    err_s=max(abs(gather(ACP.(fn).StdDeviation)-gather(LifeCycle.(fn).StdDeviation)));
    fprintf('AutoCorr Mean vs LifeCycle Mean (%-11s), should be zero: %.3e\n', fn, err_m); fail_count=fail_count+(err_m>TOL_EXACT);
    fprintf('AutoCorr StdDev vs LifeCycle StdDev (%-11s), should be zero: %.3e\n', fn, err_s); fail_count=fail_count+(err_s>TOL_EXACT);
end

% AutoCorrelation values should lie in [-1, 1] (or be NaN when StdDev=0)
for ff=1:length(FnNames)
    fn=FnNames{ff};
    ac=gather(ACP.(fn).AutoCorrelation);
    nan_mask=isnan(ac);
    in_range=all(ac(~nan_mask)>=-1-1e-10 & ac(~nan_mask)<=1+1e-10);
    fprintf('AutoCorrelation(%-11s) in [-1,1] (NaN where StdDev=0): %d\n', fn, in_range); fail_count=fail_count+(~in_range);
end

%% Analytical-truth checks
% T1 (Jnumbers): X_j = j is constant within age, so StdDev_j = 0 -> AutoCorr = NaN, AutoCov = 0 exactly
ac_J=gather(ACP.Jnumbers.AutoCovariance);
err=max(abs(ac_J));
fprintf('T1 Jnumbers AutoCovariance per-age, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
allNaN=all(isnan(gather(ACP.Jnumbers.AutoCorrelation)));
fprintf('T1 Jnumbers AutoCorrelation should be all NaN: %d\n',allNaN); fail_count=fail_count+(~allNaN);

% T2 (one): constant 1 -> AutoCov = 0 exactly, AutoCorr = NaN
ac_one=gather(ACP.one.AutoCovariance);
err=max(abs(ac_one));
fprintf('T2 one AutoCovariance per-age, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
allNaN=all(isnan(gather(ACP.one.AutoCorrelation)));
fprintf('T2 one AutoCorrelation should be all NaN: %d\n',allNaN); fail_count=fail_count+(~allNaN);

% Sanity: zvar (the markov z value) should have AutoCorrelation close to AR(1) rho (here 0.9)
% This is loose -- depends on discretization and the joint marginal under the stationary dist
ac_z=gather(ACP.zvar.AutoCorrelation);
fprintf('zvar AutoCorrelation (loose sanity, AR(0.9) -> near 0.9 on interior ages):\n')
fprintf('   median(ac_z) = %2.4f, min = %2.4f, max = %2.4f\n', median(ac_z(~isnan(ac_z))), min(ac_z(~isnan(ac_z))), max(ac_z(~isnan(ac_z))))

%% TransitionProbs checks
% one: 1 unique value per age -> cell of 1x1 matrices, each containing 1
tp_one=ACP.one.TransitionProbs;
sz1=cellfun(@(M) all(size(M)==[1,1]),tp_one);
val1=cellfun(@(M) isempty(M) || abs(M-1)<1e-10, tp_one);
fprintf('TransitionProbs(one) per-age: all 1x1 = %d, all equal 1 = %d\n', all(sz1), all(val1)); fail_count=fail_count+(~all(sz1))+(~all(val1));

% assets: smoke-test, shape sanity
tp_a=ACP.assets.TransitionProbs;
fprintf('TransitionProbs(assets): cell length %d (expected %d)\n', length(tp_a), N_j-1); fail_count=fail_count+(length(tp_a)~=N_j-1);
% Each row of each P_v should sum to 1 (with eps tolerance), where there is positive mass
worst_rowsum_err=0;
for jj=1:N_j-1
    if isempty(tp_a{jj}); continue; end
    rs=sum(tp_a{jj},2);
    err=max(abs(rs(rs>0.5)-1));
    worst_rowsum_err=max(worst_rowsum_err,err);
    if jj<=3 || jj==N_j-1
        fprintf('   age %2d: max |row-sum - 1| = %2.10f\n', jj, err)
    end
end
fprintf('TransitionProbs(assets) worst row-sum |err-1| across ages: %2.10f\n', worst_rowsum_err); fail_count=fail_count+(worst_rowsum_err>TOL_TP);

%% Quantile-binned TransitionProbs
simoptions_q=simoptions_acp;
simoptions_q.transprobs={'assets'};
simoptions_q.transprobquantiles=5;
ACP_q=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_q);
sz=size(ACP_q.assets.TransitionProbs);
fprintf('Quantile (5) TransitionProbs(assets) shape = [%d %d %d] (expected [5 5 %d])\n',sz(1),sz(2),sz(3),N_j-1); fail_count=fail_count+(sz(1)~=5)+(sz(2)~=5)+(sz(3)~=N_j-1);
% Row sums of each quintile slice should be ~1 (within mass-weighted bins)
all_ok=true;
for jj=1:N_j-1
    rs=sum(ACP_q.assets.TransitionProbs(:,:,jj),2);
    if max(abs(rs(rs>0.5)-1))>TOL_TP
        all_ok=false;
    end
end
fprintf('Quantile TransitionProbs row sums == 1 across all ages: %d\n', all_ok); fail_count=fail_count+(~all_ok);

%% Deferred-error guards
% timehorizons
simoptions_th=simoptions;
simoptions_th.timehorizons=[5,10];
try
    EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_th);
    fprintf('ERROR: timehorizons should have errored but did not\n'); fail_count=fail_count+1;
catch ME
    fprintf('Deferred-error timehorizons: %s\n', ME.message)
end

% agegroupings
simoptions_ag=simoptions;
simoptions_ag.agegroupings=1:5:N_j;
try
    EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_ag);
    fprintf('ERROR: agegroupings should have errored but did not\n'); fail_count=fail_count+1;
catch ME
    fprintf('Deferred-error agegroupings: %s\n', ME.message)
end

%% Summary for this subcode
fprintf('\nTestFnsToEvaluate_d_z_noe_nosemiz: %d tests outside tolerance.\n', fail_count)
output=struct();
output.fail_count=fail_count;
end
