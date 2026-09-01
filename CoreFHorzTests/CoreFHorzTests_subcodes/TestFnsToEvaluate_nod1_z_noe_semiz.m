function output=TestFnsToEvaluate_nod1_z_noe_semiz(n_d,n_a,n_a_big,n_z,N_j,d_grid,a_grid,a_grid_big,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,vfoptionsbaseline,simoptionsbaseline)
% Test every FHorz FnsToEvaluate consumer + cross-validations + analytical-truth tests.
% Config: no d1 (only d2), with z, without e, with semiz

fprintf('\n========== TestFnsToEvaluate_nod1_z_noe_semiz ==========\n')

%% Setup vfoptions / simoptions (semiz)
vfoptions=struct();
vfoptions.n_semiz=vfoptionsbaseline.n_semiz;
vfoptions.semiz_grid=vfoptionsbaseline.semiz_grid;
vfoptions.SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;
simoptions=struct();
simoptions.n_semiz=simoptionsbaseline.n_semiz;
simoptions.semiz_grid=simoptionsbaseline.semiz_grid;
simoptions.SemiExoStateFn=simoptionsbaseline.SemiExoStateFn;
simoptions.d_grid=d_grid;

% J1 for the at-age-J1 indicator function
Params.J1=floor(N_j/2);

jequaloneDist=zeros(n_a,vfoptions.n_semiz,n_z,'gpuArray');
jequaloneDist(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;

ReturnFn=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);

%% FnsToEvaluate (note semiz precedes z in the state-variable order)
FnsToEvaluate.assets=@(d2,aprime,a,semiz,z) a;
FnsToEvaluate.earnings=@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*z*semiz;
FnsToEvaluate.consumption=@(d2,aprime,a,semiz,z,r,w,kappa_j) (1+r)*a + w*kappa_j*z*semiz - aprime;
FnsToEvaluate.one=@(d2,aprime,a,semiz,z) 1;
FnsToEvaluate.Jnumbers=@(d2,aprime,a,semiz,z,agej) agej;
FnsToEvaluate.retired=@(d2,aprime,a,semiz,z,agej,Jr) (agej>=Jr);
FnsToEvaluate.atJ1=@(d2,aprime,a,semiz,z,agej,J1) (agej==J1);
FnNames=fieldnames(FnsToEvaluate);

% Counter incremented for every test that exceeds its tolerance (reported at end of subcode)
fail_count=0;
TOL_EXACT=1e-10;
TOL_GINI=1e-6;
TOL_SIM=0.2;

%% Solve VFI + StationaryDist (small grid, no GI)
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,simoptions);

%% Call every consumer
AggVars=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
ValuesOnGrid=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
LifeCycle=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
CovarCorr=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
AgeCondCovarCorr=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
SimPanel=SimPanelValues_FHorz_Case1(jequaloneDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions);
fprintf('-- All 7 consumers ran without error.\n')

AgeMass=Params.mewj(:)'; % 1 x N_j row of age weights (CPU)

%% ===== Section D: Consumer cross-checks =====
fprintf('\n-- Section D: consumer cross-checks --\n')

% (25) Mean agreement
for ff=1:length(FnNames)
    fn=FnNames{ff};
    m_agg=gather(AggVars.(fn).Mean);
    m_all=gather(AllStats.(fn).Mean);
    m_lcj=gather(LifeCycle.(fn).Mean);
    m_lc=sum(m_lcj.*AgeMass);
    m_cc=gather(CovarCorr.(fn).Mean);
    err=max(abs([m_all, m_lc, m_cc]-m_agg));
    fprintf('Mean(%-11s) AggVars vs AllStats vs LifeCycle vs CovarCorr, should be zero: %.3e\n', fn, err)
    fail_count=fail_count+(err>TOL_EXACT);
end

% (26) Variance/StdDev agreement
for ff=1:length(FnNames)
    fn=FnNames{ff};
    var_all=gather(AllStats.(fn).Variance);
    sd_cc=gather(CovarCorr.(fn).StdDeviation);
    var_diag=gather(CovarCorr.CovarianceMatrix(ff,ff));
    err1=abs(var_all-sd_cc^2);
    err2=abs(var_all-var_diag);
    fprintf('Var(%-11s) AllStats vs CovarCorr.StdDev^2 vs diag(CovMat), should be zero: %.3e / %.3e\n', fn, err1, err2)
    fail_count=fail_count+(err1>TOL_EXACT)+(err2>TOL_EXACT);
end

% (27) Age-conditional Mean/StdDev: LifeCycle == AgeCondCovarCorr per age
for ff=1:length(FnNames)
    fn=FnNames{ff};
    lcj_mean=gather(LifeCycle.(fn).Mean);
    acc_mean=gather(AgeCondCovarCorr.(fn).Mean);
    lcj_sd=gather(LifeCycle.(fn).StdDeviation);
    acc_sd=gather(AgeCondCovarCorr.(fn).StdDeviation);
    err_m=max(abs(lcj_mean-acc_mean));
    err_s=max(abs(lcj_sd-acc_sd));
    fprintf('LifeCycle vs AgeCondCovarCorr Mean(%-11s)/StdDev, should be zero: %.3e / %.3e\n', fn, err_m, err_s)
    fail_count=fail_count+(err_m>TOL_EXACT)+(err_s>TOL_EXACT);
end

% (28) Law of total covariance
for ff1=1:length(FnNames)
    for ff2=ff1+1:length(FnNames)
        fn1=FnNames{ff1}; fn2=FnNames{ff2};
        cov_pool=gather(CovarCorr.(fn1).CovarianceWith.(fn2));
        cov_kj=gather(reshape(AgeCondCovarCorr.CovarianceMatrix(ff1,ff2,:),1,[]));
        m1_kj=gather(LifeCycle.(fn1).Mean);
        m2_kj=gather(LifeCycle.(fn2).Mean);
        E_cov=sum(cov_kj.*AgeMass);
        E_m1=sum(m1_kj.*AgeMass); E_m2=sum(m2_kj.*AgeMass);
        cov_m=sum((m1_kj-E_m1).*(m2_kj-E_m2).*AgeMass);
        err=abs(cov_pool-(E_cov+cov_m));
        fprintf('Law-of-total-cov %s,%s, should be zero: %.3e\n', fn1, fn2, err)
        fail_count=fail_count+(err>TOL_EXACT);
    end
end

% (31) SimPanelValues mean per age ~ LifeCycle.Mean (loose)
for ff=1:length(FnNames)
    fn=FnNames{ff};
    sim_mean=mean(gather(SimPanel.(fn)),2)';
    lc_mean=gather(LifeCycle.(fn).Mean);
    err=max(abs(sim_mean-lc_mean));
    fprintf('SimPanel mean vs LifeCycle mean (%-11s), loose: %2.6f\n', fn, err)
    fail_count=fail_count+(err>TOL_SIM);
end

% (32) ValuesOnGrid consistency
SDvec=reshape(StationaryDist,[],1);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    vogvec=reshape(ValuesOnGrid.(fn),[],1);
    pooled=gather(sum(SDvec.*vogvec));
    err=abs(pooled-gather(AggVars.(fn).Mean));
    fprintf('ValuesOnGrid pooled vs AggVars.Mean (%-11s), should be zero: %.3e\n', fn, err)
    fail_count=fail_count+(err>TOL_EXACT);
end

% (33) Trivial conditional restriction
simoptions_TR=simoptions;
simoptions_TR.conditionalrestrictions.always=@(d2,aprime,a,semiz,z) 1;
AllStats_TR=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_TR);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    err=abs(gather(AllStats.(fn).Mean)-gather(AllStats_TR.always.(fn).Mean));
    fprintf('Trivial-restriction Mean(%-11s), should be zero: %.3e\n', fn, err)
    fail_count=fail_count+(err>TOL_EXACT);
end

%% ===== Section E: Analytical-truth checks =====
fprintf('\n-- Section E: analytical-truth checks --\n')

agej_vec=1:N_j;

% T1. Jnumbers
EJ_an=sum(agej_vec.*AgeMass);
VarJ_an=sum((agej_vec-EJ_an).^2.*AgeMass);
SDJ_an=sqrt(VarJ_an);
val=agej_vec.*AgeMass; totalval=sum(val);
P_cum=cumsum(AgeMass); L_cum=cumsum(val)/totalval;
P_aug=[0,P_cum]; L_aug=[0,L_cum];
area=sum(diff(P_aug).*(L_aug(1:end-1)+L_aug(2:end))/2);
Gini_an=1-2*area;
err=abs(EJ_an-gather(AllStats.Jnumbers.Mean));
fprintf('T1 Jnumbers Mean analytical vs AllStats, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(SDJ_an-gather(AllStats.Jnumbers.StdDeviation));
fprintf('T1 Jnumbers StdDev analytical vs AllStats, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(Gini_an-gather(AllStats.Jnumbers.Gini));
fprintf('T1 Jnumbers Gini analytical vs AllStats, should be near-zero: %2.6f\n',err); fail_count=fail_count+(err>TOL_GINI);
err=max(abs(gather(LifeCycle.Jnumbers.Mean)-agej_vec));
fprintf('T1 LifeCycle.Jnumbers.Mean - agej, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=max(abs(gather(LifeCycle.Jnumbers.StdDeviation)));
fprintf('T1 LifeCycle.Jnumbers.StdDev, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
idx_J=find(strcmp(FnNames,'Jnumbers'));
diag_var=gather(reshape(AgeCondCovarCorr.CovarianceMatrix(idx_J,idx_J,:),1,[]));
err=max(abs(diag_var));
fprintf('T1 AgeCond CovMat(Jnumbers,Jnumbers,kk), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

% T2. one is constant 1
err=abs(gather(AggVars.one.Mean)-1);
fprintf('T2 AggVars.one.Mean - 1, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(gather(AllStats.one.Variance));
fprintf('T2 AllStats.one.Variance, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(gather(AllStats.one.Gini));
fprintf('T2 AllStats.one.Gini, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
idx_one=find(strcmp(FnNames,'one'));
cov_one_row=gather(reshape(CovarCorr.CovarianceMatrix(idx_one,:),1,[]));
err=max(abs(cov_one_row));
fprintf('T2 CovarCorr.CovarianceMatrix(one,:), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

% T3. retired
p_ret=sum(AgeMass(Params.Jr:N_j));
err=abs(p_ret-gather(AggVars.retired.Mean));
fprintf('T3 retired Mean analytical vs AggVars, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
err=abs(p_ret*(1-p_ret)-gather(AllStats.retired.Variance));
fprintf('T3 retired Var analytical vs AllStats, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

% T4. atJ1
err=abs(AgeMass(Params.J1)-gather(AggVars.atJ1.Mean));
fprintf('T4 atJ1 Mean - AgeMass(J1), should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

% "at age J1" conditional restriction: restricted moments == LifeCycle.Mean at J1.
% Use condJ1 as the restriction name to avoid colliding with the atJ1 FnsToEvaluate field
% (the per-fn loop would otherwise overwrite AllStats.atJ1 when ff hits the atJ1 fn).
simoptions_J1=simoptions;
simoptions_J1.conditionalrestrictions.condJ1=@(d2,aprime,a,semiz,z,agej,J1) (agej==J1);
AllStats_J1=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_J1);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    if strcmp(fn,'one') || strcmp(fn,'atJ1')
        continue
    end
    m_res=gather(AllStats_J1.condJ1.(fn).Mean);
    m_lcJ1=gather(LifeCycle.(fn).Mean(Params.J1));
    err=abs(m_res-m_lcJ1);
    fprintf('T4 ConditionalRestriction at age J1 vs LifeCycle.Mean(J1) for %-11s, should be zero: %.3e\n', fn, err)
    fail_count=fail_count+(err>TOL_EXACT);
end

% T5. Law of total covariance for (Jnumbers, assets)
m_a_kj=gather(LifeCycle.assets.Mean);
EJ=sum(agej_vec.*AgeMass); Ea=sum(m_a_kj.*AgeMass);
cov_an=sum((agej_vec-EJ).*(m_a_kj-Ea).*AgeMass);
cov_num=gather(CovarCorr.Jnumbers.CovarianceWith.assets);
err=abs(cov_an-cov_num);
fprintf('T5 Cov(Jnumbers,assets) law of total cov, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% ===== Section F: agegroupings tests =====
fprintf('\n-- Section F: agegroupings --\n')

simoptions_1bin=simoptions;
simoptions_1bin.agegroupings=[1];
AgeCond_1bin=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_1bin);
err_cov=max(abs(gather(AgeCond_1bin.CovarianceMatrix(:))-gather(CovarCorr.CovarianceMatrix(:))));
err_cor=max(abs(gather(AgeCond_1bin.CorrelationMatrix(:))-gather(CovarCorr.CorrelationMatrix(:))));
fprintf('T6 single-bin AgeCond CovMat vs pooled, should be zero: %.3e\n',err_cov); fail_count=fail_count+(err_cov>TOL_EXACT);
fprintf('T6 single-bin AgeCond CorrMat vs pooled, should be zero: %.3e\n',err_cor); fail_count=fail_count+(err_cor>TOL_EXACT);

simoptions_5bin=simoptions;
simoptions_5bin.agegroupings=1:5:N_j;
AgeCond_5bin=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_5bin);
LC_5bin=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_5bin);
n_bins=length(1:5:N_j);
fprintf('5-period agegroupings AgeCond CovMat third-dim, expected %d got %d\n',n_bins,size(AgeCond_5bin.CovarianceMatrix,3)); fail_count=fail_count+(size(AgeCond_5bin.CovarianceMatrix,3)~=n_bins);
fprintf('5-period LifeCycle.assets.Mean length, expected %d got %d\n',n_bins,length(LC_5bin.assets.Mean)); fail_count=fail_count+(length(LC_5bin.assets.Mean)~=n_bins);

%% ===== Section G: nquantiles / npoints / tolerance =====
fprintf('\n-- Section G: nquantiles / npoints / tolerance --\n')

simoptions_q=simoptions;
simoptions_q.nquantiles=4;
simoptions_q.npoints=50;
simoptions_q.tolerance=1e-9;
AllStats_q=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_q);
fprintf('nquantiles=4: QuantileCutoffs length=%d (expected 5)\n',length(AllStats_q.assets.QuantileCutoffs)); fail_count=fail_count+(length(AllStats_q.assets.QuantileCutoffs)~=5);
fprintf('npoints=50: LorenzCurve length=%d (expected 50)\n',length(AllStats_q.assets.LorenzCurve)); fail_count=fail_count+(length(AllStats_q.assets.LorenzCurve)~=50);

%% ===== Section H: outputasstructure (AggVars) =====
fprintf('\n-- Section H: outputasstructure --\n')

simoptions_oas=simoptions;
simoptions_oas.outputasstructure=1;
simoptions_oas.AggVarNames=FnNames;
AggVars_oas=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_oas);
err=abs(gather(AggVars_oas.assets.Mean)-gather(AggVars.assets.Mean));
fprintf('outputasstructure=1 AggVars.assets.Mean unchanged, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);

%% ===== Section I: Gridinterplayer (big grid) =====
fprintf('\n-- Section I: gridinterplayer --\n')

vfoptions_gi=vfoptions;
vfoptions_gi.gridinterplayer=1; vfoptions_gi.ngridinterp=5;
simoptions_gi=simoptions;
simoptions_gi.gridinterplayer=1; simoptions_gi.ngridinterp=5;

jequaloneDist_big=zeros(n_a_big,vfoptions.n_semiz,n_z,'gpuArray');
jequaloneDist_big(1,ceil(vfoptions.n_semiz/2),ceil(n_z/2))=1;

[~,Policy_big_nogi]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
SD_big_nogi=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,Policy_big_nogi,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions);
[~,Policy_big_gi]=ValueFnIter_Case1_FHorz(n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions_gi);
SD_big_gi=StationaryDist_FHorz_Case1(jequaloneDist_big,AgeWeightParamNames,Policy_big_gi,n_d,n_a_big,n_z,N_j,pi_z,Params,simoptions_gi);

AggVars_big_nogi=EvalFnOnAgentDist_AggVars_FHorz_Case1(SD_big_nogi,Policy_big_nogi,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions);
AggVars_big_gi=EvalFnOnAgentDist_AggVars_FHorz_Case1(SD_big_gi,Policy_big_gi,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions_gi);
fprintf('GI big-grid assets.Mean diff (nogi vs gi), loose: %2.6f\n',abs(gather(AggVars_big_nogi.assets.Mean)-gather(AggVars_big_gi.assets.Mean)))
fprintf('GI big-grid earnings.Mean diff (nogi vs gi), loose: %2.6f\n',abs(gather(AggVars_big_nogi.earnings.Mean)-gather(AggVars_big_gi.earnings.Mean)))

EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(SD_big_gi,Policy_big_gi,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions_gi);
EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(SD_big_gi,Policy_big_gi,FnsToEvaluate,Params,[],n_d,n_a_big,n_z,N_j,d_grid,a_grid_big,z_grid,simoptions_gi);
fprintf('GI: CovarCorr and AgeCondCovarCorr ran without error\n')

%% ===== Section J: AutoCorrTransProbs_FHorz with semiz shocks =====
fprintf('\n-- Section J: AutoCorrTransProbs_FHorz (semiz shocks) --\n')
TOL_TP=1e-8;     % TransitionProbs row-sum precision

% Add the raw shock values as extra FnsToEvaluate (used only in this section)
FnsToEvaluate_acp=FnsToEvaluate;
FnsToEvaluate_acp.zvar=@(d2,aprime,a,semiz,z) z;
FnsToEvaluate_acp.semizvar=@(d2,aprime,a,semiz,z) semiz;
ACPNames=fieldnames(FnsToEvaluate_acp);

simoptions_acp=simoptions;
simoptions_acp.transprobs={'assets'};
ACP=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate_acp,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_acp);
fprintf('AutoCorrTransProbs_FHorz ran without error.\n')

% Shape checks
fprintf('Shape: Mean length=%d (expected %d)\n',length(ACP.assets.Mean),N_j); fail_count=fail_count+(length(ACP.assets.Mean)~=N_j);
fprintf('Shape: StdDeviation length=%d (expected %d)\n',length(ACP.assets.StdDeviation),N_j); fail_count=fail_count+(length(ACP.assets.StdDeviation)~=N_j);
fprintf('Shape: AutoCovariance length=%d (expected %d)\n',length(ACP.assets.AutoCovariance),N_j-1); fail_count=fail_count+(length(ACP.assets.AutoCovariance)~=N_j-1);
fprintf('Shape: AutoCorrelation length=%d (expected %d)\n',length(ACP.assets.AutoCorrelation),N_j-1); fail_count=fail_count+(length(ACP.assets.AutoCorrelation)~=N_j-1);

% Mean and StdDeviation should equal LifeCycleProfiles (original FnsToEvaluate only)
for ff=1:length(FnNames)
    fn=FnNames{ff};
    err_m=max(abs(gather(ACP.(fn).Mean)-gather(LifeCycle.(fn).Mean)));
    err_s=max(abs(gather(ACP.(fn).StdDeviation)-gather(LifeCycle.(fn).StdDeviation)));
    fprintf('AutoCorr Mean vs LifeCycle Mean (%-11s), should be zero: %.3e\n', fn, err_m); fail_count=fail_count+(err_m>TOL_EXACT);
    fprintf('AutoCorr StdDev vs LifeCycle StdDev (%-11s), should be zero: %.3e\n', fn, err_s); fail_count=fail_count+(err_s>TOL_EXACT);
end

% AutoCorrelation values should lie in [-1, 1] (or be NaN when StdDev=0)
for ff=1:length(ACPNames)
    fn=ACPNames{ff};
    ac=gather(ACP.(fn).AutoCorrelation);
    nan_mask=isnan(ac);
    in_range=all(ac(~nan_mask)>=-1-1e-10 & ac(~nan_mask)<=1+1e-10);
    fprintf('AutoCorrelation(%-11s) in [-1,1] (NaN where StdDev=0): %d\n', fn, in_range); fail_count=fail_count+(~in_range);
end

% T1 (Jnumbers): constant within age -> AutoCov = 0 exactly, AutoCorr = NaN
err=max(abs(gather(ACP.Jnumbers.AutoCovariance)));
fprintf('T1 Jnumbers AutoCovariance per-age, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
allNaN=all(isnan(gather(ACP.Jnumbers.AutoCorrelation)));
fprintf('T1 Jnumbers AutoCorrelation should be all NaN: %d\n',allNaN); fail_count=fail_count+(~allNaN);

% T2 (one): constant 1 -> AutoCov = 0 exactly, AutoCorr = NaN
err=max(abs(gather(ACP.one.AutoCovariance)));
fprintf('T2 one AutoCovariance per-age, should be zero: %.3e\n',err); fail_count=fail_count+(err>TOL_EXACT);
allNaN=all(isnan(gather(ACP.one.AutoCorrelation)));
fprintf('T2 one AutoCorrelation should be all NaN: %d\n',allNaN); fail_count=fail_count+(~allNaN);

% Sanity: zvar (the markov z value) should have AutoCorrelation close to AR(1) rho (here 0.9)
ac_z=gather(ACP.zvar.AutoCorrelation);
fprintf('zvar AutoCorrelation (loose sanity, AR(0.9) -> near 0.9 on interior ages):\n')
fprintf('   median(ac_z) = %2.4f, min = %2.4f, max = %2.4f\n', median(ac_z(~isnan(ac_z))), min(ac_z(~isnan(ac_z))), max(ac_z(~isnan(ac_z))))

% Sanity: semizvar autocorrelation (semiz transitions depend on d2, so no analytic value; loose print)
ac_s=gather(ACP.semizvar.AutoCorrelation);
fprintf('semizvar AutoCorrelation (loose sanity, persistent semi-exo employment state):\n')
fprintf('   median(ac_s) = %2.4f, min = %2.4f, max = %2.4f\n', median(ac_s(~isnan(ac_s))), min(ac_s(~isnan(ac_s))), max(ac_s(~isnan(ac_s))))

% TransitionProbs(assets): row sums should be 1 where there is positive mass
tp_a=ACP.assets.TransitionProbs;
fprintf('TransitionProbs(assets): cell length %d (expected %d)\n', length(tp_a), N_j-1); fail_count=fail_count+(length(tp_a)~=N_j-1);
worst_rowsum_err=0;
for jj=1:N_j-1
    if isempty(tp_a{jj}); continue; end
    rs=sum(tp_a{jj},2);
    worst_rowsum_err=max(worst_rowsum_err,max(abs(rs(rs>0.5)-1)));
end
fprintf('TransitionProbs(assets) worst row-sum |err-1| across ages: %2.10f\n', worst_rowsum_err); fail_count=fail_count+(worst_rowsum_err>TOL_TP);

%% Summary for this subcode
fprintf('\nTestFnsToEvaluate_nod1_z_noe_semiz: %d tests outside tolerance.\n', fail_count)
output=struct();
output.fail_count=fail_count;
end
