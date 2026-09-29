function output=CoreFHorzPType_StatsCmds_d_z_e(n_d,n_a,n_z,n_e,N_j,d_grid,a_grid,z_grid,pi_z,e_grid,pi_e,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: every PType command downstream of StationaryDist, on the config with d, with z, with e, without semiz
% (the PType counterpart of CoreFHorzTests/TestFnsToEvaluate_d_z_e_nosemiz):
%   EvalFnOnAgentDist_AggVars_FHorz_Case1_PType
%   EvalFnOnAgentDist_AllStats_FHorz_Case1_PType
%   EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType
%   LifeCycleProfiles_FHorz_Case1_PType
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType                      (new 2026-09-29)
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType  (new 2026-09-29)
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType                         (new 2026-09-29; with timehorizons and a conditional restriction)
%   PolicyInd2Val_FHorz_PType
%   SimPanelValues_FHorz_Case1_PType   (loose only: a panel is random)
% Two genuinely different types (per-type kappa_j_pt) against two solves without PType. Per-type outputs must
% equal the solo outputs exactly (every field). Both types have the same mewj, so the grouped outputs are the
% ptweights-pooling of the solo outputs, computed here by hand: means are the weighted means, variances and
% covariances add the between-type part about the pooled means, and the pair populations of the AutoCorr
% command pool with ptweights (unrestricted) or ptweights times the pair mass (restricted). Also:
%   AutoCorr Mean/StdDeviation equal LifeCycleProfiles Mean/StdDeviation (per type and grouped),
%   AgeConditional cross-section with a single age group equals the whole-population cross-section (per type and grouped),
%   the grouped matrices equal the named grouped fields,
%   groupptypesforstats=0 drops the grouped fields of the three new commands,
%   gridinterplayer=1: per-type AutoCorr equals the solo command (same small grid, ngridinterp=5).

fprintf('\n========== CoreFHorzPType_StatsCmds_d_z_e ==========\n')

Names_i={'low','high'};
N_i=length(Names_i);
ptw=Params.(PTypeDistParamNames{1});

vfoptions=struct();
vfoptions.n_e=n_e; vfoptions.e_grid=e_grid; vfoptions.pi_e=pi_e;
simoptions=struct();
simoptions.n_e=n_e; simoptions.e_grid=e_grid; simoptions.pi_e=pi_e;

ReturnFn_PT=@(d,aprime,a,z,e,r,w,kappa_j_pt,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j_pt,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_NoPT=@(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ReturnFn_d_z_e_nosemiz(d,aprime,a,z,e,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate_PT.earnings=@(d,aprime,a,z,e,w,kappa_j_pt) w*kappa_j_pt*z*e*d;
FnsToEvaluate_PT.consumption=@(d,aprime,a,z,e,r,w,kappa_j_pt) (1+r)*a+w*kappa_j_pt*z*e*d-aprime;
FnsToEvaluate_PT.hours=@(d,aprime,a,z,e) d;
FnsToEvaluate_NoPT.assets=@(d,aprime,a,z,e) a;
FnsToEvaluate_NoPT.earnings=@(d,aprime,a,z,e,w,kappa_j) w*kappa_j*z*e*d;
FnsToEvaluate_NoPT.consumption=@(d,aprime,a,z,e,r,w,kappa_j) (1+r)*a+w*kappa_j*z*e*d-aprime;
FnsToEvaluate_NoPT.hours=@(d,aprime,a,z,e) d;
FnNames=fieldnames(FnsToEvaluate_PT);
nF=length(FnNames);

jequaloneDist=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist(1,ceil(n_z/2),ceil(n_e/2))=1;

timehorizons=[2,4];
hstr={'','_k2','_k4'};
horizons=[1,2,4];
simoptions_acp=simoptions;
simoptions_acp.timehorizons=timehorizons;
simoptions_acp.transprobs={'assets'};
simoptions_acp.conditionalrestrictions.rich=@(d,aprime,a,z,e) (a>0.5);
simoptions_ag=simoptions;
simoptions_ag.agegroupings=1:5:N_j;
simoptions_1bin=simoptions;
simoptions_1bin.agegroupings=1;

%% PType
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AggVars_PT=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
ValuesOnGrid_PT=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AgeStats_PT=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
CovarCorr_PT=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AgeCondCovarCorr_PT=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AgeCondCovarCorr_PT_ag=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
AgeCondCovarCorr_PT_1bin=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_1bin);
ACP_PT=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_acp);
PolicyVals_PT=PolicyInd2Val_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);
fprintf('-- All PType consumers ran without error.\n')

% groupptypesforstats=0
simoptions_ng=simoptions; simoptions_ng.groupptypesforstats=0;
simoptions_acp_ng=simoptions_acp; simoptions_acp_ng.groupptypesforstats=0;
CovarCorr_PT_ng=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ng);
AgeCondCovarCorr_PT_ng=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ng);
ACP_PT_ng=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_acp_ng);

% Panel
simoptions_panel=simoptions;
simoptions_panel.numbersims=10^4;
SimPanel_PT=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_panel);

% Grid interpolation layer (AutoCorr only; the other commands with GI are covered by CoreFHorzTests)
vfoptions_gi=vfoptions; vfoptions_gi.gridinterplayer=1; vfoptions_gi.ngridinterp=5;
simoptions_gi=simoptions_acp; simoptions_gi.gridinterplayer=1; simoptions_gi.ngridinterp=5;
[~,Policy_PT_gi]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions_gi);
StationaryDist_PT_gi=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT_gi,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions_gi);
ACP_PT_gi=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_PT_gi,Policy_PT_gi,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_gi);

%% Two separate single-type solves
Policy_solo=struct(); Dist_solo=struct(); AggVars_solo=struct(); AllStats_solo=struct(); ValuesOnGrid_solo=struct(); AgeStats_solo=struct();
CovarCorr_solo=struct(); AgeCondCovarCorr_solo=struct(); AgeCondCovarCorr_solo_ag=struct(); AgeCondCovarCorr_solo_1bin=struct(); ACP_solo=struct(); PolicyVals_solo=struct(); ACP_solo_gi=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [~,Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    AggVars_solo.(nm)=EvalFnOnAgentDist_AggVars_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AllStats_solo.(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    ValuesOnGrid_solo.(nm)=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeStats_solo.(nm)=LifeCycleProfiles_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    CovarCorr_solo.(nm)=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeCondCovarCorr_solo.(nm)=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeCondCovarCorr_solo_ag.(nm)=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_ag);
    AgeCondCovarCorr_solo_1bin.(nm)=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_1bin);
    ACP_solo.(nm)=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_acp);
    PolicyVals_solo.(nm)=PolicyInd2Val_FHorz(Policy_solo.(nm),n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);
    [~,Policy_gi_ii]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions_gi);
    Dist_gi_ii=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_gi_ii,n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions_gi);
    ACP_solo_gi.(nm)=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(Dist_gi_ii,Policy_gi_ii,FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_gi);
end

%% Per type: Dist, PolicyInd2Val, AggVars, ValuesOnGrid
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Stats cmds d z e, Dist (type %s), this should be zero: %.3e \n',nm,max(abs(gather(StationaryDist_PT.(nm)(:))-gather(Dist_solo.(nm)(:)))))
    fprintf('Stats cmds d z e, PolicyInd2Val (type %s), this should be zero: %.3e \n',nm,max(abs(gather(PolicyVals_PT.(nm)(:))-gather(PolicyVals_solo.(nm)(:)))))
    for ff=1:nF
        fn=FnNames{ff};
        fprintf('Stats cmds d z e, AggVars %s (type %s), this should be zero: %.3e \n',fn,nm,abs(gather(AggVars_PT.(fn).(nm).Mean)-gather(AggVars_solo.(nm).(fn).Mean)))
        fprintf('Stats cmds d z e, ValuesOnGrid %s (type %s), this should be zero: %.3e \n',fn,nm,max(abs(gather(ValuesOnGrid_PT.(fn).(nm)(:))-gather(ValuesOnGrid_solo.(nm).(fn)(:)))))
    end
end

%% Grouped, by hand from the solo outputs (same mewj for both types, so the pooling weights are ptweights alone)
for ff=1:nF
    fn=FnNames{ff};
    % AggVars and AllStats Mean
    agg_mean=ptw(1)*gather(AggVars_solo.low.(fn).Mean)+ptw(2)*gather(AggVars_solo.high.(fn).Mean);
    fprintf('Stats cmds d z e, AggVars %s grouped Mean, this should be zero: %.3e \n',fn,abs(gather(AggVars_PT.(fn).Mean)-agg_mean))
    fprintf('Stats cmds d z e, AllStats %s grouped Mean, this should be zero: %.3e \n',fn,abs(gather(AllStats_PT.(fn).Mean)-agg_mean))
    % LifeCycleProfiles grouped Mean and StdDeviation
    m1=gather(AgeStats_solo.low.(fn).Mean(:)); m2=gather(AgeStats_solo.high.(fn).Mean(:));
    s1=gather(AgeStats_solo.low.(fn).StdDeviation(:)); s2=gather(AgeStats_solo.high.(fn).StdDeviation(:));
    mG=ptw(1)*m1+ptw(2)*m2;
    sG=sqrt(ptw(1)*(s1.^2+(m1-mG).^2)+ptw(2)*(s2.^2+(m2-mG).^2));
    fprintf('Stats cmds d z e, LifeCycleProfiles %s grouped Mean, this should be zero: %.3e \n',fn,max(abs(gather(AgeStats_PT.(fn).Mean(:))-mG)))
    fprintf('Stats cmds d z e, LifeCycleProfiles %s grouped StdDeviation, this should be zero: %.3e \n',fn,max(abs(gather(AgeStats_PT.(fn).StdDeviation(:))-sG)))
    % CrossSection grouped Mean and StdDeviation
    m1=gather(CovarCorr_solo.low.(fn).Mean); m2=gather(CovarCorr_solo.high.(fn).Mean);
    s1=gather(CovarCorr_solo.low.(fn).StdDeviation); s2=gather(CovarCorr_solo.high.(fn).StdDeviation);
    mG=ptw(1)*m1+ptw(2)*m2;
    sG=sqrt(ptw(1)*(s1^2+(m1-mG)^2)+ptw(2)*(s2^2+(m2-mG)^2));
    fprintf('Stats cmds d z e, CrossSection %s grouped Mean, this should be zero: %.3e \n',fn,abs(gather(CovarCorr_PT.(fn).Mean)-mG))
    fprintf('Stats cmds d z e, CrossSection %s grouped StdDeviation, this should be zero: %.3e \n',fn,abs(gather(CovarCorr_PT.(fn).StdDeviation)-sG))
    fprintf('Stats cmds d z e, CrossSection %s grouped CovarianceMatrix diagonal vs StdDeviation^2, this should be zero: %.3e \n',fn,abs(gather(CovarCorr_PT.CovarianceMatrix(ff,ff))-gather(CovarCorr_PT.(fn).StdDeviation)^2))
    % AgeConditional CrossSection grouped Mean and StdDeviation (per age, and 5-year bins)
    m1=gather(AgeCondCovarCorr_solo.low.(fn).Mean(:)); m2=gather(AgeCondCovarCorr_solo.high.(fn).Mean(:));
    s1=gather(AgeCondCovarCorr_solo.low.(fn).StdDeviation(:)); s2=gather(AgeCondCovarCorr_solo.high.(fn).StdDeviation(:));
    mG=ptw(1)*m1+ptw(2)*m2;
    sG=sqrt(ptw(1)*(s1.^2+(m1-mG).^2)+ptw(2)*(s2.^2+(m2-mG).^2));
    fprintf('Stats cmds d z e, AgeCond CrossSection %s grouped Mean, this should be zero: %.3e \n',fn,max(abs(gather(AgeCondCovarCorr_PT.(fn).Mean(:))-mG)))
    fprintf('Stats cmds d z e, AgeCond CrossSection %s grouped StdDeviation, this should be zero: %.3e \n',fn,max(abs(gather(AgeCondCovarCorr_PT.(fn).StdDeviation(:))-sG)))
    m1=gather(AgeCondCovarCorr_solo_ag.low.(fn).Mean(:)); m2=gather(AgeCondCovarCorr_solo_ag.high.(fn).Mean(:));
    s1=gather(AgeCondCovarCorr_solo_ag.low.(fn).StdDeviation(:)); s2=gather(AgeCondCovarCorr_solo_ag.high.(fn).StdDeviation(:));
    mG=ptw(1)*m1+ptw(2)*m2;
    sG=sqrt(ptw(1)*(s1.^2+(m1-mG).^2)+ptw(2)*(s2.^2+(m2-mG).^2));
    fprintf('Stats cmds d z e, AgeCond CrossSection agegroupings %s grouped Mean, this should be zero: %.3e \n',fn,max(abs(gather(AgeCondCovarCorr_PT_ag.(fn).Mean(:))-mG)))
    fprintf('Stats cmds d z e, AgeCond CrossSection agegroupings %s grouped StdDeviation, this should be zero: %.3e \n',fn,max(abs(gather(AgeCondCovarCorr_PT_ag.(fn).StdDeviation(:))-sG)))
    % AutoCorr grouped Mean/StdDeviation equal LifeCycleProfiles grouped (and per type)
    fprintf('Stats cmds d z e, AutoCorr %s grouped Mean vs LifeCycleProfiles, this should be zero: %.3e \n',fn,max(abs(gather(ACP_PT.(fn).Mean(:))-gather(AgeStats_PT.(fn).Mean(:)))))
    fprintf('Stats cmds d z e, AutoCorr %s grouped StdDeviation vs LifeCycleProfiles, this should be zero: %.3e \n',fn,max(abs(gather(ACP_PT.(fn).StdDeviation(:))-gather(AgeStats_PT.(fn).StdDeviation(:)))))
    for ii=1:N_i
        nm=Names_i{ii};
        fprintf('Stats cmds d z e, AutoCorr %s Mean (type %s) vs LifeCycleProfiles, this should be zero: %.3e \n',fn,nm,max(abs(gather(ACP_PT.(fn).(nm).Mean(:))-gather(AgeStats_PT.(fn).(nm).Mean(:)))))
        fprintf('Stats cmds d z e, AutoCorr %s StdDeviation (type %s) vs LifeCycleProfiles, this should be zero: %.3e \n',fn,nm,max(abs(gather(ACP_PT.(fn).(nm).StdDeviation(:))-gather(AgeStats_PT.(fn).(nm).StdDeviation(:)))))
    end
    % AutoCorr grouped pairs at each horizon, unrestricted: weights ptweights
    for hh=1:length(horizons)
        kk=horizons(hh);
        m1x=gather(ACP_solo.low.(fn).Mean(1:N_j-kk)'); m2x=gather(ACP_solo.high.(fn).Mean(1:N_j-kk)');
        m1y=gather(ACP_solo.low.(fn).Mean(1+kk:N_j)'); m2y=gather(ACP_solo.high.(fn).Mean(1+kk:N_j)');
        s1x=gather(ACP_solo.low.(fn).StdDeviation(1:N_j-kk)'); s2x=gather(ACP_solo.high.(fn).StdDeviation(1:N_j-kk)');
        s1y=gather(ACP_solo.low.(fn).StdDeviation(1+kk:N_j)'); s2y=gather(ACP_solo.high.(fn).StdDeviation(1+kk:N_j)');
        c1=gather(ACP_solo.low.(fn).(['AutoCovariance',hstr{hh}])(:)); c2=gather(ACP_solo.high.(fn).(['AutoCovariance',hstr{hh}])(:));
        mx=ptw(1)*m1x+ptw(2)*m2x; my=ptw(1)*m1y+ptw(2)*m2y;
        vx=ptw(1)*(s1x.^2+(m1x-mx).^2)+ptw(2)*(s2x.^2+(m2x-mx).^2);
        vy=ptw(1)*(s1y.^2+(m1y-my).^2)+ptw(2)*(s2y.^2+(m2y-my).^2);
        cG=ptw(1)*(c1+(m1x-mx).*(m1y-my))+ptw(2)*(c2+(m2x-mx).*(m2y-my));
        rG=cG./sqrt(vx.*vy); rG(sqrt(vx.*vy)<=1e-15)=NaN;
        x=gather(ACP_PT.(fn).(['AutoCovariance',hstr{hh}])(:)); dd=abs(x-cG); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr %s grouped AutoCovariance%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(cG)))
        x=gather(ACP_PT.(fn).(['AutoCorrelation',hstr{hh}])(:)); dd=abs(x-rG); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr %s grouped AutoCorrelation%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(rG)))
        % restricted: weights ptweights times the pair mass, moments from the pair fields
        p1=ptw(1)*gather(ACP_solo.low.rich.(fn).(['PairMass',hstr{hh}])(:)); p2=ptw(2)*gather(ACP_solo.high.rich.(fn).(['PairMass',hstr{hh}])(:));
        m1x=gather(ACP_solo.low.rich.(fn).(['PairMean_j',hstr{hh}])(:)); m2x=gather(ACP_solo.high.rich.(fn).(['PairMean_j',hstr{hh}])(:));
        m1y=gather(ACP_solo.low.rich.(fn).(['PairMean_jplusk',hstr{hh}])(:)); m2y=gather(ACP_solo.high.rich.(fn).(['PairMean_jplusk',hstr{hh}])(:));
        s1x=gather(ACP_solo.low.rich.(fn).(['PairStdDeviation_j',hstr{hh}])(:)); s2x=gather(ACP_solo.high.rich.(fn).(['PairStdDeviation_j',hstr{hh}])(:));
        s1y=gather(ACP_solo.low.rich.(fn).(['PairStdDeviation_jplusk',hstr{hh}])(:)); s2y=gather(ACP_solo.high.rich.(fn).(['PairStdDeviation_jplusk',hstr{hh}])(:));
        c1=gather(ACP_solo.low.rich.(fn).(['AutoCovariance',hstr{hh}])(:)); c2=gather(ACP_solo.high.rich.(fn).(['AutoCovariance',hstr{hh}])(:));
        p1(isnan(p1))=0; p2(isnan(p2))=0;
        w1=p1./(p1+p2); w2=p2./(p1+p2); % NaN where no pairs at all
        w1(p1==0)=0; w2(p2==0)=0; % a type with no pairs contributes nothing (and its moments are NaN)
        m1x(p1==0)=0; m1y(p1==0)=0; s1x(p1==0)=0; s1y(p1==0)=0; c1(p1==0)=0;
        m2x(p2==0)=0; m2y(p2==0)=0; s2x(p2==0)=0; s2y(p2==0)=0; c2(p2==0)=0;
        mx=w1.*m1x+w2.*m2x; my=w1.*m1y+w2.*m2y;
        vx=w1.*(s1x.^2+(m1x-mx).^2)+w2.*(s2x.^2+(m2x-mx).^2);
        vy=w1.*(s1y.^2+(m1y-my).^2)+w2.*(s2y.^2+(m2y-my).^2);
        cG=w1.*(c1+(m1x-mx).*(m1y-my))+w2.*(c2+(m2x-mx).*(m2y-my));
        rG=cG./sqrt(vx.*vy); rG(sqrt(vx.*vy)<=1e-15)=NaN;
        nopairs=(p1+p2==0); mx(nopairs)=NaN; vy(nopairs)=NaN; cG(nopairs)=NaN; rG(nopairs)=NaN;
        x=gather(ACP_PT.rich.(fn).(['PairMass',hstr{hh}])(:)); y=p1+p2; dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr rich %s grouped PairMass%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
        x=gather(ACP_PT.rich.(fn).(['AutoCovariance',hstr{hh}])(:)); dd=abs(x-cG); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr rich %s grouped AutoCovariance%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(cG)))
        x=gather(ACP_PT.rich.(fn).(['AutoCorrelation',hstr{hh}])(:)); dd=abs(x-rG); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr rich %s grouped AutoCorrelation%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(rG)))
        x=gather(ACP_PT.rich.(fn).(['PairMean_j',hstr{hh}])(:)); dd=abs(x-mx); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr rich %s grouped PairMean_j%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(mx)))
        x=gather(ACP_PT.rich.(fn).(['PairStdDeviation_jplusk',hstr{hh}])(:)); y=sqrt(vy); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AutoCorr rich %s grouped PairStdDeviation_jplusk%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,hstr{hh},max([dd;0]),sum(isnan(x)~=isnan(y)))
    end
    % AutoCorr grouped restricted Mean: weights ptweights times restricted mass
    r1=ptw(1)*gather(ACP_solo.low.rich.RestrictedSampleMass(:)); r2=ptw(2)*gather(ACP_solo.high.rich.RestrictedSampleMass(:));
    m1=gather(ACP_solo.low.rich.(fn).Mean(:)); m2=gather(ACP_solo.high.rich.(fn).Mean(:));
    m1(r1==0)=0; m2(r2==0)=0;
    mG=(r1.*m1+r2.*m2)./(r1+r2);
    x=gather(ACP_PT.rich.(fn).Mean(:)); dd=abs(x-mG); dd=dd(~isnan(dd));
    fprintf('Stats cmds d z e, AutoCorr rich %s grouped Mean vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',fn,max([dd;0]),sum(isnan(x)~=isnan(mG)))
    fprintf('Stats cmds d z e, AutoCorr rich RestrictedSampleMass.ByAge vs pooled solo, this should be zero: %.3e \n',max(abs(gather(ACP_PT.rich.RestrictedSampleMass.ByAge(:))-(r1+r2))))
end
% CrossSection grouped covariances and correlations: pooled about the pooled means
for ff1=1:nF
    for ff2=ff1+1:nF
        f1=FnNames{ff1}; f2=FnNames{ff2};
        m11=gather(CovarCorr_solo.low.(f1).Mean); m21=gather(CovarCorr_solo.high.(f1).Mean);
        m12=gather(CovarCorr_solo.low.(f2).Mean); m22=gather(CovarCorr_solo.high.(f2).Mean);
        s11=gather(CovarCorr_solo.low.(f1).StdDeviation); s21=gather(CovarCorr_solo.high.(f1).StdDeviation);
        s12=gather(CovarCorr_solo.low.(f2).StdDeviation); s22=gather(CovarCorr_solo.high.(f2).StdDeviation);
        c1=gather(CovarCorr_solo.low.(f1).CovarianceWith.(f2)); c2=gather(CovarCorr_solo.high.(f1).CovarianceWith.(f2));
        mu1=ptw(1)*m11+ptw(2)*m21; mu2=ptw(1)*m12+ptw(2)*m22;
        v1=ptw(1)*(s11^2+(m11-mu1)^2)+ptw(2)*(s21^2+(m21-mu1)^2);
        v2=ptw(1)*(s12^2+(m12-mu2)^2)+ptw(2)*(s22^2+(m22-mu2)^2);
        cG=ptw(1)*(c1+(m11-mu1)*(m12-mu2))+ptw(2)*(c2+(m21-mu1)*(m22-mu2));
        fprintf('Stats cmds d z e, CrossSection grouped Covariance %s-%s vs pooled solo, this should be zero: %.3e \n',f1,f2,abs(gather(CovarCorr_PT.(f1).CovarianceWith.(f2))-cG))
        fprintf('Stats cmds d z e, CrossSection grouped Correlation %s-%s vs pooled solo, this should be zero: %.3e \n',f1,f2,abs(gather(CovarCorr_PT.(f1).CorrelationWith.(f2))-cG/sqrt(v1*v2)))
        fprintf('Stats cmds d z e, CrossSection grouped CovarianceMatrix vs named field %s-%s, this should be zero: %.3e \n',f1,f2,abs(gather(CovarCorr_PT.CovarianceMatrix(ff1,ff2))-gather(CovarCorr_PT.(f1).CovarianceWith.(f2)))+abs(gather(CovarCorr_PT.CovarianceMatrix(ff2,ff1))-gather(CovarCorr_PT.(f2).CovarianceWith.(f1))))
        fprintf('Stats cmds d z e, CrossSection grouped CorrelationMatrix vs named field %s-%s, this should be zero: %.3e \n',f1,f2,abs(gather(CovarCorr_PT.CorrelationMatrix(ff1,ff2))-gather(CovarCorr_PT.(f1).CorrelationWith.(f2))))
        % age-conditional, per age
        m11=gather(AgeCondCovarCorr_solo.low.(f1).Mean(:)); m21=gather(AgeCondCovarCorr_solo.high.(f1).Mean(:));
        m12=gather(AgeCondCovarCorr_solo.low.(f2).Mean(:)); m22=gather(AgeCondCovarCorr_solo.high.(f2).Mean(:));
        s11=gather(AgeCondCovarCorr_solo.low.(f1).StdDeviation(:)); s21=gather(AgeCondCovarCorr_solo.high.(f1).StdDeviation(:));
        s12=gather(AgeCondCovarCorr_solo.low.(f2).StdDeviation(:)); s22=gather(AgeCondCovarCorr_solo.high.(f2).StdDeviation(:));
        c1=gather(AgeCondCovarCorr_solo.low.(f1).CovarianceWith.(f2)(:)); c2=gather(AgeCondCovarCorr_solo.high.(f1).CovarianceWith.(f2)(:));
        mu1=ptw(1)*m11+ptw(2)*m21; mu2=ptw(1)*m12+ptw(2)*m22;
        v1=ptw(1)*(s11.^2+(m11-mu1).^2)+ptw(2)*(s21.^2+(m21-mu1).^2);
        v2=ptw(1)*(s12.^2+(m12-mu2).^2)+ptw(2)*(s22.^2+(m22-mu2).^2);
        cG=ptw(1)*(c1+(m11-mu1).*(m12-mu2))+ptw(2)*(c2+(m21-mu1).*(m22-mu2));
        rG=cG./sqrt(v1.*v2);
        x=gather(AgeCondCovarCorr_PT.(f1).CovarianceWith.(f2)(:)); dd=abs(x-cG); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AgeCond CrossSection grouped Covariance %s-%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',f1,f2,max([dd;0]),sum(isnan(x)~=isnan(cG)))
        x=gather(AgeCondCovarCorr_PT.(f1).CorrelationWith.(f2)(:)); dd=abs(x-rG); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AgeCond CrossSection grouped Correlation %s-%s vs pooled solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',f1,f2,max([dd;0]),sum(isnan(x)~=isnan(rG)))
        x=gather(reshape(AgeCondCovarCorr_PT.CovarianceMatrix(ff1,ff2,:),[],1)); y=gather(AgeCondCovarCorr_PT.(f1).CovarianceWith.(f2)(:)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('Stats cmds d z e, AgeCond CrossSection grouped CovarianceMatrix vs named field %s-%s, this should be zero: %.3e, NaN-pattern mismatches %i \n',f1,f2,max([dd;0]),sum(isnan(x)~=isnan(y)))
        % one age group of all ages equals the whole-population cross-section
        fprintf('Stats cmds d z e, AgeCond CrossSection (one group) grouped Covariance %s-%s vs CrossSection, this should be zero: %.3e \n',f1,f2,abs(gather(AgeCondCovarCorr_PT_1bin.(f1).CovarianceWith.(f2)(1))-gather(CovarCorr_PT.(f1).CovarianceWith.(f2))))
        for ii=1:N_i
            nm=Names_i{ii};
            fprintf('Stats cmds d z e, AgeCond CrossSection (one group) Covariance %s-%s (type %s) vs CrossSection, this should be zero: %.3e \n',f1,f2,nm,abs(gather(AgeCondCovarCorr_PT_1bin.(f1).(nm).CovarianceWith.(f2)(1))-gather(CovarCorr_PT.(f1).(nm).CovarianceWith.(f2))))
        end
    end
end
for ff=1:nF
    fn=FnNames{ff};
    fprintf('Stats cmds d z e, AgeCond CrossSection (one group) grouped Mean %s vs CrossSection, this should be zero: %.3e \n',fn,abs(gather(AgeCondCovarCorr_PT_1bin.(fn).Mean(1))-gather(CovarCorr_PT.(fn).Mean)))
    fprintf('Stats cmds d z e, AgeCond CrossSection (one group) grouped StdDeviation %s vs CrossSection, this should be zero: %.3e \n',fn,abs(gather(AgeCondCovarCorr_PT_1bin.(fn).StdDeviation(1))-gather(CovarCorr_PT.(fn).StdDeviation)))
end

%% groupptypesforstats=0 must drop the grouped fields of the three new commands
for ff=1:nF
    fn=FnNames{ff};
    fprintf('Stats cmds d z e, CrossSection %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',fn,isfield(CovarCorr_PT_ng.(fn),'Mean'))
    fprintf('Stats cmds d z e, AgeCond CrossSection %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',fn,isfield(AgeCondCovarCorr_PT_ng.(fn),'Mean'))
    fprintf('Stats cmds d z e, AutoCorr %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',fn,isfield(ACP_PT_ng.(fn),'Mean'))
end
fprintf('Stats cmds d z e, CrossSection groupptypesforstats=0 has no grouped CovarianceMatrix, this should be zero: %i \n',isfield(CovarCorr_PT_ng,'CovarianceMatrix'))

%% Panel (loose): the age-conditional panel mean against the exact grouped mean
for ff=1:nF
    fn=FnNames{ff};
    fprintf('Stats cmds d z e, SimPanel %s number of NaN entries, this should be zero: %i \n',fn,sum(isnan(SimPanel_PT.(fn)(:))))
    panelmean=mean(SimPanel_PT.(fn),2);
    fprintf('Stats cmds d z e, SimPanel %s age-conditional mean vs LifeCycleProfiles, this should be close to zero: %.3e \n',fn,max(abs(gather(panelmean(:))-gather(AgeStats_PT.(fn).Mean(:)))))
end

%% Per type, every field of each (label, Test, Ref) pair
pairs=cell(0,3);
for ii=1:N_i
    nm=Names_i{ii};
    for ff=1:nF
        fn=FnNames{ff};
        pairs(end+1,:)={sprintf('Stats cmds d z e, AllStats %s (type %s)',fn,nm), AllStats_PT.(fn).(nm), AllStats_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, LifeCycleProfiles %s (type %s)',fn,nm), AgeStats_PT.(fn).(nm), AgeStats_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, CrossSection %s (type %s)',fn,nm), CovarCorr_PT.(fn).(nm), CovarCorr_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AgeCond CrossSection %s (type %s)',fn,nm), AgeCondCovarCorr_PT.(fn).(nm), AgeCondCovarCorr_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AgeCond CrossSection agegroupings %s (type %s)',fn,nm), AgeCondCovarCorr_PT_ag.(fn).(nm), AgeCondCovarCorr_solo_ag.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AutoCorr %s (type %s)',fn,nm), ACP_PT.(fn).(nm), ACP_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AutoCorr rich %s (type %s)',fn,nm), ACP_PT.rich.(fn).(nm), ACP_solo.(nm).rich.(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AutoCorr GI %s (type %s)',fn,nm), ACP_PT_gi.(fn).(nm), ACP_solo_gi.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AutoCorr GI rich %s (type %s)',fn,nm), ACP_PT_gi.rich.(fn).(nm), ACP_solo_gi.(nm).rich.(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, AutoCorr %s (type %s) groupptypesforstats=0 vs =1',fn,nm), ACP_PT_ng.(fn).(nm), ACP_PT.(fn).(nm)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Stats cmds d z e, CrossSection %s (type %s) groupptypesforstats=0 vs =1',fn,nm), CovarCorr_PT_ng.(fn).(nm), CovarCorr_PT.(fn).(nm)}; %#ok<AGROW>
    end
    fprintf('Stats cmds d z e, CrossSection CovarianceMatrix_ptype (type %s) vs solo, this should be zero: %.3e \n',nm,max(abs(gather(CovarCorr_PT.CovarianceMatrix_ptype.(nm)(:))-gather(CovarCorr_solo.(nm).CovarianceMatrix(:)))))
    x=gather(AgeCondCovarCorr_PT.CorrelationMatrix_ptype.(nm)(:)); y=gather(AgeCondCovarCorr_solo.(nm).CorrelationMatrix(:)); dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('Stats cmds d z e, AgeCond CrossSection CorrelationMatrix_ptype (type %s) vs solo, this should be zero: %.3e, NaN-pattern mismatches %i \n',nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
    fprintf('Stats cmds d z e, AutoCorr rich RestrictedSampleMass (type %s) vs solo, this should be zero: %.3e \n',nm,max(abs(gather(ACP_PT.rich.RestrictedSampleMass.(nm)(:))-gather(ACP_solo.(nm).rich.RestrictedSampleMass(:)))))
end

% Compare every field of each (label, Test, Ref) pair. Ref supplies the list of fields; per-ptype substructures and
% the CovarianceWith/CorrelationWith substructures are compared field by field, other substructures skipped (they get
% pairs of their own); cells (TransitionProbs and its labels) are compared entry by entry. NaN is compared as a
% pattern, the count of entries that are NaN in one but not the other, because a NaN printed as the value of a check
% is not read by CoreSummary at all: the check would silently vanish from the count.
for pp=1:size(pairs,1)
    Test=pairs{pp,2};
    Ref=pairs{pp,3};
    statnames=fieldnames(Ref);
    for ss=1:length(statnames)
        sn=statnames{ss};
        if ischar(Ref.(sn)) || (isstruct(Ref.(sn)) && ~any(strcmp(sn,{'MoreInequality','CovarianceWith','CorrelationWith'})))
            continue
        end
        if isstruct(Ref.(sn))
            subnames=fieldnames(Ref.(sn));
            cmpnames=cell(length(subnames),1); xs=cmpnames; ys=cmpnames;
            for s2=1:length(subnames)
                cmpnames{s2}=[sn,'.',subnames{s2}];
                ys{s2}=Ref.(sn).(subnames{s2});
                if isfield(Test,sn) && isfield(Test.(sn),subnames{s2})
                    xs{s2}=Test.(sn).(subnames{s2});
                else
                    xs{s2}=[];
                end
            end
        else
            cmpnames={sn};
            ys={Ref.(sn)};
            if isfield(Test,sn)
                xs={Test.(sn)};
            else
                xs={[]};
            end
        end
        for s2=1:length(cmpnames)
            x=xs{s2}; y=ys{s2};
            if iscell(y)
                if ~iscell(x) || numel(x)~=numel(y)
                    fprintf('%s, %s missing or wrongly sized, this should be zero: 1 \n',pairs{pp,1},cmpnames{s2})
                else
                    dcell=0; nbad=0;
                    for cc=1:numel(y)
                        if numel(x{cc})~=numel(y{cc})
                            nbad=nbad+1;
                        elseif ~isempty(y{cc})
                            dcell=max(dcell,max(abs(gather(x{cc}(:))-gather(y{cc}(:)))));
                        end
                    end
                    fprintf('%s, %s, this should be zero: %.3e, wrongly sized cells %i \n',pairs{pp,1},cmpnames{s2},dcell,nbad)
                end
            elseif numel(x)~=numel(y)
                fprintf('%s, %s missing or wrongly sized, this should be zero: 1 \n',pairs{pp,1},cmpnames{s2})
            else
                x=gather(x(:)); y=gather(y(:));
                dd=abs(x-y); dd=dd(~isnan(dd));
                fprintf('%s, %s, this should be zero: %.3e, NaN-pattern mismatches %i \n',pairs{pp,1},cmpnames{s2},max([dd;0]),sum(isnan(x)~=isnan(y)))
            end
        end
    end
end

output=struct();

end
