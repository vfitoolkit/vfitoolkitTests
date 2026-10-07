function output=CalibFHorzPType_Extensions(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% CalibFHorzPTypeTests part 7: four things CalibrateLifeCycleModel_PType supports that parts 1-6 do not exercise. Each calibration is
% started at the truth (objective zero, parameters unmoved: exact) and, where it is a recovery, away from it (close to zero). Each
% run is wrapped so that an error in the plumbing is reported as a failed check rather than ending the bank.
%   (1) the type weights calibrated: (1a) through ParametrizeParamsFn from one share (the weights then always sum to one), started at
%       the truth and away; (1b) the weight vector itself as the calibrated parameter, started at the truth: it must stay a vector over
%       the types rather than be split into per-type scalars (it was, before 2026-10-07, and the distribution command then errored)
%   (2) per-type N_j (low 20 ages, high 15; age weights, agej and kappa_j_pt per type): beta and pension from the truth and from away,
%       with per-type profile targets of different lengths; then the same with agejshifter (the high type two ages later)
%   (3) an iid shock e with EiidShockFn and sigma_e calibrated; then both shock functions with rho_z and sigma_e calibrated together
%       (either alone must trigger the per-evaluation rebuild of the grids)
%   (4) a semi-exogenous employment state with the job-finding probability calibrated: SemiExoStateFn common to the types, and as a
%       structure over the types with the low type's probability calibrated
% The grid interpolation layer is on where assets drive the targets ((1),(2),(3)); the semi-exogenous targets are employment rates, which
% move smoothly with the transition probabilities on their own.

N_i=length(Names_i);
Jr=Params.Jr;
caliboptions=struct();
caliboptions.verbose=0;
caliboptions.fminalgo=8;
vfoptions=struct();
simoptions=struct();
vfoptions.gridinterplayer=1;
vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1;
simoptions.ngridinterp=15;
simoptions_full=simoptions;
simoptions_full.whichstats=ones(1,7);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;

%% (1) the type weights
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TM1=struct(); % grouped targets: the only ones that depend on the type weights
TM1.AllStats.earnings.Mean=gather(AllFull.earnings.Mean);
t=gather(LCPFull.earnings.Mean); t(Jr:end)=NaN; TM1.AgeConditionalStats.earnings.Mean=t;
t=gather(LCPFull.assets.Mean); t(1)=NaN; TM1.AgeConditionalStats.assets.Mean=t;
% (1a) through ParametrizeParamsFn from the share of the low type
Params_1=Params;
Params_1.sharelow=Params.ptypeweights(1);
ParametrizeWeightsFn=@(Parameters) CalibFHorzPType_ParametrizeWeightsFn(Parameters);
for ss=1:2
    if ss==2
        Params_1.sharelow=0.5;
    end
    ran=0; objv=NaN; err1=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'sharelow'},TM1,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_1,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,ParametrizeWeightsFn,FnsToEvaluate,caliboptions,vfoptions,simoptions);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; err1=abs(CalibParams.sharelow-Params.ptypeweights(1));
    catch err
        fprintf('Extensions PType (1a) type weights via ParametrizeParamsFn: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if ss==1
        fprintf('Extensions PType (1a) type weights via ParametrizeParamsFn (share of the low type), started at the truth (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; share moved, this should be zero: %.3e \n',tt,1-ran,objv,err1)
    else
        fprintf('Extensions PType (1a) type weights via ParametrizeParamsFn, share started at 0.5 (truth %.2f) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; share error, this should be close to zero: %.3e \n',Params.ptypeweights(1),tt,objv,1-ran,err1)
    end
end
% (1b) the weight vector itself, at the truth
ran=0; objv=NaN; moved=NaN; formbad=NaN;
try
    [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'ptypeweights'},TM1,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions,simoptions);
    ran=1; objv=calibsummary.objvalue;
    if isnumeric(CalibParams.ptypeweights) && numel(CalibParams.ptypeweights)==N_i
        formbad=0; moved=max(abs(CalibParams.ptypeweights(:)-Params.ptypeweights(:)));
    else
        formbad=1;
    end
catch err
    fprintf('Extensions PType (1b) the type-weight vector as the calibrated parameter: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Extensions PType (1b) the type-weight vector itself calibrated, started at the truth: did not run, this should be zero: %.3e; objective, this should be zero: %.3e; result not a vector over the types, this should be zero: %.3e; weights moved, this should be zero: %.3e \n',1-ran,objv,formbad,moved)

%% (2) per-type N_j: low 20 ages, high 15
N_j2=struct(); N_j2.low=N_j; N_j2.high=15;
Params_2=Params;
Params_2.mewj=struct(); Params_2.agej=struct(); Params_2.kappa_j_pt=struct();
for ii=1:N_i
    nj=N_j2.(Names_i{ii});
    mw=0.9.^(0:nj-1); Params_2.mewj.(Names_i{ii})=mw/sum(mw);
    Params_2.agej.(Names_i{ii})=1:1:nj;
end
Params_2.kappa_j_pt.low=Params.kappa_j;
Params_2.kappa_j_pt.high=1.5*Params.kappa_j(1:15);
ran=0; tt=NaN;
try
    [~,Policy_2]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j2,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_2,DiscountFactorParamNames,vfoptions);
    StationaryDist_2=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_2,n_d,n_a,n_z,N_j2,Names_i,pi_z,Params_2,simoptions);
    AllFull_2=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_2,Policy_2,FnsToEvaluate,Params_2,n_d,n_a,n_z,N_j2,Names_i,d_grid,a_grid,z_grid,simoptions_full);
    LCPFull_2=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_2,Policy_2,FnsToEvaluate,Params_2,n_d,n_a,n_z,N_j2,Names_i,d_grid,a_grid,z_grid,simoptions_full);
    ran=1;
catch err
    fprintf('Extensions PType (2) per-type N_j: the model solve or the stats commands errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Extensions PType (2) per-type N_j (low %i, high %i): the stats commands ran, did not run, this should be zero: %.3e \n',N_j2.low,N_j2.high,1-ran)
if ran==1
    TM2=struct();
    t=gather(LCPFull_2.assets.low.Mean); t(1)=NaN; TM2.AgeConditionalStats.assets.low.Mean=t;
    t=gather(LCPFull_2.assets.high.Mean); t(1)=NaN; TM2.AgeConditionalStats.assets.high.Mean=t;
    TM2.AllStats.assets.Gini=gather(AllFull_2.assets.Gini);
    fprintf('Extensions PType (2) per-type profile targets of lengths %i and %i \n',numel(TM2.AgeConditionalStats.assets.low.Mean),numel(TM2.AgeConditionalStats.assets.high.Mean))
    truth=[Params.beta;Params.pension];
    for ss=1:3
        Params_start=Params_2; simoptions_2=simoptions; lab='started at the truth';
        if ss==2
            Params_start.beta=0.90; Params_start.pension=0.30; lab='started at beta=0.90, pension=0.30';
        elseif ss==3
            simoptions_2.agejshifter=struct(); simoptions_2.agejshifter.low=0; simoptions_2.agejshifter.high=2; lab='started at the truth, with agejshifter (high type two ages later)';
        end
        ran=0; objv=NaN; perr=NaN; tt=NaN;
        try
            t0=tic;
            [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'beta','pension'},TM2,n_d,n_a,n_z,N_j2,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,caliboptions,vfoptions,simoptions_2);
            tt=toc(t0);
            ran=1; objv=calibsummary.objvalue; perr=max(abs([CalibParams.beta;CalibParams.pension]-truth));
        catch err
            fprintf('Extensions PType (2) per-type N_j, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ss==2
            fprintf('Extensions PType (2) per-type N_j, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
        else
            fprintf('Extensions PType (2) per-type N_j, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; parameters moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
        end
    end
end

%% (3) an iid shock e with EiidShockFn
ReturnFn_e=@(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension) ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j_pt,sigma,agej,Jr,pension);
FnsToEvaluate_e.assets=@(aprime,a,z,e) a;
FnsToEvaluate_e.earnings=@(aprime,a,z,e,w,kappa_j_pt) w*kappa_j_pt*z*e;
n_e=3;
Params_3=Params; Params_3.sigma_e=0.1; Params_3.rho_z=0.9; Params_3.sigma_z=0.03;
[e_grid,pi_e]=CalibFHorzPType_EiidShockFn(Params_3.sigma_e);
vfoptions_e=vfoptions; vfoptions_e.n_e=n_e; vfoptions_e.e_grid=e_grid; vfoptions_e.pi_e=pi_e;
simoptions_e=simoptions; simoptions_e.n_e=n_e; simoptions_e.e_grid=e_grid; simoptions_e.pi_e=pi_e;
simoptions_e_full=simoptions_e; simoptions_e_full.whichstats=ones(1,7);
jequaloneDist_e=zeros(n_a,n_z,n_e,'gpuArray'); jequaloneDist_e(1,ceil(n_z/2),ceil(n_e/2))=1;
ran=0;
try
    [~,Policy_3]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_3,DiscountFactorParamNames,vfoptions_e);
    StationaryDist_3=StationaryDist_Case1_FHorz_PType(jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,Policy_3,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_3,simoptions_e);
    AllFull_3=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_3,Policy_3,FnsToEvaluate_e,Params_3,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_e_full);
    LCPFull_3=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_3,Policy_3,FnsToEvaluate_e,Params_3,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_e_full);
    ran=1;
catch err
    fprintf('Extensions PType (3) iid shock e: the model solve or the stats commands errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Extensions PType (3) iid shock e: the model with e solved and its stats ran, did not run, this should be zero: %.3e \n',1-ran)
if ran==1
    TM3=struct();
    t=gather(LCPFull_3.earnings.low.StdDeviation); t(Jr:end)=NaN; TM3.AgeConditionalStats.earnings.low.StdDeviation=t;
    t=gather(LCPFull_3.earnings.high.StdDeviation); t(Jr:end)=NaN; TM3.AgeConditionalStats.earnings.high.StdDeviation=t;
    TM3.AllStats.earnings.Gini=gather(AllFull_3.earnings.Gini);
    vfoptions_c=vfoptions; vfoptions_c.n_e=n_e; vfoptions_c.EiidShockFn=@(sigma_e) CalibFHorzPType_EiidShockFn(sigma_e);
    simoptions_c=simoptions; simoptions_c.n_e=n_e; simoptions_c.EiidShockFn=vfoptions_c.EiidShockFn;
    for ss=1:2
        Params_start=Params_3; lab='sigma_e started at the truth';
        if ss==2
            Params_start.sigma_e=0.15; lab='sigma_e started at 0.15 (truth 0.1)';
        end
        ran=0; objv=NaN; perr=NaN; tt=NaN;
        try
            t0=tic;
            [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'sigma_e'},TM3,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_start,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_e,caliboptions,vfoptions_c,simoptions_c);
            tt=toc(t0);
            ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.sigma_e-0.1);
        catch err
            fprintf('Extensions PType (3) EiidShockFn, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ss==1
            fprintf('Extensions PType (3) EiidShockFn, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; sigma_e moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
        else
            fprintf('Extensions PType (3) EiidShockFn, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; sigma_e error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
        end
    end
    % both shock functions, both parameters calibrated
    vfoptions_c.ExogShockFn=@(rho_z,sigma_z) CalibFHorzPType_ExogShockFn(rho_z,sigma_z);
    simoptions_c.ExogShockFn=vfoptions_c.ExogShockFn;
    Params_start=Params_3; Params_start.rho_z=0.8; Params_start.sigma_e=0.15;
    caliboptions_both=caliboptions; caliboptions_both.constrain0to1={'rho_z'}; % keep the autocorrelation in (0,1), as a user would: lsqnonlin's finite-difference steps took it past one
    ran=0; objv=NaN; perr=NaN; tt=NaN;
    try
        t0=tic;
        [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'rho_z','sigma_e'},TM3,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_e,Params_start,DiscountFactorParamNames,jequaloneDist_e,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_e,caliboptions_both,vfoptions_c,simoptions_c);
        tt=toc(t0);
        ran=1; objv=calibsummary.objvalue; perr=max(abs([CalibParams.rho_z;CalibParams.sigma_e]-[0.9;0.1]));
    catch err
        fprintf('Extensions PType (3) ExogShockFn and EiidShockFn together: the calibration errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('Extensions PType (3) ExogShockFn and EiidShockFn together, rho_z and sigma_e from (0.8, 0.15) (truth 0.9, 0.1) (%.1f s), objective %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,objv,1-ran,perr)
end

%% (4) a semi-exogenous employment state, the job-finding probability calibrated
n_d_s=2; d_grid_s=[0;1]; % the search-effort decision d2 only (no d1)
n_z_s=0; z_grid_s=[]; pi_z_s=[];
Params_4=Params; Params_4.uempbenefit=0.2; Params_4.searcheffortcost=0.6; Params_4.probfindjob=0.7; Params_4.problosejob=0.3;
vfoptions_s=struct(); vfoptions_s.n_semiz=2; vfoptions_s.semiz_grid=[0;1];
vfoptions_s.SemiExoStateFn=@(n,nprime,dsemiz,probfindjob,problosejob) CoreFHorzSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob,problosejob);
simoptions_s=struct(); simoptions_s.n_semiz=vfoptions_s.n_semiz; simoptions_s.semiz_grid=vfoptions_s.semiz_grid; simoptions_s.SemiExoStateFn=vfoptions_s.SemiExoStateFn; simoptions_s.d_grid=d_grid_s;
simoptions_s_full=simoptions_s; simoptions_s_full.whichstats=ones(1,7);
ReturnFn_s=@(d2,aprime,a,semiz,r,w,kappa_j_pt,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j_pt,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
FnsToEvaluate_s.assets=@(d2,aprime,a,semiz) a;
FnsToEvaluate_s.employed=@(d2,aprime,a,semiz) semiz;
jequaloneDist_s=zeros(n_a,2,'gpuArray'); jequaloneDist_s(1,2)=1; % no assets, employed
ran=0;
try
    [~,Policy_4]=ValueFnIter_Case1_FHorz_PType(n_d_s,n_a,n_z_s,N_j,Names_i,d_grid_s,a_grid,z_grid_s,pi_z_s,ReturnFn_s,Params_4,DiscountFactorParamNames,vfoptions_s);
    StationaryDist_4=StationaryDist_Case1_FHorz_PType(jequaloneDist_s,AgeWeightParamNames,PTypeDistParamNames,Policy_4,n_d_s,n_a,n_z_s,N_j,Names_i,pi_z_s,Params_4,simoptions_s);
    AllFull_4=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_4,Policy_4,FnsToEvaluate_s,Params_4,n_d_s,n_a,n_z_s,N_j,Names_i,d_grid_s,a_grid,z_grid_s,simoptions_s_full);
    LCPFull_4=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_4,Policy_4,FnsToEvaluate_s,Params_4,n_d_s,n_a,n_z_s,N_j,Names_i,d_grid_s,a_grid,z_grid_s,simoptions_s_full);
    ran=1;
catch err
    fprintf('Extensions PType (4) semi-exogenous state: the model solve or the stats commands errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('Extensions PType (4) semi-exogenous employment state: the model solved and its stats ran, did not run, this should be zero: %.3e \n',1-ran)
if ran==1
    TM4=struct();
    t=gather(LCPFull_4.employed.low.Mean); t(Jr:end)=NaN; TM4.AgeConditionalStats.employed.low.Mean=t;
    t=gather(LCPFull_4.employed.high.Mean); t(Jr:end)=NaN; TM4.AgeConditionalStats.employed.high.Mean=t;
    TM4.AllStats.assets.Gini=gather(AllFull_4.assets.Gini);
    for ss=1:2
        Params_start=Params_4; lab='probfindjob started at the truth';
        if ss==2
            Params_start.probfindjob=0.5; lab='probfindjob started at 0.5 (truth 0.7)';
        end
        ran=0; objv=NaN; perr=NaN; tt=NaN;
        try
            t0=tic;
            [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'probfindjob'},TM4,n_d_s,n_a,n_z_s,N_j,Names_i,d_grid_s,a_grid,z_grid_s,pi_z_s,ReturnFn_s,Params_start,DiscountFactorParamNames,jequaloneDist_s,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_s,caliboptions,vfoptions_s,simoptions_s);
            tt=toc(t0);
            ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.probfindjob-0.7);
        catch err
            fprintf('Extensions PType (4) common SemiExoStateFn, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ss==1
            fprintf('Extensions PType (4) common SemiExoStateFn, %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; probfindjob moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
        else
            fprintf('Extensions PType (4) common SemiExoStateFn, %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; probfindjob error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
        end
    end
    % SemiExoStateFn as a structure over the types, the low type's job-finding probability calibrated
    Params_4b=Params_4; Params_4b.probfindjob_low=0.7; Params_4b.probfindjob_high=0.7;
    vfoptions_sb=vfoptions_s; vfoptions_sb.SemiExoStateFn=struct();
    vfoptions_sb.SemiExoStateFn.low=@(n,nprime,dsemiz,probfindjob_low,problosejob) CoreFHorzSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob_low,problosejob);
    vfoptions_sb.SemiExoStateFn.high=@(n,nprime,dsemiz,probfindjob_high,problosejob) CoreFHorzSetup_SemiExoStateFn(n,nprime,dsemiz,probfindjob_high,problosejob);
    simoptions_sb=simoptions_s; simoptions_sb.SemiExoStateFn=vfoptions_sb.SemiExoStateFn;
    for ss=1:2
        Params_start=Params_4b; lab='probfindjob_low started at the truth';
        if ss==2
            Params_start.probfindjob_low=0.5; lab='probfindjob_low started at 0.5 (truth 0.7)';
        end
        ran=0; objv=NaN; perr=NaN; tt=NaN;
        try
            t0=tic;
            [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType({'probfindjob_low'},TM4,n_d_s,n_a,n_z_s,N_j,Names_i,d_grid_s,a_grid,z_grid_s,pi_z_s,ReturnFn_s,Params_start,DiscountFactorParamNames,jequaloneDist_s,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_s,caliboptions,vfoptions_sb,simoptions_sb);
            tt=toc(t0);
            ran=1; objv=calibsummary.objvalue; perr=abs(CalibParams.probfindjob_low-0.7);
        catch err
            fprintf('Extensions PType (4) per-type SemiExoStateFn, %s: the calibration errored with: %s \n',lab,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
        end
        if ss==1
            fprintf('Extensions PType (4) per-type SemiExoStateFn (structure), %s (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; probfindjob_low moved, this should be zero: %.3e \n',lab,tt,1-ran,objv,perr)
        else
            fprintf('Extensions PType (4) per-type SemiExoStateFn (structure), %s (%.1f s), objective %.3e: did not run, this should be zero: %.3e; probfindjob_low error, this should be close to zero: %.3e \n',lab,tt,objv,1-ran,perr)
        end
    end
end

output=struct();
output.TM1=TM1;

end
