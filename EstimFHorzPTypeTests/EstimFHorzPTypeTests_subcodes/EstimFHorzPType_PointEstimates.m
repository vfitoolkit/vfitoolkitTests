function output=EstimFHorzPType_PointEstimates(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,ReturnFn,jequaloneDist)
% EstimFHorzPTypeTests part 2: point estimates of EstimateLifeCycleModel_PType_MethodOfMoments. Targets generated from the model at the truth
% (grid interpolation layer on, so the moments vary smoothly with the parameters). Started AT the truth the estimator returns the
% truth with a zero objective (exact); started AWAY it gets back within the tolerances (close to zero).
%   (a),(b) beta and pension, started at the truth, lsqnonlin (fminalgo=8) and fminsearch (fminalgo=1): zero
%   (c),(d) started away (beta 0.90, pension 0.30): recovered (close); objectivefnval vs (M_d-M_m)'W(M_d-M_m) recomputed from the
%           model moments at the estimate (informational for fminsearch, whose objective is divided by the number of parameters)
%   (e) a non-identity (diagonal) W, lsqnonlin, started away: recovered (close); objectivefnval equals (M_d-M_m)'W(M_d-M_m) at the
%       estimate (the Cholesky route of lsqnonlin), exact
%   (f) constraints, started away: constrainpositive on pension (softplus), constrain0to1 on beta, constrainAtoB on beta in [0.9,0.99] (beta started at 0.92, inside the bounds):
%       recovered (close); the Jacobian at the estimate equals the test's finite differences at the reported up/down vectors (exact:
%       J is taken on the model parameters, the constraints do not enter it)
%   (g) omitestimparam: the per-type kappa_j_pt (a matrix with N_i as its first dimension) estimated at ages 1-6 for both types (one
%       omit mask serves every type), targets the per-type age-conditional mean of earnings while working; started at 0.8 times the
%       truth at those ages: the estimated ages recovered (close), the fixed ages unchanged (exact), the result a matrix of the input
%       shape; the Jacobian has twelve columns and equals the test's finite differences (exact)
%   (h) pension per type as a structure (truth low 0.5, high 0.4) with beta, from per-type targets, started away: recovered (close), the
%       result a structure
% Each driver run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
Jr=Params.Jr;

vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);

%% Targets at the truth
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllFull=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TM=struct();
t=gather(LCPFull.assets.Mean); t(1)=NaN; t(t<1e-6)=NaN; TM.AgeConditionalStats.assets.Mean=t; % (no assets at age 1; an age with none at all omitted too)
TM.AllStats.assets.Gini=gather(AllFull.assets.Gini);
nT=sum(~isnan(t))+1;
tvec=[TM.AllStats.assets.Gini; t(~isnan(t))']; % the targets in the vector's order (AllStats first), NaN dropped
W=eye(nT);
Omega=diag(0.01*(1+(1:nT)'/nT));

EstimParamNames={'beta','pension'};
truth=[Params.beta;Params.pension];
estimoptions=struct();
estimoptions.verbose=0;

%% (a),(b) started at the truth; (c),(d) started away
Params_away=Params; Params_away.beta=0.90; Params_away.pension=0.30;
for cc=1:4
    alg=8*(mod(cc,2)==1)+1*(mod(cc,2)==0);
    if cc<=2
        Params_start=Params; lab='started at the truth';
    else
        Params_start=Params_away; lab='started at beta=0.90, pension=0.30';
    end
    estimoptions_t=estimoptions; estimoptions_t.fminalgo=alg;
    ran=0; objv=NaN; perr=NaN; tt=NaN; handobj=NaN; rec=[NaN;NaN];
    try
        t0=tic;
        [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_t,vfoptions,simoptions);
        tt=toc(t0);
        ran=1; objv=estsummary.objectivefnval; rec=[EstimParams.beta;EstimParams.pension]; perr=CoreMaxAbs(rec-truth);
        m=EstimFHorzPType_ModelMoments(rec,EstimParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
        handobj=(m-tvec)'*W*(m-tvec);
    catch err
        fprintf('PointEstimates PType (%s) %s, fminalgo=%i: the estimation command errored with: %s \n',char('a'+cc-1),lab,alg,[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    if cc<=2
        fprintf('PointEstimates PType (%s) %s, fminalgo=%i (%.1f s): did not run, this should be zero: %.3e; objective, this should be zero: %.3e; max parameter move, this should be zero: %.3e \n',char('a'+cc-1),lab,alg,tt,1-ran,objv,perr)
    else
        fprintf('PointEstimates PType (%s) %s, fminalgo=%i (%.1f s), recovered beta %.6f and pension %.6f (truth %.2f, %.2f), objective %.3e, (M_d-M_m)''W(M_d-M_m) at the estimate %.3e: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',char('a'+cc-1),lab,alg,tt,rec(1),rec(2),truth(1),truth(2),objv,handobj,1-ran,perr)
    end
end

%% (e) a diagonal weighting matrix, lsqnonlin, started away
Wd=diag((1:nT)'/nT);
estimoptions_e=estimoptions; estimoptions_e.fminalgo=8;
ran=0; objv=NaN; perr=NaN; tt=NaN; handobj=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,Wd,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_away,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_e,vfoptions,simoptions);
    tt=toc(t0);
    ran=1; objv=estsummary.objectivefnval; rec=[EstimParams.beta;EstimParams.pension]; perr=CoreMaxAbs(rec-truth);
    m=EstimFHorzPType_ModelMoments(rec,EstimParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
    handobj=(m-tvec)'*Wd*(m-tvec);
catch err
    fprintf('PointEstimates PType (e): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('PointEstimates PType (e) diagonal W, lsqnonlin, started away (%.1f s): did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e; objectivefnval vs (M_d-M_m)''W(M_d-M_m) at the estimate, this should be zero: %.3e; W reported as input, this should be zero: %.3e \n',tt,1-ran,perr,abs(objv-handobj),CoreMaxAbs(estsummary.variousmatrices.W(:)-Wd(:)))

%% (f) constraints, started away, lsqnonlin: recovered, and J unaffected by the constraints
labels={'constrainpositive pension (softplus)','constrain0to1 beta','constrainAtoB beta in [0.9,0.99]'};
for cc=1:3
    estimoptions_f=estimoptions; estimoptions_f.fminalgo=8;
    estimoptions_f.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
    if cc==1
        estimoptions_f.constrainpositive={'pension'};
    elseif cc==2
        estimoptions_f.constrain0to1={'beta'};
    else
        estimoptions_f.constrainAtoB={'beta'}; estimoptions_f.constrainAtoBlimits.beta=[0.9,0.99];
    end
    Params_start=Params_away;
    if cc==3
        Params_start.beta=0.92; % (the start must be inside [A,B]: on the bound the log-odds transform is saturated and lsqnonlin cannot move, run 2 of 2026-10-07)
    end
    ran=0; perr=NaN; tt=NaN; maxd=NaN; pointd=NaN; rec=[NaN;NaN];
    try
        t0=tic;
        [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments(EstimParamNames,TM,W,Omega,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_f,vfoptions,simoptions);
        tt=toc(t0);
        ran=1; rec=[EstimParams.beta;EstimParams.pension]; perr=CoreMaxAbs(rec-truth);
        dc=estsummary.doublechecks.epsilon10000;
        x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
        Jc=zeros(nT,2);
        for pp=1:2
            xp=x0; xp(pp)=xu(pp); mu=EstimFHorzPType_ModelMoments(xp,EstimParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
            xp=x0; xp(pp)=xd(pp); md=EstimFHorzPType_ModelMoments(xp,EstimParamNames,TM,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate,vfoptions,simoptions,0,[]);
            Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
        end
        maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
        pointd=CoreMaxAbs(x0-rec); % the Jacobian must be taken at the model-parameter estimate
    catch err
        fprintf('PointEstimates PType (f) %s: the estimation command errored with: %s \n',labels{cc},[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
    end
    fprintf('PointEstimates PType (f) %s, started away (%.1f s), recovered beta %.6f and pension %.6f: did not run, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e; Jacobian taken at a point other than the estimate, this should be zero: %.3e; J_centered vs finite differences at that point, this should be zero: %.3e \n',labels{cc},tt,rec(1),rec(2),1-ran,perr,pointd,maxd)
end

%% (g) omitestimparam: the per-type kappa_j_pt at ages 1-6 (both types)
% One omit mask serves every type, and its fixed values must equal the parameter's for every type, so the truth here has the two types
% equal at the fixed ages 7:N_j (they still differ at the estimated ages 1-6); the targets are regenerated at that truth.
FnsToEvaluate_g.earnings=FnsToEvaluate.earnings;
Params_gt=Params; Params_gt.kappa_j_pt(2,7:end)=Params.kappa_j_pt(1,7:end);
[~,Policy_gt]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_gt,DiscountFactorParamNames,vfoptions);
StationaryDist_gt=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_gt,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_gt,simoptions);
LCPFull_gt=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_gt,Policy_gt,FnsToEvaluate_g,Params_gt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMg=struct();
t=gather(LCPFull_gt.earnings.low.Mean); t(Jr:end)=NaN; TMg.AgeConditionalStats.earnings.low.Mean=t; % (earnings are zero in retirement)
t=gather(LCPFull_gt.earnings.high.Mean); t(Jr:end)=NaN; TMg.AgeConditionalStats.earnings.high.Mean=t;
nTg=2*(Jr-1);
Wg=eye(nTg); Omegag=diag(0.01*(1+(1:nTg)'/nTg));
mask=Params_gt.kappa_j_pt(1,:); mask(1:6)=NaN; % NaN = estimated; one mask serves every type
estimoptions_g=estimoptions; estimoptions_g.fminalgo=8; estimoptions_g.omitestimparam.kappa_j_pt=mask;
estimoptions_g.eedefault=3; % the checks below are at the epsilon10000 step (the default is 'auto')
Params_g=Params_gt; Params_g.kappa_j_pt(:,1:6)=0.8*Params_gt.kappa_j_pt(:,1:6);
ran=0; tt=NaN; perr=NaN; fixedmove=NaN; maxd=NaN; ncols=NaN; shapebad=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments({'kappa_j_pt'},TMg,Wg,Omegag,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_g,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate_g,estimoptions_g,vfoptions,simoptions);
    tt=toc(t0);
    ran=1;
    k=EstimParams.kappa_j_pt;
    shapebad=~isequal(size(k),size(Params_gt.kappa_j_pt)); % returned as a matrix of the input shape (N_i as the first dimension)
    if shapebad==0
        perr=max(CoreMaxAbs(k(:,1:6)-Params_gt.kappa_j_pt(:,1:6))); fixedmove=max(CoreMaxAbs(k(:,7:end)-Params_gt.kappa_j_pt(:,7:end)));
    end
    dc=estsummary.doublechecks.epsilon10000;
    x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
    ncols=size(dc.J_centered,2);
    if ncols==12 && length(x0)==12 && shapebad==0
        Jc=zeros(nTg,12);
        for pp=1:12
            xp=x0; xp(pp)=xu(pp); mu=EstimFHorzPType_ModelMoments(xp,{'kappa_j_pt'},TMg,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,EstimParamsPType_to_Params(Params_gt,EstimParams),DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate_g,vfoptions,simoptions,0,estimoptions_g.omitestimparam);
            xp=x0; xp(pp)=xd(pp); md=EstimFHorzPType_ModelMoments(xp,{'kappa_j_pt'},TMg,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,EstimParamsPType_to_Params(Params_gt,EstimParams),DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,FnsToEvaluate_g,vfoptions,simoptions,0,estimoptions_g.omitestimparam);
            Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
        end
        maxd=CoreMaxAbs(dc.J_centered(:)-Jc(:));
    else
        maxd=1;
    end
catch err
    fprintf('PointEstimates PType (g): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('PointEstimates PType (g) omitestimparam, the per-type kappa_j_pt at ages 1-6 for both types from 0.8 times the truth (%.1f s), J has %i columns (hand 12): did not run, this should be zero: %.3e; result not a matrix of the input shape, this should be zero: %.3e; max error at the estimated ages, this should be close to zero: %.3e; move at the fixed ages, this should be zero: %.3e; J_centered vs finite differences, this should be zero: %.3e \n',tt,ncols,1-ran,shapebad,perr,fixedmove,maxd)

%% (h) pension per type as a structure, with beta, from per-type targets
Params_pt=Params; Params_pt.pension=struct(); Params_pt.pension.low=0.5; Params_pt.pension.high=0.4;
[~,Policy_pt]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_pt,DiscountFactorParamNames,vfoptions);
StationaryDist_pt=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_pt,n_d,n_a,n_z,N_j,Names_i,pi_z,Params_pt,simoptions);
AllFull_pt=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_pt,Policy_pt,FnsToEvaluate,Params_pt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
LCPFull_pt=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_pt,Policy_pt,FnsToEvaluate,Params_pt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_full);
TMpt=struct();
t=gather(LCPFull_pt.assets.low.Mean); t(1)=NaN; TMpt.AgeConditionalStats.assets.low.Mean=t;
t=gather(LCPFull_pt.assets.high.Mean); t(1)=NaN; TMpt.AgeConditionalStats.assets.high.Mean=t;
TMpt.AllStats.assets.Gini=gather(AllFull_pt.assets.Gini);
nTpt=2*(N_j-1)+1;
Wpt=eye(nTpt); Omegapt=diag(0.01*(1+(1:nTpt)'/nTpt));
truth_pt=[0.95;0.5;0.4];
Params_start=Params_pt; Params_start.beta=0.90; Params_start.pension.low=0.30; Params_start.pension.high=0.30;
estimoptions_h=estimoptions; estimoptions_h.fminalgo=8;
ran=0; tt=NaN; perr=NaN; formbad=NaN; rec=nan(3,1); objv=NaN;
try
    t0=tic;
    [EstimParams,~,estsummary]=EstimateLifeCycleModel_PType_MethodOfMoments({'beta','pension'},TMpt,Wpt,Omegapt,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_start,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,[],FnsToEvaluate,estimoptions_h,vfoptions,simoptions);
    tt=toc(t0);
    ran=1; objv=estsummary.objectivefnval;
    formbad=~isstruct(EstimParams.pension);
    if formbad==0
        rec=[EstimParams.beta;EstimParams.pension.low;EstimParams.pension.high];
    end
    perr=CoreMaxAbs(rec-truth_pt);
catch err
    fprintf('PointEstimates PType (h): the estimation command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('PointEstimates PType (h) pension per type as a structure (truth low 0.5, high 0.4), with beta, started away (%.1f s), recovered %.6f, %.6f, %.6f, objective %.3e: did not run, this should be zero: %.3e; result not a structure, this should be zero: %.3e; max parameter error, this should be close to zero: %.3e \n',tt,rec(1),rec(2),rec(3),objv,1-ran,formbad,perr)

output=struct();
output.ran=ran;

end
