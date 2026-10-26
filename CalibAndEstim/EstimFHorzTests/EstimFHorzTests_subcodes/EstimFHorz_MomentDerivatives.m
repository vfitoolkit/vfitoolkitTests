function output=EstimFHorz_MomentDerivatives(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,ReturnFn,jequaloneDist)
% EstimFHorzTests part 7: EstimateLifeCycleModel_MomentDerivatives, the derivatives of every AllStats and LifeCycleProfiles statistic of
% every FnsToEvaluate with respect to given parameters (centered finite differences at a relative step of 1.5e-4, grid interpolation
% layer on).
%   (a) beta and pension: the reported centered finite differences equal the test's own, taken from full AllStats and LifeCycleProfiles
%       runs at the reported up/down parameter vectors and stacked in the command's order (per function: Mean, Median, RatioMeanToMedian,
%       Variance, StdDeviation, LorenzCurve, Gini, Maximum, Minimum, QuantileCutoffs, QuantileMeans; AllStats then the age profiles)
%       (exact); the main output equals the centered ones (no constraints); the named output MomentDerivatives.wrt_beta.AllStats.assets.Mean
%       etc. equals the matrix entries
%   (b) omitestimparam: kappa_j at ages 1-6: six columns, equal to the test's finite differences (exact)
%   (c) vfoptions.ExogShockFn with rho_z as the parameter: the derivative with the grids rebuilt at each step (exact)
%   (d) SortedMomentDerivatives: per parameter a list of every moment sorted by the size of its derivative (fields names and values,
%       descending in absolute value, a permutation of the column)
% Each command run is wrapped so an error in the plumbing is reported as a failed check.

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnNames=fieldnames(FnsToEvaluate);
vfoptions=struct(); simoptions=struct();
vfoptions.gridinterplayer=1; vfoptions.ngridinterp=15;
simoptions.gridinterplayer=1; simoptions.ngridinterp=15;
simoptions.npercentiles=100; simoptions.nquantiles=20;
simoptions_full=simoptions; simoptions_full.whichstats=ones(1,7);
estimoptions=struct(); estimoptions.verbose=0;
allstatnames={'Mean','Median','RatioMeanToMedian','Variance','StdDeviation','LorenzCurve','Gini','Maximum','Minimum','QuantileCutoffs','QuantileMeans'};
nmom=length(FnNames)*(8+simoptions.npercentiles+2*simoptions.nquantiles+1)*(1+N_j);

%% (a) beta and pension
EstimParamNames={'beta','pension'};
ran=0; maxd=NaN; dmain=NaN; dnamed=NaN; ncols=NaN; nrows=NaN;
try
    [MD,SMD,mds]=EstimateLifeCycleModel_MomentDerivatives(EstimParamNames,n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions,vfoptions,simoptions);
    ran=1;
    FD=mds.epsilon10000.FiniteDifference_centered;
    [nrows,ncols]=size(FD);
    dc=mds.doublechecks.epsilon10000;
    x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
    Jc=zeros(nmom,2);
    for pp=1:2
        for side=1:2
            xp=x0;
            if side==1
                xp(pp)=xu(pp);
            else
                xp(pp)=xd(pp);
            end
            P=Params; P.beta=xp(1); P.pension=xp(2);
            [~,Pol]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,P,DiscountFactorParamNames,[],vfoptions);
            Dist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Pol,n_d,n_a,n_z,N_j,pi_z,P,simoptions);
            AllS=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist,Pol,FnsToEvaluate,P,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
            LCP=LifeCycleProfiles_FHorz_Case1(Dist,Pol,FnsToEvaluate,P,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
            stack=[];
            for ff=1:length(FnNames)
                for ss=1:length(allstatnames)
                    v=gather(AllS.(FnNames{ff}).(allstatnames{ss})); stack=[stack; v(:)];
                end
            end
            for ff=1:length(FnNames)
                for ss=1:length(allstatnames)
                    v=gather(LCP.(FnNames{ff}).(allstatnames{ss})); stack=[stack; v(:)];
                end
            end
            if side==1
                mu=stack;
            else
                md=stack;
            end
        end
        Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
    end
    if isequal(size(FD),[nmom,2])
        d=abs(FD-Jc); d=d(~isnan(d)); maxd=max([d;0]);
        dmain=max(abs(mds.epsilon10000.FiniteDifference(:)-FD(:)),[],'omitnan');
        dnamed=max([abs(MD.wrt_beta.AllStats.assets.Mean-FD(1,1)); abs(MD.wrt_pension.AllStats.earnings.Gini-FD((8+simoptions.npercentiles+2*simoptions.nquantiles+1)+7,2)); ...
            abs(MD.wrt_beta.AgeConditionalStats.earnings.Gini(:)-FD(length(FnNames)*(8+simoptions.npercentiles+2*simoptions.nquantiles+1)+(8+simoptions.npercentiles+2*simoptions.nquantiles+1)*N_j+6*N_j+(1:N_j)',1))],[],'omitnan');
    else
        maxd=1; dmain=1; dnamed=1;
    end
catch err
    fprintf('MomentDerivatives (a): the command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('MomentDerivatives (a) beta and pension, %i moments by %i parameters (hand %i by 2): did not run, this should be zero: %.3e; centered finite differences vs the test''s from full AllStats and LifeCycleProfiles runs at the reported steps, this should be zero: %.3e; main output vs centered, this should be zero: %.3e; named outputs vs the matrix, this should be zero: %.3e \n',nrows,ncols,nmom,1-ran,maxd,dmain,dnamed)
% (d) the sorted output
bad=NaN;
if ran==1
    bad=0;
    for pp=1:2
        nm=['wrt_',EstimParamNames{pp}];
        if isstruct(SMD) && isfield(SMD,nm) && isfield(SMD.(nm),'values') && isfield(SMD.(nm),'names')
            vals=SMD.(nm).values(:);
            bad=bad+(numel(vals)~=nmom)+(numel(SMD.(nm).names)~=nmom)+any(diff(abs(vals))>0);
            if numel(vals)==nmom
                bad=bad+(max(abs(sort(abs(vals),'descend')-sort(abs(FD(:,pp)),'descend')),[],'omitnan')>0);
            end
        else
            bad=bad+1;
        end
    end
end
fprintf('MomentDerivatives (d) SortedMomentDerivatives per parameter: names and values of every moment, descending in absolute value, a permutation of the column; missing or wrong, this should be zero: %.3e \n',bad)

%% (b) omitestimparam: kappa_j at ages 1-6
FnsToEvaluate_b.earnings=FnsToEvaluate.earnings;
mask=Params.kappa_j; mask(1:6)=NaN;
estimoptions_b=estimoptions; estimoptions_b.omitestimparam.kappa_j=mask;
nmomb=(8+simoptions.npercentiles+2*simoptions.nquantiles+1)*(1+N_j);
ran=0; maxd=NaN; ncols=NaN;
try
    [~,~,mds]=EstimateLifeCycleModel_MomentDerivatives({'kappa_j'},n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate_b,estimoptions_b,vfoptions,simoptions);
    ran=1;
    FD=mds.epsilon10000.FiniteDifference_centered; ncols=size(FD,2);
    dc=mds.doublechecks.epsilon10000;
    x0=dc.estimparamsvec(:); xu=dc.estimparamsvecup(:); xd=dc.estimparamsvecdown(:);
    if ncols==6 && numel(x0)==6
        Jc=zeros(nmomb,6);
        for pp=1:6
            for side=1:2
                xp=x0;
                if side==1
                    xp(pp)=xu(pp);
                else
                    xp(pp)=xd(pp);
                end
                P=Params; P.kappa_j(1:6)=xp';
                [~,Pol]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,P,DiscountFactorParamNames,[],vfoptions);
                Dist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Pol,n_d,n_a,n_z,N_j,pi_z,P,simoptions);
                AllS=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist,Pol,FnsToEvaluate_b,P,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
                LCP=LifeCycleProfiles_FHorz_Case1(Dist,Pol,FnsToEvaluate_b,P,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_full);
                stack=[];
                for ss=1:length(allstatnames)
                    v=gather(AllS.earnings.(allstatnames{ss})); stack=[stack; v(:)];
                end
                for ss=1:length(allstatnames)
                    v=gather(LCP.earnings.(allstatnames{ss})); stack=[stack; v(:)];
                end
                if side==1
                    mu=stack;
                else
                    md=stack;
                end
            end
            Jc(:,pp)=(mu-md)/(xu(pp)-xd(pp));
        end
        d=abs(FD-Jc); d=d(~isnan(d)); maxd=max([d;0]);
    else
        maxd=1;
    end
catch err
    fprintf('MomentDerivatives (b): the command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('MomentDerivatives (b) omitestimparam, kappa_j at ages 1-6, %i columns (hand 6): did not run, this should be zero: %.3e; centered finite differences vs the test''s, this should be zero: %.3e \n',ncols,1-ran,maxd)

%% (c) ExogShockFn with rho_z as the parameter
Params_c=Params; Params_c.rho_z=0.9; Params_c.sigma_z=0.03;
vfoptions_c=vfoptions; vfoptions_c.ExogShockFn=@(rho_z,sigma_z) CalibFHorz_ExogShockFn(rho_z,sigma_z);
simoptions_c=simoptions; simoptions_c.ExogShockFn=vfoptions_c.ExogShockFn;
simoptions_c_full=simoptions_c; simoptions_c_full.whichstats=ones(1,7);
ran=0; maxd=NaN;
try
    [~,~,mds]=EstimateLifeCycleModel_MomentDerivatives({'rho_z'},n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_c,DiscountFactorParamNames,jequaloneDist,AgeWeightParamNames,[],FnsToEvaluate,estimoptions,vfoptions_c,simoptions_c);
    ran=1;
    FD=mds.epsilon10000.FiniteDifference_centered;
    dc=mds.doublechecks.epsilon10000;
    for side=1:2
        if side==1
            rho=dc.estimparamsvecup(1);
        else
            rho=dc.estimparamsvecdown(1);
        end
        P=Params_c; P.rho_z=rho;
        [zg,pz]=CalibFHorz_ExogShockFn(rho,Params_c.sigma_z);
        [~,Pol]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,zg,pz,ReturnFn,P,DiscountFactorParamNames,[],vfoptions);
        Dist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Pol,n_d,n_a,n_z,N_j,pz,P,simoptions);
        AllS=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist,Pol,FnsToEvaluate,P,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zg,simoptions_full);
        LCP=LifeCycleProfiles_FHorz_Case1(Dist,Pol,FnsToEvaluate,P,[],n_d,n_a,n_z,N_j,d_grid,a_grid,zg,simoptions_full);
        stack=[];
        for ff=1:length(FnNames)
            for ss=1:length(allstatnames)
                v=gather(AllS.(FnNames{ff}).(allstatnames{ss})); stack=[stack; v(:)];
            end
        end
        for ff=1:length(FnNames)
            for ss=1:length(allstatnames)
                v=gather(LCP.(FnNames{ff}).(allstatnames{ss})); stack=[stack; v(:)];
            end
        end
        if side==1
            mu=stack;
        else
            md=stack;
        end
    end
    Jc=(mu-md)/(dc.estimparamsvecup(1)-dc.estimparamsvecdown(1));
    if isequal(size(FD),[nmom,1])
        d=abs(FD-Jc); d=d(~isnan(d)); maxd=max([d;0]);
    else
        maxd=1;
    end
catch err
    fprintf('MomentDerivatives (c): the command errored with: %s \n',[err.message,' [in ',err.stack(1).name,' line ',num2str(err.stack(1).line),']'])
end
fprintf('MomentDerivatives (c) ExogShockFn, rho_z: did not run, this should be zero: %.3e; centered finite differences vs the test''s with the grids rebuilt, this should be zero: %.3e \n',1-ran,maxd)

output=struct();
output.ran=ran;

end
