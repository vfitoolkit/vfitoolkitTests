function output=CoreFHorzPType_WhichCombos(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: simoptions.whichcombos (and per-combination simoptions.whichstats where the command has whichstats) across the PType commands:
%   LifeCycleProfiles_FHorz_Case1_PType                                whichcombos [nFns, nAgeGroups, 1+nRestr], whichstats [.., 7]
%   EvalFnOnAgentDist_AllStats_FHorz_Case1_PType                       whichcombos [nFns, 1+nRestr], whichstats [nFns, 1+nRestr, 7]
%   EvalFnOnAgentDist_AggVars_FHorz_Case1_PType                        whichcombos [nFns, 1]
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType                   whichcombos [nFns, N_j, 1+nRestr] (the start-age dimension acts in the per-type command)
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType                whichcombos [nFns, nFns] (upper triangle; diagonal = own Mean/StdDev)
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType   whichcombos [nFns, nFns, nAgeGroups]
% Pages are the unrestricted output first, then the restrictions in the fieldnames order of simoptions.conditionalrestrictions.
%
% Every check compares a selective run against the full run leaf by leaf (CoreStructLeaves lists every numeric leaf with
% its dotted path): a leaf (or, for the age-indexed commands, an age-group column of a leaf) the selection asks for must be present
% and equal to the full run, NaN pattern included; one it does not ask for must be NaN where the command preallocates its fields
% (LifeCycleProfiles, AggVars, CrossSection) or absent where it assigns whole structs (AllStats, AutoCorr). RestrictedSampleMass is
% always filled. The per-type covariance matrices of the CrossSection commands shrink to the evaluated functions and are checked
% separately. Broadcast inputs equal their expansions; invalid inputs error.
% LifeCycleProfiles specifics: an all-zero function row / restriction page keeps its (NaN) fields; a combination asking for no
% statistic is skipped like a whichcombos zero and then has no struct at all; mixed Gini 3/1 across ages keeps the Lorenz curve where
% asked. AllStats specifics: StatsFromWeightedGrid reports Minimum/Maximum whenever quantiles are requested (they are the end cutoffs).

Names_i={'low','high'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

Params.ptypeweights=[0.7; 0.3];
Params.kappa_j_pt=[Params.kappa_j; 1.5*Params.kappa_j];
Params.mewj=0.9.^(0:N_j-1); Params.mewj=Params.mewj/sum(Params.mewj);

ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate.richind=@(aprime,a,z) (a>0.5);
FnNames=fieldnames(FnsToEvaluate);
nFns=length(FnNames);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();
simoptions.whichstats=ones(1,7);
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions.conditionalrestrictions.poor=@(aprime,a,z) (a<=0.5);
RestrNames=fieldnames(simoptions.conditionalrestrictions);
nRestr=length(RestrNames);
npages=1+nRestr;

%% Solve once
[~,Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,struct());

%% ===================== LifeCycleProfiles PType =====================
Full=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreStructLeaves(Full);

% (1) explicit all-ones whichcombos equals the default
simoptions_1=simoptions; simoptions_1.whichcombos=ones(nFns,N_j,npages);
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_1);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0;
for k=1:numel(pF)
    [present,idx]=ismember(pF{k},pS);
    if present
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichCombos LifeCycleProfiles PType, explicit all-ones vs default, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pF,pS))+sum(~ismember(pS,pF)))

% (2) checkerboard over (fn, age, page): every field exists; for each leaf the age-group columns the selection asks for equal the full
% run and the rest are NaN. The rule for a leaf path: a restriction name first means that page, else page 1; the function is the next
% segment; RestrictedSampleMass leaves are always whole.
wc=zeros(nFns,N_j,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:N_j
            wc(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
simoptions_2=simoptions; simoptions_2.whichcombos=wc;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_2);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            onmask=true(1,size(vF{k},2));
        else
            onmask=reshape(wc(find(strcmp(FnNames,segs{2})),:,pp)==1,[1,N_j]);
        end
    else
        onmask=reshape(wc(find(strcmp(FnNames,segs{1})),:,1)==1,[1,N_j]);
    end
    [present,idx]=ismember(pF{k},pS);
    if ~present
        missing=missing+1;
    else
        a=vF{k}; b=vS{idx};
        if ~isequal(size(a),size(b)) || size(a,2)~=numel(onmask)
            missing=missing+1;
        else
            nanmis=nanmis+sum(sum(isnan(a(:,onmask))~=isnan(b(:,onmask)))); d=abs(a(:,onmask)-b(:,onmask)); d=d(~isnan(d)); maxd=max([maxd;d(:)]);
            notnan_off=notnan_off+sum(sum(~isnan(b(:,~onmask))));
        end
    end
end
fprintf('WhichCombos LifeCycleProfiles PType, checkerboard over (fn, age, page): selected columns vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected columns not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF)))
fprintf('WhichCombos LifeCycleProfiles PType, runtime: full run vs checkerboard (half the combinations), information only: ')
tic; LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions); t1=toc;
tic; LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_2); t2=toc;
fprintf('%.3f s vs %.3f s \n',t1,t2)

% (3) an all-zero function row (earnings) and an all-zero restriction page (poor): their fields exist and are NaN (RestrictedSampleMass
% of poor still equals the full run); everything else equals the full run
wc3=ones(nFns,N_j,npages); wc3(2,:,:)=0; wc3(:,:,3)=0;
simoptions_3=simoptions; simoptions_3.whichcombos=wc3;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_3);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            on=true;
        else
            on=(wc3(find(strcmp(FnNames,segs{2})),1,pp)==1);
        end
    else
        on=(wc3(find(strcmp(FnNames,segs{1})),1,1)==1);
    end
    [present,idx]=ismember(pF{k},pS);
    if ~present
        missing=missing+1;
    elseif on
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        notnan_off=notnan_off+sum(~isnan(vS{idx}(:)));
    end
end
fprintf('WhichCombos LifeCycleProfiles PType, all-zero row (earnings) and page (poor): the other leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; zeroed leaves not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF)))

% (4) a 2D whichcombos with restrictions equals its 3D repmat
wc2=zeros(nFns,N_j);
for ff=1:nFns
    for jj=1:N_j
        wc2(ff,jj)=mod(ff+jj,2);
    end
end
simoptions_4=simoptions; simoptions_4.whichcombos=wc2;
A=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_4);
simoptions_4.whichcombos=repmat(wc2,[1,1,npages]);
B=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_4);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    if present && isequal(size(vA{k}),size(vB{idx}))
        a=vA{k}; b=vB{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos LifeCycleProfiles PType, 2D input vs its 3D repmat, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))

% (5) age groups of five ages: the second dimension is the number of age groups
simoptions_5=simoptions; simoptions_5.agegroupings=1:5:N_j; ngroups=length(simoptions_5.agegroupings);
Full5=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_5);
wc5=zeros(nFns,ngroups,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:ngroups
            wc5(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
simoptions_5.whichcombos=wc5;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_5);
[pF5,vF5]=CoreStructLeaves(Full5); [pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    ncol=size(vF5{k},2); % grouped arrays are sized by the number of ages, per-type ones by the number of age groups; only the first ngroups columns are age groups, the rest are NaN in both
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            onmask=true(1,ncol);
        else
            onmask=false(1,ncol); onmask(1:ngroups)=(wc5(find(strcmp(FnNames,segs{2})),:,pp)==1);
        end
    else
        onmask=false(1,ncol); onmask(1:ngroups)=(wc5(find(strcmp(FnNames,segs{1})),:,1)==1);
    end
    [present,idx]=ismember(pF5{k},pS);
    if ~present
        missing=missing+1;
    else
        a=vF5{k}; b=vS{idx};
        if ~isequal(size(a),size(b))
            missing=missing+1;
        else
            nanmis=nanmis+sum(sum(isnan(a(:,onmask))~=isnan(b(:,onmask)))); d=abs(a(:,onmask)-b(:,onmask)); d=d(~isnan(d)); maxd=max([maxd;d(:)]);
            notnan_off=notnan_off+sum(sum(~isnan(b(:,~onmask))));
        end
    end
end
fprintf('WhichCombos LifeCycleProfiles PType, checkerboard with age groups of five ages: selected columns vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected columns not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF5)))

% (6) without restrictions (2D whichcombos)
simoptions_6=struct(); simoptions_6.whichstats=ones(1,7);
Full6=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_6);
simoptions_6.whichcombos=wc2;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_6);
[pF6,vF6]=CoreStructLeaves(Full6); [pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF6)
    segs=strsplit(pF6{k},'.'); segs=segs(~cellfun(@isempty,segs));
    onmask=reshape(wc2(find(strcmp(FnNames,segs{1})),:)==1,[1,N_j]);
    [present,idx]=ismember(pF6{k},pS);
    if ~present
        missing=missing+1;
    else
        a=vF6{k}; b=vS{idx};
        if ~isequal(size(a),size(b)) || size(a,2)~=N_j
            missing=missing+1;
        else
            nanmis=nanmis+sum(sum(isnan(a(:,onmask))~=isnan(b(:,onmask)))); d=abs(a(:,onmask)-b(:,onmask)); d=d(~isnan(d)); maxd=max([maxd;d(:)]);
            notnan_off=notnan_off+sum(sum(~isnan(b(:,~onmask))));
        end
    end
end
fprintf('WhichCombos LifeCycleProfiles PType, without restrictions (2D): selected columns vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected columns not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF6)))

% (7) invalid inputs error
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns,N_j+1,npages);
try
    LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns,N_j); simoptions_bad.whichcombos(1,1)=2;
try
    LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=struct(); simoptions_bad.whichstats=ones(1,7); simoptions_bad.whichcombos=ones(nFns,N_j,2); % no restrictions, but a 3D input
try
    LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos LifeCycleProfiles PType, three invalid inputs (wrong age dimension, an entry of 2, 3D without restrictions) must all error, this should be zero: %.3e \n',3-errored)

%% whichstats per combination in LifeCycleProfiles PType: [nFns, nAgeGroups, 1+nRestr, 7]
% (a) a 4D all-ones equals the full run
simoptions_w=simoptions; simoptions_w.whichstats=ones(nFns,N_j,npages,7);
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0;
for k=1:numel(pF)
    [present,idx]=ismember(pF{k},pS);
    if present
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichStats LifeCycleProfiles PType, 4D all-ones vs the full run, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pF,pS))+sum(~ismember(pS,pF)))

% (b) per-combination pattern: statistic kk is on for (fn ff, age jj, page pp) iff mod(ff+jj+pp+kk,2)==1, so every statistic is on for
% about half the ages of every (fn, page) and every field exists; a leaf's columns where its statistic is on equal the full run, the
% rest are NaN. Field -> statistic: Mean 1, Median 2, RatioMeanToMedian 1 and 2, Variance/StdDeviation 3, Gini/LorenzCurve 4,
% Minimum/Maximum 5, QuantileCutoffs/QuantileMeans 6, MoreInequality.* 7.
wsarr=zeros(nFns,N_j,npages,7);
for pp=1:npages
    for ff=1:nFns
        for jj=1:N_j
            for kk=1:7
                wsarr(ff,jj,pp,kk)=mod(ff+jj+pp+kk,2);
            end
        end
    end
end
simoptions_w.whichstats=wsarr;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        isRSM=strcmp(segs{2},'RestrictedSampleMass');
        ff=find(strcmp(FnNames,segs{2}));
    else
        pp=1; isRSM=false;
        ff=find(strcmp(FnNames,segs{1}));
    end
    if isRSM
        onmask=true(1,size(vF{k},2));
    else
        if any(strcmp(segs,'MoreInequality'))
            kk=7;
        else
            switch segs{end}
                case 'Mean'
                    kk=1;
                case 'Median'
                    kk=2;
                case 'RatioMeanToMedian'
                    kk=[1,2];
                case {'Variance','StdDeviation'}
                    kk=3;
                case {'Gini','LorenzCurve'}
                    kk=4;
                case {'Minimum','Maximum'}
                    kk=5;
                case {'QuantileCutoffs','QuantileMeans'}
                    kk=6;
            end
        end
        onmask=reshape(all(wsarr(ff,:,pp,kk)==1,4),[1,N_j]);
    end
    [present,idx]=ismember(pF{k},pS);
    if ~present
        missing=missing+1;
    else
        a=vF{k}; b=vS{idx};
        if ~isequal(size(a),size(b)) || size(a,2)~=numel(onmask)
            missing=missing+1;
        else
            nanmis=nanmis+sum(sum(isnan(a(:,onmask))~=isnan(b(:,onmask)))); d=abs(a(:,onmask)-b(:,onmask)); d=d(~isnan(d)); maxd=max([maxd;d(:)]);
            notnan_off=notnan_off+sum(sum(~isnan(b(:,~onmask))));
        end
    end
end
fprintf('WhichStats LifeCycleProfiles PType, per-combination pattern: columns with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; columns with it off not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF)))

% (c) a statistic off at every age of one (fn, page) removes its field: Median off for assets unrestricted (so RatioMeanToMedian too),
% QuantileCutoffs/QuantileMeans off for rich earnings; everything else equals the full run
wsarr3=ones(nFns,N_j,npages,7); wsarr3(1,:,1,2)=0; wsarr3(2,:,2,6)=0;
simoptions_w.whichstats=wsarr3;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    removed=false;
    if strcmp(segs{1},'assets') && any(strcmp(segs{end},{'Median','RatioMeanToMedian'}))
        removed=true;
    elseif strcmp(segs{1},'rich') && numel(segs)>=2 && strcmp(segs{2},'earnings') && any(strcmp(segs{end},{'QuantileCutoffs','QuantileMeans'}))
        removed=true;
    end
    [present,idx]=ismember(pF{k},pS);
    if removed
        extra=extra+double(present);
    elseif ~present
        missing=missing+1;
    else
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichStats LifeCycleProfiles PType, a statistic off at every age: the other leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; removed fields present + other fields missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,extra+missing+sum(~ismember(pS,pF)))

% (d) [nFns, nAgeGroups, 7] with restrictions equals its 4D repmat over the pages
simoptions_w.whichstats=reshape(wsarr(:,:,1,:),[nFns,N_j,7]);
A=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
simoptions_w.whichstats=repmat(reshape(wsarr(:,:,1,:),[nFns,N_j,1,7]),[1,1,npages,1]);
B=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    if present && isequal(size(vA{k}),size(vB{idx}))
        a=vA{k}; b=vB{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichStats LifeCycleProfiles PType, 3D input vs its 4D repmat, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))

% (e) a combination with every statistic off is skipped like a whichcombos zero: richind unrestricted all off leaves no richind struct
% (or one with no statistic fields); rich.richind is unaffected
wsarr5=ones(nFns,N_j,npages,7); wsarr5(3,:,1,:)=0;
simoptions_w.whichstats=wsarr5;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
[pS,vS]=CoreStructLeaves(Sub);
nstat=sum(cellfun(@(p) strncmp(p,'.richind.',9),pS));
fprintf('WhichStats LifeCycleProfiles PType, all-off combination (richind unrestricted): leaves under richind, this should be zero: %.3e \n',nstat)
a=gather(Full.rich.richind.Mean); b=gather(Sub.rich.richind.Mean);
fprintf('WhichStats LifeCycleProfiles PType, all-off combination: rich richind grouped Mean vs the full run, this should be zero: %.3e \n',max(abs(a-b)))

% (f) mixed Gini settings across ages: whichstats(4)=3 (Gini only) at odd ages and 1 (Gini and Lorenz curve) at even ages, for assets
% unrestricted: LorenzCurve is present, NaN at the odd ages and equal to the full run at the even ages; Gini equals the full run at all ages
wsarr6=ones(nFns,N_j,npages,7);
for jj=1:2:N_j
    wsarr6(1,jj,1,4)=3;
end
simoptions_w.whichstats=wsarr6;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
fprintf('WhichStats LifeCycleProfiles PType, mixed Gini 3/1 across ages: assets LorenzCurve field present (grouped and both ptypes), this should be zero: %.3e \n',3-isfield(Sub.assets,'LorenzCurve')-isfield(Sub.assets.low,'LorenzCurve')-isfield(Sub.assets.high,'LorenzCurve'))
LC_full=gather(Full.assets.LorenzCurve); LC=gather(Sub.assets.LorenzCurve);
fprintf('WhichStats LifeCycleProfiles PType, mixed Gini 3/1 across ages: grouped assets LorenzCurve at even ages vs the full run, this should be zero: %.3e; non-NaN entries at odd ages, this should be zero: %.3e \n',max(max(abs(LC_full(:,2:2:N_j)-LC(:,2:2:N_j)))),sum(sum(~isnan(LC(:,1:2:N_j)))))
LC_full=gather(Full.assets.low.LorenzCurve); LC=gather(Sub.assets.low.LorenzCurve);
fprintf('WhichStats LifeCycleProfiles PType, mixed Gini 3/1 across ages: ptype low assets LorenzCurve at even ages vs the full run, this should be zero: %.3e; non-NaN entries at odd ages, this should be zero: %.3e \n',max(max(abs(LC_full(:,2:2:N_j)-LC(:,2:2:N_j)))),sum(sum(~isnan(LC(:,1:2:N_j)))))
fprintf('WhichStats LifeCycleProfiles PType, mixed Gini 3/1 across ages: grouped assets Gini vs the full run, this should be zero: %.3e \n',max(abs(gather(Full.assets.Gini)-gather(Sub.assets.Gini))))

% (g) wrong sizes must error
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichstats=ones(nFns,N_j+1,npages,7);
try
    LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichstats=ones(nFns,N_j,npages,6);
try
    LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichStats LifeCycleProfiles PType, two invalid sizes must both error, this should be zero: %.3e \n',2-errored)

%% ===================== AllStats PType =====================
FullAS=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreStructLeaves(FullAS);

% (a) checkerboard over (fn, page): a skipped (fn, page) has no per-type struct and no grouped fields
wcA=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wcA(ff,pp)=mod(ff+pp,2);
    end
end
simoptions_a=simoptions; simoptions_a.whichcombos=wcA;
Sub=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_a);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            expectpresent=true;
        else
            expectpresent=(wcA(find(strcmp(FnNames,segs{2})),pp)==1);
        end
    else
        expectpresent=(wcA(find(strcmp(FnNames,segs{1})),1)==1);
    end
    [present,idx]=ismember(pF{k},pS);
    if expectpresent
        if ~present
            missing=missing+1;
        else
            a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos AllStats PType, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))

% (b) whichstats per combination: statistic kk on for (fn ff, page pp) iff mod(ff+pp+kk,2)==1; a leaf is present iff its statistic is on
% (Minimum/Maximum also whenever quantiles are on: StatsFromWeightedGrid reports them with the cutoffs)
wsA=zeros(nFns,npages,7);
for pp=1:npages
    for ff=1:nFns
        for kk=1:7
            wsA(ff,pp,kk)=mod(ff+pp+kk,2);
        end
    end
end
simoptions_b=simoptions; simoptions_b.whichstats=wsA;
Sub=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_b);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        isRSM=strcmp(segs{2},'RestrictedSampleMass');
        ff=find(strcmp(FnNames,segs{2}));
    else
        pp=1; isRSM=false;
        ff=find(strcmp(FnNames,segs{1}));
    end
    if isRSM
        expectpresent=true;
    else
        statname=segs{end};
        if any(strcmp(segs,'MoreInequality'))
            kk=7;
        else
            switch statname
                case 'Mean'
                    kk=1;
                case 'Median'
                    kk=2;
                case 'RatioMeanToMedian'
                    kk=[1,2];
                case {'Variance','StdDeviation'}
                    kk=3;
                case {'Gini','LorenzCurve'}
                    kk=4;
                case {'Minimum','Maximum'}
                    kk=5;
                case {'QuantileCutoffs','QuantileMeans'}
                    kk=6;
            end
        end
        expectpresent=all(wsA(ff,pp,kk)==1);
        if any(strcmp(statname,{'Minimum','Maximum'}))
            expectpresent=(wsA(ff,pp,5)==1 || wsA(ff,pp,6)>=1);
        end
    end
    [present,idx]=ismember(pF{k},pS);
    if expectpresent
        if ~present
            missing=missing+1;
        else
            a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichStats AllStats PType, per-combination pattern: leaves with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; on-leaves missing + off-leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))

% (c) [nFns,7] whichstats equals its repmat over pages, [nFns,1] whichcombos its repmat; invalid inputs error
simoptions_c=simoptions; simoptions_c.whichstats=reshape(wsA(:,1,:),[nFns,7]); simoptions_c.whichcombos=[1;0;1];
A=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_c);
simoptions_c.whichstats=repmat(reshape(wsA(:,1,:),[nFns,1,7]),[1,npages,1]); simoptions_c.whichcombos=repmat([1;0;1],[1,npages]);
B=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_c);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    if present
        a=vA{k}; b=vB{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichCombos/WhichStats AllStats PType, broadcast inputs vs their repmats, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pA,pB))+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns+1,npages);
try
    EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichstats=ones(nFns,npages,6);
try
    EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AllStats PType, two invalid inputs must both error, this should be zero: %.3e \n',2-errored)

%% ===================== AggVars PType =====================
simoptions_ag=struct();
FullAg=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
simoptions_ag.whichcombos=[1,0,1];
SubAg=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
[pF2,vF2]=CoreStructLeaves(FullAg); [pS2,vS2]=CoreStructLeaves(SubAg);
maxd=0; notnan_off=0; missing=0;
for k=1:numel(pF2)
    segs=strsplit(pF2{k},'.'); segs=segs(~cellfun(@isempty,segs));
    ff=find(strcmp(FnNames,segs{1}));
    [present,idx]=ismember(pF2{k},pS2);
    if ~present
        missing=missing+1;
    elseif simoptions_ag.whichcombos(ff)==1
        maxd=max(maxd,abs(vF2{k}-vS2{idx}));
    else
        notnan_off=notnan_off+double(~isnan(vS2{idx}));
    end
end
fprintf('WhichCombos AggVars PType, [1,0,1]: selected Means vs full run, this should be zero: %.3e; unselected Means that are not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,notnan_off+missing+sum(~ismember(pS2,pF2)))
simoptions_ag.whichcombos=[0;0;0];
ZAg=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
[pZ,vZ]=CoreStructLeaves(ZAg);
fprintf('WhichCombos AggVars PType, all zeros: leaves present %i (full run has %i), non-NaN entries, this should be zero: %.3e \n',numel(pZ),numel(pF2),sum(cellfun(@(x) sum(~isnan(x(:))),vZ)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=[1,0];
try
    EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AggVars PType, a wrong-length whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== AutoCorrTransProbs PType =====================
simoptions_ac=struct();
simoptions_ac.timehorizons=2;
simoptions_ac.conditionalrestrictions=simoptions.conditionalrestrictions;
FullAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pF3,vF3]=CoreStructLeaves(FullAC);
% (a) (fn, page) checkerboard, every start age
wc3=zeros(nFns,N_j,npages);
for pp=1:npages
    for ff=1:nFns
        wc3(ff,:,pp)=mod(ff+pp,2);
    end
end
simoptions_ac.whichcombos=wc3;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreStructLeaves(SubAC);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF3)
    segs=strsplit(pF3{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            expectpresent=true;
        else
            expectpresent=(wc3(find(strcmp(FnNames,segs{2})),1,pp)==1);
        end
    else
        expectpresent=(wc3(find(strcmp(FnNames,segs{1})),1,1)==1);
    end
    [present,idx]=ismember(pF3{k},pS3);
    if expectpresent
        if ~present
            missing=missing+1;
        else
            a=vF3{k}; b=vS3{idx};
            if ~isequal(size(a),size(b))
                missing=missing+1;
            else
                nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
            end
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos AutoCorr PType, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS3,pF3)))
% (b) the start-age dimension acts (in the per-type command, which the wrapper passes it to): selecting only start age 1 of (assets,
% unrestricted) gives, per type and grouped, the horizon outputs starting at age 1 equal to the full run and NaN at every other start age;
% the Mean/StdDeviation are reported at every age (computed anyway); the restricted pages and the other functions are untouched
wc3b=ones(nFns,N_j,npages); wc3b(1,2:end,1)=0;
simoptions_ac.whichcombos=wc3b;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreStructLeaves(SubAC);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF3)
    segs=strsplit(pF3{k},'.'); segs=segs(~cellfun(@isempty,segs));
    [present,idx]=ismember(pF3{k},pS3);
    if ~present || ~isequal(size(vF3{k}),size(vS3{idx}))
        missing=missing+1;
        continue
    end
    a=vF3{k}; b=vS3{idx};
    if strcmp(segs{1},'assets') && ~any(strcmp(segs{end},{'Mean','StdDeviation'})) % unrestricted assets: only start age 1 selected for the horizon outputs (the age-j Mean/StdDeviation are reported at every age regardless, they are computed anyway)
        onmask=false(size(a)); onmask(1)=true;
    else
        onmask=true(size(a));
    end
    nanmis=nanmis+sum(isnan(a(onmask))~=isnan(b(onmask))); d=abs(a(onmask)-b(onmask)); d=d(~isnan(d)); maxd=max([maxd;d(:)]);
    notnan_off=notnan_off+sum(~isnan(b(~onmask)));
end
fprintf('WhichCombos AutoCorr PType, start age 1 only for unrestricted assets: selected entries vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected start ages not NaN + leaves missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS3,pF3)))
% (c) the two 2D forms expand to the same 3D selection; invalid sizes error
simoptions_ac.whichcombos=reshape(wc3(:,1,:),[nFns,npages]);
A=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
simoptions_ac.whichcombos=wc3;
B=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    if present && isequal(size(vA{k}),size(vB{idx}))
        a=vA{k}; b=vB{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos AutoCorr PType, (fn, page) input vs its expansion over start ages, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=simoptions_ac; simoptions_bad.whichcombos=ones(nFns,N_j+1,npages);
try
    EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AutoCorr PType, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== CrossSectionCovarCorr PType =====================
simoptions_cs=struct();
FullCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
% wc: assets own stats + pair (assets,richind); earnings: nothing at all; richind: own stats + pair with assets
wcs=[1,0,1; 0,0,0; 0,0,1];
wcsym=max(wcs,wcs');
fnwanted=any(wcsym,2);
simoptions_cs.whichcombos=wcs;
SubCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
[pF4,vF4]=CoreStructLeaves(FullCS); [pS4,vS4]=CoreStructLeaves(SubCS);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF4)
    segs=strsplit(pF4{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix_ptype','CorrelationMatrix_ptype','CovarianceMatrix','CorrelationMatrix'}))
        continue % the matrices are checked below
    end
    ff=find(strcmp(FnNames,segs{1}));
    if ismember(segs{2},Names_i) % per-type leaf
        if numel(segs)==3 % .fn.ptype.Mean / .StdDeviation / .fn (self): present iff the function is evaluated, equal
            expectpresent=fnwanted(ff); on=1;
        else % .fn.ptype.CovarianceWith.fn2: present iff both evaluated, equal iff the pair is selected else NaN
            ff2=find(strcmp(FnNames,segs{4}));
            expectpresent=fnwanted(ff) && fnwanted(ff2); on=wcsym(ff,ff2);
        end
    else % grouped leaf: always present; equal iff selected (the diagonal for Mean/StdDeviation/self, the pair otherwise) else NaN
        expectpresent=true;
        if numel(segs)==2
            on=wcsym(ff,ff);
        else
            ff2=find(strcmp(FnNames,segs{3}));
            on=wcsym(ff,ff2);
        end
    end
    [present,idx]=ismember(pF4{k},pS4);
    if expectpresent
        if ~present
            missing=missing+1;
        elseif on==1
            a=vF4{k}; b=vS4{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            notnan_off=notnan_off+sum(~isnan(vS4{idx}(:)));
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos CrossSection PType, pair selection: selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS4,pF4)))
Mf=gather(FullCS.CovarianceMatrix); Ms=gather(SubCS.CovarianceMatrix); Rf=gather(FullCS.CorrelationMatrix); Rs=gather(SubCS.CorrelationMatrix);
on=(max(wcsym,diag(fnwanted))==1); % pairs where selected; the diagonal (variance) wherever the function is evaluated
dd=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; dd=dd(~isnan(dd));
fprintf('WhichCombos CrossSection PType, grouped matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected entries NaN where the full run is not (a zero-variance indicator gives NaN correlations in both), this should be zero: %.3e \n',max([dd;0]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)))+sum(isnan(Rs(on))&~isnan(Rf(on))))
dmat=0; badmat=0;
sel=find(fnwanted);
for ii=1:N_i
    Mf_ii=gather(FullCS.CovarianceMatrix_ptype.(Names_i{ii})); Ms_ii=gather(SubCS.CovarianceMatrix_ptype.(Names_i{ii}));
    if ~isequal(size(Ms_ii),[numel(sel),numel(sel)])
        badmat=badmat+1;
    else
        onsub=max(wcsym(sel,sel),eye(numel(sel)))==1;
        Mf_sub=Mf_ii(sel,sel);
        dmat=max(dmat,max(abs(Mf_sub(onsub)-Ms_ii(onsub))));
        badmat=badmat+sum(~isnan(Ms_ii(~onsub)));
    end
end
fprintf('WhichCombos CrossSection PType, per-type covariance matrices (shrunk to the evaluated functions): selected entries vs full run, this should be zero: %.3e; wrong size or unselected entries not NaN, this should be zero: %.3e \n',dmat,badmat)
simoptions_cs.whichcombos=wcsym;
SymCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
[pY,vY]=CoreStructLeaves(SymCS);
maxd=0; nanmis=0;
for k=1:numel(pS4)
    [present,idx]=ismember(pS4{k},pY);
    if present && isequal(size(vS4{k}),size(vY{idx}))
        a=vS4{k}; b=vY{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos CrossSection PType, symmetric input vs upper-triangle input, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pY,pS4)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=ones(nFns,nFns+1);
try
    EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos CrossSection PType, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== AgeConditionalStats CrossSectionCovarCorr PType =====================
simoptions_acs=struct();
simoptions_acs.agegroupings=1:5:N_j;
ngroups=length(simoptions_acs.agegroupings);
FullACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_acs);
% selection: earnings nothing at all; assets/richind own stats and their pair on a checkerboard over age groups
wca=zeros(nFns,nFns,ngroups);
for kk=1:ngroups
    wca(1,1,kk)=mod(kk,2); wca(3,3,kk)=mod(kk+1,2); wca(1,3,kk)=mod(kk,2);
end
wcasym=wca;
for kk=1:ngroups
    wcasym(:,:,kk)=max(wca(:,:,kk),wca(:,:,kk)');
end
fnwanted=any(any(wcasym,2),3);
simoptions_acs.whichcombos=wca;
SubACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_acs);
[pF5,vF5]=CoreStructLeaves(FullACS); [pS5,vS5]=CoreStructLeaves(SubACS);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix_ptype','CorrelationMatrix_ptype','CovarianceMatrix','CorrelationMatrix'}))
        continue
    end
    ff=find(strcmp(FnNames,segs{1}));
    if ismember(segs{2},Names_i)
        if numel(segs)==3
            expectpresent=fnwanted(ff); onmask=reshape(any(wcasym(ff,:,:),2),[1,ngroups]); % own stats: reported in the age groups where the function is evaluated (byproducts of the pairs)
        else
            ff2=find(strcmp(FnNames,segs{4}));
            expectpresent=fnwanted(ff) && fnwanted(ff2); onmask=reshape(wcasym(ff,ff2,:)==1,[1,ngroups]);
        end
    else
        expectpresent=true;
        if numel(segs)==2
            onmask=reshape(any(wcasym(ff,:,:),2),[1,ngroups]); % grouped own stats: reported where the function is evaluated
        else
            ff2=find(strcmp(FnNames,segs{3}));
            onmask=reshape(wcasym(ff,ff2,:)==1,[1,ngroups]);
        end
    end
    [present,idx]=ismember(pF5{k},pS5);
    if expectpresent
        if ~present
            missing=missing+1;
        else
            a=vF5{k}; b=vS5{idx};
            if ~isequal(size(a),size(b)) || numel(a)~=ngroups
                missing=missing+1;
            else
                a=reshape(a,[1,ngroups]); b=reshape(b,[1,ngroups]);
                nanmis=nanmis+sum(isnan(a(onmask))~=isnan(b(onmask))); d=abs(a(onmask)-b(onmask)); d=d(~isnan(d)); maxd=max([maxd,d]);
                notnan_off=notnan_off+sum(~isnan(b(~onmask)));
            end
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos AgeConditional CrossSection PType, pair x age-group selection: selected entries vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected entries not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS5,pF5)))
Mf=gather(FullACS.CovarianceMatrix); Ms=gather(SubACS.CovarianceMatrix); Rf=gather(FullACS.CorrelationMatrix); Rs=gather(SubACS.CorrelationMatrix);
on=(wcasym==1);
for kk=1:ngroups
    on(:,:,kk)=max(on(:,:,kk),diag(any(wcasym(:,:,kk),2))); % the diagonal (variance) wherever the function is evaluated in that age group
end
dd=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; dd=dd(~isnan(dd));
fprintf('WhichCombos AgeConditional CrossSection PType, grouped matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected entries NaN where the full run is not (a zero-variance indicator gives NaN correlations in both), this should be zero: %.3e \n',max([dd;0]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)))+sum(isnan(Rs(on))&~isnan(Rf(on))))
dmat=0; badmat=0;
sel=find(fnwanted);
for ii=1:N_i
    Mf_ii=gather(FullACS.CovarianceMatrix_ptype.(Names_i{ii})); Ms_ii=gather(SubACS.CovarianceMatrix_ptype.(Names_i{ii}));
    if ~isequal(size(Ms_ii),[numel(sel),numel(sel),ngroups])
        badmat=badmat+1;
    else
        onsub=false(numel(sel),numel(sel),ngroups);
        for kk=1:ngroups
            onsub(:,:,kk)=max(wcasym(sel,sel,kk),diag(any(wcasym(sel,:,kk),2)))==1; % pairs where selected; the diagonal (variance) where the function is evaluated in that age group
        end
        Mf_sub=Mf_ii(sel,sel,:);
        badmat=badmat+sum(isnan(Ms_ii(onsub))&~isnan(Mf_sub(onsub))); % a selected entry must not be NaN unless the full run is NaN there too
        dmat=max(dmat,max(abs(Mf_sub(onsub)-Ms_ii(onsub))));
        badmat=badmat+sum(~isnan(Ms_ii(~onsub)));
    end
end
fprintf('WhichCombos AgeConditional CrossSection PType, per-type covariance matrices (shrunk to the evaluated functions): selected entries vs full run, this should be zero: %.3e; wrong size or unselected entries not NaN, this should be zero: %.3e \n',dmat,badmat)
simoptions_acs.whichcombos=wcs;
A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_acs);
simoptions_acs.whichcombos=repmat(wcs,[1,1,ngroups]);
B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_acs);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    if present && isequal(size(vA{k}),size(vB{idx}))
        a=vA{k}; b=vB{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos AgeConditional CrossSection PType, 2D input vs its repmat over age groups, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=struct(); simoptions_bad.agegroupings=1:5:N_j; simoptions_bad.whichcombos=ones(nFns,nFns,ngroups+1);
try
    EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AgeConditional CrossSection PType, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

output=struct();
output.LifeCycleProfiles_full=Full;
output.AllStats_full=FullAS;

end
