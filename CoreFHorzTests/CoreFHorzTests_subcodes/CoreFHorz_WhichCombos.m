function output=CoreFHorz_WhichCombos(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames)
% Test: simoptions.whichcombos (and per-combination simoptions.whichstats where the command has whichstats) in the single-type
% FHorz stats commands:
%   LifeCycleProfiles_FHorz_Case1                                 whichcombos [nFns, nAgeGroups, 1+nRestr], whichstats [.., 7]
%   EvalFnOnAgentDist_AllStats_FHorz_Case1                        whichcombos [nFns, 1+nRestr], whichstats [nFns, 1+nRestr, 7]
%   EvalFnOnAgentDist_AggVars_FHorz_Case1                         whichcombos [nFns, 1]
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz                    whichcombos [nFns, N_j, 1+nRestr] (second dimension = start age)
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz                 whichcombos [nFns, nFns] (upper triangle; diagonal = own Mean/StdDev)
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz   whichcombos [nFns, nFns, nAgeGroups]
% Pages are the unrestricted output first, then the restrictions in the fieldnames order of simoptions.conditionalrestrictions.
%
% Every check compares a selective run against the full run leaf by leaf (CoreStructLeaves, in SharedSubcodes): a leaf, or an
% age/age-group column of a leaf, that the selection asks for must be present and equal to the full run, NaN pattern included; one
% it does not ask for must be NaN where the command preallocates its fields (LifeCycleProfiles, AggVars, AutoCorr, CrossSection) or
% absent where it assigns whole structs (AllStats) or never creates them (a function with nothing selected). RestrictedSampleMass
% is always filled. Broadcast inputs equal their expansions; invalid inputs error.
% AutoCorr: a one at (ff, j, page) asks for the age-j Mean/StdDeviation and every horizon output starting at age j; TransitionProbs
% (a cell, not seen by the leaf lister) are checked separately: empty at unselected start ages, equal at selected ones.
% AllStats: StatsFromWeightedGrid reports Minimum/Maximum whenever quantiles are requested (they are the end cutoffs).
% The PType counterpart is CoreFHorzPTypeTests part 21.

n_d=0;
d_grid=[];

ReturnFn=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
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
[~,Policy]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy,n_d,n_a,n_z,N_j,pi_z,Params,struct());

%% ===================== LifeCycleProfiles_FHorz_Case1 =====================
Full=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreStructLeaves(Full);

% (1) explicit all-ones equals the default
simoptions_1=simoptions; simoptions_1.whichcombos=ones(nFns,N_j,npages);
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_1);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0;
for k=1:numel(pF)
    [present,idx]=ismember(pF{k},pS);
    if present
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichCombos LifeCycleProfiles, explicit all-ones vs default, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pF,pS))+sum(~ismember(pS,pF)))

% (2) checkerboard over (fn, age, page): every field exists; selected age columns equal the full run, the rest are NaN
wc=zeros(nFns,N_j,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:N_j
            wc(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
simoptions_2=simoptions; simoptions_2.whichcombos=wc;
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_2);
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
fprintf('WhichCombos LifeCycleProfiles, checkerboard over (fn, age, page): selected columns vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected columns not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF)))
tic; LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions); t1=toc;
tic; LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_2); t2=toc;
fprintf('WhichCombos LifeCycleProfiles, runtime: full run %.3f s vs checkerboard (half the combinations) %.3f s, information only \n',t1,t2)

% (3) an all-zero function row (earnings) and an all-zero restriction page (poor): their fields exist and are NaN, everything else
% equals the full run (RestrictedSampleMass of poor included)
wc3=ones(nFns,N_j,npages); wc3(2,:,:)=0; wc3(:,:,3)=0;
simoptions_3=simoptions; simoptions_3.whichcombos=wc3;
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_3);
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
fprintf('WhichCombos LifeCycleProfiles, all-zero row (earnings) and page (poor): the other leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; zeroed leaves not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF)))

% (4) a 2D whichcombos with restrictions equals its 3D repmat; (5) age groups of five ages; (6) no restrictions; (7) invalid inputs
wc2=zeros(nFns,N_j);
for ff=1:nFns
    for jj=1:N_j
        wc2(ff,jj)=mod(ff+jj,2);
    end
end
simoptions_4=simoptions; simoptions_4.whichcombos=wc2;
A=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_4);
simoptions_4.whichcombos=repmat(wc2,[1,1,npages]);
B=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_4);
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
fprintf('WhichCombos LifeCycleProfiles, 2D input vs its 3D repmat, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))

simoptions_5=simoptions; simoptions_5.agegroupings=1:5:N_j; ngroups=length(simoptions_5.agegroupings);
Full5=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_5);
wc5=zeros(nFns,ngroups,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:ngroups
            wc5(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
simoptions_5.whichcombos=wc5;
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_5);
[pF5,vF5]=CoreStructLeaves(Full5); [pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; notnan_off=0; missing=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    ncol=size(vF5{k},2);
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
fprintf('WhichCombos LifeCycleProfiles, checkerboard with age groups of five ages: selected columns vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected columns not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF5)))

simoptions_6=struct(); simoptions_6.whichstats=ones(1,7);
Full6=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_6);
simoptions_6.whichcombos=wc2;
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_6);
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
fprintf('WhichCombos LifeCycleProfiles, without restrictions (2D): selected columns vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected columns not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF6)))

errored=0;
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns,N_j+1,npages);
try
    LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns,N_j); simoptions_bad.whichcombos(1,1)=2;
try
    LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos LifeCycleProfiles, two invalid inputs must both error, this should be zero: %.3e \n',2-errored)

%% whichstats per combination in LifeCycleProfiles_FHorz_Case1: [nFns, nAgeGroups, 1+nRestr, 7]
% (a) a 4D all-ones equals the full run
simoptions_w=simoptions; simoptions_w.whichstats=ones(nFns,N_j,npages,7);
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_w);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0;
for k=1:numel(pF)
    [present,idx]=ismember(pF{k},pS);
    if present
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichStats LifeCycleProfiles, 4D all-ones vs the full run, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pF,pS))+sum(~ismember(pS,pF)))

% (b) per-combination pattern: statistic kk on for (fn ff, age jj, page pp) iff mod(ff+jj+pp+kk,2)==1; every field exists, and a leaf's
% columns where its statistic is on equal the full run, the rest are NaN
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
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_w);
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
fprintf('WhichStats LifeCycleProfiles, per-combination pattern: columns with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; columns with it off not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS,pF)))

% (c) a statistic off at every age of one (fn, page) removes its field; a combination with every statistic off is skipped like a
% whichcombos zero (richind unrestricted: no richind struct); (d) [nFns, nAgeGroups, 7] equals its 4D repmat; (e) invalid sizes error
wsarr3=ones(nFns,N_j,npages,7); wsarr3(1,:,1,2)=0; wsarr3(2,:,2,6)=0; wsarr3(3,:,1,:)=0;
simoptions_w.whichstats=wsarr3;
Sub=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_w);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    removed=false;
    if strcmp(segs{1},'assets') && any(strcmp(segs{end},{'Median','RatioMeanToMedian'}))
        removed=true;
    elseif strcmp(segs{1},'rich') && numel(segs)>=2 && strcmp(segs{2},'earnings') && any(strcmp(segs{end},{'QuantileCutoffs','QuantileMeans'}))
        removed=true;
    elseif strcmp(segs{1},'richind')
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
fprintf('WhichStats LifeCycleProfiles, statistics off at every age and an all-off combination: the other leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; removed fields present + other fields missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,extra+missing+sum(~ismember(pS,pF)))
simoptions_w.whichstats=reshape(wsarr(:,:,1,:),[nFns,N_j,7]);
A=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_w);
simoptions_w.whichstats=repmat(reshape(wsarr(:,:,1,:),[nFns,N_j,1,7]),[1,1,npages,1]);
B=LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_w);
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
fprintf('WhichStats LifeCycleProfiles, 3D input vs its 4D repmat, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichstats=ones(nFns,N_j+1,npages,7);
try
    LifeCycleProfiles_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichStats LifeCycleProfiles, an invalid size must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_AllStats_FHorz_Case1 =====================
FullAS=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreStructLeaves(FullAS);
wcA=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wcA(ff,pp)=mod(ff+pp,2);
    end
end
simoptions_a=simoptions; simoptions_a.whichcombos=wcA;
Sub=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_a);
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
fprintf('WhichCombos AllStats, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))
wsA=zeros(nFns,npages,7);
for pp=1:npages
    for ff=1:nFns
        for kk=1:7
            wsA(ff,pp,kk)=mod(ff+pp+kk,2);
        end
    end
end
simoptions_b=simoptions; simoptions_b.whichstats=wsA;
Sub=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_b);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0; badpaths={};
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
            badpaths{end+1}=['missing: ',pF{k}]; %#ok<AGROW>
        else
            a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
        badpaths{end+1}=['present: ',pF{k}]; %#ok<AGROW>
    end
end
fprintf('WhichStats AllStats, per-combination pattern: leaves with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; on-leaves missing + off-leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))
if ~isempty(badpaths)
    fprintf('WhichStats AllStats, per-combination pattern: offending leaves (information): %s \n',strjoin(badpaths(1:min(8,end)),' | '))
end
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns+1,npages);
try
    EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AllStats, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_AggVars_FHorz_Case1 =====================
simoptions_ag=struct();
FullAg=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_ag);
simoptions_ag.whichcombos=[1,0,1];
SubAg=EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_ag);
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
fprintf('WhichCombos AggVars, [1,0,1]: selected Means vs full run, this should be zero: %.3e; unselected Means not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,notnan_off+missing+sum(~ismember(pS2,pF2)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=[1,0];
try
    EvalFnOnAgentDist_AggVars_FHorz_Case1(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AggVars, a wrong-length whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_AutoCorrTransProbs_FHorz =====================
simoptions_ac=struct();
simoptions_ac.timehorizons=[1,2]; % horizons 1 and 2 (since 2026-10-14 exactly the horizons listed are computed)
simoptions_ac.conditionalrestrictions=simoptions.conditionalrestrictions;
simoptions_ac.transprobs={'assets'};
FullAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pF3,vF3]=CoreStructLeaves(FullAC);
% (a) checkerboard over (fn, start age, page): a leaf's columns are start ages (1..N_j for Mean/StdDeviation, 1..N_j-k for the horizon
% outputs), so the selected columns equal the full run and the rest are NaN; a (fn, page) with nothing selected has no fields
wc3=zeros(nFns,N_j,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:N_j
            wc3(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
wc3(2,:,1)=0; % earnings unrestricted: nothing at all
simoptions_ac.whichcombos=wc3;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreStructLeaves(SubAC);
maxd=0; nanmis=0; notnan_off=0; missing=0; extra=0;
for k=1:numel(pF3)
    segs=strsplit(pF3{k},'.'); segs=segs(~cellfun(@isempty,segs));
    ncol=size(vF3{k},2);
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            expectpresent=true; onmask=true(1,ncol);
        else
            ff=find(strcmp(FnNames,segs{2}));
            expectpresent=any(wc3(ff,:,pp)); onmask=reshape(wc3(ff,1:ncol,pp)==1,[1,ncol]);
            if any(strcmp(segs{end},{'Mean','StdDeviation'})) % reported at every age (computed anyway)
                onmask=true(1,ncol);
            end
        end
    else
        ff=find(strcmp(FnNames,segs{1}));
        expectpresent=any(wc3(ff,:,1)); onmask=reshape(wc3(ff,1:ncol,1)==1,[1,ncol]);
        if any(strcmp(segs{end},{'Mean','StdDeviation'})) % the age-j Mean/StdDeviation are reported at every age (computed anyway)
            onmask=true(1,ncol);
        end
    end
    [present,idx]=ismember(pF3{k},pS3);
    if expectpresent
        if ~present || ~isequal(size(vF3{k}),size(vS3{idx}))
            missing=missing+1;
        else
            a=vF3{k}; b=vS3{idx};
            nanmis=nanmis+sum(isnan(a(onmask))~=isnan(b(onmask))); d=abs(a(onmask)-b(onmask)); d=d(~isnan(d)); maxd=max([maxd,d(:)']);
            notnan_off=notnan_off+sum(~isnan(b(~onmask)));
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos AutoCorr, checkerboard over (fn, start age, page): selected start ages vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected start ages not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS3,pF3)))
% TransitionProbs of assets (a cell over start ages 1..N_j-1): empty where the start age is not selected, equal where it is
badTP=0; dTP=0;
for jj=1:N_j-1
    if wc3(1,jj,1)==1
        if isempty(SubAC.assets.TransitionProbs{jj})
            badTP=badTP+1;
        else
            dTP=max(dTP,max(abs(gather(SubAC.assets.TransitionProbs{jj}(:))-gather(FullAC.assets.TransitionProbs{jj}(:)))));
            dTP=max(dTP,max(abs(gather(SubAC.assets.TransitionMass_j{jj}(:))-gather(FullAC.assets.TransitionMass_j{jj}(:)))));
        end
    else
        badTP=badTP+double(~isempty(SubAC.assets.TransitionProbs{jj}));
    end
end
fprintf('WhichCombos AutoCorr, TransitionProbs of assets: selected start ages vs full run, this should be zero: %.3e; non-empty cells at unselected start ages + empty cells at selected ones, this should be zero: %.3e \n',dTP,badTP)
% (b) the (fn, page) and (fn, start age) 2D forms expand as documented; (c) invalid sizes error
simoptions_ac.whichcombos=[1,1,0; 0,1,1; 1,0,1];
A=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
simoptions_ac.whichcombos=repmat(reshape([1,1,0; 0,1,1; 1,0,1],[nFns,1,npages]),[1,N_j,1]);
B=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
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
fprintf('WhichCombos AutoCorr, (fn, page) input vs its expansion over start ages, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=simoptions_ac; simoptions_bad.whichcombos=ones(nFns,N_j+1,npages);
try
    EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AutoCorr, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz =====================
simoptions_cs=struct();
FullCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_cs);
% assets: own stats + pair with richind; earnings: nothing; richind: pair with assets only (diagonal off)
wcs=[1,0,1; 0,0,0; 0,0,0];
wcsym=max(wcs,wcs');
fnwanted=any(wcsym,2);
simoptions_cs.whichcombos=wcs;
SubCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_cs);
[pF4,vF4]=CoreStructLeaves(FullCS); [pS4,vS4]=CoreStructLeaves(SubCS);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF4)
    segs=strsplit(pF4{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix','CorrelationMatrix'}))
        continue % checked below
    end
    ff=find(strcmp(FnNames,segs{1}));
    expectpresent=fnwanted(ff); % a function with nothing selected has no fields
    if numel(segs)==2 % .fn.Mean, .fn.StdDeviation, .fn.fn (self): reported whenever the function is evaluated (byproducts of the pairs)
        on=1;
    else % .fn.CovarianceWith.fn2 / .CorrelationWith.fn2: present iff both evaluated, equal iff the pair is on, else NaN
        ff2=find(strcmp(FnNames,segs{3}));
        expectpresent=fnwanted(ff) && fnwanted(ff2); on=wcsym(ff,ff2);
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
fprintf('WhichCombos CrossSection, pair selection: selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS4,pF4)))
Mf=gather(FullCS.CovarianceMatrix); Ms=gather(SubCS.CovarianceMatrix); Rf=gather(FullCS.CorrelationMatrix); Rs=gather(SubCS.CorrelationMatrix);
on=(max(wcsym,diag(fnwanted))==1); % pairs where selected; the diagonal (variance) wherever the function is evaluated
dd=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; dd=dd(~isnan(dd));
fprintf('WhichCombos CrossSection, matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected entries NaN where the full run is not (a zero-variance indicator gives NaN correlations in both), this should be zero: %.3e \n',max([dd;0]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)))+sum(isnan(Rs(on))&~isnan(Rf(on))))
simoptions_cs.whichcombos=wcsym;
SymCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_cs);
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
fprintf('WhichCombos CrossSection, symmetric input vs upper-triangle input, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pY,pS4)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=ones(nFns,nFns+1);
try
    EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos CrossSection, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz =====================
simoptions_acs=struct();
simoptions_acs.agegroupings=1:5:N_j;
ngroups=length(simoptions_acs.agegroupings);
FullACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_acs);
% earnings nothing at all; assets/richind own stats and their pair on a checkerboard over age groups
wca=zeros(nFns,nFns,ngroups);
for kk=1:ngroups
    wca(1,1,kk)=mod(kk,2); wca(3,3,kk)=mod(kk+1,2); wca(1,3,kk)=mod(kk,2);
end
wcasym=wca;
for kk=1:ngroups
    wcasym(:,:,kk)=max(wca(:,:,kk),wca(:,:,kk)');
end
simoptions_acs.whichcombos=wca;
SubACS=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_acs);
[pF5,vF5]=CoreStructLeaves(FullACS); [pS5,vS5]=CoreStructLeaves(SubACS);
maxd=0; nanmis=0; missing=0; notnan_off=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix','CorrelationMatrix'}))
        continue
    end
    ff=find(strcmp(FnNames,segs{1}));
    if numel(segs)==2
        onmask=reshape(any(wcasym(ff,:,:),2),[1,ngroups]); % own stats: reported in the age groups where the function is evaluated (byproducts of the pairs)
    else
        ff2=find(strcmp(FnNames,segs{3}));
        onmask=reshape(wcasym(ff,ff2,:)==1,[1,ngroups]);
    end
    [present,idx]=ismember(pF5{k},pS5);
    if ~present || ~isequal(size(vF5{k}),size(vS5{idx})) || numel(vF5{k})~=ngroups
        missing=missing+1;
    else
        a=reshape(vF5{k},[1,ngroups]); b=reshape(vS5{idx},[1,ngroups]);
        nanmis=nanmis+sum(isnan(a(onmask))~=isnan(b(onmask))); d=abs(a(onmask)-b(onmask)); d=d(~isnan(d)); maxd=max([maxd,d]);
        notnan_off=notnan_off+sum(~isnan(b(~onmask)));
    end
end
fprintf('WhichCombos AgeConditional CrossSection, pair x age-group selection: selected entries vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected entries not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+sum(~ismember(pS5,pF5)))
Mf=gather(FullACS.CovarianceMatrix); Ms=gather(SubACS.CovarianceMatrix); Rf=gather(FullACS.CorrelationMatrix); Rs=gather(SubACS.CorrelationMatrix);
on=(wcasym==1);
for kk=1:ngroups
    on(:,:,kk)=max(on(:,:,kk),diag(any(wcasym(:,:,kk),2))); % the diagonal (variance) wherever the function is evaluated in that age group
end
dd=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; dd=dd(~isnan(dd));
fprintf('WhichCombos AgeConditional CrossSection, matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected entries NaN where the full run is not (a zero-variance indicator gives NaN correlations in both), this should be zero: %.3e \n',max([dd;0]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)))+sum(isnan(Rs(on))&~isnan(Rf(on))))
simoptions_acs.whichcombos=wcs;
A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_acs);
simoptions_acs.whichcombos=repmat(wcs,[1,1,ngroups]);
B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_acs);
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
fprintf('WhichCombos AgeConditional CrossSection, 2D input vs its repmat over age groups, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=struct(); simoptions_bad.agegroupings=1:5:N_j; simoptions_bad.whichcombos=ones(nFns,nFns,ngroups+1);
try
    EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos AgeConditional CrossSection, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

output=struct();
output.LifeCycleProfiles_full=Full;
output.AutoCorr_full=FullAC;

end
