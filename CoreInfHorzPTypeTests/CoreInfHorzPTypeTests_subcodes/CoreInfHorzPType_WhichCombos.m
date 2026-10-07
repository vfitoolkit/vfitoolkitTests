function output=CoreInfHorzPType_WhichCombos(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test (CoreInfHorzPTypeTests part 15): simoptions.whichcombos (and per-combination simoptions.whichstats in AllStats) across the
% InfHorz PType stats commands:
%   EvalFnOnAgentDist_AllStats_InfHorz_PType                whichcombos [nFns, 1+nRestr], whichstats [nFns, 1+nRestr, 7]
%   EvalFnOnAgentDist_AggVars_InfHorz_PType                 whichcombos [nFns, 1]
%   EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType      whichcombos [nFns, 1+nRestr]
%   EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType   whichcombos [nFns, nFns] or [nFns, nFns, 1+nRestr] (upper triangle; diagonal = own Mean/StdDev)
% Every command but AggVars also takes a trailing type dimension of N_i+1 (the types in Names_i order, then the grouped slot; a grouped
% slot forces every type's computation of that combination; a type's output is reported wherever computed, the grouped output only where
% the grouped slot asks). Pages are the unrestricted output first, then the restrictions in the fieldnames order of
% simoptions.conditionalrestrictions. The InfHorz counterpart of CoreFHorzPTypeTests part 21 (CoreFHorzPType_WhichCombos).
%
% Every check compares a selective run against the full run leaf by leaf (CoreStructLeaves, in SharedSubcodes): a leaf the selection
% asks for must be present and equal to the full run, NaN pattern included; one it does not ask for must be NaN where the command
% preallocates it (AggVars, the grouped CrossSection leaves and the pairs of evaluated functions) or absent where it assigns whole
% structs (AllStats, AutoCorr) or never creates them (a function not evaluated). RestrictedSampleMass is always filled. The per-type
% covariance matrices of the CrossSection command shrink to the evaluated functions and are checked separately. Broadcast inputs
% equal their expansions; invalid inputs error.
% AllStats: StatsFromWeightedGrid reports Minimum/Maximum whenever quantiles are requested (they are the end cutoffs).

Names_i={'low','high'};
N_i=length(Names_i);

ReturnFn=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);

FnsToEvaluate.assets=@(d,aprime,a,z) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
FnsToEvaluate.savings=@(d,aprime,a,z) aprime-a;
FnNames=fieldnames(FnsToEvaluate);
nFns=length(FnNames);
fa=find(strcmp(FnNames,'assets')); fe=find(strcmp(FnNames,'earnings')); fs=find(strcmp(FnNames,'savings'));
ilow=find(strcmp(Names_i,'low'));

vfoptions=struct();
simoptions=struct();
simoptions.whichstats=ones(1,7);
simoptions.conditionalrestrictions.positiveassets=@(d,aprime,a,z) (a>0);
simoptions.conditionalrestrictions.highz=@(d,aprime,a,z) (z>1);
RestrNames=fieldnames(simoptions.conditionalrestrictions);
nRestr=length(RestrNames);
npages=1+nRestr;
ppos=1+find(strcmp(RestrNames,'positiveassets'));

%% Solve once
[~,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
StationaryDist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);

%% ===================== AllStats PType =====================
FullAS=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreStructLeaves(FullAS);

% (a) checkerboard over (fn, page): a skipped (fn, page) has no per-type struct and no grouped fields
wcA=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wcA(ff,pp)=mod(ff+pp,2);
    end
end
simoptions_a=simoptions; simoptions_a.whichcombos=wcA;
Sub=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_a);
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
fprintf('WhichCombos InfHorz AllStats PType, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))

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
Sub=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_b);
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
fprintf('WhichStats InfHorz AllStats PType, per-combination pattern: leaves with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; on-leaves missing + off-leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))
if ~isempty(badpaths)
    fprintf('WhichStats InfHorz AllStats PType, per-combination pattern: offending leaves (information): %s \n',strjoin(badpaths(1:min(8,end)),' | '))
end

% (c) [nFns,7] whichstats equals its repmat over pages, [nFns,1] whichcombos its repmat; the 2D whichcombos equals its 3D broadcast over
% the ptypes and the grouped slot, the 3D whichstats its 4D broadcast; invalid inputs error
simoptions_c=simoptions; simoptions_c.whichstats=reshape(wsA(:,1,:),[nFns,7]); simoptions_c.whichcombos=[1;0;1];
A=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_c);
simoptions_c.whichstats=repmat(reshape(wsA(:,1,:),[nFns,1,7]),[1,npages,1]); simoptions_c.whichcombos=repmat([1;0;1],[1,npages]);
B=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_c);
simoptions_c.whichstats=repmat(reshape(wsA(:,1,:),[nFns,1,1,7]),[1,npages,N_i+1,1]); simoptions_c.whichcombos=repmat([1;0;1],[1,npages,N_i+1]);
C=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_c);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B); [pC,vC]=CoreStructLeaves(C);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    [presentC,idxC]=ismember(pA{k},pC);
    if present && presentC && isequal(size(vA{k}),size(vB{idx})) && isequal(size(vA{k}),size(vC{idxC}))
        a=vA{k}; b=vB{idx}; c=vC{idxC}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:)))+sum(isnan(a(:))~=isnan(c(:)));
        d=[abs(a(:)-b(:));abs(a(:)-c(:))]; d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos/WhichStats InfHorz AllStats PType, broadcast inputs (over pages, and over the ptypes and the grouped slot) vs their repmats, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA))+sum(~ismember(pC,pA)))
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns+1,npages);
try
    EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns,npages,N_i);
try
    EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichstats=ones(nFns,npages,6);
try
    EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz AllStats PType, three invalid inputs (wrong fn count, type dimension N_i instead of N_i+1, whichstats of 6) must all error, this should be zero: %.3e \n',3-errored)

% (d) type-specific selection: type low gets assets (page 1), type high nothing of its own, the grouped slot gets earnings on page 1 and on
% the positiveassets page. A grouped slot forces the per-type stats of that combination for every type, so earnings must be present for
% both types too. AllStats assigns whole structs, so a per-type or grouped leaf is present iff selected (own slot, or forced);
% RestrictedSampleMass leaves are always present.
wcAT=zeros(nFns,npages,N_i+1);
wcAT(fa,1,ilow)=1;
wcAT(fe,1,N_i+1)=1;
wcAT(fe,ppos,N_i+1)=1;
simoptions_t=simoptions; simoptions_t.whichcombos=wcAT;
Sub=EvalFnOnAgentDist_AllStats_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_t);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    pp=1;
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1})); segs=segs(2:end);
    end
    if strcmp(segs{1},'RestrictedSampleMass')
        expectpresent=true;
    else
        ff=find(strcmp(FnNames,segs{1}));
        if numel(segs)>=2 && ismember(segs{2},Names_i)
            slot=find(strcmp(Names_i,segs{2}));
            expectpresent=(wcAT(ff,pp,slot)==1 || wcAT(ff,pp,N_i+1)==1); % own slot, or forced by the grouped slot
        else
            expectpresent=(wcAT(ff,pp,N_i+1)==1);
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
fprintf('WhichCombos InfHorz AllStats PType, type-specific selection (low: assets; grouped: earnings, unrestricted and positiveassets): selected leaves, incl. the per-type stats a grouped slot forces, vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))

%% ===================== AggVars PType =====================
simoptions_ag=struct();
FullAg=EvalFnOnAgentDist_AggVars_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
simoptions_ag.whichcombos=[1,0,1];
SubAg=EvalFnOnAgentDist_AggVars_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
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
fprintf('WhichCombos InfHorz AggVars PType, [1,0,1]: selected Means (per type and grouped) vs full run, this should be zero: %.3e; unselected Means that are not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,notnan_off+missing+sum(~ismember(pS2,pF2)))
simoptions_ag.whichcombos=[0;0;0];
ZAg=EvalFnOnAgentDist_AggVars_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
[pZ,vZ]=CoreStructLeaves(ZAg);
fprintf('WhichCombos InfHorz AggVars PType, all zeros: leaves present %i (full run has %i), non-NaN entries + missing leaves, this should be zero: %.3e \n',numel(pZ),numel(pF2),sum(cellfun(@(x) sum(~isnan(x(:))),vZ))+sum(~ismember(pF2,pZ)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=[1,0];
try
    EvalFnOnAgentDist_AggVars_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz AggVars PType, a wrong-length whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== AutoCorrTransProbs PType =====================
simoptions_ac=struct();
simoptions_ac.timehorizons=2;
simoptions_ac.conditionalrestrictions=simoptions.conditionalrestrictions;
simoptions_ac.transprobs={'assets'};
FullAC=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pF3,vF3]=CoreStructLeaves(FullAC);
% (a) checkerboard over (fn, page), and earnings unrestricted off too: a selected (fn, page) has every per-type and grouped output (every
% horizon, the per-type TransitionProbs of assets on the unrestricted page) equal to the full run, an unselected one has none
wc3=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wc3(ff,pp)=mod(ff+pp,2);
    end
end
wc3(fe,1)=0;
simoptions_ac.whichcombos=wc3;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreStructLeaves(SubAC);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF3)
    segs=strsplit(pF3{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            expectpresent=true;
        else
            expectpresent=(wc3(find(strcmp(FnNames,segs{2})),pp)==1);
        end
    else
        expectpresent=(wc3(find(strcmp(FnNames,segs{1})),1)==1);
    end
    [present,idx]=ismember(pF3{k},pS3);
    if expectpresent
        if ~present || ~isequal(size(vF3{k}),size(vS3{idx}))
            missing=missing+1;
        else
            a=vF3{k}; b=vS3{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos InfHorz AutoCorr PType, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present + extra, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS3,pF3)))
% (b) a vector equals its repmat over pages, a 2D input its 3D broadcast over the ptypes and the grouped slot; invalid sizes error
simoptions_ac.whichcombos=[1,0,1];
A=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
simoptions_ac.whichcombos=repmat([1;0;1],[1,npages]);
B=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
simoptions_ac.whichcombos=repmat([1;0;1],[1,npages,N_i+1]);
C=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B); [pC,vC]=CoreStructLeaves(C);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    [presentC,idxC]=ismember(pA{k},pC);
    if present && presentC && isequal(size(vA{k}),size(vB{idx})) && isequal(size(vA{k}),size(vC{idxC}))
        a=vA{k}; b=vB{idx}; c=vC{idxC}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:)))+sum(isnan(a(:))~=isnan(c(:)));
        d=[abs(a(:)-b(:));abs(a(:)-c(:))]; d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos InfHorz AutoCorr PType, vector input vs its repmat over pages and its broadcast over the ptypes and the grouped slot, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA))+sum(~ismember(pC,pA)))
errored=0;
simoptions_bad=simoptions_ac; simoptions_bad.whichcombos=ones(nFns,npages+1);
try
    EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions_ac; simoptions_bad.whichcombos=ones(nFns,npages,N_i);
try
    EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz AutoCorr PType, two wrong-size whichcombos must both error, this should be zero: %.3e \n',2-errored)
% (c) type-specific selection: type low gets unrestricted assets, type high nothing of its own, the grouped slot gets earnings,
% unrestricted and on the positiveassets page. A grouped slot forces that computation for every type (the grouped outputs are pooled
% from every type's), so earnings must be present for both types too; a (fn, page) with nothing selected for a slot has no leaf there;
% RestrictedSampleMass leaves are always present.
wcT=zeros(nFns,npages,N_i+1);
wcT(fa,1,ilow)=1;
wcT(fe,1,N_i+1)=1;
wcT(fe,ppos,N_i+1)=1;
simoptions_ac.whichcombos=wcT;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreStructLeaves(SubAC);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF3)
    segs=strsplit(pF3{k},'.'); segs=segs(~cellfun(@isempty,segs));
    pp=1;
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1})); segs=segs(2:end);
    end
    if strcmp(segs{1},'RestrictedSampleMass')
        expectpresent=true;
    else
        ff=find(strcmp(FnNames,segs{1}));
        if numel(segs)>=2 && ismember(segs{2},Names_i)
            slot=find(strcmp(Names_i,segs{2}));
            expectpresent=(wcT(ff,pp,slot)==1 || wcT(ff,pp,N_i+1)==1); % own slot, or forced by the grouped slot
        else
            expectpresent=(wcT(ff,pp,N_i+1)==1);
        end
    end
    [present,idx]=ismember(pF3{k},pS3);
    if expectpresent
        if ~present || ~isequal(size(vF3{k}),size(vS3{idx}))
            missing=missing+1;
        else
            a=vF3{k}; b=vS3{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos InfHorz AutoCorr PType, type-specific selection (low: assets; grouped: earnings, unrestricted and positiveassets): selected leaves, incl. the per-type outputs a grouped slot forces, vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present + extra, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS3,pF3)))

%% ===================== CrossSectionCovarCorr PType =====================
% (a) with restrictions, a different upper-triangle selection on every page (no type dimension, so every type and the grouped stats get
% the same selection). Per page: a per-type function leaf (.fn.type.Mean etc) is present iff the function is evaluated on that page, a
% per-type pair leaf iff both are (equal iff the pair is on, else NaN); a grouped leaf is always present, equal iff selected (its own
% stats wherever any grouped entry of its row is on, a pair where the pair is on) else NaN; RestrictedSampleMass always filled
simoptions_cs=struct(); simoptions_cs.conditionalrestrictions=simoptions.conditionalrestrictions;
FullCS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
[pF4,vF4]=CoreStructLeaves(FullCS);
wcr=zeros(nFns,nFns,npages);
wcr(:,:,1)=[1,0,1; 0,0,0; 0,0,1];
wcr(:,:,2)=[0,1,0; 0,0,0; 0,0,1];
wcr(:,:,3)=[0,0,0; 0,1,1; 0,0,0];
wcrsym=wcr;
for pp=1:npages
    wcrsym(:,:,pp)=max(wcr(:,:,pp),wcr(:,:,pp)');
end
fnwantedpage=reshape(any(wcrsym,2),[nFns,npages]);
simoptions_cs.whichcombos=wcr;
SubCS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
[pS4,vS4]=CoreStructLeaves(SubCS);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF4)
    segs=strsplit(pF4{k},'.'); segs=segs(~cellfun(@isempty,segs));
    pp=1;
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1})); segs=segs(2:end);
    end
    if any(strcmp(segs{1},{'CovarianceMatrix_ptype','CorrelationMatrix_ptype','CovarianceMatrix','CorrelationMatrix'}))
        continue % the matrices are checked below
    end
    if strcmp(segs{1},'RestrictedSampleMass')
        expectpresent=true; on=1;
    else
        ff=find(strcmp(FnNames,segs{1}));
        if ismember(segs{2},Names_i) % per-type leaf
            if numel(segs)==3 % .fn.type.Mean / .StdDeviation / .fn (self)
                expectpresent=fnwantedpage(ff,pp); on=1;
            else % .fn.type.CovarianceWith.fn2
                ff2=find(strcmp(FnNames,segs{4}));
                expectpresent=fnwantedpage(ff,pp) && fnwantedpage(ff2,pp); on=wcrsym(ff,ff2,pp);
            end
        else % grouped leaf
            expectpresent=true;
            if numel(segs)==2
                on=any(wcrsym(ff,:,pp));
            else
                ff2=find(strcmp(FnNames,segs{3}));
                on=wcrsym(ff,ff2,pp);
            end
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
fprintf('WhichCombos InfHorz CrossSection PType, per-page pair selection with restrictions: selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS4,pF4)))
dd=[]; bad=0; dmat=0;
selall=find(any(fnwantedpage,2)); % the functions each type evaluates (on some page): its matrices shrink to them
for pp=1:npages
    if pp==1
        PF=FullCS; PS=SubCS;
    else
        PF=FullCS.(RestrNames{pp-1}); PS=SubCS.(RestrNames{pp-1});
    end
    Mf=gather(PF.CovarianceMatrix); Ms=gather(PS.CovarianceMatrix); Rf=gather(PF.CorrelationMatrix); Rs=gather(PS.CorrelationMatrix);
    on=(max(wcrsym(:,:,pp),diag(any(wcrsym(:,:,pp),2)))==1); % grouped pairs where selected; the diagonal wherever any grouped entry of the row is on
    dd=[dd;abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; %#ok<AGROW>
    bad=bad+sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)));
    for ii=1:N_i
        Mf_ii=gather(PF.CovarianceMatrix_ptype.(Names_i{ii})); Ms_ii=gather(PS.CovarianceMatrix_ptype.(Names_i{ii}));
        if ~isequal(size(Ms_ii),[numel(selall),numel(selall)])
            bad=bad+1;
        else
            onsub=max(wcrsym(selall,selall,pp),diag(fnwantedpage(selall,pp)))==1;
            Mf_sub=Mf_ii(selall,selall);
            dmat=max([dmat;abs(Mf_sub(onsub)-Ms_ii(onsub))]);
            bad=bad+sum(~isnan(Ms_ii(~onsub)));
        end
    end
end
dd=dd(~isnan(dd));
fprintf('WhichCombos InfHorz CrossSection PType, grouped and per-type (shrunk) covariance matrices on every page: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected grouped entries NaN where the full run is not + wrong sizes, this should be zero: %.3e \n',max([dd;dmat;0]),bad)

% (b) without restrictions: symmetric input vs upper triangle; a 2D input equals its 3D broadcast over the ptypes and the grouped slot;
% with restrictions a 3D input equals its 4D broadcast; invalid sizes error
wcs=[1,0,1; 0,0,0; 0,0,1];
wcsym=max(wcs,wcs');
simoptions_cs0=struct(); simoptions_cs0.whichcombos=wcs;
A=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs0);
simoptions_cs0.whichcombos=wcsym;
B=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs0);
simoptions_cs0.whichcombos=repmat(wcs,[1,1,N_i+1]);
C=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs0);
simoptions_cs.whichcombos=wcr;
D=SubCS;
simoptions_cs.whichcombos=repmat(wcr,[1,1,1,N_i+1]);
E=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
[pA,vA]=CoreStructLeaves(A); [pB,vB]=CoreStructLeaves(B); [pC,vC]=CoreStructLeaves(C); [pD,vD]=CoreStructLeaves(D); [pE,vE]=CoreStructLeaves(E);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    [presentC,idxC]=ismember(pA{k},pC);
    if present && presentC && isequal(size(vA{k}),size(vB{idx})) && isequal(size(vA{k}),size(vC{idxC}))
        a=vA{k}; b=vB{idx}; c=vC{idxC}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:)))+sum(isnan(a(:))~=isnan(c(:)));
        d=[abs(a(:)-b(:));abs(a(:)-c(:))]; d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
for k=1:numel(pD)
    [present,idx]=ismember(pD{k},pE);
    if present && isequal(size(vD{k}),size(vE{idx}))
        a=vD{k}; b=vE{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos InfHorz CrossSection PType, symmetric vs upper-triangle input, and inputs without the type dimension vs their broadcast over the ptypes and the grouped slot (with and without restrictions), every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA))+sum(~ismember(pC,pA))+sum(~ismember(pE,pD)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=ones(nFns,nFns+1);
try
    EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=struct(); simoptions_bad.conditionalrestrictions=simoptions.conditionalrestrictions; simoptions_bad.whichcombos=ones(nFns,nFns,npages,N_i);
try
    EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz CrossSection PType, two wrong-size whichcombos must both error, this should be zero: %.3e \n',2-errored)

% (c) type-specific selection without restrictions: type low gets the own stats of assets and the pair (assets, savings), type high
% nothing of its own, the grouped slot gets the pair (earnings, savings) only. The grouped pair forces its two functions and their pair
% for every type, so for type high earnings and savings are evaluated with that pair only; the own stats of an evaluated function are
% reported for the type (byproducts); a per-type pair of two evaluated functions is NaN unless selected for that type (own or forced);
% the grouped leaves are always present and NaN unless the grouped slot asks (the grouped own stats of a function wherever any grouped
% entry of its row is on). The matrices are checked separately.
simoptions_cs0=struct();
FullCS0=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs0);
[pF5,vF5]=CoreStructLeaves(FullCS0);
wcT=zeros(nFns,nFns,N_i+1);
wcT(fa,fa,ilow)=1; wcT(fa,fs,ilow)=1;
wcT(fe,fs,N_i+1)=1;
effG=max(wcT(:,:,N_i+1),wcT(:,:,N_i+1)');
eff=cell(N_i,1);
for ii=1:N_i
    e=max(wcT(:,:,ii),wcT(:,:,N_i+1)); eff{ii}=max(e,e');
end
simoptions_cs0.whichcombos=wcT;
SubCS0=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,simoptions_cs0);
[pS5,vS5]=CoreStructLeaves(SubCS0);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix_ptype','CorrelationMatrix_ptype','CovarianceMatrix','CorrelationMatrix'}))
        continue % the matrices are checked below
    end
    ff=find(strcmp(FnNames,segs{1}));
    if ismember(segs{2},Names_i) % per-type leaf
        Ei=eff{find(strcmp(Names_i,segs{2}))}; evaluated=any(Ei,2);
        if numel(segs)==3
            expectpresent=evaluated(ff); on=1;
        else
            ff2=find(strcmp(FnNames,segs{4}));
            expectpresent=evaluated(ff) && evaluated(ff2); on=Ei(ff,ff2);
        end
    else % grouped leaf
        expectpresent=true;
        if numel(segs)==2
            on=any(effG(ff,:));
        else
            ff2=find(strcmp(FnNames,segs{3}));
            on=effG(ff,ff2);
        end
    end
    [present,idx]=ismember(pF5{k},pS5);
    if expectpresent
        if ~present
            missing=missing+1;
        elseif on==1
            a=vF5{k}; b=vS5{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            notnan_off=notnan_off+sum(~isnan(vS5{idx}(:)));
        end
    elseif present
        extra=extra+1;
    end
end
fprintf('WhichCombos InfHorz CrossSection PType, type-specific selection (low: assets own + assets-savings; grouped: earnings-savings): selected leaves, incl. the per-type outputs a grouped slot forces, vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS5,pF5)))
Mf=gather(FullCS0.CovarianceMatrix); Ms=gather(SubCS0.CovarianceMatrix); Rf=gather(FullCS0.CorrelationMatrix); Rs=gather(SubCS0.CorrelationMatrix);
on=(max(effG,diag(any(effG,2)))==1);
dd=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; dd=dd(~isnan(dd));
dmat=0; badmat=sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)));
for ii=1:N_i
    Ei=eff{ii}; sel=find(any(Ei,2));
    Mf_ii=gather(FullCS0.CovarianceMatrix_ptype.(Names_i{ii})); Ms_ii=gather(SubCS0.CovarianceMatrix_ptype.(Names_i{ii}));
    if ~isequal(size(Ms_ii),[numel(sel),numel(sel)])
        badmat=badmat+1;
    else
        onsub=max(Ei(sel,sel),eye(numel(sel)))==1;
        Mf_sub=Mf_ii(sel,sel);
        dmat=max([dmat;abs(Mf_sub(onsub)-Ms_ii(onsub))]);
        badmat=badmat+sum(~isnan(Ms_ii(~onsub)));
    end
end
fprintf('WhichCombos InfHorz CrossSection PType, type-specific selection, grouped and per-type covariance matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + wrong sizes, this should be zero: %.3e \n',max([dd;dmat;0]),badmat)

output=struct();
output.AllStats_full=FullAS;
output.CrossSection_full=FullCS;

end
