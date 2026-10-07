function output=CoreInfHorz_WhichCombos(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames)
% Test: simoptions.whichcombos (and per-combination simoptions.whichstats in AllStats) in the single-type InfHorz stats commands:
%   EvalFnOnAgentDist_AllStats_InfHorz                 whichcombos [nFns, 1+nRestr], whichstats [nFns, 1+nRestr, 7]
%   EvalFnOnAgentDist_AggVars_InfHorz                  whichcombos [nFns, 1]
%   EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz       whichcombos [nFns, 1+nRestr] (no start-age dimension)
%   EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz    whichcombos [nFns, nFns] or [nFns, nFns, 1+nRestr] (upper triangle; diagonal = own Mean/StdDev)
% Pages are the unrestricted output first, then the restrictions in the fieldnames order of simoptions.conditionalrestrictions.
% The InfHorz counterpart of CoreFHorzTests part 13 (CoreFHorz_WhichCombos).
%
% Every check compares a selective run against the full run leaf by leaf (CoreStructLeaves, in SharedSubcodes): a leaf the
% selection asks for must be present and equal to the full run, NaN pattern included; one it does not ask for must be NaN where
% the command preallocates it (AggVars, the CrossSection pairs and matrices) or absent where it assigns whole structs (AllStats,
% AutoCorr) or never creates them (a function with nothing selected). RestrictedSampleMass is always filled. Broadcast inputs equal
% their expansions; invalid inputs error.
% AllStats: StatsFromWeightedGrid reports Minimum/Maximum whenever quantiles are requested (they are the end cutoffs).

ReturnFn=@(d,aprime,a,z,r,w,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,sigma,eta,varphi);

FnsToEvaluate.assets=@(d,aprime,a,z) a;
FnsToEvaluate.earnings=@(d,aprime,a,z,w) w*d*z;
FnsToEvaluate.savings=@(d,aprime,a,z) aprime-a;
FnNames=fieldnames(FnsToEvaluate);
nFns=length(FnNames);

vfoptions=struct();
simoptions=struct();
simoptions.whichstats=ones(1,7);
simoptions.conditionalrestrictions.positiveassets=@(d,aprime,a,z) (a>0);
simoptions.conditionalrestrictions.highz=@(d,aprime,a,z) (z>1);
RestrNames=fieldnames(simoptions.conditionalrestrictions);
nRestr=length(RestrNames);
npages=1+nRestr;

%% Solve once
[~,Policy]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist=StationaryDist_InfHorz(Policy,n_d,n_a,n_z,pi_z,struct(),Params,[]);

%% ===================== EvalFnOnAgentDist_AllStats_InfHorz =====================
FullAS=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreStructLeaves(FullAS);

% (1) explicit all-ones equals the default
simoptions_1=simoptions; simoptions_1.whichcombos=ones(nFns,npages);
Sub=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_1);
[pS,vS]=CoreStructLeaves(Sub);
maxd=0; nanmis=0;
for k=1:numel(pF)
    [present,idx]=ismember(pF{k},pS);
    if present
        a=vF{k}; b=vS{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichCombos InfHorz AllStats, explicit all-ones vs default, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pF,pS))+sum(~ismember(pS,pF)))

% (2) checkerboard over (fn, page): selected leaves equal the full run, unselected ones are absent
wcA=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wcA(ff,pp)=mod(ff+pp,2);
    end
end
simoptions_a=simoptions; simoptions_a.whichcombos=wcA;
Sub=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_a);
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
fprintf('WhichCombos InfHorz AllStats, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))

% (3) per-combination whichstats: statistic kk on for (fn ff, page pp) iff mod(ff+pp+kk,2)==1
wsA=zeros(nFns,npages,7);
for pp=1:npages
    for ff=1:nFns
        for kk=1:7
            wsA(ff,pp,kk)=mod(ff+pp+kk,2);
        end
    end
end
simoptions_b=simoptions; simoptions_b.whichstats=wsA;
Sub=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_b);
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
fprintf('WhichStats InfHorz AllStats, per-combination pattern: leaves with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; on-leaves missing + off-leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS,pF)))
if ~isempty(badpaths)
    fprintf('WhichStats InfHorz AllStats, per-combination pattern: offending leaves (information): %s \n',strjoin(badpaths(1:min(8,end)),' | '))
end

% (4) a vector whichcombos equals its repmat over pages; [nFns,7] whichstats equals its repmat over pages
simoptions_c=simoptions; simoptions_c.whichcombos=[1;0;1];
A=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_c);
simoptions_c.whichcombos=repmat([1;0;1],[1,npages]);
simoptions_c.whichstats=repmat(reshape(ones(1,7),[1,1,7]),[nFns,npages,1]);
B=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_c);
simoptions_c.whichcombos=[1;0;1];
simoptions_c.whichstats=ones(nFns,7);
C=EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_c);
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
fprintf('WhichCombos InfHorz AllStats, vector whichcombos and [nFns,7] whichstats vs their repmats over pages, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA))+sum(~ismember(pC,pA)))

% (5) invalid inputs error
errored=0;
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns+1,npages);
try
    EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichcombos=ones(nFns,npages); simoptions_bad.whichcombos(1,1)=2;
try
    EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=simoptions; simoptions_bad.whichstats=ones(nFns,npages+1,7);
try
    EvalFnOnAgentDist_AllStats_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz AllStats, three invalid inputs must all error, this should be zero: %.3e \n',3-errored)

%% ===================== EvalFnOnAgentDist_AggVars_InfHorz =====================
simoptions_ag=struct();
FullAg=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_ag);
simoptions_ag.whichcombos=[1,0,1];
SubAg=EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_ag);
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
fprintf('WhichCombos InfHorz AggVars, [1,0,1]: selected Means vs full run, this should be zero: %.3e; unselected Means not NaN + leaves missing + extra leaves, this should be zero: %.3e \n',maxd,notnan_off+missing+sum(~ismember(pS2,pF2)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=[1,0];
try
    EvalFnOnAgentDist_AggVars_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz AggVars, a wrong-length whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz =====================
simoptions_ac=struct();
simoptions_ac.timehorizons=2;
simoptions_ac.conditionalrestrictions=simoptions.conditionalrestrictions;
simoptions_ac.transprobs={'assets'};
FullAC=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pF3,vF3]=CoreStructLeaves(FullAC);
% (a) checkerboard over (fn, page), and earnings unrestricted off too: a selected (fn, page) has every output (every horizon, and
% the TransitionProbs of assets on the unrestricted page) equal to the full run, an unselected one has no fields
wc3=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wc3(ff,pp)=mod(ff+pp,2);
    end
end
wc3(2,1)=0;
simoptions_ac.whichcombos=wc3;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
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
fprintf('WhichCombos InfHorz AutoCorr, checkerboard over (fn, page): selected leaves (every horizon, TransitionProbs included) vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present + extra leaves, this should be zero: %.3e \n',maxd,nanmis,missing+extra+sum(~ismember(pS3,pF3)))
% (b) a vector equals its repmat over pages; (c) invalid sizes error
simoptions_ac.whichcombos=[1,0,1];
A=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
simoptions_ac.whichcombos=repmat([1;0;1],[1,npages]);
B=EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
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
fprintf('WhichCombos InfHorz AutoCorr, vector input vs its repmat over pages, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=simoptions_ac; simoptions_bad.whichcombos=ones(nFns,npages+1);
try
    EvalFnOnAgentDist_AutoCorrTransProbs_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz AutoCorr, a wrong-size whichcombos must error, this should be zero: %.3e \n',1-errored)

%% ===================== EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz =====================
% (a) without restrictions. assets: own stats + pair with savings; earnings: nothing; savings: pair with assets only (diagonal off)
simoptions_cs=struct();
FullCS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_cs);
wcs=[1,0,1; 0,0,0; 0,0,0];
wcsym=max(wcs,wcs');
fnwanted=any(wcsym,2);
simoptions_cs.whichcombos=wcs;
SubCS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_cs);
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
fprintf('WhichCombos InfHorz CrossSection, pair selection: selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS4,pF4)))
Mf=gather(FullCS.CovarianceMatrix); Ms=gather(SubCS.CovarianceMatrix); Rf=gather(FullCS.CorrelationMatrix); Rs=gather(SubCS.CorrelationMatrix);
on=(max(wcsym,diag(fnwanted))==1); % pairs where selected; the diagonal (variance) wherever the function is evaluated
dd=[abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; dd=dd(~isnan(dd));
fprintf('WhichCombos InfHorz CrossSection, matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected entries NaN where the full run is not, this should be zero: %.3e \n',max([dd;0]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)))+sum(isnan(Rs(on))&~isnan(Rf(on))))
simoptions_cs.whichcombos=wcsym;
SymCS=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_cs);
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
fprintf('WhichCombos InfHorz CrossSection, symmetric input vs upper-triangle input, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pY,pS4)))

% (b) with restrictions: a different upper-triangle selection on every page; per page, a function's own fields are present iff it is
% evaluated on that page, a pair's fields iff both are (equal iff the pair is on, else NaN); RestrictedSampleMass always filled
simoptions_csr=struct(); simoptions_csr.conditionalrestrictions=simoptions.conditionalrestrictions;
FullCSR=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_csr);
wcr=zeros(nFns,nFns,npages);
wcr(:,:,1)=[1,0,1; 0,0,0; 0,0,0];
wcr(:,:,2)=[0,1,0; 0,0,0; 0,0,1];
wcr(:,:,3)=[0,0,0; 0,1,1; 0,0,0];
wcrsym=wcr;
for pp=1:npages
    wcrsym(:,:,pp)=max(wcr(:,:,pp),wcr(:,:,pp)');
end
fnwantedpage=reshape(any(wcrsym,2),[nFns,npages]);
simoptions_csr.whichcombos=wcr;
SubCSR=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_csr);
[pF5,vF5]=CoreStructLeaves(FullCSR); [pS5,vS5]=CoreStructLeaves(SubCSR);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1})); segs=segs(2:end);
    else
        pp=1;
    end
    if any(strcmp(segs{1},{'CovarianceMatrix','CorrelationMatrix'}))
        continue % checked below
    end
    [present,idx]=ismember(pF5{k},pS5);
    if strcmp(segs{1},'RestrictedSampleMass')
        expectpresent=true; on=1;
    else
        ff=find(strcmp(FnNames,segs{1}));
        expectpresent=fnwantedpage(ff,pp);
        if numel(segs)==2
            on=1;
        else
            ff2=find(strcmp(FnNames,segs{3}));
            expectpresent=fnwantedpage(ff,pp) && fnwantedpage(ff2,pp); on=wcrsym(ff,ff2,pp);
        end
    end
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
fprintf('WhichCombos InfHorz CrossSection, per-page pair selection with restrictions: selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra+sum(~ismember(pS5,pF5)))
dd=[]; bad=0;
for pp=1:npages
    if pp==1
        PF=FullCSR; PS=SubCSR;
    else
        PF=FullCSR.(RestrNames{pp-1}); PS=SubCSR.(RestrNames{pp-1});
    end
    Mf=gather(PF.CovarianceMatrix); Ms=gather(PS.CovarianceMatrix); Rf=gather(PF.CorrelationMatrix); Rs=gather(PS.CorrelationMatrix);
    on=(max(wcrsym(:,:,pp),diag(fnwantedpage(:,pp)))==1);
    dd=[dd;abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]; %#ok<AGROW>
    bad=bad+sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on)))+sum(isnan(Ms(on))&~isnan(Mf(on)))+sum(isnan(Rs(on))&~isnan(Rf(on)));
end
dd=dd(~isnan(dd));
fprintf('WhichCombos InfHorz CrossSection, matrices on every page: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN + selected entries NaN where the full run is not, this should be zero: %.3e \n',max([dd;0]),bad)

% (c) with restrictions, a 2D input equals its repmat over pages; (d) invalid sizes error
simoptions_csr.whichcombos=wcs;
A=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_csr);
simoptions_csr.whichcombos=repmat(wcs,[1,1,npages]);
B=EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_csr);
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
fprintf('WhichCombos InfHorz CrossSection, 2D input with restrictions vs its repmat over pages, every leaf, this should be zero: %.3e, mismatches %i \n',maxd,nanmis+sum(~ismember(pB,pA)))
errored=0;
simoptions_bad=struct(); simoptions_bad.whichcombos=ones(nFns,nFns+1);
try
    EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
simoptions_bad=struct(); simoptions_bad.conditionalrestrictions=simoptions.conditionalrestrictions; simoptions_bad.whichcombos=ones(nFns,nFns,npages+1);
try
    EvalFnOnAgentDist_CrossSectionCovarCorr_InfHorz(StationaryDist,Policy,FnsToEvaluate,Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos InfHorz CrossSection, two wrong-size whichcombos must both error, this should be zero: %.3e \n',2-errored)

output=struct();
output.AllStats_full=FullAS;
output.AutoCorr_full=FullAC;

end
