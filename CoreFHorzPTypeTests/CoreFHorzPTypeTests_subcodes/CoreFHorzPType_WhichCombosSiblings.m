function output=CoreFHorzPType_WhichCombosSiblings(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: simoptions.whichcombos (and, for AllStats, per-combination simoptions.whichstats) in the other PType commands:
%   EvalFnOnAgentDist_AllStats_FHorz_Case1_PType                       whichcombos [nFns, 1+nRestr], whichstats [nFns, 1+nRestr, 7]
%   EvalFnOnAgentDist_AggVars_FHorz_Case1_PType                        whichcombos [nFns, 1]
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType                   whichcombos [nFns, N_j, 1+nRestr] (age dimension validated, not acted on)
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType                whichcombos [nFns, nFns] (upper triangle; diagonal = own Mean/StdDev)
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType   whichcombos [nFns, nFns, nAgeGroups]
% Every check compares a selective run against the full run leaf by leaf (CoreFHorzPType_StructLeaves): a leaf the selection asks
% for must be present and equal to the full run (NaN pattern included); a leaf it does not ask for must be absent (AllStats,
% AutoCorr: whole (fn, page) blocks are not produced) or NaN (AggVars, CrossSection: the fields exist and are NaN). The
% per-type covariance matrices of the CrossSection commands shrink to the evaluated functions, so they are checked separately.
% Invalid inputs must error.

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

%% ===================== AllStats PType =====================
Full=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
[pF,vF]=CoreFHorzPType_StructLeaves(Full);

% (a) checkerboard over (fn, page)
wc=zeros(nFns,npages);
for pp=1:npages
    for ff=1:nFns
        wc(ff,pp)=mod(ff+pp,2);
    end
end
simoptions_a=simoptions; simoptions_a.whichcombos=wc;
Sub=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_a);
[pS,vS]=CoreFHorzPType_StructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            expectpresent=true;
        else
            expectpresent=(wc(find(strcmp(FnNames,segs{2})),pp)==1);
        end
    else
        expectpresent=(wc(find(strcmp(FnNames,segs{1})),1)==1);
    end
    [present,idx]=ismember(pF{k},pS);
    if expectpresent
        if ~present
            missing=missing+1;
        else
            a=vF{k}; b=vS{idx};
            nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
    end
end
extra=extra+sum(~ismember(pS,pF));
fprintf('WhichCombos AllStats PType, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra)

% (b) whichstats per combination: statistic kk on for (fn ff, page pp) iff mod(ff+pp+kk,2)==1; a leaf is present iff its statistic is on
ws=zeros(nFns,npages,7);
for pp=1:npages
    for ff=1:nFns
        for kk=1:7
            ws(ff,pp,kk)=mod(ff+pp+kk,2);
        end
    end
end
simoptions_b=simoptions; simoptions_b.whichstats=ws;
Sub=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_b);
[pS,vS]=CoreFHorzPType_StructLeaves(Sub);
maxd=0; nanmis=0; missing=0; extra=0;
for k=1:numel(pF)
    segs=strsplit(pF{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if ismember(segs{1},RestrNames)
        pp=1+find(strcmp(RestrNames,segs{1}));
        if strcmp(segs{2},'RestrictedSampleMass')
            ff=0;
        else
            ff=find(strcmp(FnNames,segs{2}));
        end
    else
        pp=1;
        ff=find(strcmp(FnNames,segs{1}));
    end
    if ff==0
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
                    kk=5; % (see below: StatsFromWeightedGrid also reports Minimum/Maximum whenever quantiles are requested)
                case {'QuantileCutoffs','QuantileMeans'}
                    kk=6;
            end
        end
        expectpresent=all(ws(ff,pp,kk)==1);
        if any(strcmp(statname,{'Minimum','Maximum'})) % StatsFromWeightedGrid reports Minimum and Maximum whenever quantiles are requested too (they are the end cutoffs)
            expectpresent=(ws(ff,pp,5)==1 || ws(ff,pp,6)>=1);
        end
    end
    [present,idx]=ismember(pF{k},pS);
    if expectpresent
        if ~present
            missing=missing+1;
        else
            a=vF{k}; b=vS{idx};
            nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        end
    elseif present
        extra=extra+1;
    end
end
extra=extra+sum(~ismember(pS,pF));
fprintf('WhichStats AllStats PType, per-combination pattern: leaves with the statistic on vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; on-leaves missing + off-leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra)

% (c) [nFns,7] whichstats equals its repmat over pages; [nFns,1] whichcombos equals its repmat; invalid inputs error
simoptions_c=simoptions; simoptions_c.whichstats=reshape(ws(:,1,:),[nFns,7]); simoptions_c.whichcombos=[1;0;1];
A=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_c);
simoptions_c.whichstats=repmat(reshape(ws(:,1,:),[nFns,1,7]),[1,npages,1]); simoptions_c.whichcombos=repmat([1;0;1],[1,npages]);
B=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_c);
[pA,vA]=CoreFHorzPType_StructLeaves(A); [pB,vB]=CoreFHorzPType_StructLeaves(B);
maxd=0; nanmis=0;
for k=1:numel(pA)
    [present,idx]=ismember(pA{k},pB);
    if present
        a=vA{k}; b=vB{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    end
end
fprintf('WhichCombos/WhichStats AllStats PType, broadcast inputs vs their repmats, every leaf, this should be zero: %.3e, NaN-pattern mismatches + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pA,pB))+sum(~ismember(pB,pA)))
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
[pF2,vF2]=CoreFHorzPType_StructLeaves(FullAg); [pS2,vS2]=CoreFHorzPType_StructLeaves(SubAg);
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
[pZ,vZ]=CoreFHorzPType_StructLeaves(ZAg);
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
[pF3,vF3]=CoreFHorzPType_StructLeaves(FullAC);
% (a) (fn, page) checkerboard, every start age
wc3=zeros(nFns,N_j,npages);
for pp=1:npages
    for ff=1:nFns
        wc3(ff,:,pp)=mod(ff+pp,2);
    end
end
simoptions_ac.whichcombos=wc3;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreFHorzPType_StructLeaves(SubAC);
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
extra=extra+sum(~ismember(pS3,pF3));
fprintf('WhichCombos AutoCorr PType, checkerboard over (fn, page): selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; selected leaves missing + unselected leaves present, this should be zero: %.3e \n',maxd,nanmis,missing+extra)
% (b) the age dimension is validated but not acted on: selecting only the first start age of (assets, unrestricted) gives the same output for it as the full run
wc3b=ones(nFns,N_j,npages); wc3b(1,2:end,1)=0;
simoptions_ac.whichcombos=wc3b;
SubAC=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pS3,vS3]=CoreFHorzPType_StructLeaves(SubAC);
maxd=0; nanmis=0;
for k=1:numel(pF3)
    [present,idx]=ismember(pF3{k},pS3);
    if present && isequal(size(vF3{k}),size(vS3{idx}))
        a=vF3{k}; b=vS3{idx}; nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
    else
        nanmis=nanmis+1;
    end
end
fprintf('WhichCombos AutoCorr PType, a partial start-age selection (not acted on yet) vs full run, every leaf, this should be zero: %.3e, NaN-pattern + leaf-set mismatches %i \n',maxd,nanmis+sum(~ismember(pS3,pF3)))
% (c) the two 2D forms expand to the same 3D selection; invalid sizes error
simoptions_ac.whichcombos=reshape(wc3(:,1,:),[nFns,npages]);
A=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
simoptions_ac.whichcombos=wc3;
B=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_ac);
[pA,vA]=CoreFHorzPType_StructLeaves(A); [pB,vB]=CoreFHorzPType_StructLeaves(B);
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
[pF4,vF4]=CoreFHorzPType_StructLeaves(FullCS); [pS4,vS4]=CoreFHorzPType_StructLeaves(SubCS);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF4)
    segs=strsplit(pF4{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix_ptype','CorrelationMatrix_ptype','CovarianceMatrix','CorrelationMatrix'}))
        continue % the matrices are checked below
    end
    ff=find(strcmp(FnNames,segs{1}));
    if ismember(segs{2},Names_i) % per-type leaf
        if numel(segs)==3 % .fn.ptype.Mean / .StdDeviation: present iff the function is evaluated, equal
            expectpresent=fnwanted(ff); on=1;
        else % .fn.ptype.CovarianceWith.fn2: present iff both evaluated, equal iff the pair is selected else NaN
            ff2=find(strcmp(FnNames,segs{4}));
            expectpresent=fnwanted(ff) && fnwanted(ff2); on=wcsym(ff,ff2);
        end
    else % grouped leaf: always present; equal iff selected (diagonal for Mean/StdDeviation/self, pair otherwise) else NaN
        expectpresent=true;
        if numel(segs)==2 % .fn.Mean, .fn.StdDeviation and the self-correlation field .fn.fn all follow the diagonal
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
extra=extra+sum(~ismember(pS4,pF4));
fprintf('WhichCombos CrossSection PType, pair selection: selected leaves vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected leaves not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra)
% grouped matrices: equal where selected (incl. the diagonal), NaN elsewhere; per-type matrices: over the evaluated functions, equal where selected or on the diagonal, NaN elsewhere
Mf=gather(FullCS.CovarianceMatrix); Ms=gather(SubCS.CovarianceMatrix); Rf=gather(FullCS.CorrelationMatrix); Rs=gather(SubCS.CorrelationMatrix);
on=(wcsym==1);
fprintf('WhichCombos CrossSection PType, grouped matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN, this should be zero: %.3e \n',max([abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on))))
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
% symmetric input equals upper-triangle input; invalid size errors
simoptions_cs.whichcombos=wcsym;
SymCS=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cs);
[pY,vY]=CoreFHorzPType_StructLeaves(SymCS);
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
[pF5,vF5]=CoreFHorzPType_StructLeaves(FullACS); [pS5,vS5]=CoreFHorzPType_StructLeaves(SubACS);
maxd=0; nanmis=0; missing=0; extra=0; notnan_off=0;
for k=1:numel(pF5)
    segs=strsplit(pF5{k},'.'); segs=segs(~cellfun(@isempty,segs));
    if any(strcmp(segs{1},{'CovarianceMatrix_ptype','CorrelationMatrix_ptype','CovarianceMatrix','CorrelationMatrix'}))
        continue
    end
    ff=find(strcmp(FnNames,segs{1}));
    if ismember(segs{2},Names_i)
        if numel(segs)==3
            expectpresent=fnwanted(ff); onmask=true(1,ngroups);
        else
            ff2=find(strcmp(FnNames,segs{4}));
            expectpresent=fnwanted(ff) && fnwanted(ff2); onmask=reshape(wcasym(ff,ff2,:)==1,[1,ngroups]);
        end
    else
        expectpresent=true;
        if numel(segs)==2
            onmask=reshape(wcasym(ff,ff,:)==1,[1,ngroups]);
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
extra=extra+sum(~ismember(pS5,pF5));
fprintf('WhichCombos AgeConditional CrossSection PType, pair x age-group selection: selected entries vs full run, this should be zero: %.3e, NaN-pattern mismatches %i; unselected entries not NaN + missing + extra, this should be zero: %.3e \n',maxd,nanmis,notnan_off+missing+extra)
Mf=gather(FullACS.CovarianceMatrix); Ms=gather(SubACS.CovarianceMatrix); Rf=gather(FullACS.CorrelationMatrix); Rs=gather(SubACS.CorrelationMatrix);
on=(wcasym==1);
fprintf('WhichCombos AgeConditional CrossSection PType, grouped matrices: selected entries vs full run, this should be zero: %.3e; unselected entries not NaN, this should be zero: %.3e \n',max([abs(Mf(on)-Ms(on));abs(Rf(on)-Rs(on))]),sum(~isnan(Ms(~on)))+sum(~isnan(Rs(~on))))
dmat=0; badmat=0;
sel=find(fnwanted);
for ii=1:N_i
    Mf_ii=gather(FullACS.CovarianceMatrix_ptype.(Names_i{ii})); Ms_ii=gather(SubACS.CovarianceMatrix_ptype.(Names_i{ii}));
    if ~isequal(size(Ms_ii),[numel(sel),numel(sel),ngroups])
        badmat=badmat+1;
    else
        onsub=false(numel(sel),numel(sel),ngroups);
        for kk=1:ngroups
            onsub(:,:,kk)=max(wcasym(sel,sel,kk),eye(numel(sel)))==1;
        end
        Mf_sub=Mf_ii(sel,sel,:);
        dmat=max(dmat,max(abs(Mf_sub(onsub)-Ms_ii(onsub))));
        badmat=badmat+sum(~isnan(Ms_ii(~onsub)));
    end
end
fprintf('WhichCombos AgeConditional CrossSection PType, per-type covariance matrices (shrunk to the evaluated functions): selected entries vs full run, this should be zero: %.3e; wrong size or unselected entries not NaN, this should be zero: %.3e \n',dmat,badmat)
% a 2D input equals its repmat over age groups; invalid size errors
simoptions_acs.whichcombos=wcs;
A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_acs);
simoptions_acs.whichcombos=repmat(wcs,[1,1,ngroups]);
B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_acs);
[pA,vA]=CoreFHorzPType_StructLeaves(A); [pB,vB]=CoreFHorzPType_StructLeaves(B);
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
output.AllStats_full=Full;
output.CrossSection_full=FullCS;

end
