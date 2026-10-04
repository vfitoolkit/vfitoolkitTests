function output=CoreFHorzPType_WhichCombos(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: simoptions.whichcombos, and simoptions.whichstats per combination, in LifeCycleProfiles_FHorz_Case1_PType.
%
% whichcombos is [nFns, nAgeGroups, 1+nRestr] (page 1 unrestricted, pages 2:end the restrictions in fieldnames order; the third
% dimension is dropped when there are no restrictions). Ones are computed, zeros skipped; skipped entries stay NaN; the output has
% the same fields either way; RestrictedSampleMass is always filled.
%
% Checks (all against a full run with the default, which must itself equal an explicit all-ones whichcombos):
%  1. Explicit all-ones equals the default, every field.
%  2. Checkerboard selection: every selected (fn, age group, page) entry, per-ptype and grouped, every statistic, equals the full
%     run exactly (same code path); every unselected entry is NaN.
%  3. An all-zero function row and an all-zero restriction page: the fields are still present, all NaN; RestrictedSampleMass of the
%     zeroed restriction equals the full run.
%  4. A 2D whichcombos with restrictions present equals the 3D repmat of it, NaN pattern included.
%  5. Age groups of five ages: the second dimension is the number of age groups, checkerboard vs full run as in 2.
%  6. Wrong size, and entries other than 0/1, error.
%  7. Runtime of the full run vs a one-row selection is printed (information, not a check).
% Without restrictions (2D whichcombos) check 2 is repeated.
% whichstats per combination ([nFns, nAgeGroups, 1+nRestr, 7]): (a) all-ones equals the full run; (b) a per-combination pattern:
% entries with a statistic on equal the full run, off ones are NaN, every field present; (c) a statistic off at every age of one
% (fn, page) removes its field, the rest equals the full run; (d) a [nFns, nAgeGroups, 7] input equals its 4D repmat; (e) a
% combination with every statistic off is skipped like a whichcombos zero; (f) mixed Gini 3/1 across ages keeps the Lorenz curve
% where asked; (g) wrong sizes error.

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

%% Full run (default whichcombos) and explicit all-ones
tic;
Full=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
time_full=toc;
simoptions_ones=simoptions;
simoptions_ones.whichcombos=ones(nFns,N_j,npages);
Ones=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ones);

% Walk every leaf of the two outputs (page level -> fn -> grouped stats, per-ptype stats, MoreInequality) and compare
maxd=0; nanmis=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=Full.(FnNames{ff}); B=Ones.(FnNames{ff});
        else
            A=Full.(RestrNames{pp-1}).(FnNames{ff}); B=Ones.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i]; % '' is the grouped level, then each ptype
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    for mm=1:length(minames)
                        a=gather(SA.MoreInequality.(minames{mm})); b=gather(SB.MoreInequality.(minames{mm}));
                        nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:)));
                        d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    end
                else
                    a=gather(SA.(statnames{ss})); b=gather(SB.(statnames{ss}));
                    nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:)));
                    d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                end
            end
        end
    end
end
fprintf('WhichCombos, explicit all-ones vs default, every field, this should be zero: %.3e, NaN-pattern mismatches %i \n',maxd,nanmis)

%% Checkerboard selection, each age its own group (nAgeGroups=N_j)
whichcombos=zeros(nFns,N_j,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:N_j
            whichcombos(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
simoptions_cb=simoptions;
simoptions_cb.whichcombos=whichcombos;
tic;
Sub=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_cb);
time_cb=toc;

% For every (page, fn): selected age groups equal the full run exactly, unselected ones are NaN; per-ptype and grouped alike
for pp=1:npages
    if pp==1
        pagename='unrestricted';
    else
        pagename=RestrNames{pp-1};
    end
    for ff=1:nFns
        if pp==1
            A=Full.(FnNames{ff}); B=Sub.(FnNames{ff});
        else
            A=Full.(RestrNames{pp-1}).(FnNames{ff}); B=Sub.(RestrNames{pp-1}).(FnNames{ff});
        end
        maxd_sel=0; nonnan_unsel=0; nleaves=0;
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    leavesA=cell(length(minames),1); leavesB=cell(length(minames),1);
                    for mm=1:length(minames)
                        leavesA{mm}=gather(SA.MoreInequality.(minames{mm})); leavesB{mm}=gather(SB.MoreInequality.(minames{mm}));
                    end
                else
                    leavesA={gather(SA.(statnames{ss}))}; leavesB={gather(SB.(statnames{ss}))};
                end
                for kk=1:length(leavesA)
                    a=leavesA{kk}; b=leavesB{kk}; % columns are age groups (a row vector or an npoints/nquantiles-by-ngroups matrix)
                    nleaves=nleaves+1;
                    for jj=1:N_j
                        if whichcombos(ff,jj,pp)==1
                            d=abs(a(:,jj)-b(:,jj)); d=d(~isnan(d)); maxd_sel=max([maxd_sel;d]);
                            nonnan_unsel=nonnan_unsel+sum(isnan(a(:,jj))~=isnan(b(:,jj))); % a selected entry must have the full run's NaN pattern
                        else
                            nonnan_unsel=nonnan_unsel+sum(~isnan(b(:,jj)));
                        end
                    end
                end
            end
        end
        fprintf('WhichCombos checkerboard, %s %s (%i leaves): selected entries vs full run, this should be zero: %.3e; unselected entries that are not NaN (and selected NaN-pattern mismatches), this should be zero: %.3e \n',pagename,FnNames{ff},nleaves,maxd_sel,nonnan_unsel)
    end
end
% RestrictedSampleMass is always filled
for rr=1:nRestr
    a=gather(Full.(RestrNames{rr}).RestrictedSampleMass.ByAge); b=gather(Sub.(RestrNames{rr}).RestrictedSampleMass.ByAge);
    fprintf('WhichCombos checkerboard, RestrictedSampleMass.ByAge of %s vs full run, this should be zero: %.3e \n',RestrNames{rr},max(abs(a-b)))
end
fprintf('WhichCombos, runtime: full run %.3f s, checkerboard (half the combinations) %.3f s \n',time_full,time_cb)

%% All-zero function row (earnings) and all-zero restriction page (poor)
whichcombos=ones(nFns,N_j,npages);
whichcombos(2,:,:)=0; % earnings: nothing
whichcombos(:,:,1+2)=0; % poor: nothing
simoptions_z=simoptions;
simoptions_z.whichcombos=whichcombos;
tic;
Zrow=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_z);
time_z=toc;
fprintf('WhichCombos all-zero row/page, earnings field present at top level and under each restriction, this should be zero: %.3e \n',double(~isfield(Zrow,'earnings'))+double(~isfield(Zrow.rich,'earnings'))+double(~isfield(Zrow.poor,'earnings')))
% everything under Zrow.earnings and everything under Zrow.poor (except RestrictedSampleMass) must be NaN
nonnan=0;
checkstructs={Zrow.earnings,Zrow.rich.earnings,Zrow.poor.assets,Zrow.poor.earnings,Zrow.poor.richind};
for cc=1:length(checkstructs)
    S=checkstructs{cc};
    levelnames=[{''},Names_i];
    for ll=1:length(levelnames)
        if isempty(levelnames{ll})
            SA=S;
        else
            SA=S.(levelnames{ll});
        end
        statnames=fieldnames(SA);
        for ss=1:length(statnames)
            if ismember(statnames{ss},Names_i)
                continue
            elseif strcmp(statnames{ss},'MoreInequality')
                minames=fieldnames(SA.MoreInequality);
                for mm=1:length(minames)
                    a=gather(SA.MoreInequality.(minames{mm})); nonnan=nonnan+sum(~isnan(a(:)));
                end
            else
                a=gather(SA.(statnames{ss})); nonnan=nonnan+sum(~isnan(a(:)));
            end
        end
    end
end
fprintf('WhichCombos all-zero row/page, non-NaN entries under earnings and under poor, this should be zero: %.3e \n',nonnan)
a=gather(Full.poor.RestrictedSampleMass.ByAge); b=gather(Zrow.poor.RestrictedSampleMass.ByAge);
fprintf('WhichCombos all-zero page, RestrictedSampleMass.ByAge of poor vs full run, this should be zero: %.3e \n',max(abs(a-b)))
fprintf('WhichCombos all-zero page, RestrictedSampleMass.Total of poor vs full run, this should be zero: %.3e \n',abs(gather(Full.poor.RestrictedSampleMass.Total)-gather(Zrow.poor.RestrictedSampleMass.Total)))
% the remaining (assets, richind; unrestricted and rich) equal the full run
maxd=0; nanmis=0;
for pp=1:2 % unrestricted, rich
    for ff=[1,3]
        if pp==1
            A=Full.(FnNames{ff}); B=Zrow.(FnNames{ff});
        else
            A=Full.rich.(FnNames{ff}); B=Zrow.rich.(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    for mm=1:length(minames)
                        a=gather(SA.MoreInequality.(minames{mm})); b=gather(SB.MoreInequality.(minames{mm}));
                        nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    end
                else
                    a=gather(SA.(statnames{ss})); b=gather(SB.(statnames{ss}));
                    nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                end
            end
        end
    end
end
fprintf('WhichCombos all-zero row/page, the other combinations vs full run, this should be zero: %.3e, NaN-pattern mismatches %i \n',maxd,nanmis)
fprintf('WhichCombos, runtime: full run %.3f s, one fn and one restriction dropped %.3f s \n',time_full,time_z)

%% 2D whichcombos with restrictions present broadcasts over the pages
wc2=zeros(nFns,N_j);
for ff=1:nFns
    for jj=1:N_j
        wc2(ff,jj)=mod(ff+jj,2);
    end
end
simoptions_2d=simoptions; simoptions_2d.whichcombos=wc2;
simoptions_3d=simoptions; simoptions_3d.whichcombos=repmat(wc2,[1,1,npages]);
Two=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_2d);
Three=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_3d);
maxd=0; nanmis=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=Two.(FnNames{ff}); B=Three.(FnNames{ff});
        else
            A=Two.(RestrNames{pp-1}).(FnNames{ff}); B=Three.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    for mm=1:length(minames)
                        a=gather(SA.MoreInequality.(minames{mm})); b=gather(SB.MoreInequality.(minames{mm}));
                        nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    end
                else
                    a=gather(SA.(statnames{ss})); b=gather(SB.(statnames{ss}));
                    nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                end
            end
        end
    end
end
fprintf('WhichCombos, 2D input vs its 3D repmat, every field, this should be zero: %.3e, NaN-pattern mismatches %i \n',maxd,nanmis)

%% Age groups of five ages: the second dimension is the number of age groups
simoptions_ag=simoptions;
simoptions_ag.agegroupings=1:5:N_j;
ngroups=length(simoptions_ag.agegroupings);
Full_ag=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
whichcombos=zeros(nFns,ngroups,npages);
for pp=1:npages
    for ff=1:nFns
        for jj=1:ngroups
            whichcombos(ff,jj,pp)=mod(ff+jj+pp,2);
        end
    end
end
simoptions_ag.whichcombos=whichcombos;
Sub_ag=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_ag);
maxd_sel=0; nonnan_unsel=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=Full_ag.(FnNames{ff}); B=Sub_ag.(FnNames{ff});
        else
            A=Full_ag.(RestrNames{pp-1}).(FnNames{ff}); B=Sub_ag.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    leavesA=cell(length(minames),1); leavesB=cell(length(minames),1);
                    for mm=1:length(minames)
                        leavesA{mm}=gather(SA.MoreInequality.(minames{mm})); leavesB{mm}=gather(SB.MoreInequality.(minames{mm}));
                    end
                else
                    leavesA={gather(SA.(statnames{ss}))}; leavesB={gather(SB.(statnames{ss}))};
                end
                for kk=1:length(leavesA)
                    a=leavesA{kk}; b=leavesB{kk};
                    for jj=1:ngroups % (grouped arrays are sized by the number of ages, only the first ngroups columns are filled; the rest are NaN in both)
                        if whichcombos(ff,jj,pp)==1
                            d=abs(a(:,jj)-b(:,jj)); d=d(~isnan(d)); maxd_sel=max([maxd_sel;d]);
                            nonnan_unsel=nonnan_unsel+sum(isnan(a(:,jj))~=isnan(b(:,jj)));
                        else
                            nonnan_unsel=nonnan_unsel+sum(~isnan(b(:,jj)));
                        end
                    end
                end
            end
        end
    end
end
fprintf('WhichCombos checkerboard with age groups of five ages: selected entries vs full run, this should be zero: %.3e; unselected entries that are not NaN, this should be zero: %.3e \n',maxd_sel,nonnan_unsel)

%% Without restrictions: 2D whichcombos
simoptions_nr=struct();
simoptions_nr.whichstats=ones(1,7);
Full_nr=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_nr);
simoptions_nr.whichcombos=wc2;
Sub_nr=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_nr);
maxd_sel=0; nonnan_unsel=0;
for ff=1:nFns
    A=Full_nr.(FnNames{ff}); B=Sub_nr.(FnNames{ff});
    levelnames=[{''},Names_i];
    for ll=1:length(levelnames)
        if isempty(levelnames{ll})
            SA=A; SB=B;
        else
            SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
        end
        statnames=fieldnames(SA);
        for ss=1:length(statnames)
            if ismember(statnames{ss},Names_i)
                continue
            elseif strcmp(statnames{ss},'MoreInequality')
                minames=fieldnames(SA.MoreInequality);
                leavesA=cell(length(minames),1); leavesB=cell(length(minames),1);
                for mm=1:length(minames)
                    leavesA{mm}=gather(SA.MoreInequality.(minames{mm})); leavesB{mm}=gather(SB.MoreInequality.(minames{mm}));
                end
            else
                leavesA={gather(SA.(statnames{ss}))}; leavesB={gather(SB.(statnames{ss}))};
            end
            for kk=1:length(leavesA)
                a=leavesA{kk}; b=leavesB{kk};
                for jj=1:N_j
                    if wc2(ff,jj)==1
                        d=abs(a(:,jj)-b(:,jj)); d=d(~isnan(d)); maxd_sel=max([maxd_sel;d]);
                        nonnan_unsel=nonnan_unsel+sum(isnan(a(:,jj))~=isnan(b(:,jj)));
                    else
                        nonnan_unsel=nonnan_unsel+sum(~isnan(b(:,jj)));
                    end
                end
            end
        end
    end
end
fprintf('WhichCombos without restrictions (2D): selected entries vs full run, this should be zero: %.3e; unselected entries that are not NaN, this should be zero: %.3e \n',maxd_sel,nonnan_unsel)

%% Wrong size and non-0/1 entries must error
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
simoptions_bad=simoptions_nr; simoptions_bad.whichcombos=ones(nFns,N_j,2); % no restrictions, but a 3D input
try
    LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_bad);
catch
    errored=errored+1;
end
fprintf('WhichCombos, three invalid inputs (wrong age dimension, an entry of 2, 3D without restrictions) must all error, this should be zero: %.3e \n',3-errored)

%% whichstats per combination: [nFns, nAgeGroups, 1+nRestr, 7]
% (a) a 4D whichstats of all ones equals the full run (which used the 1x7 vector of ones)
simoptions_w=simoptions;
simoptions_w.whichstats=ones(nFns,N_j,npages,7);
W1=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
maxd=0; nanmis=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=Full.(FnNames{ff}); B=W1.(FnNames{ff});
        else
            A=Full.(RestrNames{pp-1}).(FnNames{ff}); B=W1.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    for mm=1:length(minames)
                        a=gather(SA.MoreInequality.(minames{mm})); b=gather(SB.MoreInequality.(minames{mm}));
                        nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    end
                else
                    a=gather(SA.(statnames{ss})); b=gather(SB.(statnames{ss}));
                    nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                end
            end
        end
    end
end
fprintf('WhichStats 4D all-ones vs the full run, every field, this should be zero: %.3e, NaN-pattern mismatches %i \n',maxd,nanmis)

% (b) per-combination pattern: statistic kk is on for (fn ff, age jj, page pp) iff mod(ff+jj+pp+kk,2)==1, so every statistic is on for
% about half the ages of every (fn, page) and every field is present. Entries at ages where the statistic is on equal the full run,
% the rest are NaN. Field -> statistic: Mean 1, Median 2, RatioMeanToMedian 1 and 2, Variance/StdDeviation 3, Gini/LorenzCurve 4,
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
W2=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
maxd_on=0; bad_off=0; missingfields=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=Full.(FnNames{ff}); B=W2.(FnNames{ff});
        else
            A=Full.(RestrNames{pp-1}).(FnNames{ff}); B=W2.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i)
                    continue
                end
                switch statnames{ss}
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
                    case 'MoreInequality'
                        kk=7;
                end
                if ~isfield(SB,statnames{ss})
                    missingfields=missingfields+1;
                    continue
                end
                if strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    leavesA=cell(length(minames),1); leavesB=cell(length(minames),1);
                    for mm=1:length(minames)
                        leavesA{mm}=gather(SA.MoreInequality.(minames{mm})); leavesB{mm}=gather(SB.MoreInequality.(minames{mm}));
                    end
                else
                    leavesA={gather(SA.(statnames{ss}))}; leavesB={gather(SB.(statnames{ss}))};
                end
                for ll2=1:length(leavesA)
                    a=leavesA{ll2}; b=leavesB{ll2};
                    for jj=1:N_j
                        if all(wsarr(ff,jj,pp,kk)==1)
                            d=abs(a(:,jj)-b(:,jj)); d=d(~isnan(d)); maxd_on=max([maxd_on;d]);
                            bad_off=bad_off+sum(isnan(a(:,jj))~=isnan(b(:,jj)));
                        else
                            bad_off=bad_off+sum(~isnan(b(:,jj)));
                        end
                    end
                end
            end
        end
    end
end
fprintf('WhichStats per-combination pattern: entries with the statistic on vs the full run, this should be zero: %.3e; entries with it off that are not NaN (and on-entry NaN-pattern mismatches), this should be zero: %.3e; fields missing, this should be zero: %.3e \n',maxd_on,bad_off,missingfields)

% (c) a statistic off for every age of one (fn, page) removes its field: Median off for assets unrestricted (so RatioMeanToMedian too),
% QuantileCutoffs/QuantileMeans off for rich earnings; everything else equals the full run
wsarr3=ones(nFns,N_j,npages,7);
wsarr3(1,:,1,2)=0;
wsarr3(2,:,2,6)=0;
simoptions_w.whichstats=wsarr3;
W3=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
present=isfield(W3.assets,'Median')+isfield(W3.assets,'RatioMeanToMedian')+isfield(W3.assets.low,'Median')+isfield(W3.assets.high,'RatioMeanToMedian')+isfield(W3.rich.earnings,'QuantileCutoffs')+isfield(W3.rich.earnings.low,'QuantileMeans')+isfield(W3.rich.earnings.high,'QuantileCutoffs');
fprintf('WhichStats all-ages-off statistic: fields that must be absent (7 looked for), this should be zero: %.3e \n',present)
maxd=0; nanmis=0; nleaves=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=Full.(FnNames{ff}); B=W3.(FnNames{ff});
        else
            A=Full.(RestrNames{pp-1}).(FnNames{ff}); B=W3.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i) || ~isfield(SB,statnames{ss})
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    for mm=1:length(minames)
                        a=gather(SA.MoreInequality.(minames{mm})); b=gather(SB.MoreInequality.(minames{mm}));
                        nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]); nleaves=nleaves+1;
                    end
                else
                    a=gather(SA.(statnames{ss})); b=gather(SB.(statnames{ss}));
                    nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]); nleaves=nleaves+1;
                end
            end
        end
    end
end
fprintf('WhichStats all-ages-off statistic: the other fields (%i leaves) vs the full run, this should be zero: %.3e, NaN-pattern mismatches %i \n',nleaves,maxd,nanmis)

% (d) [nFns, nAgeGroups, 7] with restrictions present equals its 4D repmat over the pages
simoptions_w.whichstats=reshape(wsarr(:,:,1,:),[nFns,N_j,7]);
W4=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
simoptions_w.whichstats=repmat(reshape(wsarr(:,:,1,:),[nFns,N_j,1,7]),[1,1,npages,1]);
W5=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
maxd=0; nanmis=0; fieldmis=0;
for pp=1:npages
    for ff=1:nFns
        if pp==1
            A=W4.(FnNames{ff}); B=W5.(FnNames{ff});
        else
            A=W4.(RestrNames{pp-1}).(FnNames{ff}); B=W5.(RestrNames{pp-1}).(FnNames{ff});
        end
        levelnames=[{''},Names_i];
        for ll=1:length(levelnames)
            if isempty(levelnames{ll})
                SA=A; SB=B;
            else
                SA=A.(levelnames{ll}); SB=B.(levelnames{ll});
            end
            fieldmis=fieldmis+numel(setxor(fieldnames(SA),fieldnames(SB)));
            statnames=fieldnames(SA);
            for ss=1:length(statnames)
                if ismember(statnames{ss},Names_i) || ~isfield(SB,statnames{ss})
                    continue
                elseif strcmp(statnames{ss},'MoreInequality')
                    minames=fieldnames(SA.MoreInequality);
                    for mm=1:length(minames)
                        a=gather(SA.MoreInequality.(minames{mm})); b=gather(SB.MoreInequality.(minames{mm}));
                        nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                    end
                else
                    a=gather(SA.(statnames{ss})); b=gather(SB.(statnames{ss}));
                    nanmis=nanmis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                end
            end
        end
    end
end
fprintf('WhichStats 3D input vs its 4D repmat, every field, this should be zero: %.3e, NaN-pattern mismatches + field-set mismatches %i \n',maxd,nanmis+fieldmis)

% (e) a combination with every statistic off is skipped like a whichcombos zero: richind unrestricted all off leaves no statistic
% fields under richind (grouped) or under richind.(ptype); rich.richind is unaffected
wsarr5=ones(nFns,N_j,npages,7);
wsarr5(3,:,1,:)=0;
simoptions_w.whichstats=wsarr5;
W6=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
if isfield(W6,'richind') % with nothing asked for, the richind struct need not exist at all; if it does it must hold no statistic fields
    nstatfields=numel(setdiff(fieldnames(W6.richind),Names_i));
    for ii=1:N_i
        if isfield(W6.richind,Names_i{ii})
            nstatfields=nstatfields+numel(fieldnames(W6.richind.(Names_i{ii})));
        end
    end
else
    nstatfields=0;
end
fprintf('WhichStats all-off combination (richind unrestricted): statistic fields present, this should be zero: %.3e \n',nstatfields)
a=gather(Full.rich.richind.Mean); b=gather(W6.rich.richind.Mean);
fprintf('WhichStats all-off combination: rich richind grouped Mean vs the full run, this should be zero: %.3e \n',max(abs(a-b)))

% (f) mixed Gini settings across ages: whichstats(4)=3 (Gini only) at odd ages and 1 (Gini and Lorenz curve) at even ages, for assets
% unrestricted: LorenzCurve is present, NaN at the odd ages and equal to the full run at the even ages; Gini equals the full run at all ages
wsarr6=ones(nFns,N_j,npages,7);
for jj=1:2:N_j
    wsarr6(1,jj,1,4)=3;
end
simoptions_w.whichstats=wsarr6;
W7=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_w);
fprintf('WhichStats mixed Gini 3/1 across ages: assets LorenzCurve field present (grouped and both ptypes), this should be zero: %.3e \n',3-isfield(W7.assets,'LorenzCurve')-isfield(W7.assets.low,'LorenzCurve')-isfield(W7.assets.high,'LorenzCurve'))
LC_full=gather(Full.assets.LorenzCurve); LC=gather(W7.assets.LorenzCurve);
fprintf('WhichStats mixed Gini 3/1 across ages: grouped assets LorenzCurve at even ages vs the full run, this should be zero: %.3e; non-NaN entries at odd ages, this should be zero: %.3e \n',max(max(abs(LC_full(:,2:2:N_j)-LC(:,2:2:N_j)))),sum(sum(~isnan(LC(:,1:2:N_j)))))
LC_full=gather(Full.assets.low.LorenzCurve); LC=gather(W7.assets.low.LorenzCurve);
fprintf('WhichStats mixed Gini 3/1 across ages: ptype low assets LorenzCurve at even ages vs the full run, this should be zero: %.3e; non-NaN entries at odd ages, this should be zero: %.3e \n',max(max(abs(LC_full(:,2:2:N_j)-LC(:,2:2:N_j)))),sum(sum(~isnan(LC(:,1:2:N_j)))))
fprintf('WhichStats mixed Gini 3/1 across ages: grouped assets Gini vs the full run, this should be zero: %.3e \n',max(abs(gather(Full.assets.Gini)-gather(W7.assets.Gini))))

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
fprintf('WhichStats, two invalid sizes must both error, this should be zero: %.3e \n',2-errored)

output=struct();
output.Full=Full;
output.Sub=Sub;

end
