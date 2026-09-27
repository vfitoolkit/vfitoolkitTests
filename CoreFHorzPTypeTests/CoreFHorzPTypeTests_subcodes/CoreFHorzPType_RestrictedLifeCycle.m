function output=CoreFHorzPType_RestrictedLifeCycle(n_a,n_z,N_j,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: LifeCycleProfiles_FHorz_Case1_PType with conditional restrictions, for agegroupings of one age and of several ages.
%
% Two defects were fixed in the restricted stats (2026-09-28):
%  A. Per-ptype, agegrouping of several ages: RestrictedStationaryDistVec is normalized to mass one at each age, and the
%     ages of the agegrouping were then pooled as they were, so each age counted equally rather than by its restricted mass.
%  B. Grouped (all ptypes together): the weight of each ptype was sum(AgeMasses.*restrictedsamplemass), which omits
%     ptweights (so a rare ptype counted as much as a common one) and counts the age weights twice. It is now
%     ptweights(ii)*sum(restrictedsamplemass) (the rule 2ab58e30 applied to AllStats PType).
% To catch both, the setup has unequal ptweights, unequal age weights (mewj), and a restriction ('rich') whose mass differs
% by age and by ptype.
%
% The checks are exact identities, not regression values:
%  1. Per-ptype restricted stats (every stat) equal LifeCycleProfiles_FHorz_Case1 of a solo solve of that type with the
%     same restriction and agegroupings (Case1 had defect A fixed in 617cf5a1 and is checked by CoreFHorzTests T6b).
%  2. Grouped restricted Mean equals E[X*1_R]/E[1_R], computed from the unrestricted grouped Means of X*1_R and 1_R
%     (every type has the same age weights, so the ptypes enter the unrestricted grouped stats by ptweights alone).
%  3. With a single agegrouping over all ages, the grouped restricted stats (every stat) equal the grouped restricted
%     AllStats of EvalFnOnAgentDist_AllStats_FHorz_Case1_PType (whose grouped restricted stats were fixed in 2ab58e30).
% It also checks the RestrictedSampleMass outputs, which were fixed at the same time (C): each ptype's StationaryDist puts
% mass mewj(j) on age j, so restrictedsamplemass is a mass within the ptype (it already includes mewj), and ByPType and
% Total used to multiply it by AgeMasses again. ByPType must equal the within-ptype shares of AllStats PType, Total its
% TotalAllPTypes, and ByAge (a share of the whole population, as in Case1) sum_ii ptweights(ii)*mewj(j)*P(R|ii,j).

Names_i={'low','high'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

Params.ptypeweights=[0.7; 0.3];
Params.kappa_j_pt=[Params.kappa_j; 1.5*Params.kappa_j]; % high type has the higher earnings, so more of it is 'rich'
Params.mewj=0.9.^(0:N_j-1); Params.mewj=Params.mewj/sum(Params.mewj); % unequal age weights

ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(aprime,a,z) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate_PT.richind=@(aprime,a,z) (a>0.5);
FnsToEvaluate_PT.assetsXrich=@(aprime,a,z) a*(a>0.5);
FnsToEvaluate_PT.earningsXrich=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z*(a>0.5);
FnsToEvaluate_NoPT.assets=@(aprime,a,z) a;
FnsToEvaluate_NoPT.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate_NoPT.richind=@(aprime,a,z) (a>0.5);
FnsToEvaluate_NoPT.assetsXrich=@(aprime,a,z) a*(a>0.5);
FnsToEvaluate_NoPT.earningsXrich=@(aprime,a,z,w,kappa_j) w*kappa_j*z*(a>0.5);
FnNames={'assets','earnings'}; % the stats that are checked (the X*1_R and 1_R functions are only there for identity 2)

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();
simoptions.whichstats=ones(1,7);
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% PType solve
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,struct());
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);

%% Solo solves of each type
Policy_solo=struct(); Dist_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [~,Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z,Params_ii,struct());
end

%% Premise of RestrictedSampleMass: each ptype's StationaryDist puts mass mewj(j) on age j
for ii=1:N_i
    nm=Names_i{ii};
    agemass=sum(reshape(gather(StationaryDist_PT.(nm)),[n_a*n_z,N_j]),1);
    fprintf('Restricted LifeCycle PType, StationaryDist.%s mass by age minus mewj, this should be zero: %.3e \n',nm,max(abs(agemass-Params.mewj)))
end

%% RestrictedSampleMass.ByPType, .Total and .ByAge
% restrictedsamplemass(ii,j) is the restricted mass of ptype ii at age j, as a share of ptype ii (it includes mewj(j)).
% So the share of ptype ii that is in the restriction is sum_j restrictedsamplemass(ii,j) (= AllStats RestrictedSampleMass.(ii)),
% and the share of the population is sum_ii ptweights(ii)*that (= AllStats RestrictedSampleMass.TotalAllPTypes).
% ByAge is checked against the unrestricted per-ptype Mean of the indicator, which is P(R|ii,j) (default agegroupings: one age each).
AgeStats_all=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
RSM=AgeStats_all.rich.RestrictedSampleMass;
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Restricted LifeCycle PType, RestrictedSampleMass.ByPType(%s) vs AllStats PType within-ptype share, this should be zero: %.3e \n',nm,abs(gather(RSM.ByPType(ii))-gather(AllStats_PT.rich.RestrictedSampleMass.(nm))))
end
fprintf('Restricted LifeCycle PType, RestrictedSampleMass.Total vs AllStats PType TotalAllPTypes, this should be zero: %.3e \n',abs(gather(RSM.Total)-gather(AllStats_PT.rich.RestrictedSampleMass.TotalAllPTypes)))
ByAge_direct=zeros(1,N_j);
for ii=1:N_i
    nm=Names_i{ii};
    ByAge_direct=ByAge_direct+Params.ptypeweights(ii)*Params.mewj.*gather(AgeStats_all.richind.(nm).Mean(1:N_j));
end
fprintf('Restricted LifeCycle PType, RestrictedSampleMass.ByAge vs sum_ii ptweights*mewj*P(R|ii,j), this should be zero: %.3e \n',max(abs(gather(RSM.ByAge(1:N_j))-ByAge_direct)))

%% For each agegroupings: one age each, five ages each, and one group of all ages
AgeGroupingsList={1:N_j, 1:5:N_j, 1};
AgeGroupingsNames={'one age each','five ages each','all ages in one group'};
for gg=1:length(AgeGroupingsList)
    simoptions_gg=simoptions;
    simoptions_gg.agegroupings=AgeGroupingsList{gg};
    ngroups=length(AgeGroupingsList{gg});
    AgeStats_PT=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_gg);
    AgeStats_solo=struct();
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params;
        Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
        AgeStats_solo.(nm)=LifeCycleProfiles_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions_gg);
    end

    % Identity 2: grouped restricted Mean = E[X*1_R]/E[1_R] from the unrestricted grouped Means (grouped arrays are sized by
    % the number of ages, so only the first ngroups entries are filled)
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        MeanXR=gather(AgeStats_PT.([fn,'Xrich']).Mean(1:ngroups));
        MeanR=gather(AgeStats_PT.richind.Mean(1:ngroups));
        MeanRestricted=gather(AgeStats_PT.rich.(fn).Mean(1:ngroups));
        dd=abs(MeanRestricted-MeanXR./MeanR); dd=dd(~isnan(dd));
        fprintf('Restricted LifeCycle PType (%s), grouped rich %s Mean vs E[X*1_R]/E[1_R], this should be zero: %.3e, NaN-pattern mismatches %i \n',AgeGroupingsNames{gg},fn,max([dd,0]),sum(isnan(MeanRestricted)~=isnan(MeanXR./MeanR)))
        for ii=1:N_i
            nm=Names_i{ii};
            MeanXR=gather(AgeStats_PT.([fn,'Xrich']).(nm).Mean(1:ngroups));
            MeanR=gather(AgeStats_PT.richind.(nm).Mean(1:ngroups));
            MeanRestricted=gather(AgeStats_PT.rich.(fn).(nm).Mean(1:ngroups));
            dd=abs(MeanRestricted-MeanXR./MeanR); dd=dd(~isnan(dd));
            fprintf('Restricted LifeCycle PType (%s), rich %s Mean (type %s) vs E[X*1_R]/E[1_R], this should be zero: %.3e, NaN-pattern mismatches %i \n',AgeGroupingsNames{gg},fn,nm,max([dd,0]),sum(isnan(MeanRestricted)~=isnan(MeanXR./MeanR)))
        end
    end

    % Identity 1 (per-ptype restricted = solo Case1 restricted) and, for one group of all ages, identity 3 (grouped restricted = AllStats PType grouped restricted)
    pairs=cell(0,3);
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        for ii=1:N_i
            nm=Names_i{ii};
            pairs(end+1,:)={sprintf('Restricted LifeCycle PType (%s), rich %s (type %s) vs solo Case1',AgeGroupingsNames{gg},fn,nm), AgeStats_PT.rich.(fn).(nm), AgeStats_solo.(nm).rich.(fn)}; %#ok<AGROW>
        end
        if ngroups==1
            pairs(end+1,:)={sprintf('Restricted LifeCycle PType (%s), grouped rich %s vs AllStats PType',AgeGroupingsNames{gg},fn), AgeStats_PT.rich.(fn), AllStats_PT.rich.(fn)}; %#ok<AGROW>
        end
    end
    % Compare every statistic of each (label, Test, Ref) pair; Ref supplies the list of stats, and only the first ngroups
    % columns of Test are used (the grouped arrays have a column per age). Per-ptype substructures of Ref are skipped.
    % NaN is compared as a pattern, as a NaN printed as the value of a check is not read by CoreSummary.
    for pp=1:size(pairs,1)
        Test=pairs{pp,2};
        Ref=pairs{pp,3};
        statnames=fieldnames(Ref);
        for ss=1:length(statnames)
            sn=statnames{ss};
            if ischar(Ref.(sn)) || iscell(Ref.(sn)) || (isstruct(Ref.(sn)) && ~strcmp(sn,'MoreInequality'))
                continue
            end
            if strcmp(sn,'MoreInequality')
                subnames=fieldnames(Ref.MoreInequality);
                cmpnames=cell(length(subnames),1); xs=cmpnames; ys=cmpnames;
                for s2=1:length(subnames)
                    cmpnames{s2}=['MoreInequality.',subnames{s2}];
                    ys{s2}=Ref.MoreInequality.(subnames{s2});
                    if isfield(Test,'MoreInequality') && isfield(Test.MoreInequality,subnames{s2})
                        xs{s2}=Test.MoreInequality.(subnames{s2});
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
                if ~isempty(x) && size(x,2)>=ngroups
                    x=x(:,1:ngroups);
                end
                if numel(x)~=numel(y)
                    fprintf('%s, %s missing or wrongly sized, this should be zero: 1 \n',pairs{pp,1},cmpnames{s2})
                else
                    x=gather(x(:)); y=gather(y(:));
                    dd=abs(x-y); dd=dd(~isnan(dd));
                    fprintf('%s, %s, this should be zero: %.3e, NaN-pattern mismatches %i \n',pairs{pp,1},cmpnames{s2},max([dd;0]),sum(isnan(x)~=isnan(y)))
                end
            end
        end
    end
end

output=struct();

end
