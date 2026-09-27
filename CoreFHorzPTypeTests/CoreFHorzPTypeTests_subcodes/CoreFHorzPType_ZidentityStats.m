function output=CoreFHorzPType_ZidentityStats(n_a,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test: every statistic of the grouped (pooled-across-ptypes) AllStats and LifeCycleProfiles.
% Uses the construction of part 2: a PType model with N_i types is the same model as one without
% PType where z encodes the type and pi_z is the identity. Grouping across ptypes is then exactly
% pooling across z, so the grouped PType stats must equal the plain no-PType stats, for every
% statistic (Median, quantiles, Lorenz curve, Gini, StdDeviation, MoreInequality, ...), not just
% the Mean. That is the only exact reference for the grouped stats beyond the Mean.
%
% The no-PType model also gives exact references for the per-type stats and for the conditional
% restrictions, via conditional restrictions on z:
%   type<ii>      z==sigma_pt(ii)                  -> the PType per-type stats of type ii
%   rich          a>0.5                            -> the PType grouped restricted stats
%   richtype<ii>  a>0.5 and z==sigma_pt(ii)        -> the PType per-type restricted stats
% and the restricted sample masses: RestrictedSampleMass of richtype<ii> (a population share) is
% ptweight(ii) times the PType per-type RestrictedSampleMass (a within-type share).
%
% Also, all against the default-options PType run:
%   LifeCycleProfiles with agegroupings (5-year bins)
%   LifeCycleProfiles with lowmemory=1 (outer loop over FnsToEvaluate rather than ptypes)
%   LifeCycleProfiles and AllStats with groupptypesforstats=0 (per-type unchanged; LifeCycleProfiles
%     must drop the grouped fields; AllStats always reports both, so it must be unchanged)
%   AllStats with a subset of whichstats (the stats that are computed must be unchanged)
%
% FnsToEvaluate.earnings does not depend on a, so every type has many tied values at each age,
% which is what the pooling (unique + accumarray) has to get right.

Params.sigma_pt=[2; 3];
ptw=Params.(PTypeDistParamNames{1});

n_d=0;
d_grid=[];

vfoptions=struct();

%% (A) PType, no z
ReturnFn_A=@(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension);
n_z_A=0;
z_grid_A=[];
pi_z_A=[];

jequaloneDist_A=zeros(n_a,1,'gpuArray');
jequaloneDist_A(1)=1; % no assets

FnsToEvaluate_A.assets=@(aprime,a) a;
FnsToEvaluate_A.earnings=@(aprime,a,w,kappa_j) w*kappa_j;
FnNames=fieldnames(FnsToEvaluate_A);

simoptions_A=struct();
simoptions_A.conditionalrestrictions.rich=@(aprime,a) (a>0.5);

[~,Policy_A]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,ReturnFn_A,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_A=StationaryDist_Case1_FHorz_PType(jequaloneDist_A,AgeWeightParamNames,PTypeDistParamNames,Policy_A,n_d,n_a,n_z_A,N_j,N_i,pi_z_A,Params,struct());
names_A=fieldnames(Policy_A);

AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeStats_A=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);

simoptions_A_ag=simoptions_A;
simoptions_A_ag.agegroupings=1:5:N_j;
AgeStats_A_ag=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ag);

% lowmemory=1 cannot be combined with conditionalrestrictions, so compare it against a run without them
simoptions_A_plain=struct();
AgeStats_A_plain=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_plain);
simoptions_A_lowmem=struct();
simoptions_A_lowmem.lowmemory=1;
AgeStats_A_lowmem=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_lowmem);

simoptions_A_ng=simoptions_A;
simoptions_A_ng.groupptypesforstats=0;
AllStats_A_ng=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ng);
AgeStats_A_ng=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ng);

simoptions_A_ws=simoptions_A;
simoptions_A_ws.whichstats=[1,0,1,0,1,0,0]; % mean, std dev/variance, min/max only
AllStats_A_ws=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ws);

%% (B) No PType, z encodes the type with identity transitions
n_z_B=N_i;
z_grid_B=Params.sigma_pt; % [2; 3]
pi_z_B=eye(N_i,'gpuArray');

ReturnFn_B=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) ...
    ReturnFn_zAsSigma(aprime,a,z,r,w,kappa_j,agej,Jr,pension);

jequaloneDist_B=zeros(n_a,n_z_B,'gpuArray');
jequaloneDist_B(1,:)=reshape(ptw,1,n_z_B);

FnsToEvaluate_B.assets=@(aprime,a,z) a;
FnsToEvaluate_B.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j;

% type<ii> and richtype<ii>: z is sigma_pt, which is 2 or 3, so 2.5 separates the two types
simoptions_B=struct();
simoptions_B.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions_B.conditionalrestrictions.type1=@(aprime,a,z) (z<2.5);
simoptions_B.conditionalrestrictions.type2=@(aprime,a,z) (z>2.5);
simoptions_B.conditionalrestrictions.richtype1=@(aprime,a,z) (a>0.5)*(z<2.5);
simoptions_B.conditionalrestrictions.richtype2=@(aprime,a,z) (a>0.5)*(z>2.5);

[~,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,n_z_B,N_j,pi_z_B,Params,struct());

AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B);
AgeStats_B=LifeCycleProfiles_FHorz_Case1(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B);

simoptions_B_ag=simoptions_B;
simoptions_B_ag.agegroupings=1:5:N_j;
AgeStats_B_ag=LifeCycleProfiles_FHorz_Case1(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B_ag);

%% Restricted sample masses
fprintf('z-identity stats, AllStats rich RestrictedSampleMass TotalAllPTypes, this should be zero: %.3e \n',abs(AllStats_A.rich.RestrictedSampleMass.TotalAllPTypes-AllStats_B.rich.RestrictedSampleMass))
for ii=1:N_i
    fprintf('z-identity stats, AllStats rich RestrictedSampleMass (type %s) times ptweight, this should be zero: %.3e \n',names_A{ii},abs(ptw(ii)*AllStats_A.rich.RestrictedSampleMass.(names_A{ii})-AllStats_B.(['richtype',num2str(ii)]).RestrictedSampleMass))
end
x=gather(AgeStats_A.rich.RestrictedSampleMass.ByAge(:)); y=gather(AgeStats_B.rich.RestrictedSampleMass(:));
dd=abs(x-y); dd=dd(~isnan(dd));
fprintf('z-identity stats, LifeCycleProfiles rich RestrictedSampleMass ByAge, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([dd;0]),sum(isnan(x)~=isnan(y)))
for ii=1:N_i
    x=gather(ptw(ii)*AgeStats_A.rich.RestrictedSampleMass.(names_A{ii})(:)); y=gather(AgeStats_B.(['richtype',num2str(ii)]).RestrictedSampleMass(:));
    dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('z-identity stats, LifeCycleProfiles rich RestrictedSampleMass (type %s) times ptweight, this should be zero: %.3e, NaN-pattern mismatches %i \n',names_A{ii},max([dd;0]),sum(isnan(x)~=isnan(y)))
end

%% groupptypesforstats=0 must drop the grouped LifeCycleProfiles fields
for ff=1:length(FnNames)
    fprintf('z-identity stats, LifeCycleProfiles %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',FnNames{ff},isfield(AgeStats_A_ng.(FnNames{ff}),'Mean'))
end

%% Every statistic of each (label, Test, Ref) pair
pairs=cell(0,3);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    % AllStats
    pairs(end+1,:)={sprintf('z-identity stats, AllStats %s grouped',fn), AllStats_A.(fn), AllStats_B.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity stats, AllStats rich %s grouped',fn), AllStats_A.rich.(fn), AllStats_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity stats, AllStats %s grouped, groupptypesforstats=0 vs =1',fn), AllStats_A_ng.(fn), AllStats_A.(fn)}; %#ok<AGROW>
    % whichstats subset: Ref is the subset run, so only the stats it computed are compared
    pairs(end+1,:)={sprintf('z-identity stats, AllStats %s grouped, whichstats subset vs all',fn), AllStats_A.(fn), AllStats_A_ws.(fn)}; %#ok<AGROW>
    % LifeCycleProfiles
    pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles %s grouped',fn), AgeStats_A.(fn), AgeStats_B.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles rich %s grouped',fn), AgeStats_A.rich.(fn), AgeStats_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles agegroupings %s grouped',fn), AgeStats_A_ag.(fn), AgeStats_B_ag.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles agegroupings rich %s grouped',fn), AgeStats_A_ag.rich.(fn), AgeStats_B_ag.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles %s grouped, lowmemory=1 vs =0',fn), AgeStats_A_lowmem.(fn), AgeStats_A_plain.(fn)}; %#ok<AGROW>
    for ii=1:N_i
        nm=names_A{ii};
        typestr=['type',num2str(ii)];
        richtypestr=['richtype',num2str(ii)];
        pairs(end+1,:)={sprintf('z-identity stats, AllStats %s (type %s)',fn,nm), AllStats_A.(fn).(nm), AllStats_B.(typestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, AllStats rich %s (type %s)',fn,nm), AllStats_A.rich.(fn).(nm), AllStats_B.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, AllStats %s (type %s), groupptypesforstats=0 vs =1',fn,nm), AllStats_A_ng.(fn).(nm), AllStats_A.(fn).(nm)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, AllStats %s (type %s), whichstats subset vs all',fn,nm), AllStats_A.(fn).(nm), AllStats_A_ws.(fn).(nm)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles %s (type %s)',fn,nm), AgeStats_A.(fn).(nm), AgeStats_B.(typestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles rich %s (type %s)',fn,nm), AgeStats_A.rich.(fn).(nm), AgeStats_B.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles agegroupings %s (type %s)',fn,nm), AgeStats_A_ag.(fn).(nm), AgeStats_B_ag.(typestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles %s (type %s), lowmemory=1 vs =0',fn,nm), AgeStats_A_lowmem.(fn).(nm), AgeStats_A_plain.(fn).(nm)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity stats, LifeCycleProfiles %s (type %s), groupptypesforstats=0 vs =1',fn,nm), AgeStats_A_ng.(fn).(nm), AgeStats_A.(fn).(nm)}; %#ok<AGROW>
    end
end

% Compare every statistic of each (label, Test, Ref) pair. Ref supplies the list of stats; per-ptype
% substructures are skipped (they get pairs of their own). NaN is compared as a pattern, the count of
% entries that are NaN in one but not the other, because a NaN printed as the value of a check is
% not read by CoreSummary at all: the check would silently vanish from the count.
for pp=1:size(pairs,1)
    Test=pairs{pp,2};
    Ref=pairs{pp,3};
    statnames=fieldnames(Ref);
    for ss=1:length(statnames)
        sn=statnames{ss};
        if ischar(Ref.(sn)) || (isstruct(Ref.(sn)) && ~strcmp(sn,'MoreInequality'))
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

output=struct();

end
