function output=CoreFHorzPType_CondlRest_Zidentity(n_a,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test: the conditional restrictions of every PType stats command, via the construction of part 2: a PType model with N_i types
% is the same model as one without PType where z encodes the type and pi_z is the identity. Grouping across ptypes is then
% exactly pooling across z, and restrictions on z pick out the types, so the no-PType model gives exact references:
%   rich          a>0.5                        -> the PType grouped restricted outputs (every field)
%   richtype<ii>  a>0.5 and z==sigma_pt(ii)    -> the PType per-type restricted outputs (every field)
% and the restricted sample masses: those of richtype<ii> (population shares) are ptweight(ii) times the PType per-type ones
% (within-type shares). Commands: AllStats, LifeCycleProfiles (per age and in age groups of five), CrossSectionCovarCorr,
% AgeConditionalStats_CrossSectionCovarCorr (per age and in age groups of five) and AutoCorrTransProbs (horizons 1,2,3), each
% _PType. For AutoCorr this is the only exact reference for the grouped pair statistics of a restriction that does not come
% from a by-hand pooling. (Moved here from parts 10 and 20, whose unrestricted z-identity checks stay there; the CrossSection
% commands are new, as their restrictions are new 2026-10-07.)

Params.sigma_pt=[2; 3];
ptw=Params.(PTypeDistParamNames{1});

n_d=0;
d_grid=[];

vfoptions=struct();
timehorizons=[1,2,3]; % exactly the horizons computed (since 2026-10-14 the command no longer adds horizon 1)
hstr={'','_horizon2','_horizon3'};
horizons=[1,2,3];

%% (A) PType, no z
ReturnFn_A=@(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma_pt,agej,Jr,pension);
n_z_A=0;
z_grid_A=[];
pi_z_A=[];

jequaloneDist_A=zeros(n_a,1,'gpuArray');
jequaloneDist_A(1)=1; % no assets

FnsToEvaluate_A.assets=@(aprime,a) a;
FnsToEvaluate_A.consumption=@(aprime,a,r,w,kappa_j,agej,Jr,pension) (agej<Jr)*((1+r)*a+w*kappa_j-aprime)+(agej>=Jr)*((1+r)*a+pension-aprime);
FnsToEvaluate_A.earnings=@(aprime,a,w,kappa_j) w*kappa_j;
FnNames=fieldnames(FnsToEvaluate_A);
nF=length(FnNames);

simoptions_A=struct();
simoptions_A.conditionalrestrictions.rich=@(aprime,a) (a>0.5);
simoptions_A_ag=simoptions_A;
simoptions_A_ag.agegroupings=1:5:N_j;
simoptions_A_acp=simoptions_A;
simoptions_A_acp.timehorizons=timehorizons;

[~,Policy_A]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,ReturnFn_A,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_A=StationaryDist_Case1_FHorz_PType(jequaloneDist_A,AgeWeightParamNames,PTypeDistParamNames,Policy_A,n_d,n_a,n_z_A,N_j,N_i,pi_z_A,Params,struct());
names_A=fieldnames(Policy_A);

AllStats_A=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeStats_A=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeStats_A_ag=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ag);
CovarCorr_A=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeCond_A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeCond_A_ag=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ag);
ACP_A=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,simoptions_A_acp);

%% (B) No PType, z encodes the type with identity transitions
n_z_B=N_i;
z_grid_B=Params.sigma_pt; % [2; 3]
pi_z_B=eye(N_i,'gpuArray');

ReturnFn_B=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) ReturnFn_zAsSigma(aprime,a,z,r,w,kappa_j,agej,Jr,pension);

jequaloneDist_B=zeros(n_a,n_z_B,'gpuArray');
jequaloneDist_B(1,:)=reshape(ptw,1,n_z_B);

FnsToEvaluate_B.assets=@(aprime,a,z) a;
FnsToEvaluate_B.consumption=@(aprime,a,z,r,w,kappa_j,agej,Jr,pension) (agej<Jr)*((1+r)*a+w*kappa_j-aprime)+(agej>=Jr)*((1+r)*a+pension-aprime);
FnsToEvaluate_B.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j;

% richtype<ii>: z is sigma_pt, which is 2 or 3, so 2.5 separates the two types
simoptions_B=struct();
simoptions_B.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions_B.conditionalrestrictions.richtype1=@(aprime,a,z) (a>0.5)*(z<2.5);
simoptions_B.conditionalrestrictions.richtype2=@(aprime,a,z) (a>0.5)*(z>2.5);
simoptions_B_ag=simoptions_B;
simoptions_B_ag.agegroupings=1:5:N_j;
simoptions_B_acp=simoptions_B;
simoptions_B_acp.timehorizons=timehorizons;

[~,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,n_z_B,N_j,pi_z_B,Params,struct());

AllStats_B=EvalFnOnAgentDist_AllStats_FHorz_Case1(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B);
AgeStats_B=LifeCycleProfiles_FHorz_Case1(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B);
AgeStats_B_ag=LifeCycleProfiles_FHorz_Case1(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B_ag);
CovarCorr_B=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B);
AgeCond_B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B);
AgeCond_B_ag=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B_ag);
ACP_B=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,simoptions_B_acp);

%% Restricted sample masses: grouped (population shares) equal B's rich; per type (within-type shares) times ptweight equal B's richtype<ii>
fprintf('z-identity restricted, AllStats rich RestrictedSampleMass TotalAllPTypes, this should be zero: %.3e \n',abs(AllStats_A.rich.RestrictedSampleMass.TotalAllPTypes-AllStats_B.rich.RestrictedSampleMass))
fprintf('z-identity restricted, CrossSection rich RestrictedSampleMass Total, this should be zero: %.3e \n',abs(gather(CovarCorr_A.rich.RestrictedSampleMass.Total)-gather(CovarCorr_B.rich.RestrictedSampleMass)))
byage={'LifeCycleProfiles',AgeStats_A.rich.RestrictedSampleMass.ByAge,AgeStats_B.rich.RestrictedSampleMass; ...
    'AgeCond CrossSection',AgeCond_A.rich.RestrictedSampleMass.ByAge,AgeCond_B.rich.RestrictedSampleMass; ...
    'AutoCorr',ACP_A.rich.RestrictedSampleMass.ByAge,ACP_B.rich.RestrictedSampleMass};
for cc=1:size(byage,1)
    x=gather(byage{cc,2}(:)); y=gather(byage{cc,3}(:)); dd=abs(x-y); dd=dd(~isnan(dd));
    fprintf('z-identity restricted, %s rich RestrictedSampleMass ByAge, this should be zero: %.3e, NaN-pattern mismatches %i \n',byage{cc,1},max([dd;0]),sum(isnan(x)~=isnan(y)))
end
for ii=1:N_i
    nm=names_A{ii};
    richtypestr=['richtype',num2str(ii)];
    fprintf('z-identity restricted, AllStats rich RestrictedSampleMass (type %s) times ptweight, this should be zero: %.3e \n',nm,abs(ptw(ii)*AllStats_A.rich.RestrictedSampleMass.(nm)-AllStats_B.(richtypestr).RestrictedSampleMass))
    fprintf('z-identity restricted, CrossSection rich RestrictedSampleMass (type %s) times ptweight, this should be zero: %.3e \n',nm,abs(ptw(ii)*gather(CovarCorr_A.rich.RestrictedSampleMass.(nm))-gather(CovarCorr_B.(richtypestr).RestrictedSampleMass)))
    pertype={'LifeCycleProfiles',AgeStats_A.rich.RestrictedSampleMass.(nm),AgeStats_B.(richtypestr).RestrictedSampleMass; ...
        'AgeCond CrossSection',AgeCond_A.rich.RestrictedSampleMass.(nm),AgeCond_B.(richtypestr).RestrictedSampleMass; ...
        'AutoCorr',ACP_A.rich.RestrictedSampleMass.(nm),ACP_B.(richtypestr).RestrictedSampleMass};
    for cc=1:size(pertype,1)
        x=ptw(ii)*gather(pertype{cc,2}(:)); y=gather(pertype{cc,3}(:)); dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('z-identity restricted, %s rich RestrictedSampleMass (type %s) times ptweight, this should be zero: %.3e, NaN-pattern mismatches %i \n',pertype{cc,1},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
    end
end

%% The matrices: grouped vs rich, per type vs richtype<ii>
mats={'CrossSection grouped',CovarCorr_A.rich.CovarianceMatrix,CovarCorr_B.rich.CovarianceMatrix,CovarCorr_A.rich.CorrelationMatrix,CovarCorr_B.rich.CorrelationMatrix; ...
    'AgeCond CrossSection grouped',AgeCond_A.rich.CovarianceMatrix,AgeCond_B.rich.CovarianceMatrix,AgeCond_A.rich.CorrelationMatrix,AgeCond_B.rich.CorrelationMatrix; ...
    'AgeCond CrossSection agegroupings grouped',AgeCond_A_ag.rich.CovarianceMatrix,AgeCond_B_ag.rich.CovarianceMatrix,AgeCond_A_ag.rich.CorrelationMatrix,AgeCond_B_ag.rich.CorrelationMatrix};
for ii=1:N_i
    nm=names_A{ii};
    richtypestr=['richtype',num2str(ii)];
    mats(end+1,:)={sprintf('CrossSection (type %s)',nm),CovarCorr_A.rich.CovarianceMatrix_ptype.(nm),CovarCorr_B.(richtypestr).CovarianceMatrix,CovarCorr_A.rich.CorrelationMatrix_ptype.(nm),CovarCorr_B.(richtypestr).CorrelationMatrix}; %#ok<AGROW>
    mats(end+1,:)={sprintf('AgeCond CrossSection (type %s)',nm),AgeCond_A.rich.CovarianceMatrix_ptype.(nm),AgeCond_B.(richtypestr).CovarianceMatrix,AgeCond_A.rich.CorrelationMatrix_ptype.(nm),AgeCond_B.(richtypestr).CorrelationMatrix}; %#ok<AGROW>
end
for cc=1:size(mats,1)
    x=gather(mats{cc,2}(:)); y=gather(mats{cc,3}(:));
    if numel(x)~=numel(y)
        fprintf('z-identity restricted, %s rich CovarianceMatrix wrongly sized, this should be zero: 1 \n',mats{cc,1})
    else
        dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('z-identity restricted, %s rich CovarianceMatrix, this should be zero: %.3e, NaN-pattern mismatches %i \n',mats{cc,1},max([dd;0]),sum(isnan(x)~=isnan(y)))
    end
    x=gather(mats{cc,4}(:)); y=gather(mats{cc,5}(:));
    if numel(x)~=numel(y)
        fprintf('z-identity restricted, %s rich CorrelationMatrix wrongly sized, this should be zero: 1 \n',mats{cc,1})
    else
        dd=abs(x-y); dd=dd(~isnan(dd));
        fprintf('z-identity restricted, %s rich CorrelationMatrix, this should be zero: %.3e, NaN-pattern mismatches %i \n',mats{cc,1},max([dd;0]),sum(isnan(x)~=isnan(y)))
    end
end

%% Every field of each (label, Test, Ref) pair
pairs=cell(0,3);
for ff=1:nF
    fn=FnNames{ff};
    pairs(end+1,:)={sprintf('z-identity restricted, AllStats rich %s grouped',fn), AllStats_A.rich.(fn), AllStats_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity restricted, LifeCycleProfiles rich %s grouped',fn), AgeStats_A.rich.(fn), AgeStats_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity restricted, LifeCycleProfiles agegroupings rich %s grouped',fn), AgeStats_A_ag.rich.(fn), AgeStats_B_ag.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity restricted, CrossSection rich %s grouped',fn), CovarCorr_A.rich.(fn), CovarCorr_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity restricted, AgeCond CrossSection rich %s grouped',fn), AgeCond_A.rich.(fn), AgeCond_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity restricted, AgeCond CrossSection agegroupings rich %s grouped',fn), AgeCond_A_ag.rich.(fn), AgeCond_B_ag.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity restricted, AutoCorr rich %s grouped',fn), ACP_A.rich.(fn), ACP_B.rich.(fn)}; %#ok<AGROW>
    for ii=1:N_i
        nm=names_A{ii};
        richtypestr=['richtype',num2str(ii)];
        pairs(end+1,:)={sprintf('z-identity restricted, AllStats rich %s (type %s)',fn,nm), AllStats_A.rich.(fn).(nm), AllStats_B.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity restricted, LifeCycleProfiles rich %s (type %s)',fn,nm), AgeStats_A.rich.(fn).(nm), AgeStats_B.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity restricted, LifeCycleProfiles agegroupings rich %s (type %s)',fn,nm), AgeStats_A_ag.rich.(fn).(nm), AgeStats_B_ag.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity restricted, CrossSection rich %s (type %s)',fn,nm), CovarCorr_A.rich.(fn).(nm), CovarCorr_B.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity restricted, AgeCond CrossSection rich %s (type %s)',fn,nm), AgeCond_A.rich.(fn).(nm), AgeCond_B.(richtypestr).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity restricted, AgeCond CrossSection agegroupings rich %s (type %s)',fn,nm), AgeCond_A_ag.rich.(fn).(nm), AgeCond_B_ag.(richtypestr).(fn)}; %#ok<AGROW>
        % the PairMass fields of the PType per-type output are within-type shares, those of the richtype<ii> restriction are
        % population shares (ptweight(ii) times the within-type share), so the reference is rescaled before the comparison
        Ref_richtype=ACP_B.(richtypestr).(fn);
        for hh=1:length(horizons)
            Ref_richtype.(['PairMass',hstr{hh}])=Ref_richtype.(['PairMass',hstr{hh}])/ptw(ii);
        end
        pairs(end+1,:)={sprintf('z-identity restricted, AutoCorr rich %s (type %s)',fn,nm), ACP_A.rich.(fn).(nm), Ref_richtype}; %#ok<AGROW>
    end
end

% Compare every field of each (label, Test, Ref) pair. Ref supplies the list of fields; the MoreInequality, CovarianceWith and
% CorrelationWith substructures are compared field by field, other substructures (per-type ones) skipped. NaN is compared as a
% pattern, the count of entries that are NaN in one but not the other, because a NaN printed as the value of a check is not
% read by CoreSummary at all: the check would silently vanish from the count.
for pp=1:size(pairs,1)
    Test=pairs{pp,2};
    Ref=pairs{pp,3};
    statnames=fieldnames(Ref);
    for ss=1:length(statnames)
        sn=statnames{ss};
        if ischar(Ref.(sn)) || iscell(Ref.(sn)) || (isstruct(Ref.(sn)) && ~any(strcmp(sn,{'MoreInequality','CovarianceWith','CorrelationWith'})))
            continue
        end
        if isstruct(Ref.(sn))
            subnames=fieldnames(Ref.(sn));
            cmpnames=cell(length(subnames),1); xs=cmpnames; ys=cmpnames;
            for s2=1:length(subnames)
                cmpnames{s2}=[sn,'.',subnames{s2}];
                ys{s2}=Ref.(sn).(subnames{s2});
                if isfield(Test,sn) && isfield(Test.(sn),subnames{s2})
                    xs{s2}=Test.(sn).(subnames{s2});
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
