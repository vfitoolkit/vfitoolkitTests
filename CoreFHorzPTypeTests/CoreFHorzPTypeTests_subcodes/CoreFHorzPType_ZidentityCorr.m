function output=CoreFHorzPType_ZidentityCorr(n_a,N_j,a_grid,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,N_i)
% Test: the grouped (pooled-across-ptypes) outputs of the three PType correlation commands
%   EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType
%   EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType
%   EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType
% via the construction of part 2: a PType model with N_i types is the same model as one without PType where z
% encodes the type and pi_z is the identity. Grouping across ptypes is then exactly pooling across z, so the grouped
% PType outputs must equal the plain no-PType outputs of the single-type commands, for every field.
%
% For the AutoCorr command the no-PType model also gives exact references via conditional restrictions on z:
%   type<ii>      z==sigma_pt(ii)              -> the PType per-type outputs of type ii (Mean, StdDeviation, AutoCovariance,
%                                                  AutoCorrelation at every horizon; and the pair fields of type<ii> are the
%                                                  age-j and age-j+k moments of type ii, since every agent of a type is a pair)
%   rich          a>0.5                        -> the PType grouped restricted outputs (every field, pair fields included)
%   richtype<ii>  a>0.5 and z==sigma_pt(ii)    -> the PType per-type restricted outputs (every field)
% and the restricted sample masses. This is the only exact reference for the grouped pair statistics of a restriction
% that does not come from a by-hand pooling.

Params.sigma_pt=[2; 3];
ptw=Params.(PTypeDistParamNames{1});

n_d=0;
d_grid=[];

vfoptions=struct();
timehorizons=[2,3];
hstr={'','_k2','_k3'};
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
simoptions_A.timehorizons=timehorizons;
simoptions_A_ag=simoptions_A;
simoptions_A_ag.agegroupings=1:5:N_j;

[~,Policy_A]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,ReturnFn_A,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_A=StationaryDist_Case1_FHorz_PType(jequaloneDist_A,AgeWeightParamNames,PTypeDistParamNames,Policy_A,n_d,n_a,n_z_A,N_j,N_i,pi_z_A,Params,struct());
names_A=fieldnames(Policy_A);

CovarCorr_A=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeCond_A=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A);
AgeCond_A_ag=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,simoptions_A_ag);
ACP_A=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,simoptions_A);
simoptions_A_lowmem=simoptions_A;
simoptions_A_lowmem.lowmemory=1;
ACP_A_lowmem=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz_PType(StationaryDist_A,Policy_A,FnsToEvaluate_A,Params,n_d,n_a,n_z_A,N_j,N_i,d_grid,a_grid,z_grid_A,pi_z_A,simoptions_A_lowmem);

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

% type<ii> and richtype<ii>: z is sigma_pt, which is 2 or 3, so 2.5 separates the two types
simoptions_B=struct();
simoptions_B.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);
simoptions_B.conditionalrestrictions.type1=@(aprime,a,z) (z<2.5);
simoptions_B.conditionalrestrictions.type2=@(aprime,a,z) (z>2.5);
simoptions_B.conditionalrestrictions.richtype1=@(aprime,a,z) (a>0.5)*(z<2.5);
simoptions_B.conditionalrestrictions.richtype2=@(aprime,a,z) (a>0.5)*(z>2.5);
simoptions_B.timehorizons=timehorizons;
simoptions_B_ag=simoptions_B;
simoptions_B_ag.agegroupings=1:5:N_j;

[~,Policy_B]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,ReturnFn_B,Params,DiscountFactorParamNames,[],vfoptions);
StationaryDist_B=StationaryDist_FHorz_Case1(jequaloneDist_B,AgeWeightParamNames,Policy_B,n_d,n_a,n_z_B,N_j,pi_z_B,Params,struct());

CovarCorr_B=EvalFnOnAgentDist_CrossSectionCovarCorr_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,struct());
AgeCond_B=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,struct());
AgeCond_B_ag=EvalFnOnAgentDist_AgeConditionalStats_CrossSectionCovarCorr_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,simoptions_B_ag);
ACP_B=EvalFnOnAgentDist_AutoCorrTransProbs_FHorz(StationaryDist_B,Policy_B,FnsToEvaluate_B,Params,[],n_d,n_a,n_z_B,N_j,d_grid,a_grid,z_grid_B,pi_z_B,simoptions_B);

%% Restricted sample masses
x=gather(ACP_A.rich.RestrictedSampleMass.ByAge(:)); y=gather(ACP_B.rich.RestrictedSampleMass(:)); dd=abs(x-y); dd=dd(~isnan(dd));
fprintf('z-identity corr, AutoCorr rich RestrictedSampleMass.ByAge, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([dd;0]),sum(isnan(x)~=isnan(y)))
for ii=1:N_i
    fprintf('z-identity corr, AutoCorr rich RestrictedSampleMass (type %s) times ptweight, this should be zero: %.3e \n',names_A{ii},max(abs(ptw(ii)*gather(ACP_A.rich.RestrictedSampleMass.(names_A{ii})(:))-gather(ACP_B.(['richtype',num2str(ii)]).RestrictedSampleMass(:)))))
    fprintf('z-identity corr, AutoCorr type%i RestrictedSampleMass vs ptweight times mewj, this should be zero: %.3e \n',ii,max(abs(gather(ACP_B.(['type',num2str(ii)]).RestrictedSampleMass(:))-ptw(ii)*Params.mewj(:))))
end

%% The pair fields of type<ii> in B are the age-j and age-j+k moments of type ii (every agent of the type is a pair)
for ii=1:N_i
    nm=names_A{ii};
    typestr=['type',num2str(ii)];
    for ff=1:nF
        fn=FnNames{ff};
        for hh=1:length(horizons)
            kk=horizons(hh);
            x=gather(ACP_B.(typestr).(fn).(['PairMean_j',hstr{hh}])(:)); y=gather(ACP_A.(fn).(nm).Mean(1:N_j-kk)'); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('z-identity corr, AutoCorr %s PairMean_j%s vs Mean (type %s), this should be zero: %.3e, NaN-pattern mismatches %i \n',typestr,hstr{hh},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP_B.(typestr).(fn).(['PairMean_jplusk',hstr{hh}])(:)); y=gather(ACP_A.(fn).(nm).Mean(1+kk:N_j)'); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('z-identity corr, AutoCorr %s PairMean_jplusk%s vs Mean (type %s), this should be zero: %.3e, NaN-pattern mismatches %i \n',typestr,hstr{hh},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP_B.(typestr).(fn).(['PairStdDeviation_j',hstr{hh}])(:)); y=gather(ACP_A.(fn).(nm).StdDeviation(1:N_j-kk)'); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('z-identity corr, AutoCorr %s PairStdDeviation_j%s vs StdDeviation (type %s), this should be zero: %.3e, NaN-pattern mismatches %i \n',typestr,hstr{hh},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP_B.(typestr).(fn).(['PairStdDeviation_jplusk',hstr{hh}])(:)); y=gather(ACP_A.(fn).(nm).StdDeviation(1+kk:N_j)'); dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('z-identity corr, AutoCorr %s PairStdDeviation_jplusk%s vs StdDeviation (type %s), this should be zero: %.3e, NaN-pattern mismatches %i \n',typestr,hstr{hh},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
            x=gather(ACP_B.(typestr).(fn).(['PairMass',hstr{hh}])(:)); y=ptw(ii)*Params.mewj(1:N_j-kk)'; dd=abs(x-y); dd=dd(~isnan(dd));
            fprintf('z-identity corr, AutoCorr %s PairMass%s vs ptweight times mewj (type %s), this should be zero: %.3e, NaN-pattern mismatches %i \n',typestr,hstr{hh},nm,max([dd;0]),sum(isnan(x)~=isnan(y)))
        end
    end
end

%% Every field of each (label, Test, Ref) pair
pairs=cell(0,3);
for ff=1:nF
    fn=FnNames{ff};
    pairs(end+1,:)={sprintf('z-identity corr, CrossSection %s grouped',fn), CovarCorr_A.(fn), CovarCorr_B.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AgeCond CrossSection %s grouped',fn), AgeCond_A.(fn), AgeCond_B.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AgeCond CrossSection agegroupings %s grouped',fn), AgeCond_A_ag.(fn), AgeCond_B_ag.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AutoCorr %s grouped',fn), ACP_A.(fn), ACP_B.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AutoCorr rich %s grouped',fn), ACP_A.rich.(fn), ACP_B.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AutoCorr %s grouped, lowmemory=1 vs =0',fn), ACP_A_lowmem.(fn), ACP_A.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('z-identity corr, AutoCorr rich %s grouped, lowmemory=1 vs =0',fn), ACP_A_lowmem.rich.(fn), ACP_A.rich.(fn)}; %#ok<AGROW>
    for ii=1:N_i
        nm=names_A{ii};
        typestr=['type',num2str(ii)];
        richtypestr=['richtype',num2str(ii)];
        % the per-type unrestricted output has no pair fields, so compare the fields it does have (Ref is the PType output)
        pairs(end+1,:)={sprintf('z-identity corr, AutoCorr %s (type %s) vs %s restriction',fn,nm,typestr), ACP_B.(typestr).(fn), ACP_A.(fn).(nm)}; %#ok<AGROW>
        % the PairMass fields of the PType per-type output are within-type shares, those of the richtype<ii> restriction are
        % population shares (ptweight(ii) times the within-type share), so the reference is rescaled before the comparison
        Ref_richtype=ACP_B.(richtypestr).(fn);
        for hh=1:length(horizons)
            Ref_richtype.(['PairMass',hstr{hh}])=Ref_richtype.(['PairMass',hstr{hh}])/ptw(ii);
        end
        pairs(end+1,:)={sprintf('z-identity corr, AutoCorr rich %s (type %s) vs %s restriction',fn,nm,richtypestr), ACP_A.rich.(fn).(nm), Ref_richtype}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity corr, AutoCorr %s (type %s), lowmemory=1 vs =0',fn,nm), ACP_A_lowmem.(fn).(nm), ACP_A.(fn).(nm)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('z-identity corr, AutoCorr rich %s (type %s), lowmemory=1 vs =0',fn,nm), ACP_A_lowmem.rich.(fn).(nm), ACP_A.rich.(fn).(nm)}; %#ok<AGROW>
    end
end
% the matrices
fprintf('z-identity corr, CrossSection grouped CovarianceMatrix, this should be zero: %.3e \n',max(abs(gather(CovarCorr_A.CovarianceMatrix(:))-gather(CovarCorr_B.CovarianceMatrix(:)))))
fprintf('z-identity corr, CrossSection grouped CorrelationMatrix, this should be zero: %.3e \n',max(abs(gather(CovarCorr_A.CorrelationMatrix(:))-gather(CovarCorr_B.CorrelationMatrix(:)))))
x=gather(AgeCond_A.CovarianceMatrix(:)); y=gather(AgeCond_B.CovarianceMatrix(:)); dd=abs(x-y); dd=dd(~isnan(dd));
fprintf('z-identity corr, AgeCond CrossSection grouped CovarianceMatrix, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([dd;0]),sum(isnan(x)~=isnan(y)))
x=gather(AgeCond_A_ag.CorrelationMatrix(:)); y=gather(AgeCond_B_ag.CorrelationMatrix(:)); dd=abs(x-y); dd=dd(~isnan(dd));
fprintf('z-identity corr, AgeCond CrossSection agegroupings grouped CorrelationMatrix, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([dd;0]),sum(isnan(x)~=isnan(y)))

% Compare every field of each (label, Test, Ref) pair. Ref supplies the list of fields; per-ptype substructures and
% the CovarianceWith/CorrelationWith substructures are compared field by field, other substructures skipped. NaN is
% compared as a pattern, the count of entries that are NaN in one but not the other, because a NaN printed as the
% value of a check is not read by CoreSummary at all: the check would silently vanish from the count.
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
