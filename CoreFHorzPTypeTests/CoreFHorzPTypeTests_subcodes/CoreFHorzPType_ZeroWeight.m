function output=CoreFHorzPType_ZeroWeight(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: a ptype with zero mass (ptypeweight 0). It is the FIRST type, as the panel gives any
% leftover sims to the first types, and it is the one with the higher earnings (1.2*kappa_j), so
% its values reach beyond those of the type that actually exists.
%
% A type of zero mass is not in the population, so every grouped statistic must equal the
% statistics of the other type alone (a solo solve of it): all of AllStats and of
% LifeCycleProfiles, the restricted ones too, AggVars, and the restricted sample masses. In
% particular the grouped Minimum/Maximum and the ends of the QuantileCutoffs must not reach into
% the values that only the zero-mass type holds. (Until 2026-09-27 the grouped Minimum/Maximum
% were overwritten by the min/max over every ptype's own Minimum/Maximum, the zero-mass one
% included; they now skip ptypes with ptweights==0.)
% The per-type stats of the zero-mass type are still its within-type (conditional) stats, so they
% must equal a solo solve of it. The panel must contain no sims of it and no NaN.

Names_i={'ghost','real'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

Params.ptypeweights=[0; 1];
ptw=Params.(PTypeDistParamNames{1});
Params.kappa_j_pt=[1.2*Params.kappa_j; Params.kappa_j]; % ghost has the higher earnings

ReturnFn_PT=@(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j_pt,sigma,agej,Jr,pension);
ReturnFn_NoPT=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);

FnsToEvaluate_PT.assets=@(aprime,a,z) a;
FnsToEvaluate_PT.earnings=@(aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z;
FnsToEvaluate_NoPT.assets=@(aprime,a,z) a;
FnsToEvaluate_NoPT.earnings=@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnNames=fieldnames(FnsToEvaluate_PT);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();
simoptions.conditionalrestrictions.rich=@(aprime,a,z) (a>0.5);

%% PType
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,struct());
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AgeStats_PT=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AggVars_PT=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,struct());
simoptions_panel=struct();
simoptions_panel.numbersims=10^4;
SimPanel_PT=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_panel);

%% Solo solves of each type
Policy_solo=struct(); Dist_solo=struct(); AllStats_solo=struct(); AgeStats_solo=struct(); AggVars_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [~,Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z,Params_ii,struct());
    AllStats_solo.(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeStats_solo.(nm)=LifeCycleProfiles_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AggVars_solo.(nm)=EvalFnOnAgentDist_AggVars_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,struct());
end

%% Scalars
fprintf('Zero-weight ptype, ptweights, this should be zero: %.3e \n',max(abs(StationaryDist_PT.ptweights-ptw)))
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Zero-weight ptype, Pol  (type %s), this should be zero: %.3e \n',nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
    fprintf('Zero-weight ptype, Dist (type %s), this should be zero: %.3e \n',nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
end
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('Zero-weight ptype, AggVars %s grouped Mean, this should be zero: %.3e \n',fn,abs(AggVars_PT.(fn).Mean-AggVars_solo.real.(fn).Mean))
    for ii=1:N_i
        nm=Names_i{ii};
        fprintf('Zero-weight ptype, AggVars %s (type %s), this should be zero: %.3e \n',fn,nm,abs(AggVars_PT.(fn).(nm).Mean-AggVars_solo.(nm).(fn).Mean))
    end
end
fprintf('Zero-weight ptype, AllStats rich RestrictedSampleMass TotalAllPTypes, this should be zero: %.3e \n',abs(AllStats_PT.rich.RestrictedSampleMass.TotalAllPTypes-AllStats_solo.real.rich.RestrictedSampleMass))
x=gather(AgeStats_PT.rich.RestrictedSampleMass.ByAge(:)); y=gather(AgeStats_solo.real.rich.RestrictedSampleMass(:));
dd=abs(x-y); dd=dd(~isnan(dd));
fprintf('Zero-weight ptype, LifeCycleProfiles rich RestrictedSampleMass ByAge, this should be zero: %.3e, NaN-pattern mismatches %i \n',max([dd;0]),sum(isnan(x)~=isnan(y)))

%% Panel: every sim is of the real type, and none is NaN (loose: the age-conditional mean)
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('Zero-weight ptype, SimPanel %s number of sims minus numbersims, this should be zero: %i \n',fn,size(SimPanel_PT.(fn),2)-simoptions_panel.numbersims)
    fprintf('Zero-weight ptype, SimPanel %s number of NaN entries, this should be zero: %i \n',fn,sum(isnan(SimPanel_PT.(fn)(:))))
    panelmean=mean(SimPanel_PT.(fn),2);
    fprintf('Zero-weight ptype, SimPanel %s age-conditional mean vs solo LifeCycleProfiles of the real type, this should be close to zero: %.3e \n',fn,max(abs(gather(panelmean(:))-gather(AgeStats_solo.real.(fn).Mean(:)))))
end

%% Every statistic of each (label, Test, Ref) pair
pairs=cell(0,3);
for ff=1:length(FnNames)
    fn=FnNames{ff};
    pairs(end+1,:)={sprintf('Zero-weight ptype, AllStats %s grouped',fn), AllStats_PT.(fn), AllStats_solo.real.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('Zero-weight ptype, AllStats rich %s grouped',fn), AllStats_PT.rich.(fn), AllStats_solo.real.rich.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('Zero-weight ptype, LifeCycleProfiles %s grouped',fn), AgeStats_PT.(fn), AgeStats_solo.real.(fn)}; %#ok<AGROW>
    pairs(end+1,:)={sprintf('Zero-weight ptype, LifeCycleProfiles rich %s grouped',fn), AgeStats_PT.rich.(fn), AgeStats_solo.real.rich.(fn)}; %#ok<AGROW>
    for ii=1:N_i
        nm=Names_i{ii};
        pairs(end+1,:)={sprintf('Zero-weight ptype, AllStats %s (type %s)',fn,nm), AllStats_PT.(fn).(nm), AllStats_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Zero-weight ptype, LifeCycleProfiles %s (type %s)',fn,nm), AgeStats_PT.(fn).(nm), AgeStats_solo.(nm).(fn)}; %#ok<AGROW>
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
