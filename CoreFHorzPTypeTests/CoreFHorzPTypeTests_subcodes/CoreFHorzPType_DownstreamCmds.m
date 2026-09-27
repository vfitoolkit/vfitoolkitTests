function output=CoreFHorzPType_DownstreamCmds(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: the PType commands that the other parts of this bank never call.
%   EvalFnOnAgentDist_AggVars_FHorz_Case1_PType
%   EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType
%   LifeCycleProfiles_FHorz_Case1_PType
%   PolicyInd2Val_FHorz_PType
%   SimPanelValues_FHorz_Case1_PType   (loose only: a panel with z is random, see part 9 for an exact one)
% Two genuinely different types (per-type kappa_j_pt) with z, against two solves without PType,
% as in part 4. Per-type outputs must equal the solo outputs exactly. Both types have the same
% mewj, so the per-age masses are the same across types and the grouped means (overall and
% age-conditional) are the plain ptypeweights-weighted sums of the solo means; the grouped
% Minimum/Maximum are the min/max over the solo ones.
% Also simoptions.groupptypesforstats=0, which must drop the grouped fields of AggVars and
% LifeCycleProfiles while leaving the per-type fields as they were.

Names_i={'low','high'};
N_i=length(Names_i);

n_d=0;
d_grid=[];

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

ptw=Params.(PTypeDistParamNames{1});

%% PType
[~,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AggVars_PT=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
ValuesOnGrid_PT=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AgeStats_PT=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
PolicyVals_PT=PolicyInd2Val_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);

% groupptypesforstats=0
simoptions_nogroup=simoptions;
simoptions_nogroup.groupptypesforstats=0;
AggVars_PT_ng=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_nogroup);
AgeStats_PT_ng=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions_nogroup);

% Panel
simoptions_panel=simoptions;
simoptions_panel.numbersims=10^4;
SimPanel_PT=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,simoptions_panel);

%% Two separate single-type solves
Policy_solo=struct(); Dist_solo=struct(); AggVars_solo=struct(); ValuesOnGrid_solo=struct(); AgeStats_solo=struct(); PolicyVals_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [~,Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    AggVars_solo.(nm)=EvalFnOnAgentDist_AggVars_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    ValuesOnGrid_solo.(nm)=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeStats_solo.(nm)=LifeCycleProfiles_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    PolicyVals_solo.(nm)=PolicyInd2Val_FHorz(Policy_solo.(nm),n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);
end

%% Per type
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('Downstream cmds, Dist (type %s), this should be zero: %.3e \n',nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
    fprintf('Downstream cmds, PolicyInd2Val (type %s), this should be zero: %.3e \n',nm,max(abs(PolicyVals_PT.(nm)(:)-PolicyVals_solo.(nm)(:))))
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        fprintf('Downstream cmds, AggVars %s (type %s), this should be zero: %.3e \n',fn,nm,abs(AggVars_PT.(fn).(nm).Mean-AggVars_solo.(nm).(fn).Mean))
        fprintf('Downstream cmds, ValuesOnGrid %s (type %s), this should be zero: %.3e \n',fn,nm,max(abs(ValuesOnGrid_PT.(fn).(nm)(:)-ValuesOnGrid_solo.(nm).(fn)(:))))
        fprintf('Downstream cmds, AggVars %s (type %s) groupptypesforstats=0 vs =1, this should be zero: %.3e \n',fn,nm,abs(AggVars_PT_ng.(fn).(nm).Mean-AggVars_PT.(fn).(nm).Mean))
    end
end

%% Grouped
for ff=1:length(FnNames)
    fn=FnNames{ff};
    agg_mean=ptw(1)*AggVars_solo.(Names_i{1}).(fn).Mean+ptw(2)*AggVars_solo.(Names_i{2}).(fn).Mean;
    fprintf('Downstream cmds, AggVars %s grouped Mean, this should be zero: %.3e \n',fn,abs(AggVars_PT.(fn).Mean-agg_mean))
    agg_agemean=ptw(1)*AgeStats_solo.(Names_i{1}).(fn).Mean+ptw(2)*AgeStats_solo.(Names_i{2}).(fn).Mean;
    agg_agemin=min(AgeStats_solo.(Names_i{1}).(fn).Minimum,AgeStats_solo.(Names_i{2}).(fn).Minimum);
    agg_agemax=max(AgeStats_solo.(Names_i{1}).(fn).Maximum,AgeStats_solo.(Names_i{2}).(fn).Maximum);
    fprintf('Downstream cmds, LifeCycleProfiles %s grouped Mean,    this should be zero: %.3e \n',fn,max(abs(AgeStats_PT.(fn).Mean(:)-agg_agemean(:))))
    fprintf('Downstream cmds, LifeCycleProfiles %s grouped Minimum, this should be zero: %.3e \n',fn,max(abs(AgeStats_PT.(fn).Minimum(:)-agg_agemin(:))))
    fprintf('Downstream cmds, LifeCycleProfiles %s grouped Maximum, this should be zero: %.3e \n',fn,max(abs(AgeStats_PT.(fn).Maximum(:)-agg_agemax(:))))
    % groupptypesforstats=0 must not produce the grouped fields
    fprintf('Downstream cmds, AggVars %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',fn,isfield(AggVars_PT_ng.(fn),'Mean'))
    fprintf('Downstream cmds, LifeCycleProfiles %s groupptypesforstats=0 has no grouped Mean, this should be zero: %i \n',fn,isfield(AgeStats_PT_ng.(fn),'Mean'))
end

%% Panel (loose): the age-conditional panel mean against the exact grouped mean
% 10^4 draws, so the sampling error is of order 1e-2; this catches a panel that is wrong
% (types mislabelled, weights ignored, a type missing), not a small inaccuracy.
for ff=1:length(FnNames)
    fn=FnNames{ff};
    fprintf('Downstream cmds, SimPanel %s number of sims minus numbersims, this should be zero: %i \n',fn,size(SimPanel_PT.(fn),2)-simoptions_panel.numbersims)
    fprintf('Downstream cmds, SimPanel %s number of NaN entries, this should be zero: %i \n',fn,sum(isnan(SimPanel_PT.(fn)(:))))
    panelmean=mean(SimPanel_PT.(fn),2);
    fprintf('Downstream cmds, SimPanel %s age-conditional mean vs LifeCycleProfiles, this should be close to zero: %.3e \n',fn,max(abs(gather(panelmean(:))-gather(AgeStats_PT.(fn).Mean(:)))))
end

%% LifeCycleProfiles per type: every statistic against the solo solve; and groupptypesforstats=0 against =1
pairs=cell(0,3);
for ii=1:N_i
    nm=Names_i{ii};
    for ff=1:length(FnNames)
        pairs(end+1,:)={sprintf('Downstream cmds, LifeCycleProfiles %s (type %s)',FnNames{ff},nm), AgeStats_PT.(FnNames{ff}).(nm), AgeStats_solo.(nm).(FnNames{ff})}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('Downstream cmds, LifeCycleProfiles %s (type %s) groupptypesforstats=0 vs =1',FnNames{ff},nm), AgeStats_PT_ng.(FnNames{ff}).(nm), AgeStats_PT.(FnNames{ff}).(nm)}; %#ok<AGROW>
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
