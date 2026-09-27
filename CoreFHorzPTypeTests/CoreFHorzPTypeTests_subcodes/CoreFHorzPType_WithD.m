function output=CoreFHorzPType_WithD(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames)
% Test: PType with a standard decision variable d (labour supply). Everywhere else in this bank
% n_d=0, bar the binary d2 of the semiz types in the ShockTests. Two genuinely different types
% (per-type kappa_j_pt) against two solves without PType, as in part 4, running every PType
% command: ValueFnIter, ValueFnFromPolicy, StationaryDist, AllStats (every statistic),
% AggVars, LifeCycleProfiles (every statistic), ValuesOnGrid and PolicyInd2Val.
% FnsToEvaluate use d, so the Policy-to-values step for d is exercised in each.

Names_i={'low','high'};
N_i=length(Names_i);

ReturnFn_PT=@(d,aprime,a,z,r,w,kappa_j_pt,sigma,eta,varphi,agej,Jr,pension) ...
    ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j_pt,sigma,eta,varphi,agej,Jr,pension);
ReturnFn_NoPT=@(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension) ...
    ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w,kappa_j,sigma,eta,varphi,agej,Jr,pension);

FnsToEvaluate_PT.hours=@(d,aprime,a,z) d;
FnsToEvaluate_PT.assets=@(d,aprime,a,z) a;
FnsToEvaluate_PT.earnings=@(d,aprime,a,z,w,kappa_j_pt) w*kappa_j_pt*z*d;
FnsToEvaluate_NoPT.hours=@(d,aprime,a,z) d;
FnsToEvaluate_NoPT.assets=@(d,aprime,a,z) a;
FnsToEvaluate_NoPT.earnings=@(d,aprime,a,z,w,kappa_j) w*kappa_j*z*d;
FnNames=fieldnames(FnsToEvaluate_PT);

jequaloneDist=zeros(n_a,n_z,'gpuArray');
jequaloneDist(1,ceil(n_z/2))=1;

vfoptions=struct();
simoptions=struct();

ptw=Params.(PTypeDistParamNames{1});

%% PType
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
V_PT_vfp=ValueFnFromPolicy_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn_PT,Params,DiscountFactorParamNames,vfoptions);
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,N_j,Names_i,pi_z,Params,simoptions);
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AggVars_PT=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
AgeStats_PT=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
ValuesOnGrid_PT=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1_PType(Policy_PT,FnsToEvaluate_PT,Params,n_d,n_a,n_z,N_j,Names_i,d_grid,a_grid,z_grid,simoptions);
PolicyVals_PT=PolicyInd2Val_FHorz_PType(Policy_PT,n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);

%% Two separate single-type solves
V_solo=struct(); Policy_solo=struct(); V_solo_vfp=struct(); Dist_solo=struct(); AllStats_solo=struct(); AggVars_solo=struct(); AgeStats_solo=struct(); ValuesOnGrid_solo=struct(); PolicyVals_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};
    Params_ii=Params;
    Params_ii.kappa_j=Params.kappa_j_pt(ii,:);
    [V_solo.(nm),Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,[],vfoptions);
    V_solo_vfp.(nm)=ValueFnFromPolicy_FHorz(Policy_solo.(nm),n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,pi_z,ReturnFn_NoPT,Params_ii,DiscountFactorParamNames,vfoptions);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist,AgeWeightParamNames,Policy_solo.(nm),n_d,n_a,n_z,N_j,pi_z,Params_ii,simoptions);
    AllStats_solo.(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AggVars_solo.(nm)=EvalFnOnAgentDist_AggVars_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    AgeStats_solo.(nm)=LifeCycleProfiles_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    ValuesOnGrid_solo.(nm)=EvalFnOnAgentDist_ValuesOnGrid_FHorz_Case1(Policy_solo.(nm),FnsToEvaluate_NoPT,Params_ii,[],n_d,n_a,n_z,N_j,d_grid,a_grid,z_grid,simoptions);
    PolicyVals_solo.(nm)=PolicyInd2Val_FHorz(Policy_solo.(nm),n_d,n_a,n_z,N_j,d_grid,a_grid,vfoptions);
end

%% Per type
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('With d, V    (type %s), this should be zero: %.3e \n',nm,max(abs(V_PT.(nm)(:)-V_solo.(nm)(:))))
    fprintf('With d, Pol  (type %s), this should be zero: %.3e \n',nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
    fprintf('With d, Dist (type %s), this should be zero: %.3e \n',nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
    fprintf('With d, VFP  vs V (PType,   type %s), this should be zero: %.3e \n',nm,max(abs(V_PT_vfp.(nm)(:)-V_PT.(nm)(:))))
    fprintf('With d, VFP  vs V (noPType, type %s), this should be zero: %.3e \n',nm,max(abs(V_solo_vfp.(nm)(:)-V_solo.(nm)(:))))
    fprintf('With d, PolicyInd2Val (type %s), this should be zero: %.3e \n',nm,max(abs(PolicyVals_PT.(nm)(:)-PolicyVals_solo.(nm)(:))))
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        fprintf('With d, AggVars %s (type %s), this should be zero: %.3e \n',fn,nm,abs(AggVars_PT.(fn).(nm).Mean-AggVars_solo.(nm).(fn).Mean))
        fprintf('With d, ValuesOnGrid %s (type %s), this should be zero: %.3e \n',fn,nm,max(abs(ValuesOnGrid_PT.(fn).(nm)(:)-ValuesOnGrid_solo.(nm).(fn)(:))))
    end
end

%% Grouped (same mewj for both types, so ptw-weighted sums of the solo means)
for ff=1:length(FnNames)
    fn=FnNames{ff};
    agg_mean=ptw(1)*AllStats_solo.(Names_i{1}).(fn).Mean+ptw(2)*AllStats_solo.(Names_i{2}).(fn).Mean;
    fprintf('With d, AllStats %s grouped Mean, this should be zero: %.3e \n',fn,abs(AllStats_PT.(fn).Mean-agg_mean))
    fprintf('With d, AggVars  %s grouped Mean, this should be zero: %.3e \n',fn,abs(AggVars_PT.(fn).Mean-agg_mean))
    agg_agemean=ptw(1)*AgeStats_solo.(Names_i{1}).(fn).Mean+ptw(2)*AgeStats_solo.(Names_i{2}).(fn).Mean;
    fprintf('With d, LifeCycleProfiles %s grouped Mean, this should be zero: %.3e \n',fn,max(abs(AgeStats_PT.(fn).Mean(:)-agg_agemean(:))))
end

%% Every statistic of each (label, Test, Ref) pair
pairs=cell(0,3);
for ii=1:N_i
    nm=Names_i{ii};
    for ff=1:length(FnNames)
        fn=FnNames{ff};
        pairs(end+1,:)={sprintf('With d, AllStats %s (type %s)',fn,nm), AllStats_PT.(fn).(nm), AllStats_solo.(nm).(fn)}; %#ok<AGROW>
        pairs(end+1,:)={sprintf('With d, LifeCycleProfiles %s (type %s)',fn,nm), AgeStats_PT.(fn).(nm), AgeStats_solo.(nm).(fn)}; %#ok<AGROW>
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
