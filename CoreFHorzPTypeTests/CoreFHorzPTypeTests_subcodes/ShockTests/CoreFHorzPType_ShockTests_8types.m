function output=CoreFHorzPType_ShockTests_8types(n_d,n_a,n_z,n_d_semiz,d_grid_semiz,n_d2_semiz,d2_grid_semiz,N_j,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,AgeWeightParamNames,PTypeDistParamNames,vfoptionsbaseline,simoptionsbaseline)
% Shock-combination cross test for PType. One PType solve with 8 types covering
% every (z, e, semiz) combination should give the same per-type V, Policy and
% StationaryDist as 8 separate non-PType solves. Downstream of that, the per-type
% AllStats (every statistic), AggVars and PolicyInd2Val must equal the solo ones, and
% the grouped AllStats/AggVars must pool all eight types, whose state spaces all
% differ, which is the hardest case for the grouping. The FnsToEvaluate are given per
% type (struct of Names_i) as their leading inputs differ with the state space.
%
% Per-type n_d, d_grid, n_z, z_grid, pi_z, ReturnFn, jequaloneDist and the
% e/semiz pieces of vfoptions/simoptions are all passed as Names_i-keyed
% structures (required when these differ across types).
%
% All semiz types use the same nod1 + binary d2 layout used in CoreFHorzTests.

Names_i={'none','z','e','semiz','ze','zsemiz','esemiz','zesemiz'};
N_i=length(Names_i);

% Arbitrary unequal type weights
Params.ptypeweights=(1:N_i)'/sum(1:N_i);

% Pull the e and semiz pieces out of vfoptionsbaseline for convenience
n_e=vfoptionsbaseline.n_e;
e_grid=vfoptionsbaseline.e_grid;
pi_e=vfoptionsbaseline.pi_e;
n_semiz=vfoptionsbaseline.n_semiz;
semiz_grid=vfoptionsbaseline.semiz_grid;
SemiExoStateFn=vfoptionsbaseline.SemiExoStateFn;

%% Per-type structs for the single PType call
n_d_PT=struct();           d_grid_PT=struct();
n_z_PT=struct();           z_grid_PT=struct();           pi_z_PT=struct();
ReturnFn_PT=struct();      jequaloneDist_PT=struct();

vfopt_PT=struct();         simopt_PT=struct();
vfopt_PT.n_e=struct();     vfopt_PT.e_grid=struct();     vfopt_PT.pi_e=struct();
vfopt_PT.n_semiz=struct(); vfopt_PT.semiz_grid=struct(); vfopt_PT.SemiExoStateFn=struct();
simopt_PT.n_e=struct();    simopt_PT.e_grid=struct();    simopt_PT.pi_e=struct();
simopt_PT.n_semiz=struct();simopt_PT.semiz_grid=struct();simopt_PT.SemiExoStateFn=struct();

% type 1: none
nm='none';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_PT.(nm)=zeros(n_a,1,'gpuArray');
jequaloneDist_PT.(nm)(1)=1;

% type 2: z only
nm='z';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w,kappa_j,sigma,agej,Jr,pension);
jequaloneDist_PT.(nm)=zeros(n_a,n_z,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2))=1;

% type 3: e only
nm='e';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_noz_e_nosemiz(aprime,a,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_e/2))=1;

% type 4: semiz only
nm='semiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_noz_noe_semiz(d2,aprime,a,semiz,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2))=1;

% type 5: z + e
nm='ze';
n_d_PT.(nm)=0; d_grid_PT.(nm)=[];
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension) ...
    ReturnFn_nod_z_e_nosemiz(aprime,a,z,e,r,w,kappa_j,sigma,agej,Jr,pension);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
jequaloneDist_PT.(nm)=zeros(n_a,n_z,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_z/2),ceil(n_e/2))=1;

% type 6: z + semiz
nm='zsemiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_z_noe_semiz(d2,aprime,a,semiz,z,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,n_z,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2),ceil(n_z/2))=1;

% type 7: e + semiz
nm='esemiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=0; z_grid_PT.(nm)=[]; pi_z_PT.(nm)=[];
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_noz_e_semiz(d2,aprime,a,semiz,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2),ceil(n_e/2))=1;

% type 8: z + e + semiz
nm='zesemiz';
n_d_PT.(nm)=n_d2_semiz; d_grid_PT.(nm)=d2_grid_semiz;
n_z_PT.(nm)=n_z; z_grid_PT.(nm)=z_grid; pi_z_PT.(nm)=pi_z;
ReturnFn_PT.(nm)=@(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost) ...
    ReturnFn_nod1_z_e_semiz(d2,aprime,a,semiz,z,e,r,w,kappa_j,sigma,agej,Jr,pension,uempbenefit,searcheffortcost);
vfopt_PT.n_e.(nm)=n_e; vfopt_PT.e_grid.(nm)=e_grid; vfopt_PT.pi_e.(nm)=pi_e;
simopt_PT.n_e.(nm)=n_e;simopt_PT.e_grid.(nm)=e_grid;simopt_PT.pi_e.(nm)=pi_e;
vfopt_PT.n_semiz.(nm)=n_semiz; vfopt_PT.semiz_grid.(nm)=semiz_grid; vfopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
simopt_PT.n_semiz.(nm)=n_semiz;simopt_PT.semiz_grid.(nm)=semiz_grid;simopt_PT.SemiExoStateFn.(nm)=SemiExoStateFn;
jequaloneDist_PT.(nm)=zeros(n_a,n_semiz,n_z,n_e,'gpuArray');
jequaloneDist_PT.(nm)(1,ceil(n_semiz/2),ceil(n_z/2),ceil(n_e/2))=1;

%% FnsToEvaluate, one set per type (the leading inputs differ with the type's state space)
% The same two field names for every type, so the grouped stats pool all eight. income uses the
% value of every shock the type has, so a shock read from the wrong position shows up.
FnsToEvaluate_PT=struct();
FnsToEvaluate_PT.assets.none   =@(aprime,a) a;
FnsToEvaluate_PT.assets.z      =@(aprime,a,z) a;
FnsToEvaluate_PT.assets.e      =@(aprime,a,e) a;
FnsToEvaluate_PT.assets.semiz  =@(d2,aprime,a,semiz) a;
FnsToEvaluate_PT.assets.ze     =@(aprime,a,z,e) a;
FnsToEvaluate_PT.assets.zsemiz =@(d2,aprime,a,semiz,z) a;
FnsToEvaluate_PT.assets.esemiz =@(d2,aprime,a,semiz,e) a;
FnsToEvaluate_PT.assets.zesemiz=@(d2,aprime,a,semiz,z,e) a;
FnsToEvaluate_PT.income.none   =@(aprime,a,w,kappa_j) w*kappa_j;
FnsToEvaluate_PT.income.z      =@(aprime,a,z,w,kappa_j) w*kappa_j*z;
FnsToEvaluate_PT.income.e      =@(aprime,a,e,w,kappa_j) w*kappa_j*e;
FnsToEvaluate_PT.income.semiz  =@(d2,aprime,a,semiz,w,kappa_j) w*kappa_j*semiz;
FnsToEvaluate_PT.income.ze     =@(aprime,a,z,e,w,kappa_j) w*kappa_j*z*e;
FnsToEvaluate_PT.income.zsemiz =@(d2,aprime,a,semiz,z,w,kappa_j) w*kappa_j*semiz*z;
FnsToEvaluate_PT.income.esemiz =@(d2,aprime,a,semiz,e,w,kappa_j) w*kappa_j*semiz*e;
FnsToEvaluate_PT.income.zesemiz=@(d2,aprime,a,semiz,z,e,w,kappa_j) w*kappa_j*semiz*z*e;
FnNames=fieldnames(FnsToEvaluate_PT);

%% PType solve (one call, eight types)
[V_PT,Policy_PT]=ValueFnIter_Case1_FHorz_PType(n_d_PT,n_a,n_z_PT,N_j,Names_i,d_grid_PT,a_grid,z_grid_PT,pi_z_PT,ReturnFn_PT,Params,DiscountFactorParamNames,vfopt_PT);
simopt_PT.d_grid=d_grid_PT;
StationaryDist_PT=StationaryDist_Case1_FHorz_PType(jequaloneDist_PT,AgeWeightParamNames,PTypeDistParamNames,Policy_PT,n_d_PT,n_a,n_z_PT,N_j,Names_i,pi_z_PT,Params,simopt_PT);
AllStats_PT=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d_PT,n_a,n_z_PT,N_j,Names_i,d_grid_PT,a_grid,z_grid_PT,simopt_PT);
AggVars_PT=EvalFnOnAgentDist_AggVars_FHorz_Case1_PType(StationaryDist_PT,Policy_PT,FnsToEvaluate_PT,Params,n_d_PT,n_a,n_z_PT,N_j,Names_i,d_grid_PT,a_grid,z_grid_PT,simopt_PT);
PolicyVals_PT=PolicyInd2Val_FHorz_PType(Policy_PT,n_d_PT,n_a,n_z_PT,N_j,d_grid_PT,a_grid,vfopt_PT);

%% Eight independent single-type solves
V_solo=struct(); Policy_solo=struct(); Dist_solo=struct();
AllStats_solo=struct(); AggVars_solo=struct(); PolicyVals_solo=struct();
for ii=1:N_i
    nm=Names_i{ii};

    % Per-type vfoptions/simoptions for the non-PType call
    vfopt_solo=struct();
    simopt_solo=struct();
    if isfield(vfopt_PT.n_e,nm)
        vfopt_solo.n_e=vfopt_PT.n_e.(nm);
        vfopt_solo.e_grid=vfopt_PT.e_grid.(nm);
        vfopt_solo.pi_e=vfopt_PT.pi_e.(nm);
        simopt_solo.n_e=vfopt_solo.n_e;
        simopt_solo.e_grid=vfopt_solo.e_grid;
        simopt_solo.pi_e=vfopt_solo.pi_e;
    end
    if isfield(vfopt_PT.n_semiz,nm)
        vfopt_solo.n_semiz=vfopt_PT.n_semiz.(nm);
        vfopt_solo.semiz_grid=vfopt_PT.semiz_grid.(nm);
        vfopt_solo.SemiExoStateFn=vfopt_PT.SemiExoStateFn.(nm);
        simopt_solo.n_semiz=vfopt_solo.n_semiz;
        simopt_solo.semiz_grid=vfopt_solo.semiz_grid;
        simopt_solo.SemiExoStateFn=vfopt_solo.SemiExoStateFn;
        simopt_solo.d_grid=d_grid_PT.(nm);
    end

    [V_solo.(nm),Policy_solo.(nm)]=ValueFnIter_Case1_FHorz(n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid,z_grid_PT.(nm),pi_z_PT.(nm),ReturnFn_PT.(nm),Params,DiscountFactorParamNames,[],vfopt_solo);
    Dist_solo.(nm)=StationaryDist_FHorz_Case1(jequaloneDist_PT.(nm),AgeWeightParamNames,Policy_solo.(nm),n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,pi_z_PT.(nm),Params,simopt_solo);

    FnsToEvaluate_solo=struct();
    for ff=1:length(FnNames)
        FnsToEvaluate_solo.(FnNames{ff})=FnsToEvaluate_PT.(FnNames{ff}).(nm);
    end
    AllStats_solo.(nm)=EvalFnOnAgentDist_AllStats_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_solo,Params,[],n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid,z_grid_PT.(nm),simopt_solo);
    AggVars_solo.(nm)=EvalFnOnAgentDist_AggVars_FHorz_Case1(Dist_solo.(nm),Policy_solo.(nm),FnsToEvaluate_solo,Params,[],n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid,z_grid_PT.(nm),simopt_solo);
    PolicyVals_solo.(nm)=PolicyInd2Val_FHorz(Policy_solo.(nm),n_d_PT.(nm),n_a,n_z_PT.(nm),N_j,d_grid_PT.(nm),a_grid,vfopt_solo);
end

%% Compare per type
for ii=1:N_i
    nm=Names_i{ii};
    fprintf('ShockTests 8 types, V    (type %s), this should be zero: %.3e \n',nm,max(abs(V_PT.(nm)(:)-V_solo.(nm)(:))))
    fprintf('ShockTests 8 types, Pol  (type %s), this should be zero: %.3e \n',nm,max(abs(Policy_PT.(nm)(:)-Policy_solo.(nm)(:))))
    fprintf('ShockTests 8 types, Dist (type %s), this should be zero: %.3e \n',nm,max(abs(StationaryDist_PT.(nm)(:)-Dist_solo.(nm)(:))))
    fprintf('ShockTests 8 types, PolicyInd2Val (type %s), this should be zero: %.3e \n',nm,max(abs(PolicyVals_PT.(nm)(:)-PolicyVals_solo.(nm)(:))))
    for ff=1:length(FnNames)
        fprintf('ShockTests 8 types, AggVars %s (type %s), this should be zero: %.3e \n',FnNames{ff},nm,abs(AggVars_PT.(FnNames{ff}).(nm).Mean-AggVars_solo.(nm).(FnNames{ff}).Mean))
    end
end

%% Grouped stats: every type has the same mewj and every FnsToEvaluate is relevant to every type,
% so the grouped Mean is the ptypeweights-weighted sum of the solo means, and the grouped
% Minimum/Maximum are the min/max over the solo ones (all weights are positive).
ptw=Params.ptypeweights;
for ff=1:length(FnNames)
    fn=FnNames{ff};
    agg_mean=0; agg_min=Inf; agg_max=-Inf;
    for ii=1:N_i
        nm=Names_i{ii};
        agg_mean=agg_mean+ptw(ii)*AllStats_solo.(nm).(fn).Mean;
        agg_min=min(agg_min,AllStats_solo.(nm).(fn).Minimum);
        agg_max=max(agg_max,AllStats_solo.(nm).(fn).Maximum);
    end
    fprintf('ShockTests 8 types, AllStats %s grouped Mean,    this should be zero: %.3e \n',fn,abs(AllStats_PT.(fn).Mean-agg_mean))
    fprintf('ShockTests 8 types, AllStats %s grouped Minimum, this should be zero: %.3e \n',fn,abs(AllStats_PT.(fn).Minimum-agg_min))
    fprintf('ShockTests 8 types, AllStats %s grouped Maximum, this should be zero: %.3e \n',fn,abs(AllStats_PT.(fn).Maximum-agg_max))
    fprintf('ShockTests 8 types, AggVars  %s grouped Mean,    this should be zero: %.3e \n',fn,abs(AggVars_PT.(fn).Mean-agg_mean))
end

%% Per-type AllStats: every statistic against the solo solve
pairs=cell(0,3);
for ii=1:N_i
    nm=Names_i{ii};
    for ff=1:length(FnNames)
        pairs(end+1,:)={sprintf('ShockTests 8 types, AllStats %s (type %s)',FnNames{ff},nm), AllStats_PT.(FnNames{ff}).(nm), AllStats_solo.(nm).(FnNames{ff})}; %#ok<AGROW>
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
