function output=CoreInfHorzPType_ZeroWeight(n_a,n_z,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: a ptype with zero mass (ptweight 0). It is the FIRST type, and the one with the higher wage, so its values reach beyond
% those of the type that actually exists. A type of zero mass is not in the population, so every grouped statistic of AllStats
% (Minimum/Maximum and the ends of the QuantileCutoffs included) and the grouped AggVars must equal those of the other type alone
% (a single-type solve of it); the per-type stats of the zero-mass type are still its within-type stats, so they must equal a
% single-type solve of it. (InfHorz counterpart of CoreFHorzPTypeTests part 11.)

Names_i={'ghost','real'};
N_i=2;
n_d=0; d_grid=[];
Params.(PTypeDistParamNames{1})=[0;1];
Params.w_pt=[1.4;1];
ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnNames=fieldnames(FnsToEvaluate);
vfoptions=struct();

[~,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
AggVars_PT=EvalFnOnAgentDist_AggVars_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
AllStats_PT=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
fprintf('Zero-weight ptype, ptweights, this should be zero: %.3e \n',max(abs(gather(Dist_PT.ptweights(:))-[0;1])))

AllStats_solo=cell(N_i,1); AggVars_solo=cell(N_i,1);
for ii=1:N_i
    Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
    [~,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
    Dist_ii=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptionsPT,Params_ii,[]);
    AllStats_solo{ii}=EvalFnOnAgentDist_AllStats_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
    AggVars_solo{ii}=EvalFnOnAgentDist_AggVars_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
end
for tt=0:N_i
    maxd=0; mis=0;
    for ff=1:numel(FnNames)
        fn=FnNames{ff};
        if tt==0 % grouped vs the real type alone
            [p0,v0]=CoreStructLeaves(AllStats_solo{2}.(fn)); [p1,v1]=CoreStructLeaves(AllStats_PT.(fn)); ustr='grouped vs the real type alone';
            fprintf('Zero-weight ptype, AggVars %s grouped vs the real type alone, this should be zero: %.3e \n',fn,abs(gather(AggVars_PT.(fn).Mean)-gather(AggVars_solo{2}.(fn).Mean)))
        else
            [p0,v0]=CoreStructLeaves(AllStats_solo{tt}.(fn)); [p1,v1]=CoreStructLeaves(AllStats_PT.(fn).(Names_i{tt})); ustr=['type ',Names_i{tt},' vs single-type'];
        end
        for k=1:numel(p0)
            [present,idx]=ismember(p0{k},p1);
            if present && isequal(size(v0{k}),size(v1{idx}))
                x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
            else
                mis=mis+1;
            end
        end
    end
    fprintf('Zero-weight ptype, AllStats %s, every leaf, this should be zero: %.3e, mismatches %i \n',ustr,maxd,mis)
end
output=struct();
end
