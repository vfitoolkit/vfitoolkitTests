function output=CoreInfHorzPType_NoZ(n_d,n_a,d_grid,a_grid,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: PType without an exogenous shock (n_z=0), with and without d, and FnsToEvaluate that take parameters (a per-type one,
% w_pt, and a common one, r): each type's V, Policy, StationaryDist, AggVars and AllStats (every leaf) against a single-type solve
% of it, and the grouped AggVars/AllStats Mean against the ptweights-pooling. With no shock the distribution of each type is a
% point mass, so the statistics are degenerate, but the parameters of the FnsToEvaluate must still be found: their names follow
% the (d,aprime,a) inputs, and there is no z input to skip.

Names_i={'low','high'};
N_i=2;
ptw=Params.(PTypeDistParamNames{1});
n_z=0; z_grid=[]; pi_z=[];
n_d_withd=n_d; d_grid_withd=d_grid;
vfoptions=struct();
for cc=1:2
    FnsToEvaluate=struct();
    if cc==1
        mstr='without d';
        n_d=0; d_grid=[];
        ReturnFn=@(aprime,a,r,w_pt,sigma) ReturnFn_nod_noz_noe_nosemiz(aprime,a,r,w_pt,sigma);
        FnsToEvaluate.assets=@(aprime,a) a;
        FnsToEvaluate.income=@(aprime,a,r,w_pt) r*a+w_pt;
    else
        mstr='with d';
        n_d=n_d_withd; d_grid=d_grid_withd;
        ReturnFn=@(d,aprime,a,r,w_pt,sigma,eta,varphi) ReturnFn_d_noz_noe_nosemiz(d,aprime,a,r,w_pt,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a) a;
        FnsToEvaluate.income=@(d,aprime,a,r,w_pt) r*a+w_pt*d;
    end
    FnNames=fieldnames(FnsToEvaluate);
    [V_PT,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
    AggVars_PT=EvalFnOnAgentDist_AggVars_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
    AllStats_PT=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
    m=zeros(numel(FnNames),N_i);
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
        [V_ii,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
        Dist_ii=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptionsPT,Params_ii,[]);
        AggVars_ii=EvalFnOnAgentDist_AggVars_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
        AllStats_ii=EvalFnOnAgentDist_AllStats_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
        fprintf('No z %s, V (type %s) vs single-type, this should be zero: %.3e \n',mstr,nm,max(abs(gather(V_PT.(nm)(:))-gather(V_ii(:)))))
        fprintf('No z %s, Policy (type %s) vs single-type, this should be zero: %.3e \n',mstr,nm,max(abs(gather(Policy_PT.(nm)(:))-gather(Policy_ii(:)))))
        fprintf('No z %s, StationaryDist (type %s) vs single-type, this should be zero: %.3e \n',mstr,nm,max(abs(gather(Dist_PT.(nm)(:))-gather(Dist_ii(:)))))
        maxd=0; mis=0;
        for ff=1:numel(FnNames)
            fn=FnNames{ff};
            fprintf('No z %s, AggVars %s (type %s) vs single-type, this should be zero: %.3e \n',mstr,fn,nm,abs(gather(AggVars_PT.(fn).(nm).Mean)-gather(AggVars_ii.(fn).Mean)))
            [p0,v0]=CoreStructLeaves(AllStats_ii.(fn)); [p1,v1]=CoreStructLeaves(AllStats_PT.(fn).(nm));
            for k=1:numel(p0)
                [present,idx]=ismember(p0{k},p1);
                if present && isequal(size(v0{k}),size(v1{idx}))
                    x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
            m(ff,ii)=gather(AllStats_ii.(fn).Mean);
        end
        fprintf('No z %s, AllStats (type %s) vs single-type, every leaf, this should be zero: %.3e, mismatches %i \n',mstr,nm,maxd,mis)
    end
    for ff=1:numel(FnNames)
        fn=FnNames{ff};
        fprintf('No z %s, AggVars and AllStats %s grouped Mean vs pooled single-type, this should be zero: %.3e \n',mstr,fn,max(abs(gather(AggVars_PT.(fn).Mean)-sum(ptw'.*m(ff,:))),abs(gather(AllStats_PT.(fn).Mean)-sum(ptw'.*m(ff,:)))))
    end
end
output=struct();
end
