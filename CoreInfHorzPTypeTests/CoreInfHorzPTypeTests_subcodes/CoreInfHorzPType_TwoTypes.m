function output=CoreInfHorzPType_TwoTypes(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT,identical)
% Test: a PType model against separate single-type solves of each type, without d and with d. Per type, V, Policy,
% StationaryDist, PolicyInd2Val, AggVars and every leaf of AllStats must equal the single-type commands exactly; grouped, the
% AggVars Mean and the AllStats Mean/StdDeviation must equal the ptweights-pooling of the single-type outputs (the types each
% have mass one, so the pooled mean is sum_i ptweights(i)*mean_i, the pooled variance adds the between-type part).
% (InfHorz counterpart of CoreFHorzPTypeTests parts 4 and 12.)
% With identical=1 the two types are identical (the same w): then every grouped statistic of AllStats must equal the single-type
% AllStats of the one solve, every leaf. (InfHorz counterpart of CoreFHorzPTypeTests part 3.)

if identical==1
    Params.w_pt=[1;1];
    tstr='identical types';
else
    tstr='two types';
end
Names_i={'low','high'};
N_i=2;
ptw=Params.(PTypeDistParamNames{1});
n_d_withd=n_d; d_grid_withd=d_grid;

output=struct();
for cc=1:2
    FnsToEvaluate=struct();
    if cc==1
        mstr='without d';
        n_d=0; d_grid=[];
        ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
        FnsToEvaluate.assets=@(aprime,a,z) a;
        FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
        FnsToEvaluate.consumption=@(aprime,a,z,r,w_pt) (1+r)*a+w_pt*z-aprime;
    else
        mstr='with d';
        n_d=n_d_withd; d_grid=d_grid_withd;
        ReturnFn=@(d,aprime,a,z,r,w_pt,sigma,eta,varphi) ReturnFn_d_z_noe_nosemiz(d,aprime,a,z,r,w_pt,sigma,eta,varphi);
        FnsToEvaluate.assets=@(d,aprime,a,z) a;
        FnsToEvaluate.earnings=@(d,aprime,a,z,w_pt) w_pt*d*z;
        FnsToEvaluate.hours=@(d,aprime,a,z) d;
    end
    FnNames=fieldnames(FnsToEvaluate);
    vfoptions=struct();

    [V_PT,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
    AggVars_PT=EvalFnOnAgentDist_AggVars_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
    AllStats_PT=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
    PolicyVals_PT=PolicyInd2Val_InfHorz_PType(Policy_PT,n_d,n_a,n_z,d_grid,a_grid,vfoptions);

    m=zeros(numel(FnNames),N_i); s=zeros(numel(FnNames),N_i); a_solo=zeros(numel(FnNames),N_i);
    AllStats_solo=cell(N_i,1);
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii);
        [V_ii,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,DiscountFactorParamNames,[],vfoptions);
        Dist_ii=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptionsPT,Params_ii,[]);
        AggVars_ii=EvalFnOnAgentDist_AggVars_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
        AllStats_solo{ii}=EvalFnOnAgentDist_AllStats_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
        PolicyVals_ii=PolicyInd2Val_InfHorz(Policy_ii,n_d,n_a,n_z,d_grid,a_grid,vfoptions);
        fprintf('%s %s, V (type %s) vs single-type, this should be zero: %.3e \n',tstr,mstr,nm,max(abs(gather(V_PT.(nm)(:))-gather(V_ii(:)))))
        fprintf('%s %s, Policy (type %s) vs single-type, this should be zero: %.3e \n',tstr,mstr,nm,max(abs(gather(Policy_PT.(nm)(:))-gather(Policy_ii(:)))))
        fprintf('%s %s, StationaryDist (type %s) vs single-type, this should be zero: %.3e \n',tstr,mstr,nm,max(abs(gather(Dist_PT.(nm)(:))-gather(Dist_ii(:)))))
        fprintf('%s %s, PolicyInd2Val (type %s) vs single-type, this should be zero: %.3e \n',tstr,mstr,nm,max(abs(gather(PolicyVals_PT.(nm)(:))-gather(PolicyVals_ii(:)))))
        maxd=0; mis=0;
        for ff=1:numel(FnNames)
            fprintf('%s %s, AggVars %s (type %s) vs single-type, this should be zero: %.3e \n',tstr,mstr,FnNames{ff},nm,abs(gather(AggVars_PT.(FnNames{ff}).(nm).Mean)-gather(AggVars_ii.(FnNames{ff}).Mean)))
            [p0,v0]=CoreStructLeaves(AllStats_solo{ii}.(FnNames{ff})); [p1,v1]=CoreStructLeaves(AllStats_PT.(FnNames{ff}).(nm));
            for k=1:numel(p0)
                [present,idx]=ismember(p0{k},p1);
                if present && isequal(size(v0{k}),size(v1{idx}))
                    x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
            mis=mis+sum(~ismember(p1,p0));
            m(ff,ii)=gather(AllStats_solo{ii}.(FnNames{ff}).Mean);
            s(ff,ii)=gather(AllStats_solo{ii}.(FnNames{ff}).StdDeviation);
            a_solo(ff,ii)=gather(AggVars_ii.(FnNames{ff}).Mean);
        end
        fprintf('%s %s, AllStats (type %s) vs single-type, every leaf, this should be zero: %.3e, mismatches %i \n',tstr,mstr,nm,maxd,mis)
    end
    % grouped: the ptweights-pooling of the single-type outputs
    for ff=1:numel(FnNames)
        mG=sum(ptw'.*m(ff,:));
        sG=sqrt(sum(ptw'.*(s(ff,:).^2+(m(ff,:)-mG).^2)));
        fprintf('%s %s, AggVars %s grouped Mean vs pooled single-type, this should be zero: %.3e \n',tstr,mstr,FnNames{ff},abs(gather(AggVars_PT.(FnNames{ff}).Mean)-sum(ptw'.*a_solo(ff,:))))
        fprintf('%s %s, AllStats %s grouped Mean vs pooled single-type, this should be zero: %.3e \n',tstr,mstr,FnNames{ff},abs(gather(AllStats_PT.(FnNames{ff}).Mean)-mG))
        fprintf('%s %s, AllStats %s grouped StdDeviation vs pooled single-type, this should be zero: %.3e \n',tstr,mstr,FnNames{ff},abs(gather(AllStats_PT.(FnNames{ff}).StdDeviation)-sG))
        fprintf('%s %s, AllStats %s grouped Variance vs pooled single-type, this should be zero: %.3e \n',tstr,mstr,FnNames{ff},abs(gather(AllStats_PT.(FnNames{ff}).Variance)-sG^2))
    end
    if identical==1 % identical types: the grouped stats are those of the one type, every leaf
        maxd=0; mis=0;
        for ff=1:numel(FnNames)
            [p0,v0]=CoreStructLeaves(AllStats_solo{1}.(FnNames{ff})); [p1,v1]=CoreStructLeaves(AllStats_PT.(FnNames{ff}));
            for k=1:numel(p0)
                [present,idx]=ismember(p0{k},p1);
                if present && isequal(size(v0{k}),size(v1{idx}))
                    x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
        end
        fprintf('%s %s, AllStats grouped vs the single-type AllStats of the one type, every leaf, this should be zero: %.3e, mismatches %i \n',tstr,mstr,maxd,mis)
    end
    output.(strrep(mstr,' ','')).AllStats=AllStats_PT;
end
end
