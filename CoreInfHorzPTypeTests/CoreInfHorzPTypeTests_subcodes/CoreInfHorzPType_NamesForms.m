function output=CoreInfHorzPType_NamesForms(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,Params,DiscountFactorParamNames,PTypeDistParamNames,simoptionsPT)
% Test: Names_i given as the number of types (N_i=2) or as a cell of names ({'ptype001','ptype002'}, the names the numeric form
% creates) gives identical V, Policy, StationaryDist, AggVars and AllStats, per type and grouped. (InfHorz counterpart of
% CoreFHorzPTypeTests part 1.)

n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
vfoptions=struct();

out=cell(2,1);
for form=1:2
    if form==1
        Names_i=2;
    else
        Names_i={'ptype001','ptype002'};
    end
    [V,Policy]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    Dist=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
    out{form}.V=V; out{form}.Policy=Policy; out{form}.Dist=Dist;
    out{form}.AggVars=EvalFnOnAgentDist_AggVars_InfHorz_PType(Dist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
    out{form}.AllStats=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist,Policy,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
end
objs={'V','Policy','Dist','AggVars','AllStats'};
for oo=1:numel(objs)
    [pa,va]=CoreStructLeaves(out{1}.(objs{oo})); [pb,vb]=CoreStructLeaves(out{2}.(objs{oo}));
    maxd=0; mis=0;
    for k=1:numel(pa)
        [present,idx]=ismember(pa{k},pb);
        if present && isequal(size(va{k}),size(vb{idx}))
            a=va{k}; b=vb{idx}; mis=mis+sum(isnan(a(:))~=isnan(b(:))); d=abs(a(:)-b(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
        else
            mis=mis+1;
        end
    end
    mis=mis+sum(~ismember(pb,pa));
    fprintf('Names forms, %s with N_i=2 vs Names_i={ptype001,ptype002}, every leaf, this should be zero: %.3e, mismatches %i \n',objs{oo},maxd,mis)
end
output=out;
end
