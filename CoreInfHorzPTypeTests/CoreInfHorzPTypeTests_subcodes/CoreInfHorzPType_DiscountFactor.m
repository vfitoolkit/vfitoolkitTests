function output=CoreInfHorzPType_DiscountFactor(n_a,n_z,a_grid,z_grid,pi_z,Params,PTypeDistParamNames,simoptionsPT)
% Test: the discount factor differs by type, two ways, against separate single-type solves of each type (V, Policy,
% StationaryDist, and AllStats every leaf):
%   (a) one name, a parameter with one value per type (beta_pt)
%   (b) DiscountFactorParamNames a structure with a field per type (beta_low, beta_high)
% (InfHorz counterpart of CoreFHorzPTypeTests part 15.)

Names_i={'low','high'};
N_i=2;
n_d=0; d_grid=[];
ReturnFn=@(aprime,a,z,r,w_pt,sigma) ReturnFn_nod_z_noe_nosemiz(aprime,a,z,r,w_pt,sigma);
FnsToEvaluate.assets=@(aprime,a,z) a;
FnsToEvaluate.earnings=@(aprime,a,z,w_pt) w_pt*z;
FnNames=fieldnames(FnsToEvaluate);
Params.beta_pt=[0.94;0.96];
Params.beta_low=0.94; Params.beta_high=0.96;
vfoptions=struct();

for form=1:2
    if form==1
        fstr='(a) beta_pt';
        DiscountFactorParamNames={'beta_pt'};
    else
        fstr='(b) names per type';
        DiscountFactorParamNames=struct();
        DiscountFactorParamNames.low={'beta_low'};
        DiscountFactorParamNames.high={'beta_high'};
    end
    [V_PT,Policy_PT]=ValueFnIter_InfHorz_PType(n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params,DiscountFactorParamNames,vfoptions);
    Dist_PT=StationaryDist_InfHorz_PType(PTypeDistParamNames,Policy_PT,n_d,n_a,n_z,Names_i,pi_z,Params,simoptionsPT);
    AllStats_PT=EvalFnOnAgentDist_AllStats_InfHorz_PType(Dist_PT,Policy_PT,FnsToEvaluate,Params,n_d,n_a,n_z,Names_i,d_grid,a_grid,z_grid,struct());
    for ii=1:N_i
        nm=Names_i{ii};
        Params_ii=Params; Params_ii.w_pt=Params.w_pt(ii); Params_ii.beta=Params.beta_pt(ii);
        [V_ii,Policy_ii]=ValueFnIter_InfHorz(n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,ReturnFn,Params_ii,{'beta'},[],vfoptions);
        Dist_ii=StationaryDist_InfHorz(Policy_ii,n_d,n_a,n_z,pi_z,simoptionsPT,Params_ii,[]);
        AllStats_ii=EvalFnOnAgentDist_AllStats_InfHorz(Dist_ii,Policy_ii,FnsToEvaluate,Params_ii,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,struct());
        fprintf('Discount factor per type %s, V (type %s) vs single-type, this should be zero: %.3e \n',fstr,nm,max(abs(gather(V_PT.(nm)(:))-gather(V_ii(:)))))
        fprintf('Discount factor per type %s, Policy (type %s) vs single-type, this should be zero: %.3e \n',fstr,nm,max(abs(gather(Policy_PT.(nm)(:))-gather(Policy_ii(:)))))
        fprintf('Discount factor per type %s, StationaryDist (type %s) vs single-type, this should be zero: %.3e \n',fstr,nm,max(abs(gather(Dist_PT.(nm)(:))-gather(Dist_ii(:)))))
        maxd=0; mis=0;
        for ff=1:numel(FnNames)
            [p0,v0]=CoreStructLeaves(AllStats_ii.(FnNames{ff})); [p1,v1]=CoreStructLeaves(AllStats_PT.(FnNames{ff}).(nm));
            for k=1:numel(p0)
                [present,idx]=ismember(p0{k},p1);
                if present && isequal(size(v0{k}),size(v1{idx}))
                    x=v0{k}; y=v1{idx}; mis=mis+sum(isnan(x(:))~=isnan(y(:))); d=abs(x(:)-y(:)); d=d(~isnan(d)); maxd=max([maxd;d]);
                else
                    mis=mis+1;
                end
            end
        end
        fprintf('Discount factor per type %s, AllStats (type %s) vs single-type, every leaf, this should be zero: %.3e, mismatches %i \n',fstr,nm,maxd,mis)
    end
end
output=struct();
end
